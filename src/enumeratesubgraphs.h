#ifndef HELIUM_ENUMERATESUBGRAPHS_H
#define HELIUM_ENUMERATESUBGRAPHS_H

#include "molecule.h"
#include "components.h"
#include "substructure.h"
#include "tie.h"
#include "util.h"
//#include "timeout.h"

#include <set>

//#include <boost/functional/hash.hpp>

namespace Helium {

  struct Subgraph
  {
    Subgraph()
    {
    }

    Subgraph(unsigned int numAtoms, unsigned int numBonds) : atoms(numAtoms), bonds(numBonds)
    {
    }

    Subgraph(const std::vector<bool> &atoms_, const std::vector<bool> &bonds_) : atoms(atoms_), bonds(bonds_)
    {
    }

    std::vector<bool> hashable() const
    {
      std::vector<bool> h(atoms);
      std::copy(bonds.begin(), bonds.end(), std::back_inserter(h));
      return h;
    }

    std::vector<bool> atoms;
    std::vector<bool> bonds;
  };



  template<typename MoleculeType>
  bool is_cyclic(MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom,
      std::vector<bool> &visitedAtoms, std::vector<bool> &visitedBonds)
  {
    typedef typename molecule_traits<MoleculeType>::incident_iter incident_iter;

    visitedAtoms[get_index(mol, atom)] = true;

    incident_iter bond, end_bonds;
    tie(bond, end_bonds) = get_bonds(mol, atom);
    for (; bond != end_bonds; ++bond) {
      if (visitedBonds[get_index(mol, *bond)])
        continue;
      visitedBonds[get_index(mol, *bond)] = true;

      if (visitedAtoms[get_index(mol, get_other(mol, *bond, atom))])
        return true;
      if (is_cyclic(mol, get_other(mol, *bond, atom), visitedAtoms, visitedBonds))
        return true;
    }
    return false;
  }


  // replace with check for cyclomatic number????
  template<typename MoleculeType>
  bool is_cyclic(MoleculeType &mol)
  {
    typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;

    std::vector<bool> visitedAtoms(num_atoms(mol));
    std::vector<bool> visitedBonds(num_bonds(mol));

    atom_iter atom, end_atoms;
    tie(atom, end_atoms) = get_atoms(mol);
    for (; atom != end_atoms; ++atom) {
      if (visitedAtoms[get_index(mol, *atom)])
        continue;
      if (is_cyclic(mol, *atom, visitedAtoms, visitedBonds))
        return true;
    }

    return false;
  }


  /**
   * Brute force subgraph enumeration for testing purposes.
   */
  template<typename MoleculeType, typename CallbackType>
  void enumerate_subgraphs_correct(MoleculeType &mol, CallbackType &callback, int maxSize, bool trees = false)
  {
    typedef typename molecule_traits<MoleculeType>::bond_iter bond_iter;

    assert(maxSize >= 0);
    if (maxSize == 0)
      return;

    // generate single atom subgraphs
    for (unsigned int i = 0; i < num_atoms(mol); ++i) {
      std::vector<bool> atoms(num_atoms(mol));
      atoms[i] = true;
      callback(Subgraph(atoms, std::vector<bool>(num_bonds(mol))));
    }

    if (maxSize == 1 || !num_bonds(mol))
      return;

    std::vector<unsigned int> components = connected_bond_components(mol);

    int size = 1;
    bool foundSubgraph;
    while (true) {
      foundSubgraph = false;

      if (size > num_bonds(mol))
        break;

      std::vector<unsigned int> comb;
      for (unsigned int i = 0; i < num_bonds(mol); ++i)
        comb.push_back(i);

      do {
        //std::cout << "comb: " << comb << std::endl;
        // make sure all bonds are in the same component
        unsigned int component = components[comb[0]];
        bool sameComponent = true;
        for (int i = 1; i < size; ++i)
          if (components[comb[i]] != component) {
            sameComponent = false;
            break;
          }

        //std::cout << "same component: " << sameComponent << std::endl;

        if (!sameComponent)
          continue;
        
        /*
        std::cout << "combination: ";
        for (int i = 0; i < size; ++i)
          std::cout << comb[i] << " ";
        std::cout << std::endl;
        */

        std::vector<bool> atoms(num_atoms(mol));
        std::vector<bool> bonds(num_bonds(mol));

        for (int i = 0; i < size; ++i) {
          bonds[comb[i]] = true;
          unsigned int source = get_index(mol, get_source(mol, get_bond(mol, comb[i])));
          unsigned int target = get_index(mol, get_target(mol, get_bond(mol, comb[i])));
          atoms[source] = true;
          atoms[target] = true;
        }
        
        if (std::count(atoms.begin(), atoms.end(), true) > maxSize)
          continue;

        if (size > 1) {
          Substructure<MoleculeType> substruct(mol, atoms, bonds);
          // don't allow cyclic subgraphs when enumerating trees
          if (trees && is_cyclic(substruct))
            continue;
          // make sure the all subgraph bonds are connected
          std::vector<unsigned int> subcomponents = connected_bond_components(substruct);
          if (unique_elements(subcomponents) > 1)
            continue;
        }
 

        callback(Subgraph(atoms, bonds));
        foundSubgraph = true;

      } while (next_combination(comb.begin(), comb.begin() + size, comb.end()));

      if (!foundSubgraph)
        break;
      ++size;
    }
  }


  /*

  // CONTAINS A BUG!
  //
  // this is not a problem though, enumerate_subgraphs_correct is used for
  // testing and the faster enumerate_subgraphs is used in production.

  template<typename T>
  class bloom_filter
  {
    public:
      bloom_filter(std::size_t prime1 = 571, std::size_t prime2 = 1039, std::size_t prime3 = 2011)
          : m_prime1(prime1), m_prime2(prime2), m_prime3(prime3), m_filter(prime3)
      {
      }

      void add(const T &value)
      {
        std::size_t hash = m_hash(value);
        m_filter[hash % m_prime1] = true;
        m_filter[hash % m_prime2] = true;
        m_filter[hash % m_prime3] = true;
        m_set.insert(value);
      }

      bool contains(const T &value)
      {
        std::size_t hash = m_hash(value);
        if (m_filter[hash % m_prime1] && m_filter[hash % m_prime2] && m_filter[hash % m_prime3])
          return true;
        return m_set.find(value) != m_set.end();
      }

      const std::set<T>& set() const
      {
        return m_set;
      }

    private:
      std::size_t m_prime1;
      std::size_t m_prime2;
      std::size_t m_prime3;
      boost::hash<T> m_hash;
      std::vector<bool> m_filter;
      std::set<T> m_set;
  };

  namespace impl {

    // http://dalkescientific.com/writings/diary/archive/2011/01/10/subgraph_enumeration.html
    struct EnumerateSubgraphsSlow
    {
      // Add a new bond to the subgraph when both atoms are already in the subgraph.
      Subgraph new_subgraph_with_bond(const Subgraph &subgraph, std::size_t bond)
      {
        //std::cout << "new_subgraph_with_bond(" << bond << ")" << std::endl;
        assert(!subgraph.bonds[bond]);
        Subgraph newSubgraph(subgraph.atoms, subgraph.bonds);
        newSubgraph.bonds[bond] = true;
        return newSubgraph;
      }

      // Add a new atom to the subgraph, and the bond which connects it to the subgraph.
      Subgraph new_subgraph_with_atom_and_bond(const Subgraph &subgraph, std::size_t atom, std::size_t bond)
      {
        //std::cout << "new_subgraph_with_atom_and_bond(" << atom << ", " << bond << ")" << std::endl;
        assert(!subgraph.atoms[atom]);
        assert(!subgraph.bonds[bond]);
        Subgraph newSubgraph(subgraph.atoms, subgraph.bonds);
        newSubgraph.atoms[atom] = true;
        newSubgraph.bonds[bond] = true;
        return newSubgraph;
      }

      template<typename MoleculeType, typename CallbackType>
      void find_subgraphs(
          MoleculeType &mol,
          CallbackType &callback,
          int maxSize,
          bool trees)
      {
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;
        typedef typename molecule_traits<MoleculeType>::bond_iter bond_iter;

        std::vector<Subgraph> subgraphs;
        bloom_filter<std::vector<bool> > filter;

        // set up the initial set of subgraphs
        atom_iter atom, end_atoms;
        tie(atom, end_atoms) = get_atoms(mol);
        for (; atom != end_atoms; ++atom) {
          // make an initial subgraph with the atom and no bonds
          std::vector<bool> atoms(num_atoms(mol));
          atoms[get_index(mol, *atom)] = true;
          // add the subgraph to the collection
          subgraphs.push_back(Subgraph(atoms, std::vector<bool>(num_bonds(mol))));
          callback(subgraphs.back());
        }

        while (!subgraphs.empty()) {
          // select and remove a subgraph from the collection
          Subgraph subgraph = subgraphs.back();
          subgraphs.pop_back();

          // for each bond in the molecule
          bond_iter bond, end_bonds;
          tie(bond, end_bonds) = get_bonds(mol);
          for (; bond != end_bonds; ++bond) {
            // skip the bond if it is already in the subgraph
            if (subgraph.bonds[get_index(mol, *bond)])
              continue;

            atom_type newAtom = 0; // TODO molecule_traits<MoleculeType>::null_atom();
            // if neither atom is in the subgraph, skip the bond
            if (!subgraph.atoms[get_index(mol, get_source(mol, *bond))]) {
              if (!subgraph.atoms[get_index(mol, get_target(mol, *bond))])
                continue;
              newAtom = get_source(mol, *bond);
            } else
              if (!subgraph.atoms[get_index(mol, get_target(mol, *bond))])
                newAtom = get_target(mol, *bond);

            // if one of the atoms is not in the subgraph and the subgraph is at maximum size, skip the bond
            if (newAtom != 0) { // TODO null_atom();
              if (std::count(subgraph.atoms.begin(), subgraph.atoms.end(), true) == maxSize)
                continue;
              // create the new subgraph with this atom and bond
              subgraphs.push_back(new_subgraph_with_atom_and_bond(subgraph, get_index(mol, newAtom), get_index(mol, *bond)));
              std::vector<bool> hashable(subgraphs.back().hashable());
              if (!filter.contains(hashable)) {
                callback(subgraphs.back());
                filter.add(hashable);
              }
            } else {
              subgraphs.push_back(new_subgraph_with_bond(subgraph, get_index(mol, *bond)));
              std::vector<bool> hashable(subgraphs.back().hashable());
              if (!filter.contains(hashable)) {
                callback(subgraphs.back());
                filter.add(hashable);
              }
            }


          }
        }
      }

    };

  } // namespace impl

  template<typename MoleculeType, typename CallbackType>
  void enumerate_subgraphs_slow(MoleculeType &mol, CallbackType &callback, int maxSize, bool trees = false)
  {
    impl::EnumerateSubgraphsSlow es;
    es.find_subgraphs(mol, callback, maxSize, trees);
  }
  */

  namespace impl {

    typedef std::pair<unsigned int, unsigned int> SubgraphExtension;
    typedef std::vector<SubgraphExtension> SubgraphExtensions;

    struct SubgraphSeed
    {
      SubgraphSeed(const std::vector<bool> &visited_, const Subgraph &subgraph_, const SubgraphExtensions &extensions_)
          : visited(visited_), subgraph(subgraph_), extensions(extensions_)
      {
      }

      std::vector<bool> visited; // visited bonds
      Subgraph subgraph;
      SubgraphExtensions extensions;
    };

    // [(bond index, target atom index)]
    template<typename MoleculeType>
    std::vector<std::pair<unsigned int, unsigned int> > find_extensions(
        MoleculeType &mol,
        const std::vector<bool> &visited,
        const std::vector<bool> &newAtoms, const std::vector<bool> &allAtoms, bool trees)
    {
      typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
      typedef typename molecule_traits<MoleculeType>::incident_iter incident_iter;

      /*
      std::cout << "find_extensions(" << std::endl;
      std::cout << "    visited: " << visited << std::endl;
      std::cout << "    new atoms: " << newAtoms << std::endl;
      std::cout << "    all atoms: " << allAtoms << std::endl;
      std::cout << ")" << std::endl;
      */

      std::vector<std::pair<unsigned int, unsigned int> > extensions; // [(bond index, target atom index)]
      std::set<unsigned int> internalExtensions; // [bond index]

      // find the extensions going out of the new atoms
      for (std::size_t i = 0; i < newAtoms.size(); ++i) {
        if (!newAtoms[i])
          continue;

        incident_iter bond, end_bonds;
        tie(bond, end_bonds) = get_bonds(mol, get_atom(mol, i));
        for (; bond != end_bonds; ++bond) {
          if (visited[get_index(mol, *bond)])
            continue;

          atom_type other = get_other(mol, *bond, get_atom(mol, i));
          if (allAtoms[get_index(mol, other)])
            // this is an unconsidered bond going to another atom in the same graph
            // this bond will appear twice so prevent duplicates
            internalExtensions.insert(get_index(mol, *bond));
          else
            extensions.push_back(std::make_pair(get_index(mol, *bond), get_index(mol, other)));
        }
      }

      //std::cout << "internal: " << internalExtensions << std::endl;
      // add the *unique* internal extensions to the list of extensions
      if (!trees)
        for (std::set<unsigned int>::const_iterator i = internalExtensions.begin(); i != internalExtensions.end(); ++i)
          extensions.push_back(std::make_pair(*i, molecule_traits<MoleculeType>::null_index()));

      return extensions;      
    }

    /**
     * generate all combinations:
     *
     * [A, B, C] -> [A], [B], [C], [A, B], [A, C], [B, C], [A, B, C]
     */
    template<typename T>
    void _all_combinations(const std::vector<T> &extensions, int last, int i,
        std::vector<std::vector<T> > &result)
    {
      //std::cout << "_all_combinations(last: " << last << ", i: " << i << ")" << std::endl;
      if (i == last) {
        result.push_back(std::vector<T>());
        result.push_back(std::vector<T>(1, extensions[i]));
      } else {
        std::vector<std::vector<T> > subcombinations;
        _all_combinations(extensions, last, i + 1, subcombinations);
        for (std::size_t j = 0; j < subcombinations.size(); ++j) {
          result.push_back(subcombinations[j]);
          subcombinations[j].push_back(extensions[i]);
          result.push_back(subcombinations[j]);
        }
      }
    }

    // result = [(set(new atom indices), [extensions])]
    template<typename MoleculeType>
    void all_combinations(const SubgraphExtensions &extensions, int limit, std::vector<std::pair<std::set<unsigned int>, SubgraphExtensions> > &result)
    {
      // generate all 2^(num_extensions)-1 ways to combine the extensions such that
      // there is at least one extension in the combination and no combination has
      // more than limit atoms
      int n = extensions.size();
      assert(n >= 1);
      std::vector<SubgraphExtensions> combinations;
      _all_combinations(extensions, n - 1, 0, combinations);
      if (combinations.empty())
        return;
      for (std::size_t i = 1; i < combinations.size(); ++i) {
        std::set<unsigned int> atoms;
        for (std::size_t j = 0; j < combinations[i].size(); ++j)
          if (combinations[i][j].second != molecule_traits<MoleculeType>::null_index())
            atoms.insert(combinations[i][j].second);
        if (atoms.size() > limit)
          continue;
        result.push_back(std::make_pair(atoms, combinations[i]));      
      }
    }

  }


  template<typename MoleculeType, typename CallbackType>
  void enumerate_subgraphs(MoleculeType &mol, CallbackType &callback, int maxSize, bool trees = false)
  {
    typedef typename molecule_traits<MoleculeType>::bond_iter bond_iter;

    //Timeout timeout(500000); // 5s timeout

    assert(maxSize >= 0);
    if (maxSize == 0)
      return;

    // generate single atom subgraphs
    for (unsigned int i = 0; i < num_atoms(mol); ++i) {
      std::vector<bool> atoms(num_atoms(mol));
      atoms[i] = true;
      callback(Subgraph(atoms, std::vector<bool>(num_bonds(mol))));
    }

    if (maxSize == 1)
      return;

    std::vector<impl::SubgraphSeed> seeds;
    std::vector<bool> visited(num_bonds(mol)); // visited bonds

    // generate the initial seeds
    // seeds[i] starts with bond i and bonds 0-i will not be used to extend the seed
    // for each seed, we also keep track of all possible ways to extend it
    bond_iter bond, end_bonds;
    tie(bond, end_bonds) = get_bonds(mol);
    for (; bond != end_bonds; ++bond) {
      visited[get_index(mol, *bond)] = true;

      Subgraph subgraph(num_atoms(mol), num_bonds(mol));
      subgraph.atoms[get_index(mol, get_source(mol, *bond))] = true;
      subgraph.atoms[get_index(mol, get_target(mol, *bond))] = true;
      subgraph.bonds[get_index(mol, *bond)] = true;

      callback(subgraph);

      impl::SubgraphExtensions extensions = impl::find_extensions(mol, visited, subgraph.atoms, subgraph.atoms, trees);
      //std::cout << "extensions: " << extensions << std::endl;
      if (!extensions.empty())
        seeds.push_back(impl::SubgraphSeed(visited, subgraph, extensions));
    }

    if (maxSize == 2)
      return;

    while (!seeds.empty()) {
      // check timeout
      //timeout.check();

      impl::SubgraphSeed seed = seeds.back();
      seeds.pop_back();

      //std::cout << "Seed:" << std::endl;
      //std::cout << "    atoms: " << seed.subgraph.atoms << std::endl;
      //std::cout << "    bonds: " << seed.subgraph.bonds << std::endl;
      //std::cout << "    visited: " << seed.visited << std::endl;

      std::vector<bool> newVisited(seed.visited);
      // handle all 2^(n-1) ways to expand using these sets of bonds, so there
      // is no need to consider them during any of the future expansions
      for (std::size_t i = 0; i < seed.extensions.size(); ++i)
        newVisited[seed.extensions[i].first] = true;

      int newMaxSize = maxSize - std::count(seed.subgraph.atoms.begin(), seed.subgraph.atoms.end(), true);

      // for each possible extension which is small enough
      std::vector<std::pair<std::set<unsigned int>, impl::SubgraphExtensions> > combinations;
      impl::all_combinations<MoleculeType>(seed.extensions, newMaxSize, combinations);

      //std::cout << "    combinations: " << combinations << std::endl;

      for (std::size_t i = 0; i < combinations.size(); ++i) {
        std::vector<bool> atoms(seed.subgraph.atoms);
        std::vector<bool> newAtoms(seed.subgraph.atoms.size());

        //std::cout << "    new atoms: " << combinations[i].first << std::endl;

        for (std::set<unsigned int>::const_iterator index = combinations[i].first.begin(); index != combinations[i].first.end(); ++index) {
          atoms[*index] = true;
          newAtoms[*index] = true;
        }

        
        std::vector<bool> bonds(seed.subgraph.bonds);
        for (std::size_t j = 0; j < combinations[i].second.size(); ++j)
          bonds[combinations[i].second[j].first] = true;

        Subgraph subgraph(atoms, bonds);
        assert(std::count(subgraph.atoms.begin(), subgraph.atoms.end(), true) <= maxSize);
        

        Substructure<MoleculeType> substruct(mol, atoms, bonds);
        if (trees && is_cyclic(substruct))
          continue;

        callback(subgraph);
        
        // if no new atoms were added, and all ways to expand from the old atoms
        // has been explorered, than there is no other way to expand this seed
        if (combinations[i].first.empty())
          continue;
      
        // start from the new atoms to find additional bonds for further expansion
        impl::SubgraphExtensions extensions = impl::find_extensions(mol, newVisited, newAtoms, atoms, trees);
        //std::cout << "    extensions: " << extensions << std::endl;

        if (!extensions.empty())
          seeds.push_back(impl::SubgraphSeed(newVisited, subgraph, extensions));
      }

    }

  }



}

#endif
