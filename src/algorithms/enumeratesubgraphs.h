/*
 * Copyright (c) 2013, Tim Vandermeersch
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef HELIUM_ENUMERATESUBGRAPHS_H
#define HELIUM_ENUMERATESUBGRAPHS_H

#include <Helium/molecule.h>
#include <Helium/algorithms/components.h>
#include <Helium/util.h>
//#include "timeout.h"

#include <set>
#include <iterator>

//#include <boost/functional/hash.hpp>

namespace Helium {

  /**
   * @file algorithms/enumeratesubgraphs.h
   * @brief Subgraph enumeration.
   */

  /**
   * @struct Subgraph algorithms/enumeratesubgraphs.h <Helium/algorithms/enumeratesubgraph.h>
   * @brief Return type for subgraph enumeration.
   *
   * This object is passed as a parameter to the callback functor.
   */
  struct Subgraph
  {
    /**
     * @brief Default constructor.
     */
    Subgraph()
    {
    }

    /**
     * @brief Constructor.
     *
     * @param numAtoms The number of atoms in the molecule.
     * @param numBonds The number of bonds in the molecule.
     */
    Subgraph(unsigned int numAtoms, unsigned int numBonds) : atoms(numAtoms), bonds(numBonds)
    {
    }

    /**
     * @brief Constructor.
     *
     * @param atoms_ Bitvec specifiying the atoms.
     * @param bonds_ Bitvec specifiying the bonds.
     */
    Subgraph(const std::vector<bool> &atoms_, const std::vector<bool> &bonds_) : atoms(atoms_), bonds(bonds_)
    {
    }

    /**
     * @brief Get a hashable object for the subgraph.
     *
     * @return A hashable object.
     */
    std::vector<bool> hashable() const
    {
      std::vector<bool> h(atoms);
      std::copy(bonds.begin(), bonds.end(), std::back_inserter(h));
      return h;
    }

    std::vector<bool> atoms; //!< The atoms.
    std::vector<bool> bonds; //!< The bonds.
  };

  namespace impl {

    template<typename MoleculeType>
    bool is_cyclic(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom,
        std::vector<bool> &visitedAtoms, std::vector<bool> &visitedBonds)
    {
      typedef typename molecule_traits<MoleculeType>::incident_iter incident_iter;

      visitedAtoms[get_index(mol, atom)] = true;

      FOREACH_INCIDENT (bond, atom, mol) {
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
    bool is_cyclic(const MoleculeType &mol)
    {
      typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;

      std::vector<bool> visitedAtoms(num_atoms(mol));
      std::vector<bool> visitedBonds(num_bonds(mol));

      FOREACH_ATOM (atom, mol) {
        if (visitedAtoms[get_index(mol, *atom)])
          continue;
        if (is_cyclic(mol, *atom, visitedAtoms, visitedBonds))
          return true;
      }

      return false;
    }

  }

  /**
   * Brute force subgraph enumeration for testing purposes.
   */
  template<typename MoleculeType, typename CallbackType>
  void enumerate_subgraphs_correct(const MoleculeType &mol, CallbackType &callback, int maxSize, bool trees = false)
  {
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
          MoleculeType substruct;
          make_substructure(substruct, mol, atoms, bonds);
          // don't allow cyclic subgraphs when enumerating trees
          if (trees && impl::is_cyclic(substruct))
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

        FOREACH_INCIDENT (bond, get_atom(mol, i), mol) {
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


  /**
   * @brief Enumerate subgraphs of a molecule.
   *
   * @param mol The molecule.
   * @param callback The callback functor for reporting the subgraphs.
   * @param maxSize The maximum size of the subgraphs.
   * @param trees When true, only trees (i.e. non-cyclic subgraphs) are enumerated.
   */
  template<typename MoleculeType, typename CallbackType>
  void enumerate_subgraphs(const MoleculeType &mol, CallbackType &callback, int maxSize, bool trees = false)
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
    FOREACH_BOND (bond, mol) {
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


        MoleculeType substruct;
        make_substructure(substruct, mol, atoms, bonds);
        if (trees && impl::is_cyclic(substruct))
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
