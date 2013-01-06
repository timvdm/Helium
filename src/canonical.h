#ifndef HELIUM_CANONICAL_H
#define HELIUM_CANONICAL_H

#include "invariants.h"

#define DEBUG_CANON 0

namespace Helium {

  namespace impl {

    template<typename MoleculeType, typename T>
    std::vector<unsigned long> canonical_path_code(MoleculeType *mol, const std::vector<T> &path)
    {
      std::vector<unsigned long> result;
      if (path.empty())
        return result;

      for (std::size_t i = 1; i < path.size(); ++i) {
        result.push_back(atom_invariant(mol, get_atom(mol, path[i-1])));
        typename molecule_traits<MoleculeType>::bond_type bond = get_bond(mol, get_atom(mol, path[i-1]), get_atom(mol, path[i]));
        result.push_back(bond_invariant(mol, bond));
      }
      result.push_back(atom_invariant(mol, get_atom(mol, path.back())));

      return result;
    }

  }

  template<typename MoleculeType, typename T>
  std::pair<std::vector<unsigned int>, std::vector<unsigned long> > canonicalize_path(MoleculeType *mol, const std::vector<T> &path)
  {
    std::vector<unsigned long> forwardCode = impl::canonical_path_code(mol, path);
    std::vector<T> pathCopy(path.size());
    for (std::size_t i = 0; i < path.size(); ++i)
      pathCopy[i] = path[path.size() - i - 1];
    std::vector<unsigned long> backwardCode = impl::canonical_path_code(mol, pathCopy);

    if (forwardCode < backwardCode)
      return std::make_pair(path, forwardCode);
    return std::make_pair(pathCopy, backwardCode);
  }


  namespace impl {
  
    template<typename MoleculeType>
    class Canonicalize
    {
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
        typedef typename molecule_traits<MoleculeType>::mol_atom_iter mol_atom_iter;
        typedef typename molecule_traits<MoleculeType>::atom_bond_iter atom_bond_iter;
       public:
        Canonicalize(MoleculeType *mol, const std::vector<unsigned long> &symmetry)
            : m_mol(mol), m_symmetry(symmetry), m_visited(num_bonds(mol))
        {
        }

        void canonicalize()
        {
          // select atom(s) with lowest symmetry class
          mol_atom_iter atom, end_atoms;
          tie(atom, end_atoms) = get_atoms(m_mol);
          for (; atom != end_atoms; ++atom) {
            if (m_symmetry[get_index(m_mol, *atom)])
              continue;
            std::vector<std::pair<bond_type, atom_type> > stack(1, std::make_pair(molecule_traits<MoleculeType>::null_bond(), *atom));
            next(stack);
          }
        }

        const std::vector<unsigned int>& labels() const
        {
          return m_labels;
        }

        const std::vector<unsigned long>& code() const
        {
          return m_code;
        }

      private:
        void createCode()
        {
          std::vector<unsigned long> code;
          //
          // encode graph
          //

          // FROM atoms (encodes spanning tree)
          std::copy(m_from.begin(), m_from.end(), std::back_inserter(code));

          // RING-CLOSURES (encodes ring closures)
          unsigned int numClosures = 0;
          for (std::size_t i = 0; i < m_atoms.size(); ++i) {
            atom_type atom = get_atom(m_mol, m_atoms[i]);
            // still need to sort [1 3] and [1 4]
            std::vector<std::pair<unsigned int, unsigned int> > closures; // [(bond index, other atom index)]

            atom_bond_iter bond, end_bonds;
            tie(bond, end_bonds) = get_bonds(m_mol, atom);
            for (; bond != end_bonds; ++bond)
              // a closure bond is a bond not found while generating the FROM spanning tree.
              if (!m_visited[get_index(m_mol, *bond)]) {
                closures.push_back(std::make_pair(get_index(m_mol, *bond), get_index(m_mol, get_other(m_mol, *bond, atom))));
                m_visited[get_index(m_mol, *bond)] = true;
              }
            
            // do the sorting: [1 3] < [1 4]
            std::sort(closures.begin(), closures.end(), compare_second<unsigned int, unsigned int>());
            numClosures += closures.size();
        
            for (std::size_t j = 0; j < closures.size(); ++j) {
              // add the closure bond to the code
              code.push_back(get_index(m_mol, get_other(m_mol, get_bond(m_mol, closures[j].first), get_atom(m_mol, closures[j].second)))); // canonical source
              code.push_back(closures[j].second); // canonical target
              // add the bond to the list (needed for BOND-TYPES below)
              m_bonds.push_back(closures[j].first);
            }
          }

          //
          // encode ATOM attributes
          //
          for (std::size_t i = 0; i < m_atoms.size(); ++i)
            code.push_back(atom_invariant(m_mol, get_atom(m_mol, m_atoms[i])));

          //
          // encode BOND attributes
          //
          for (std::size_t i = 0; i < m_bonds.size(); ++i)
            code.push_back(bond_invariant(m_mol, get_bond(m_mol, m_bonds[i])));

          // backtrack closure bonds
          for (unsigned int i = 0; i < numClosures; ++i) {
            m_visited[m_bonds.back()] = false;
            m_bonds.pop_back();
          }

          if (DEBUG_CANON)
            std::cout << "code: " << code << std::endl;

          if (m_code.empty() || code < m_code) {
            m_labels = m_atoms;
            m_code = code;
          }
        }

        void next(std::vector<std::pair<bond_type, atom_type> > &stack)
        {
          if (DEBUG_CANON)
            std::cout << "stack: " << stack << std::endl;

          // pop next atom from stack
          atom_type atom;
          bond_type fromBond;
          tie(fromBond, atom) = stack.back();
          stack.pop_back();

          if (DEBUG_CANON)
            std::cout << "next(" << get_index(m_mol, atom) << ")" << std::endl;

          // map the atom
          m_atoms.push_back(get_index(m_mol, atom));
          if (fromBond != molecule_traits<MoleculeType>::null_bond()) {
            m_bonds.push_back(get_index(m_mol, fromBond));
            // add from atom
            m_from.push_back(get_index(m_mol, get_other(m_mol, fromBond, atom)));
            // mark bond as visited
            m_visited[get_index(m_mol, fromBond)] = true;
          }


          if (m_atoms.size() == num_atoms(m_mol)) {
            // found a mapping
            if (DEBUG_CANON)
              std::cout << "mapping: " << m_atoms << ", from: " << m_from << std::endl;
            createCode();
          } else {

            std::vector<std::pair<bond_type, atom_type> > stackCopy(stack);
            std::vector<std::pair<bond_type, atom_type> > bonds;
          

            // append unvisited bonds around atom to stack
            atom_bond_iter bond, end_bonds;
            tie(bond, end_bonds) = get_bonds(m_mol, atom);
            for (; bond != end_bonds; ++bond)
              if (!m_visited[get_index(m_mol, *bond)])
                bonds.push_back(std::make_pair(*bond, get_other(m_mol, *bond, atom)));
          

            // sort the new bonds
            std::sort(bonds.begin(), bonds.end(), compare_first<bond_type, atom_type>());

            // recursive call for each permutation of the new bonds
            bool last = false;
            do {
              stack = stackCopy;
              std::copy(bonds.begin(), bonds.end(), std::back_inserter(stack));
              next(stack);
            } while (!last && (last = std::next_permutation(bonds.begin(), bonds.end(), compare_first<bond_type, atom_type>())));

            // process stack 
            while (!stack.empty())
              next(stack);

          }

          if (DEBUG_CANON)
            std::cout << "backtrack..." << std::endl;

          // backtrack
          m_atoms.pop_back();
          if (fromBond != molecule_traits<MoleculeType>::null_bond()) {
            m_bonds.pop_back();
            m_from.pop_back();
            m_visited[get_index(m_mol, fromBond)] = false;
          }
        }

        MoleculeType *m_mol;
        const std::vector<unsigned long> &m_symmetry;
        std::vector<unsigned int> m_atoms; // canonical atom order
        std::vector<unsigned int> m_bonds; // canonical bond order
        std::vector<unsigned int> m_from; // from atoms
        std::vector<unsigned int> m_visited; // visited bonds

        std::vector<unsigned int> m_labels;
        std::vector<unsigned long> m_code;
    };

  }



  template<typename MoleculeType, typename T>
  std::pair<std::vector<unsigned int>, std::vector<unsigned long> > canonicalize(MoleculeType *mol, const std::vector<T> &symmetry)
  {
    impl::Canonicalize<MoleculeType> can(mol, symmetry);
    can.canonicalize();
    if (DEBUG_CANON)
      std::cout << "labels: " << can.labels() << ", code: " << can.code() << std::endl;
    return std::make_pair(can.labels(), can.code());
  }



}

#endif
