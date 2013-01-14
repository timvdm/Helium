#ifndef HELIUM_ENUMERATEPATHS_H
#define HELIUM_ENUMERATEPATHS_H

#include "molecule.h"
#include "tie.h"

#include <vector>
#include <algorithm>

namespace Helium {

  namespace impl {

    template<typename MoleculeType>
    class EnumeratePaths
    {
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::mol_atom_iter mol_atom_iter;
        typedef typename molecule_traits<MoleculeType>::atom_atom_iter atom_atom_iter;

      public:
        EnumeratePaths(MoleculeType *mol, int size) : m_mol(mol), m_size(size)
        {
          m_paths.resize(size);
        }

        std::vector<std::vector<unsigned int> > paths()
        {
          mol_atom_iter atom, end_atoms;
          tie(atom, end_atoms) = get_atoms(m_mol);
          for (; atom != end_atoms; ++atom) {
            std::vector<unsigned int> path;
            enumerate(*atom, path);
          }

          std::vector<std::vector<unsigned int> > result;
          for (std::size_t i = 0; i < m_paths.size(); ++i)
            for (std::size_t j = 0; j < m_paths[i].size(); ++j)
              result.push_back(m_paths[i][j]);

          return result;
        }

      private:
        void enumerate(atom_type atom, std::vector<unsigned int> &path)
        {
          path.push_back(get_index(m_mol, atom));

          addPath(path);

          if (path.size() < m_size) {
            atom_atom_iter nbr, end_nbrs;
            tie(nbr, end_nbrs) = get_nbrs(m_mol, atom);
            for (; nbr != end_nbrs; ++nbr) {
              if (std::find(path.begin(), path.end(), get_index(m_mol, *nbr)) != path.end())
                continue;

              enumerate(*nbr, path);
            }
          }

          path.pop_back();
        }

        void addPath(const std::vector<unsigned int> &path)
        {
          std::vector<std::vector<unsigned int> > &paths = m_paths[path.size() - 1];
          for (std::size_t i = 0; i < paths.size(); ++i) {
            bool forward = true, backward = true;
            for (std::size_t j = 0; j < paths[i].size(); ++j) {
              if (paths[i][j] != path[j])
                forward = false;
              if (paths[i][path.size() - j - 1] != path[j])
                backward = false;
              if (!forward && !backward)
                continue;
            }
            if (forward || backward)
              return;        
          }

          paths.push_back(path);
        }

        MoleculeType *m_mol;
        std::vector<std::vector<std::vector<unsigned int> > > m_paths;
        int m_size;
    };

  } // namespace impl

  template<typename MoleculeType>
  std::vector<std::vector<unsigned int> > enumerate_paths(MoleculeType *mol, int size)
  {
    impl::EnumeratePaths<MoleculeType> ep(mol, size);
    return ep.paths();
  }

}

#endif
