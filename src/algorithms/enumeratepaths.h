/**
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
#ifndef HELIUM_ENUMERATEPATHS_H
#define HELIUM_ENUMERATEPATHS_H

#include <Helium/molecule.h>
#include <Helium/tie.h>

#include <vector>
#include <algorithm>

namespace Helium {

  namespace impl {

    /**
     * Internal class for enumerating paths in a molecule.
     */
    template<typename MoleculeType>
    class EnumeratePaths
    {
        typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<MoleculeType>::atom_iter atom_iter;
        typedef typename molecule_traits<MoleculeType>::nbr_iter nbr_iter;

      public:
        /**
         * Constructor.
         *
         * @param mol The molecule.
         * @param size The maximum path size (i.e. number of atoms in the path).
         */
        EnumeratePaths(MoleculeType &mol, int size) : m_mol(mol), m_size(size)
        {
          m_paths.resize(size);
        }

        /**
         * Get the enumerated paths.
         */
        std::vector<std::vector<unsigned int> > paths()
        {
          // enumerate the paths starting from each atom
          atom_iter atom, end_atoms;
          tie(atom, end_atoms) = get_atoms(m_mol);
          for (; atom != end_atoms; ++atom) {
            std::vector<unsigned int> path;
            enumerate(*atom, path);
          }

          // convert the internal stored paths to the right datastructure
          std::vector<std::vector<unsigned int> > result;
          for (std::size_t i = 0; i < m_paths.size(); ++i)
            for (std::size_t j = 0; j < m_paths[i].size(); ++j)
              result.push_back(m_paths[i][j]);

          return result;
        }

      private:
        /**
         * The enumeration function that does the work.
         *
         * @param atom The atom to considernext.
         * @param path The current path.
         */
        void enumerate(atom_type atom, std::vector<unsigned int> &path)
        {
          // add the new atom
          path.push_back(get_index(m_mol, atom));

          // add the path if it is unique
          addPath(path);

          // if the maximum path size is reached the path isn't extended
          if (path.size() < m_size) {
            // create new paths by adding neighbors of atom
            nbr_iter nbr, end_nbrs;
            tie(nbr, end_nbrs) = get_nbrs(m_mol, atom);
            for (; nbr != end_nbrs; ++nbr) {
              if (std::find(path.begin(), path.end(), get_index(m_mol, *nbr)) != path.end())
                continue;

              // recursive call
              enumerate(*nbr, path);
            }
          }

          // backtrack
          path.pop_back();
        }

        /**
         * Add a path to the list of found paths if it is unique.
         */
        void addPath(const std::vector<unsigned int> &path)
        {
          // check if the newly found path is unique, all found paths must be
          // considered and the new path is not unique when it matches an
          // already found path in forward/backward direction
          std::vector<std::vector<unsigned int> > &paths = m_paths[path.size() - 1];
          for (std::size_t i = 0; i < paths.size(); ++i) {
            bool forward = true, backward = true;
            for (std::size_t j = 0; j < paths[i].size(); ++j) {
              // path still matches in forward direction?
              if (paths[i][j] != path[j])
                forward = false;
              // path still matches in backward direction?
              if (paths[i][path.size() - j - 1] != path[j])
                backward = false;
              // if the newly found path does not match the considered path in
              // any direction, continue to the next path
              if (!forward && !backward)
                continue;
            }
            // if the path is not unique, don't add it
            if (forward || backward)
              return;
          }

          // the newly found path is unique, add it to the list
          paths.push_back(path);
        }

        MoleculeType &m_mol; //!< The molecule
        std::vector<std::vector<std::vector<unsigned int> > > m_paths; //!< List of unique found paths, ordered by path size
        int m_size; //!< Maximum path size
    };

  } // namespace impl

  /**
   * Enumerate all paths in a molecule upto a given @p size.
   *
   * @param mol The molecule.
   * @param size The maximum size of the paths (i.e. number of atoms).
   *
   * @return A list containing the paths in the molecule. The paths are also
   *         lists containing the atom indices for the path. The bonds in the
   *         path are between the (i,i+1) atom pairs from the path list.
   */
  template<typename MoleculeType>
  std::vector<std::vector<unsigned int> > enumerate_paths(MoleculeType &mol, int size)
  {
    // enumerate the paths
    impl::EnumeratePaths<MoleculeType> ep(mol, size);
    // return the found paths
    return ep.paths();
  }

}

#endif
