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
#ifndef HELIUM_SMARTCOMPONENTS_H
#define HELIUM_SMARTCOMPONENTS_H

#include <Helium/smartmol.h>
#include <Helium/algorithms/components.h>
#include <Helium/algorithms/dfs.h>
#include <Helium/algorithms/cycles.h>

//@cond DEV

namespace Helium {

  class LazyBondComponents : public SmartAttribute
  {
    public:
      LazyBondComponents(const SmartMol &mol) : SmartAttribute(mol),
          m_valid(false)
      {
      }

      std::string name() const
      {
        return "Helium::bond_components";
      }

      void addAtom(Index index)
      {
        m_valid = false;
      }

      void addBond(Index index)
      {
        m_valid = false;
      }

      void removeAtom(Index index)
      {
        m_valid = false;
      }

      void removeBond(Index index)
      {
        m_valid = false;
      }

      void clear()
      {
        m_valid = false;
      }

      void perceive()
      {
        if (!m_valid) {
          m_components = connected_bond_components(molecule());
          m_valid = true;
        }
      }

      const std::vector<unsigned int>& components()
      {
        perceive();
        return m_components;
      }

      Size numComponents()
      {
        perceive();
        return unique_elements(m_components);
      }

    private:
      bool m_valid;
      std::vector<unsigned int> m_components;
  };

  class LazyAtomComponents : public SmartAttribute
  {
    public:
      LazyAtomComponents(const SmartMol &mol) : SmartAttribute(mol),
          m_valid(false)
      {
      }

      std::string name() const
      {
        return "Helium::atom_components";
      }

      void addAtom(Index index)
      {
        m_valid = false;
      }

      void addBond(Index index)
      {
        m_valid = false;
      }

      void removeAtom(Index index)
      {
        m_valid = false;
      }

      void removeBond(Index index)
      {
        m_valid = false;
      }

      void clear()
      {
        m_valid = false;
      }

      void perceive()
      {
        if (!m_valid) {
          m_components = connected_atom_components(molecule());
          m_valid = true;
        }
      }

      const std::vector<unsigned int>& components()
      {
        perceive();
        return m_components;
      }

      Size numComponents()
      {
        perceive();
        return unique_elements(m_components);
      }

    private:
      bool m_valid;
      std::vector<unsigned int> m_components;
  };

  class DynamicAtomComponents : public SmartAttribute
  {
    public:
      DynamicAtomComponents(const SmartMol &mol) : SmartAttribute(mol),
          m_numComponents(0)
      {
      }

      std::string name() const
      {
        return "Helium::atom_components";
      }

      void addAtom(Index index)
      {
        // adding an atom creates a new component
        unsigned int max = m_components.empty() ? -1 : *std::max_element(m_components.begin(), m_components.end());
        m_components.push_back(max + 1);
        // increment the number of components
        ++m_numComponents;
      }

      void addBond(Index index)
      {
        // get the bond and it's source and target atoms
        SmartMol::bond_type bond = get_bond(molecule(), index);
        SmartMol::atom_type source = get_source(molecule(), bond);
        SmartMol::atom_type target = get_target(molecule(), bond);

        // get the components for the atoms
        unsigned int sourceComponent = m_components[get_index(molecule(), source)];
        unsigned int targetComponent = m_components[get_index(molecule(), target)];

        // if both atoms belong to the same component, nothing has to be done
        if (sourceComponent == targetComponent)
          return;

        // find the minimum and maximum label for the atoms
        unsigned int minComponent = std::min(sourceComponent, targetComponent);
        unsigned int maxComponent = std::max(sourceComponent, targetComponent);

        // join the two components that are connected by the bond:
        // assign the minimum label to all atoms in the component with the
        // maximum label, also decrement all labels greater than the maximum
        // to ensure there are no gaps
        for (std::size_t i = 0; i < m_components.size(); ++i) {
          if (m_components[i] == maxComponent)
            m_components[i] = minComponent;
          else if (m_components[i] > maxComponent)
            --m_components[i];
        }

        // decrement the number of components
        --m_numComponents;
      }

      void removeAtom(Index index)
      {
        SmartMol::atom_type atom = get_atom(molecule(), index);
        unsigned int component = m_components[get_index(molecule(), atom)];

        // get the highest component label
        unsigned int maxComponent = *std::max_element(m_components.begin(), m_components.end());

        // create list of neighbors
        std::set<Index> nbrs;
        for (auto &nbr : get_nbrs(molecule(), u))
          nbrs.insert(get_index(molecule(), nbr));

        // keep track of visited neighbors
        std::set<Index> visitedNbrs;

        // decrement all component labels that are higher than component
        // the new components will be assigned labels greater than maxComponent
        for (std::size_t i = 0; i < m_components.size(); ++i)
          if (m_components[i] > component)
            --m_components[i];

        --m_numComponents;

        // the loop below labels the new components
        FOREACH_INCIDENT (bond, atom, molecule()) {
          SmartMol::atom_type nbr = get_other(molecule(), *bond, atom);
          // continue if the neighbor is already visited, this can happen if
          // the atom is cyclic
          if (visitedNbrs.find(get_index(molecule(), nbr)) != visitedNbrs.end())
            continue;

          // increment the number of components
          ++m_numComponents;

          // mark the atom and it's bonds as visited
          std::vector<bool> visited(num_atoms(molecule()) + num_bonds(molecule()));
          visited[get_index(molecule(), atom)] = true;
          FOREACH_INCIDENT (bond2, atom, molecule())
            visited[num_atoms(molecule()) + get_index(molecule(), *bond2)] = true;

          DFSAtomOrderVisitor<SmartMol> visitor;
          impl::dfs_visit(molecule(), nbr, visitor, visited);

          for (std::size_t i = 0 ; i < visitor.atoms.size(); ++i) {
            Index index2 = visitor.atoms[i];
            // assign new label
            m_components[index2] = maxComponent;
            // if the atom is one of the neighbors, add it to visitedNbrs
            if (nbrs.find(index2) != nbrs.end())
              visitedNbrs.insert(index2);
          }
          ++maxComponent;

          // mark nbr as visited
          visitedNbrs.insert(get_index(molecule(), nbr));
        }

        // finally remove the atom from m_components
        m_components.erase(m_components.begin() + get_index(molecule(), atom));
      }

      void removeBond(Index index)
      {
      }

      void clear()
      {
        m_components.clear();
      }

      const std::vector<unsigned int>& components()
      {
        return m_components;
      }

      Size numComponents()
      {
        return unique_elements(m_components);
      }

      /*
      std::vector<Subgraph<SmartMol> > subgraphs()
      {
      }
      */

    private:
      std::vector<bool> m_cyclicAtoms;
      std::vector<bool> m_cyclicBonds;
      std::vector<unsigned int> m_components;
      Size m_numComponents;
  };

}

//@endcond

#endif
