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
#ifndef HELIUM_RPPATH_H
#define HELIUM_RPPATH_H

#include <Helium/config.h>
#include <Helium/molecule.h>
#include <Helium/tie.h>
#include <Helium/distancematrix.h>
#include <Helium/util.h>

#include <vector>
#include <iostream>
#include <algorithm>
#include <cassert>

namespace Helium {

  namespace impl {

    class PIDMatrix
    {
      public:
        typedef std::vector<Index> Path;
        typedef std::vector<Path> Paths;

        PIDMatrix(Size n) : m_paths(n * n), m_dim(n)
        {
        }

        Size dim() const
        {
          return m_dim;
        }

        const std::vector<std::vector<Index> >& operator()(Index i, Index j) const
        {
          return m_paths[index(i, j)];
        }

        std::vector<std::vector<Index> >& operator()(Index i, Index j)
        {
          return m_paths[index(i, j)];
        }

        void changePath(Index i, Index j, Index k)
        {
          Paths paths;
          mergePaths(i, j, k, paths);
          this->operator()(i, j) = paths;
        }

        void appendPath(Index i, Index j, Index k)
        {
          Paths paths;
          mergePaths(i, j, k, paths);
          std::copy(paths.begin(), paths.end(), std::back_inserter(this->operator()(i, j)));
        }

        void appendPath(Index i, Index j, const Paths &paths)
        {
          std::copy(paths.begin(), paths.end(), std::back_inserter(this->operator()(i, j)));
        }

        void mergePaths(Index i, Index j, Index k, Paths &p) const
        {
          const Paths &P_ik = this->operator()(i, k);
          const Paths &P_kj = this->operator()(k, j);
          if (P_ik.empty() || P_kj.empty())
            return;

          for (Paths::const_iterator p1 = P_ik.begin(); p1 != P_ik.end(); ++p1)
            for (Paths::const_iterator p2 = P_kj.begin(); p2 != P_kj.end(); ++p2) {
              p.push_back(*p1);
              std::copy(p2->begin(), p2->end(), std::back_inserter(p.back()));
            }
        }

      private:
        Index index(Index i, Index j) const
        {
          return i * m_dim + j;
        }

        std::vector<std::vector<std::vector<Index> > > m_paths;
        Size m_dim;
    };

    std::ostream& operator<<(std::ostream &os, const PIDMatrix &pid)
    {
      for (Size i = 0; i < pid.dim(); ++i) {
        std::cout << "[ ";
        for (Size j = 0; j < pid.dim(); ++j) {
          std::cout << "[ ";
          for (std::size_t k = 0; k < pid(i, j).size(); ++k) {
            std::cout << "{ ";
            for (std::size_t l = 0; l < pid(i, j)[k].size(); ++l)
              std::cout << pid(i, j)[k][l] << " ";
            std::cout << "} ";
          }
          std::cout << "] ";
        }
        std::cout << "] " << std::endl;
      }

      return os;
    }

    template<typename MoleculeType>
    void makePIDmatrix(const MoleculeType &mol, DistanceMatrix &D,
        PIDMatrix &P1, PIDMatrix &P2)
    {
      // initialize matrices
      FOREACH_BOND (bond, mol, MoleculeType) {
        D(get_index(mol, get_source(mol, *bond)), get_index(mol, get_target(mol, *bond))) = 1;
        D(get_index(mol, get_target(mol, *bond)), get_index(mol, get_source(mol, *bond))) = 1;

        P1(get_index(mol, get_source(mol, *bond)), get_index(mol, get_target(mol, *bond))) =
          PIDMatrix::Paths(1, PIDMatrix::Path(1, get_index(mol, *bond)));
        P1(get_index(mol, get_target(mol, *bond)), get_index(mol, get_source(mol, *bond))) =
          PIDMatrix::Paths(1, PIDMatrix::Path(1, get_index(mol, *bond)));
      }


      for (Size k = 0; k < num_atoms(mol); ++k) {
        for (Size i = 0; i < num_atoms(mol); ++i)
          for (Size j = 0; j < num_atoms(mol); ++j) {
            if (D(i, k) == DistanceMatrix::infinity() || D(k, j) == DistanceMatrix::infinity())
              continue;

            if (D(i, j) > D(i, k) + D(k ,j)) {
              // a new shortest path has been found
              if (D(i, j) == D(i, k) + D(k, j) + 1) {
                // the new shortest path = old shortest path - 1
                P2(i, j) = P1(i, j);
              } else {
                P2(i, j).clear();
              }
              // update distance matrix
              D(i, j) = D(i, k) + D(k, j);
              P1.changePath(i, j, k);
            } else if (D(i, j) == D(i, k) + D(k, j)) {
              // another shortest path found, append it
              P1.appendPath(i, j, k);
            } else if (D(i, j) == D(i, k) + D(k, j) - 1) {
              // newly found path = shortest path + 1, append it
              PIDMatrix::Paths paths;
              P1.mergePaths(i, j, k, paths);
              P2.appendPath(i, j, paths);
            }
          }

      }

    }

    struct CycleCandidate
    {
      CycleCandidate() : Cnum(0), i(0), j(0)
      {
      }

      CycleCandidate(Size Cnum_, Index i_, Index j_)
        : Cnum(Cnum_), i(i_), j(j_)
      {
      }

      bool operator<(const CycleCandidate &other) const
      {
        return Cnum < other.Cnum;
      }

      Size Cnum;
      Index i, j;
    };

    void makeCset(const DistanceMatrix &D, const PIDMatrix &P1, const PIDMatrix &P2,
        std::vector<CycleCandidate> &Cset)
    {
      for (Size i = 0; i < D.dim(); ++i)
        for (Size j = 0; j < D.dim(); ++j) {
          if (D(i, j) == 0 || (P1(i, j).size() == 1 && P2(i, j).size() == 0))
            continue;
          if (P2(i, j).size())
            Cset.push_back(CycleCandidate(2 * (D(i, j) + 0.5), i, j)); // odd cycles
          if (P1(i, j).size() > 1)
            Cset.push_back(CycleCandidate(2 * D(i, j), i, j)); // even cycles
        }

      std::sort(Cset.begin(), Cset.end());

      std::cout << "############### Cset #################" << std::endl;
      for (std::size_t i = 0; i < Cset.size(); ++i) {
        std::cout << "Cnum: " << Cset[i].Cnum << " (" << Cset[i].i << ", " << Cset[i].j << ")" << std::endl;
        if (Cset[i].Cnum % 2)
          std::cout << "    # cycles: " << P1(Cset[i].i, Cset[j].j).size() * P2(Cset[i].i, Cset[j].j).size() << std::endl;
        else
          std::cout << "    # cycles: " << P1(Cset[i].i, Cset[j].j).size() * P2(Cset[i].i, Cset[j].j).size() << std::endl;
      }

    }



    template<typename CycleType>
    bool isRelevant(const std::vector<CycleType> &cycles, const CycleType &cycle)
    {
      std::size_t n = 0;

      // check for identical
      for (std::size_t i = 0; i < cycles.size(); ++i) {
        if (cycle == cycles[i])
          return false;
        if (cycles[i].edges().size() < cycle.edges().size())
          ++n;
      }

      // check if cycle is sum of smaller cycles
      std::vector<std::size_t> indices;
      for (std::size_t i = 0; i < n; ++i)
        indices.push_back(i);

      for (std::size_t size = 2; size <= indices.size(); ++size)
        do {
          CycleType sum;
          for (std::size_t i = 0; i < size; ++i)
            sum += cycles[indices[i]];
          if (cycle == sum)
            return false;
        } while (next_combination(indices.begin(), indices.begin() + size, indices.end()));

      return true;
    }

    bool has_empty_intersection(const PIDMatrix::Path &p1, const PIDMatrix::Path &p2)
    {
      for (std::size_t i = 0; i < p1.size(); ++i)
        for (std::size_t j = 0; j < p2.size(); ++j)
          if (p1[i] == p2[j])
            return false;
      return true;
    }

    template<typename CycleType>
    void findRelevant(const DistanceMatrix &D, const PIDMatrix &P1, const PIDMatrix &P2,
        const std::vector<CycleCandidate> &Cset, std::vector<CycleType> &cycles)
    {

      for (std::vector<CycleCandidate>::const_iterator C = Cset.begin(); C != Cset.end(); ++C) {
        if (C->Cnum % 2) {
          // odd cycle
          for (std::size_t i = 0; i < P2(C->i, C->j).size(); ++i) {
            if (!has_empty_intersection(P1(C->i, C->j)[0], P2(C->i, C->j)[i]))
              continue;

            CycleType cycle;
            for (std::size_t j = 0; j < P1(C->i, C->j)[0].size(); ++j)
              cycle.addEdge(P1(C->i, C->j)[0][j]);
            for (std::size_t j = 0; j < P2(C->i, C->j)[i].size(); ++j)
              cycle.addEdge(P2(C->i, C->j)[i][j]);

            if (isRelevant(cycles, cycle))
              cycles.push_back(cycle);
          }
        } else {
          // even cycle
          for (std::size_t i = 0; i < P1(C->i, C->j).size() - 1; ++i) {
            if (!has_empty_intersection(P1(C->i, C->j)[i], P1(C->i, C->j)[i + 1]))
              continue;

            CycleType cycle;
            for (std::size_t j = 0; j < P1(C->i, C->j)[i].size(); ++j)
              cycle.addEdge(P1(C->i, C->j)[i][j]);
            for (std::size_t j = 0; j < P1(C->i, C->j)[i + 1].size(); ++j)
              cycle.addEdge(P1(C->i, C->j)[i + 1][j]);

            if (isRelevant(cycles, cycle))
              cycles.push_back(cycle);
          }
        }
      }

      std::cout << "############### CYCLES #################" << std::endl;
      for (std::size_t i = 0; i < cycles.size(); ++i)
        std::cout << cycles[i] << std::endl;


    }

  }

  class IncidentVectorCycle
  {
    public:
      IncidentVectorCycle(Size n = 0)
      {
      }

      void addEdge(Index e)
      {
        m_edges.insert(e);
      }

      bool operator==(const IncidentVectorCycle &other) const
      {
        return m_edges == other.m_edges;
      }

      IncidentVectorCycle operator+(const IncidentVectorCycle &other) const
      {
        IncidentVectorCycle result;
        std::set_symmetric_difference(m_edges.begin(), m_edges.end(),
                                      other.m_edges.begin(), other.m_edges.end(),
                                      std::inserter(result.m_edges, result.m_edges.begin()));
        return result;
      }

      IncidentVectorCycle& operator+=(const IncidentVectorCycle &other)
      {
        std::set<Index> result;
        std::set_symmetric_difference(m_edges.begin(), m_edges.end(),
                                      other.m_edges.begin(), other.m_edges.end(),
                                      std::inserter(result, result.begin()));
        m_edges = result;
        return *this;
      }

      const std::set<Index>& edges() const
      {
        return m_edges;
      }

    private:
      friend std::ostream& operator<<(std::ostream &os, const IncidentVectorCycle &cycle);

      std::set<Index> m_edges;
  };

  std::ostream& operator<<(std::ostream &os, const IncidentVectorCycle &cycle)
  {
    for (std::set<Index>::const_iterator i = cycle.m_edges.begin(); i != cycle.m_edges.end(); ++i)
      std::cout << *i << " ";
    return os;
  }

}

#endif
