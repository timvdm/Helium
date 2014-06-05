/*
 * Copyright (c) 2014, Tim Vandermeersch
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
#ifndef HELIUM_GTD_H
#define HELIUM_GTD_H

#include <Helium/molecule.h>
#include <Helium/algorithms/floydwarshall.h>

namespace Helium {

  /**
   * @brief Compute the Graph Theoretical Distance.
   *
   * The Graph Theoretical Distance is, for each atom, the number of bonds plus
   * one to the most distant heavy atom (i.e. not hydrogen).
   *
   * @param mol The molecule
   *
   * @return The Graph Theoretical Distance.
   */
  template<typename MoleculeType>
  std::vector<unsigned int> graph_theoretical_distance(const MoleculeType &mol)
  {
    DistanceMatrix D = floyd_warshall(mol);

    std::vector<unsigned int> gtd;
    for (std::size_t i = 0; i < num_atoms(mol); ++i) {
      unsigned int maxDist = 0;
      for (std::size_t j = 0; j < num_atoms(mol); ++j) {
        // do not count distances to hydrogen atoms
        if (is_hydrogen(mol, get_atom(mol, j)))
          continue;

        unsigned int d = D(i, j);
        // skip atoms in other components
        if (d == DistanceMatrix::infinity())
          continue;

        // if distance is greater than current max, update max
        if (d > maxDist)
          maxDist = d;
      }

      gtd.push_back(maxDist + 1);
    }

    return gtd;
  }

}

#endif
