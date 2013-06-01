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

namespace Helium {

/**
 * @mainpage
 *
 * @section Helium
 *
 * @li @ref components
 * @li @ref fingerprints
 * @li @ref canonical
 * @li @ref enumerate
 * @li @ref morganec
 *
 * @section mol_concept The Molecule Concept
 *
 * Helium is C++ header only library and uses generic programming to achieve
 * maximal compatibility with other cheminformatics libraries. Almost all
 * data structures that represent molecules can be used as input for the
 * algorithms in Helium. In order for this to work, the molecule concept is
 * used as an interface between the data structure storing the actual
 * molecule and Helium's algorithms. The molecule concept is similar to the
 * STL iterator concept that acts as an interface between containers and the
 * STL's algorithms. However, since a molecule is complexer than a container of
 * arbitrary objects, the molecule concept is also complexer.
 *
 * @subsection mol_concept_traits The molecule_traits Structure
 *
 *
 *
 * @subsection mol_concept_substruct Substructures
 *
 * The Substructure class can be used to create a view on a molecule's substructure.
 * An instance of the Substructure class can be used anywhere a molecule argument is
 * expected.
 *
 * @section components Connected Components
 *
 * The connected components functions are defined in the components.h header file.
 *
 * @li num_connected_components()
 * @li connected_bond_components()
 * @li connected_atom_components()
 *
 * @section fingerprints Fingerprints
 *
 * The fingerprint functions are defined in the fingerprints.h header file.
 *
 * @li path_fingerprint()
 * @li tree_fingerprint()
 * @li subgraph_fingerprint()
 *
 * There are also various classes and functions for @ref fingerprints_page "storing and indexing fingerprints".
 *
 * @section canonical Canonicalization
 *
 * The canonicalization functions are defined in the canonical.h header file.
 *
 * @li canonicalize()
 *
 * @section enumerate Path, Tree and Subgraph Enumeration
 *
 * Paths in a molecule can be enumerated using the enumerate_paths() function
 * defined in enumeratepaths.h header file.
 *
 * Trees and subgraphs in a molecule can both be enumerated using the
 * enumerate_subgraphs() function defined in the enumeratesubgraphs.h header
 * file.
 *
 * @section morganec Morgan's Extended Connectivities
 *
 * The Morgan extended connectivities can be computed using the
 * extended_connectivities() function defined in the
 * extendedconnectivities.h header file.
 */

}
