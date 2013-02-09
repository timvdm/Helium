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
