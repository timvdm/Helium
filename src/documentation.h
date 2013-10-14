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
 * @li @ref mol_concept
 * @li @ref depthfirst
 * @li @ref components
 * @li @ref fingerprints
 * @li @ref canonical
 * @li @ref enumerate
 * @li @ref morganec
 * @li @ref cycles
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
 * The molecule_traits struct contains typedefs and simple functions that define
 * the molecule traits. Again, this is similar to the STL iterator_traits struct.
 * When implementing a new model, the default molecule_traits implementation can
 * be used by ensuring the right typedefs and functions are present in the new
 * molecule class. For making an existing molecule datastructure compatible with
 * without changing it, the molecule_traits struct can be reimplemented for the
 * molecule's type by template specialization.
 *
 * @li molecule_traits<MoleculeType>::atom_type: The atom type.
 * @li molecule_traits<MoleculeType>::bond_type: The bond type.
 * @li molecule_traits<MoleculeType>::atom_iter: The atom iterator type (i.e. return type of get_atoms(const MoleculeType &mol)).
 * @li molecule_traits<MoleculeType>::bond_iter: The bond iterator type (i.e. return type of get_bonds(const MoleculeType &mol)).
 * @li molecule_traits<MoleculeType>::nbr_iter: The neighbor iterator type (i.e. return type of get_nbrs(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)).
 * @li molecule_traits<MoleculeType>::incident_iter: The incident bond iterator type (i.e. return type of get_bonds(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)).
 * @li molecule_traits<MoleculeType>::null_index(): An invalid index for atoms/bonds (used for comparison).
 * @li molecule_traits<MoleculeType>::null_atom(): An invalid atom (used for comparison).
 * @li molecule_traits<MoleculeType>::null_bond(): An invalid bond (used for comparison).
 *
 * To write generic functions, the following examples can be used as guidelines.
 * First a simple function with only a molecule as parameter is shown:
 *
 * @code
 * template<typename MoleculeType>
 * void print_num_atoms(const MoleculeType &mol)
 * {
 *   std::cout << num_atoms(mol) << std::endl;
 * }
 * @endcode
 *
 * When additional types are required inside the function, the molecule_traits
 * struct can be used to add typedefs:
 *
 * @code
 * template<typename MoleculeType>
 * void print_degree_of_first_atom(const MoleculeType &mol)
 * {
 *   // the "typename" is needed since molecule_traits is dependent on the
 *   // MoleculeType template parameter...
 *   typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
 *
 *   atom_type atom = get_atom(mol, 0);
 *
 *   std::cout << get_degree(mol, atom) << std::endl;
 * }
 * @endcode
 *
 * If a bond or atom needs to be passed as argument, there are two ways this can
 * be achieved. The first requires more typing so the second way is recommended:
 *
 * @code
 * // the first and long way
 * template<typename MoleculeType>
 * void print_atom_index(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
 * {
 *   std::cout << get_index(mol, atom) << std::endl;
 * }
 *
 * // the second and easier way
 * template<typename MoleculeType, typename AtomType>
 * void print_atom_index(const MoleculeType &mol, AtomType atom)
 * {
 *   std::cout << get_index(mol, atom) << std::endl;
 * }
 * @endcode
 *
 *
 * @subsection mol_concept_api Molecule API
 *
 * The API for working with molecules consists of free functions to ensure any
 * molecule datastructure can be made compatiable with Helium without modifying
 * the original class. It is recommended to work with these free functions
 * instead of using the molecule datastructures (e.g. HeMol) directly. However,
 * for personal projects this might be easier for less experienced programmers.
 *
 * The functions in the list below operate on the @b molecule alone:
 *
 * @li num_atoms(const MoleculeType &mol): Get the number of atoms.
 * @li num_bonds(const MoleculeType &mol): Get the number of bonds.
 * @li get_atom(const MoleculeType &mol, Index index): Get the atom with the specified index.
 * @li get_bond(const MoleculeType &mol, Index index): Get the bond with the specified index.
 * @li get_atoms(const MoleculeType &mol): Get an iterator pair over all atoms in the molecule.
 * @li get_bonds(const MoleculeType &mol): Get an iterator pair over all atoms in the molecule.
 * @li get_bond(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type source, typename molecule_traits<MoleculeType>::atom_type target): Get the bond between source and target.
 *
 * The last two functions are used for iterating over the atoms/bonds. Although
 * these can be used directly, Helium also provides macros (i.e. FOREACH_ATOM()
 * and FOREACH_BOND() to make iteration easier. The usage of all the functions
 * above is illustrated in the simple example below:
 *
 * @code
 * // examples/molecule1.cpp
 * #include <Helium/molecule.h>
 * #include <Helium/hemol.h> // for HeMol and hemol_from_smiles()
 *
 * #include <iostream>
 *
 * using namespace Helium;
 *
 * template<typename MoleculeType>
 * void print_stuff(const MoleculeType &mol)
 * {
 *   // num_atoms() and num_bonds()
 *   std::cout << "# atoms: " << num_atoms(mol) << std::endl;
 *   std::cout << "# bonds: " << num_bonds(mol) << std::endl;
 *
 *   // get_atom() and get_bond()
 *   std::cout << "degree of 2nd atom: " << get_degree(mol, get_atom(mol, 1)) << std::endl;
 *   std::cout << "order of 3th bond: " << get_order(mol, get_bond(mol, 2)) << std::endl;
 *   std::cout << "the bond C=O has order: " << get_order(mol, get_bond(mol, get_atom(mol, 1), get_atom(mol, 2))) << std::endl;
 *
 *   // iterate over the atoms
 *   FOREACH_ATOM (atom, mol, MoleculeType) {
 *     std::cout << "atom " << get_index(mol, *atom) << " has element " << get_element(mol, *atom) << std::endl;
 *   }
 *
 *   // iterate over the bonds
 *   FOREACH_BOND (bond, mol, MoleculeType) {
 *     std::cout << "bond " << get_index(mol, *bond) << " has order " << get_order(mol, *bond) << std::endl;
 *   }
 * }
 *
 * int main()
 * {
 *   HeMol mol = hemol_from_smiles("CC(=O)C");
 *   print_stuff(mol);
 * }
 * @endcode
 *
 * Once an @b atom has been obtained (e.g. using get_atom(), FOREACH_ATOM(),
 * FOREACH_NBR(), ...), the following set of functions can be used to retrieve
 * information about the atom:
 *
 * @li get_index(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom): Get the atom's index.
 * @li get_element(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom): Get the atom's element.
 * @li get_mass(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom): Get the atom's mass.
 * @li get_degree(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom): Get the atom's degree.
 * @li get_charge(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom): Get the atom's charge.
 * @li num_hydrogens(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom): Get the number of hydrogens attached to this atom.
 * @li get_bonds(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom): Get an iterator pair over the atom's incident bonds.
 * @li get_nbrs(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom): Get an iterator pair over the atom's neighbor atoms.
 *
 * Again, the last two functions return iterator pairs. For these the FOREACH_NBR()
 * and FOREACH_INCIDENT() marcos provided.
 *
 * @code
 * // examples/molecule2.cpp
 * #include <Helium/molecule.h>
 * #include <Helium/hemol.h> // for HeMol and hemol_from_smiles()
 *
 * #include <iostream>
 *
 * using namespace Helium;
 *
 * template<typename MoleculeType>
 * void print_stuff(const MoleculeType &mol)
 * {
 *   // iterate over the atoms
 *   FOREACH_ATOM (atom, mol, MoleculeType) {
 *     std::cout << "atom " << get_index(mol, *atom) << ":" << std::endl;
 *     std::cout << "    element: " << get_element(mol, *atom) << std::endl;
 *     std::cout << "    mass: " << get_mass(mol, *atom) << std::endl;
 *     std::cout << "    charge: " << get_charge(mol, *atom) << std::endl;
 *     std::cout << "    degree: " << get_degree(mol, *atom) << std::endl;
 *     std::cout << "    number of hydrogens: " << num_hydrogens(mol, *atom) << std::endl;
 *   }
 *
 *   // print neighbor indices for atom 1
 *   std::cout << "neighbor indices for atom 1: ";
 *   FOREACH_NBR (nbr, get_atom(mol, 1), mol, MoleculeType)
 *     std::cout << get_index(mol, *nbr) << " ";
 *   std::cout << std::endl;
 *
 *   // print incident bond indices for atom 1
 *   std::cout << "incident bond indices for atom 1: ";
 *   FOREACH_INCIDENT (bond, get_atom(mol, 1), mol, MoleculeType)
 *     std::cout << get_index(mol, *bond) << " ";
 *   std::cout << std::endl;
 * }
 *
 * int main()
 * {
 *   HeMol mol = hemol_from_smiles("CC(=O)C");
 *   print_stuff(mol);
 * }
 * @endcode
 *
 * For @b bonds the following functions are provided.
 *
 * @li get_index(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond): Get the bond's index.
 * @li get_source(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond): Get the bond's source atom.
 * @li get_target(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond): Get the bond's target atom.
 * @li get_order(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond): Get the bond's order.
 * @li get_other(const MoleculeType &mol, typename molecule_traits<MoleculeType>::bond_type bond, typename molecule_traits<MoleculeType>::atom_type atom): Get the other bond atom.
 *
 * Simple example showing usage of the various bond functions:
 *
 * @code
 * // examples/molecule3.cpp
 * #include <Helium/molecule.h>
 * #include <Helium/hemol.h> // for HeMol and hemol_from_smiles()
 *
 * #include <iostream>
 *
 * using namespace Helium;
 *
 * template<typename MoleculeType>
 * void print_stuff(const MoleculeType &mol)
 * {
 *   // iterate over the bonds
 *   FOREACH_BOND (bond, mol, MoleculeType) {
 *     std::cout << "bond " << get_index(mol, *bond) << ":" << std::endl;
 *     std::cout << "    source atom index: " << get_index(mol, get_source(mol, *bond)) << std::endl;
 *     std::cout << "    target atom index: " << get_index(mol, get_target(mol, *bond)) << std::endl;
 *     std::cout << "    source atom index retrieved using get_other: "
 *               << get_index(mol, get_other(mol, *bond, get_target(mol, *bond))) << std::endl;
 *     std::cout << "    order: " << get_order(mol, *bond) << std::endl;
 *   }
 * }
 *
 * int main()
 * {
 *   HeMol mol = hemol_from_smiles("CC(=O)C");
 *   print_stuff(mol);
 * }
 * @endcode
 * 
 * All of the functions above are primitive functions in the sense that each model
 * of the molecule concept has to implement them. Using these primitive functions,
 * a number of useful functions are defined.
 * 
 * @li is_carbon(): Check if an atom is a carbon atom.
 * @li is_nitrogen(): Check if an atom is a nitrogen atom.
 * @li is_oxygen(): Check if an atom is an oxygen atom.
 * @li is_phosphorus(): Check if an atom is an phosphorus atom.
 * @li is_sulfur(): Check if an atom is an sulfur atom.
 * @li get_heavy_degree(): Get the number of heavy atoms (i.e. not hydrogen) attached to an atom.
 * @li get_valence(): Get the valence (i.e. number of attached heavy atoms + implicit/explicit hydrogens) of an atom.
 *
 * @subsection mol_predicates Atom & Bond Predicates
 * 
 * Helium provides atom and bond predicates to make development easier.
 * A predicate in this context is a boolean values function that checks some
 * attributes of an atom or bond and returns the result of the comparison.
 *
 * All atom predicates are derived from the AtomPredicates class. Although all
 * predicate classes can be used directly, this would require a lot of typing
 * since template parameters have to be specified etc. To aid developers in
 * writing code faster, functions are provided for creating each  type of
 * predicate with various comparison methods.
 *  
 * @li ElementPredicate: Compare an atom's element (element_eq_predicate(),
 *     element_lt_predicate(), element_gt_predicate(), element_leq_predicate(),
 *     element_geq_predicate()).
 * @li MassPredicate: Compare an atom's mass (mass_eq_predicate(),
 *     mass_lt_predicate(), mass_gt_predicate(), mass_leq_predicate(),
 *     mass_geq_predicate()).
 * @li ChargePredicate: Compare an atom's charge (charge_eq_predicate(),
 *     charge_lt_predicate(), charge_gt_predicate(), charge_leq_predicate(),
 *     charge_geq_predicate()).
 * @li NumHydrogensPredicate: Compare an atom's number of hydrogens
 *     (num_hydrogens_eq_predicate(), num_hydrogens_lt_predicate(),
 *     num_hydrogens_gt_predicate(), num_hydrogens_leq_predicate(),
 *     num_hydrogens_geq_predicate()).
 * @li DegreePredicate: Compare an atom's degree (degree_eq_predicate(),
 *     degree_lt_predicate(), degree_gt_predicate(), degree_leq_predicate(),
 *     degree_geq_predicate()).
 * @li AromaticAtomPredicate : Check an atom's aromaticity
 *     (aromatic_atom_predicate(), not_aromatic_atom_predicate()). 
 *  
 * Atom predicates can be logically combined using the AtomAndPredicates and
 * AtomOrPredicates classes. These also have functions to make construction
 * easier (atom_and_predicates() and atom_or_predicates()).
 *
 * The case for bonds is similar and bond predicates can also be combined using
 * BondAndPredicates (bond_and_predicates()) and BondOrPredicates
 * (bond_or_predicates()).
 *
 * @li OrderPredicate: Compare a bond's order (order_eq_predicate(),
 *     order_lt_predicate(), order_gt_predicate(), order_leq_predicate(),
 *     order_geq_predicate()).
 * @li AromaticBondPredicate : Check a bond's aromaticity
 *     (aromatic_bond_predicate(), not_aromatic_bond_predicate()). 
 *
 * 
 *
 * @subsection mol_concept_substruct Substructures
 *
 * The Substructure class can be used to create a view on a molecule's substructure.
 * An instance of the Substructure class can be used anywhere a molecule argument is
 * expected.
 *
 * @section depthfirst Depth-First Search
 *
 * Helium provides the depth_first_search() and exhaustive_depth_first_search()
 * functions for performing depth-first searches (dfs). Both functions make use
 * of visitor functors to perform various actions at each of the steps involved.
 * There are a number of visitors (e.g. DFSAtomOrderVisitor, DFSBondOrderVisitor,
 * DFSClosureRecorderVisitor, DFSDebugVisitor, ...) available that can be used
 * directly. Additional visitors can be added by subclassing the DFSVisitor class
 * and reimplementing the required functions.
 *
 * The depth_first_search() is the simplest function and considers all connected
 * components in the molecule. In other words, a single spanning tree is created
 * for each connected component. The exhaustive_depth_first_search() considers
 * only one connected component (i.e. the one which contains the specified atom)
 * but when there are multiple neighbors to visit next, each permutation of these
 * neighbors is tried in succession. As a result, all possible spanning trees
 * having the specified atom as root are found. Both functions are defined in the
 * algorithms/dfs.h header.
 *
 * The example below reads a SMILES string and invokes the DFSDebugVisitor to
 * print the actions of the dfs algorithm to standard output.
 *
 * @code
 * // examples/dfs1.cpp
 * #include <Helium/algorithms/dfs.h>
 * #include <Helium/hemol.h>
 *
 * using namespace Helium;
 *
 * int main()
 * {
 *   HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");
 *
 *   DFSDebugVisitor<HeMol> visitor;
 *   depth_first_search(mol, visitor);
 * }
 * @endcode
 *
 * The second example below is similar but uses the exhaustive_depth_first_search().
 *
 * @code
 * // examples/dfs2.cpp
 * #include <Helium/algorithms/dfs.h>
 * #include <Helium/hemol.h>
 *
 * using namespace Helium;
 *
 * int main()
 * {
 *   HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");
 *
 *   DFSDebugVisitor<HeMol> visitor;
 *   exhaustive_depth_first_search(mol, get_atom(mol, 0), visitor);
 * }
 * @endcode
 *
 * @section breadthfirst Breadth-First Search
 *
 * Similar to the depth-first search, Helium also provides a function for
 * performing breadth-first searches. The breadth_first_search() function
 * is very similar in usage to depth_first_search() and defined in the
 * algorithms/bfs.h header.
 *
 * Example using breadth-first search:
 *
 * @code
 * // examples/bfs1.cpp
 * #include <Helium/algorithms/bfs.h>
 * #include <Helium/hemol.h>
 *
 * using namespace Helium;
 *
 * int main()
 * {
 *   HeMol mol = hemol_from_smiles("C1CCC(CC)CC1");
 *
 *   BFSDebugVisitor<HeMol> visitor;
 *   breadth_first_search(mol, visitor);
 * }
 * @endcode
 *
 * @section components Connected Components
 *
 * The connected components functions are defined in the components.h header file.
 *
 * @li num_connected_components()
 * @li connected_bond_components()
 * @li connected_atom_components()
 *
 * The example below illustrates the use of these functions:
 *
 * @code
 * // examples/components.cpp
 * #include <Helium/algorithms/components.h>
 * #include <Helium/hemol.h>
 *
 * using namespace Helium;
 *
 * int main()
 * {
 *   HeMol mol = hemol_from_smiles("CCC.CC.C");
 *
 *   Size numComponents = num_connected_components(mol);
 *
 *   std::cout << "# components: " << numComponents << std::endl;
 *
 *   std::vector<unsigned int> atomComponentIndices = connected_atom_components(mol);
 *   std::vector<unsigned int> bondComponentIndices = connected_bond_components(mol);
 *
 *   for (Size c = 0; c < numComponents; ++c) {
 *     std::cout << "component " << c << ":" << std::endl;
 *
 *     // print component atoms
 *     std::cout << "  atom indices: ";
 *     for (Index i = 0; i < num_atoms(mol); ++i)
 *       if (atomComponentIndices[i] == c)
 *         std::cout << i << " ";
 *     std::cout << std::endl;
 *
 *     // print component bonds
 *     std::cout << "  bond indices: ";
 *     for (Index i = 0; i < num_bonds(mol); ++i)
 *       if (bondComponentIndices[i] == c)
 *         std::cout << i << " ";
 *     std::cout << std::endl;
 *   }
 * }
 * @endcode
 *
 * The output of the program above:
 *
 * @verbatim
   # components: 3
   component 0:
     atom indices: 0 1 2
     bond indices: 0 1
   component 1:
     atom indices: 3 4
     bond indices: 2
   component 2:
     atom indices: 5
     bond indices:
   @endverbatim
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
 *
 * @section cycles Cycles
 *
 * There are several algorithms available related to cycles:
 *
 * @li cyclomatic_number()
 * @li cycle_membership()
 * @li relevant_cycles()
 */

}
