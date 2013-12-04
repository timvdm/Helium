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

/**
 * @brief Helium namespace.
 */
namespace Helium {

/**
 * @mainpage
 *
 * @section Helium
 *
 * @li @ref design
 * @li @ref mol_concept
 * @li @ref depthfirst
 * @li @ref breadthfirst
 * @li @ref components
 * @li @ref cycles
 * @li @ref smarts
 * @li @ref smirks
 * @li @ref fingerprints
 * @li @ref morganec
 * @li @ref canonical
 * @li @ref enumerate
 *
 * @section design Design Philosophy
 *
 * @note This section can be skipped when in a hurry. It gives more details
 * on my personal experiences with working on cheminformatics projects.
 *
 * For several years, I have been working on OpenBabel [1] and I still find it
 * a very useful cheminformatics toolkit. Especially the capability to convert
 * between a large number of file formats is impressive. However, while working
 * on OpenBabel several limitations arose that I try to address with the Helium
 * project.
 *
 * OpenBabel was developed by multiple people over many years. This
 * often resulted in suboptimal design choices for which API/ABI compatibility
 * has to be maintained. For some of these I am myself responsible. For example,
 * the forcefield code's design isn't great since it was my first large piece
 * of code I wrote. The same is true for the OBMol class, it contains far too
 * many functions that should not be there. Both the OBForceField and OBMol
 * class do their work but they are far from elegant. The remaining sections
 * will address more issues in detail.
 *
 * @subsection design_modular Modularity
 *
 * This probably doesn't require much explanation since most developers will
 * know that this is a good practise.
 * Classes in Helium are as small as possible. For example, the molecule concept
 * explained in detail below contains as few functions as possible. This makes
 * the code more modular and also results in easier to maintain code.
 *
 * @subsection design_just_ask Data and Behaviour: Just Ask...
 *
 * OpenBabel heavily relies on lazy-perception for various properties (e.g.
 * ring perception, aromaticity, ...) which is great if you just want to get
 * work done quickly (i.e. you don't have to worry about all of this). On the
 * other hand, adding a new ring perception algorithm to OpenBabel and
 * benchmarking it's performance is far from straightforward. For users this
 * might not be so important but better performance or quality should also
 * benefit users who only use the programs based on a cheminformatics toolkit.
 *
 * The main idea used in Helium was already mentioned in this section's title:
 * "Just Ask...". To clarify, when a function or a class needs some @b data, it
 * should simply accept this as a parameter. Functions can be overloaded to
 * provide some sensible defaults to make their use easier but hard coded
 * defaults should be avoided as much as possible. For example the functions
 * to perform a SMARTS search on a molecule are given below:
 *
 @code
 // provide previously computed ring set (needed for queries such as [R], [r6], ...)
 template<typename MoleculeType>
 bool search(const MoleculeType &mol, const RingSet<MoleculeType> &rings)

 // compute ring set in function, slower but faster for prototyping
 template<typename MoleculeType>
 bool search(const MoleculeType &mol)
 @endcode
 *
 * This idea is used everywhere in Helium to allow maximum flexibility. For
 * example, SMARTS originally used the SSSR ring set which doesn't return
 * unique results (i.e. these depend on the SSSR, a minimum cycle basis,
 * chosen). Helium mainly uses the relevant cycles but users might prefer
 * to still perform the search using a SSSR ring set.
 *
 * The "Just Ask..." idea can be further extended to specify @b behaviour. This
 * technique is used very successfully in the STL. Many of the STL algorithms
 * accept functors to change the behaviour of the underlying algorithm (e.g.
 * Compare, ...). An example can be found in the code for doing isomorphism
 * searches on molecules (i.e. substructure search). The code for matching
 * query atoms against atoms in the queried target is a functor that is
 * provided by the user (default are provided). This separates the logic for
 * performing the isomorphism search from the query language (e.g. SMARTS)
 * used for matching atoms and bonds. This is illustrated below where the
 * AtomMatcher and BondMatcher are used for matching atoms and bond.
 *
 @code
 template<typename MoleculeType, typename AtomType, typename QueryType, typename MappingType, typename AtomMatcher, typename BondMatcher>
 bool isomorphism_search(const MoleculeType &mol, AtomType atom, QueryType &query, MappingType &mapping,
     const AtomMatcher &atomMatcher, const BondMatcher &bondMatcher) // <------
 @endcode
 *
 * @subsection design_toolkits One Toolkit to Rule Them All?
 *
 * Helium does not try (yet :) ) to be the only toolkit. Instead, it provides
 * building blocks that other toolkits can use. This is achieved without a
 * significant loss of performance (i.e. making copies of molecules each time)
 * by the use of @b generic @b programming. Although @b concepts are not part
 * of the C++11 standard, the idea can still be used albeit many beginning C++
 * programmers might suffer from panic attacks while debugging their code.
 * The isomorphism_search function with all the template parameters already
 * showed that @b templates are used extensively to program generically in C++.
 *
 * To motivate new users, the example/openbabel.cpp (@ref mol_toolkits) file
 * illustrates how Helium's SMIRKS implementation (a feature which isn't
 * really working in OpenBabel) to be used with OBMol objects that are read
 * using OpenBabel.
 *
 * @subsection design_license License
 *
 * The 3-clause BSD license (http://en.wikipedia.org/wiki/BSD_licenses) is used
 * since it allows the code to be used, modified and distributed without to much
 * restrictions (i.e. also for commercial use). Due to the origins of OpenBabel
 * it is currently stuck with the GPLv2 license. For most users (e.g. academia)
 * this isn't a big issue but for others this might hold them back from spending
 * their precious time on a project which they might not be able to use eventually
 * due to licensing issues.
 *
 @verbatim
 [1] N. M. O'Boyle, M. Banck, C. A. James, C. Morley, T. Vandermeersch, G. Hutchison,
     OpenBabel: An open chemical toolbox, Journal of Cheminformatics, 2001, 3:33
     http://www.jcheminf.com/content/3/1/33
 @endverbatim
 *
 * @section mol_concept The Molecule Concept
 *
 * Helium is C++ library and uses generic programming to achieve
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
 @code
 template<typename MoleculeType>
 void print_num_atoms(const MoleculeType &mol)
 {
   std::cout << num_atoms(mol) << std::endl;
 }
 @endcode
 *
 * When additional types are required inside the function, the molecule_traits
 * struct can be used to add typedefs:
 *
 @code
 template<typename MoleculeType>
 void print_degree_of_first_atom(const MoleculeType &mol)
 {
   // the "typename" is needed since molecule_traits is dependent on the
   // MoleculeType template parameter...
   typedef typename molecule_traits<MoleculeType>::atom_type atom_type;

   atom_type atom = get_atom(mol, 0);

   std::cout << get_degree(mol, atom) << std::endl;
 }
 @endcode
 *
 * If a bond or atom needs to be passed as argument, there are two ways this can
 * be achieved. The first requires more typing so the second way is recommended:
 *
 @code
 // the first and long way
 template<typename MoleculeType>
 void print_atom_index(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
 {
   std::cout << get_index(mol, atom) << std::endl;
 }

 // the second and easier way
 template<typename MoleculeType, typename AtomType>
 void print_atom_index(const MoleculeType &mol, AtomType atom)
 {
   std::cout << get_index(mol, atom) << std::endl;
 }
 @endcode
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
 *
 * @li @b num_atoms(mol): Get the number of atoms.
 * @li @b num_bonds(mol): Get the number of bonds.
 * @li @b get_atom(mol, index): Get the atom with the specified index.
 * @li @b get_bond(mol, index): Get the bond with the specified index.
 * @li @b get_atoms(mol): Get an iterator pair over all atoms in the molecule.
 * @li @b get_bonds(mol): Get an iterator pair over all atoms in the molecule.
 * @li @b get_bond(mol, source, target): Get the bond between source and target.
 *
 * The last two functions are used for iterating over the atoms/bonds. Although
 * these can be used directly, Helium also provides macros (i.e. FOREACH_ATOM()
 * and FOREACH_BOND() to make iteration easier. The usage of all the functions
 * above is illustrated in the simple example below:
 *
 * @include molecule1.cpp
 *
 * Once an @b atom has been obtained (e.g. using get_atom(), FOREACH_ATOM(),
 * FOREACH_NBR(), ...), the following set of functions can be used to retrieve
 * information about the atom:
 *
 * @li @b get_index(mol, atom): Get the atom's index.
 * @li @b get_element(mol, atom): Get the atom's element.
 * @li @b get_mass(mol, atom): Get the atom's mass.
 * @li @b get_degree(mol, atom): Get the atom's degree.
 * @li @b get_charge(mol, atom): Get the atom's charge.
 * @li @b num_hydrogens(mol, atom): Get the number of implicit hydrogens attached to this atom.
 * @li @b get_bonds(mol, atom): Get an iterator pair over the atom's incident bonds.
 * @li @b get_nbrs(mol, atom): Get an iterator pair over the atom's neighbor atoms.
 * @li @b get_nbrs(mol, atom): Get an iterator pair over the atom's neighbor atoms.
 * @li @b is_aromatic(mol, atom): Get the atom's aromaticity.
 *
 * Again, the last two functions return iterator pairs. For these the FOREACH_NBR()
 * and FOREACH_INCIDENT() marcos provided.
 *
 * @include molecule2.cpp
 *
 * For @b bonds the following functions are provided.
 *
 * @li @b get_index(mol, bond): Get the bond's index.
 * @li @b get_source(mol, bond): Get the bond's source atom.
 * @li @b get_target(mol, bond): Get the bond's target atom.
 * @li @b get_order(mol, bond): Get the bond's order.
 * @li @b get_other(mol, bond, atom): Get the other bond atom.
 * @li @b is_aromatic(mol, bond): Get the bond's aromaticity.
 *
 * Simple example showing usage of the various bond functions:
 *
 * @include molecule3.cpp
 *
 * All of the functions above are primitive functions in the sense that each model
 * of the molecule concept has to implement them. Using these primitive functions,
 * a number of useful functions are defined.
 *
 * @li is_hydrogen(): Check if an atom is a hydrogen atom.
 * @li is_carbon(): Check if an atom is a carbon atom.
 * @li is_nitrogen(): Check if an atom is a nitrogen atom.
 * @li is_oxygen(): Check if an atom is an oxygen atom.
 * @li is_phosphorus(): Check if an atom is an phosphorus atom.
 * @li is_sulfur(): Check if an atom is an sulfur atom.
 * @li get_heavy_degree(): Get the number of heavy atoms (i.e. not hydrogen) attached to an atom.
 * @li get_valence(): Get the valence (i.e. the bond order sum) of an atom.
 * @li get_connectivity(): Get the connectivity (i.e. the degree + the number of implicit hydrogens) of an atom.
 *
 * @subsection mol_concept_predicates Atom & Bond Predicates
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
 * In the example below, the SMILES for benzoic acid is read into a molecule and
 * a SMARTS search is done for a benzene ring. Next, the found mapping is used to
 * create a substructure from the molecule resuling in a benzene ring view of the
 * molecule.
 *
 * @include substructure.cpp
 *
 * @subsection mol_concept_editable Editable Molecules
 *
 * The Molecule concept is further extended to the EditableMolecule concept to
 * allow molecules to be edited. The following functions operate on the molecule
 * to add/remove atoms and bonds or clear the molecule:
 *
 * @li @b add_atom(mol): Add an atom to the molecule.
 * @li @b add_bond(mol, source, target): Add a bond to the molecule.
 * @li @b clear_molecule(mol): Clear the molecule (i.e. remove all atoms/bonds).
 *
 * Once an atom is created (or obtained another way) several set functions can
 * be used to change the atom properties:
 *
 * @li @b set_element(mol, atom, value): Set the atom's element.
 * @li @b set_mass(mol, atom, value): Set the atom's mass.
 * @li @b set_charge(mol, atom, value): Set the atom's charge.
 * @li @b set_hydrogens(mol, atom, value): Set the number of implicit hydrogens attached to this atom.
 * @li @b set_aromatic(mol, bond, value): Set the atom's aromaticity.
 *
 * Similar functions are provided for bonds:
 *
 * @li @b set_order(mol, bond, value): Set the bond's order.
 * @li @b set_aromatic(mol, bond, value): Set the bond's aromaticity.
 *
 * The Element class provides several functions to obtain default properties for
 * atoms that can be used to set the various atom properties. The example below
 * illustrates how acetone can be created:
 *
 * @include molecule4.cpp
 *
 * @subsection mol_toolkits Other Toolkits
 *
 * The toolkits/openbabel.h header can be included to make an OBMol a model
 * of the Molecule concept. This is illustrated in the example below which
 * is identical to the SMIRKS example apart from the molecule class used and
 * the way this molecule is first constructed.
 *
 * @include openbabel.cpp
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
 * @include dfs1.cpp
 *
 * The second example below is similar but uses the exhaustive_depth_first_search().
 *
 * @include dfs2.cpp
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
 * @include bfs1.cpp
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
 * @include components.cpp
 *
 * The output of the program above:
 *
 @verbatim
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
 * @section cycles Cycles
 *
 * There are several algorithms available related to cycles:
 *
 * @li cyclomatic_number()
 * @li cycle_membership()
 * @li relevant_cycles()
 *
 * Although most cheminformatics provide the Smallest Set of Smallest Rings
 * (SSSR) as the main ring set, Helium uses the set of relevant cycles.
 * Vismara [1] described this set as the union of all the minimum cycles
 * bases of a graph. Definition 1 and lemma 1 from the paper further describe
 * what the relevant cycles are:
 *
 * @li @b Definition @b 1: A cycle is said to be @b relevant if it belongs to at
 *     leat one minimum cycle basis.
 * @li @b Lemma @b 1: @f$C \in C_{all}@f$ is relevant @f$\iff@f$ no elementary cycles
 *     @f$C_1, ..., C_k@f$ exist such that @f$C = C_1 + ... + C_k@f$ and
 *     @f$\forall i \in [1,k], w(C_i) < w(C)@f$. Here, @f$C@f$ is the candidate ring,
 *     @f$C_{all}@f$ is the set of all cycles and @f$w(C)@f$ is the weight (i.e. size)
 *     of the ring @f$C@f$.
 *
 * The reason(s) why the relevant cycles were chosen as ring set is best illustrated
 * in the Berger paper [2]. The relevant cylces have the following properties:
 *
 * @li Works for general graphs (i.e. no restricted to planar graphs).
 * @li The relevant cycles are @b unique.
 * @li The relevant cycles contain the minimum cycle bases (@b MCB).
 *
 * Although the size of the relevant cycles ring set is exponential in theory,
 * this is not a problem since this behavior is not observed for molecular
 * graphs.
 *
 * The example below shows how these functions can be used.
 *
 * @include cycles.cpp
 *
 @verbatim
 [1] P. Vismara, Union of all the minimum cycle bases of a graph, The electronic
     Journal of Combinatorics, 1997, Volume 4, Issue 1.
 [2] F. Berger, C. Flamm, P. M. Gleiss, J. Leydold, P. F. Stadler, Counterexamples
     in Chemical Ring Perception, J. Chem. Inf. Comput. Sci., 2004, Volume 44,
     Issue 2: 323-331.
 @endverbatim
 *
 * @section smarts SMARTS
 *
 * Using the Smarts class SMARTS (http://www.daylight.com/dayhtml/doc/theory/theory.smarts.html) queries can be perofmed on molecules. There
 * are various types of mappings that can be obtained:
 *
 * @li NoMapping: Just check if the SMARTS matches the molecule.
 * @li CountMapping: Count the number of matches, do not retrieve the actual mappings.
 * @li SingleMapping: Get a single mapping.
 * @li MappingList: Get all the (unique) mappings of the SMARTS query to the molecule.
 *
 * A simple example using a single mapping is given below:
 *
 * @include smarts.cpp
 *
 * The example also illustrates the usage of the Smarts::requiresCycles()
 * function. This function returns true if the SMARTS query contains atom or
 * bond primitives tht requires knowledge of cylces (e.g. [R], [r6], *@*, ...).
 * If this is not the case, an empty ring set can be used to obtain better
 * performance.
 *
 * @section smirks SMIRKS
 *
 * Using the Smirks class it is possible to apply SMIRKS (http://www.daylight.com/dayhtml/doc/theory/theory.smirks.html) transformations to
 * editable molecules. In the example below an amine is reacted with an
 * acylchloride to form an amide.
 *
 * @include smirks.cpp
 *
 * The output is:
 *
 @verbatim
 CCN + CCC(=O)Cl -> CCNC(CC)=O.Cl
 CCN + c1ccccc1CC(=O)Cl -> CCNC(Cc1ccccc1)=O.Cl
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
 * @section morganec Morgan's Extended Connectivities
 *
 * The Morgan extended connectivities [1] can be computed using the
 * extended_connectivities() function defined in the
 * extendedconnectivities.h header file.
 *
 @verbatim
 [1] Morgan, H. L. The Generation of a Unique Machine Description for Chemical
     Structures - A Technique Developed at Chemical Abstracts Service. J. Chem.
     Doc. 1965, 5: 107-112.
 @endverbatim
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
 * The example below enumerates all subgraphs up to size 7 and prints out the
 * unique subgraphs with their associated counts:
 *
 * @include enumeratesubgraphs.cpp
 */

 /**
  * @defgroup Production Code in Production Phase
  * @brief Code in production phase.
  */

 /**
  * @defgroup Alpha Code in Alpha Phase
  * @brief Code in aplha phase.
  */

 /**
  * @defgroup Beta Code in Beta Phase
  * @brief Code in beta phase.
  */

}
