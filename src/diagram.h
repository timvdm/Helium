/* This file is based on the public domain MCDL program:
 *
 * Copyright (C) 2007,2008 by Sergei V. Trepalin sergey_trepalin@chemical-block.com
 * Copyright (C) 2007,2008 by Andrei Gakh andrei.gakh@nnsa.doe.gov
 *
 * http://sourceforge.net/projects/mcdl/
 */

/*
  Diagram is generated using templates, which are stored in SD file templates.sdf
  The SD file is usual SD file, which contain chemical structures and might contain data.
  Only chemical structures are used. Subgraph isomorphisme search is executed and coordinates
  of atoms are determined from templates. See Molecules, 11, 129-141 (2006) for algorithm decription.
  Structures in SD file are converted in next manner:
  1. All atoms, except explicit hydrogens, are replaced with generic ANY_ATOM (matched with any atom in subgraph isomorphisme search)
  2. All bonds are replaces with generic ANY_BOND, which can be matched with any bond in molecule
  3. All hydrogen are removed, but they are used for search-query and structure atom matching is believed fo be
     sucessfukk if chemical structure contains more or equal number of hydrogens, than query. Using explicitly-defined hydrogens
	 on query enables ones to remove substitutors attachment for atom, which are sterically hidden on templates
  if the file will not be found, predefined templates will be used
*/

#ifndef HELIUM_DIAGRAM_H
#define HELIUM_DIAGRAM_H

#include <Helium/molecule.h>
#include <Helium/algorithms/gtd.h>

namespace Helium {
  
  namespace impl {
    namespace diagram {
      std::vector<std::pair<double, double> > generate_diagram(const std::vector<int> &elements,
          const std::vector<int> &charges, const std::vector<unsigned int> &gtd,
          const std::vector<std::pair<int, int> > &bonds, const std::vector<int> &orders);
    }
  }

  template<typename MoleculeType>
  std::vector<std::pair<double, double> > generate_diagram(const MoleculeType &mol)
  {
    typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
    typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
    
    // Stereo perception should not be triggered if this function is called
    // from MDLFormat::ReadMolecule->Alias::Expand->Alias::FromNameLookup->(MCDL)groupRedraw
    // as this triggers a failure in the test suite for InChI conversion
    // (specifically the ferrocene).
    /* FIXME: stereo
    bool perceive_stereo = true;
    if (pmol->GetMod() == 1)
      perceive_stereo = false;
    OBStereoFacade facade(pmol);
    */

    // When laying out, we may need to iterate over the most central bonds
    // first - that's why we remember the GTD: the GTD is the distance from
    // atom i to every other atom j. Atoms on the "inside" of the molecule
    // will have a lower GTD value than atoms on the "outside".
    std::vector<unsigned int> gtd = graph_theoretical_distance(mol);

    std::vector<int> elements;
    std::vector<int> charges;
    std::vector<std::pair<int, int> > bonds;
    std::vector<int> orders;

    for (Size i = 0; i < num_atoms(mol); i++) {
      atom_type atom = get_atom(mol, i);
      elements.push_back(get_element(mol, atom));
      charges.push_back(get_charge(mol, atom));
    }

    for (Size i = 0; i < num_bonds(mol); i++) {
      bond_type bond = get_bond(mol, i);
      bonds.push_back(std::make_pair(get_index(mol, get_source(mol, bond)),
                                     get_index(mol, get_target(mol, bond))));
      orders.push_back(get_order(mol, bond));
    }

    return impl::diagram::generate_diagram(elements, charges, gtd, bonds, orders);
  }

}

#endif
