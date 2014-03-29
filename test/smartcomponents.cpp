#include <Helium/algorithms/smartcomponents.h>
#include <Helium/smartmol.h>
#include <Helium/smiles.h>
#include <Helium/fileio/moleculefile.h>

#include "test.h"

using namespace Helium;

void test_lazy_bond_components(const std::string &smiles)
{
  SmartMol mol;
  std::vector<unsigned int> components;

  LazyBondComponents *lazyBondComponents = new LazyBondComponents(mol);
  COMPARE("Helium::bond_components", lazyBondComponents->name());
  mol.addAttribute(lazyBondComponents);

  try {
    parse_smiles(smiles, mol);
    components = connected_bond_components(mol);
    COMPARE(components, lazyBondComponents->components());
    COMPARE(unique_elements(components), lazyBondComponents->numComponents());
  } catch (Smiley::Exception &e) {
    std::cerr << e.what();
  }
}

void test_lazy_atom_components(const std::string &smiles)
{
  SmartMol mol;
  std::vector<unsigned int> components;

  LazyAtomComponents *lazyAtomComponents = new LazyAtomComponents(mol);
  COMPARE("Helium::atom_components", lazyAtomComponents->name());
  mol.addAttribute(lazyAtomComponents);

  try {
    parse_smiles(smiles, mol);
    components = connected_atom_components(mol);
    COMPARE(components, lazyAtomComponents->components());
    COMPARE(unique_elements(components), lazyAtomComponents->numComponents());
  } catch (Smiley::Exception &e) {
    std::cerr << e.what();
  }
}

bool compare_components(const std::vector<unsigned int> &componentsRef, const std::vector<unsigned int> &components)
{
  //std::cout << "components ref: " << componentsRef << std::endl;
  //std::cout << "components:     " << components << std::endl;

  Size numComponentsRef = unique_elements(componentsRef);
  Size numComponents = unique_elements(components);
  COMPARE(numComponentsRef, numComponents);
  if (numComponentsRef != numComponents)
    return false;

  std::vector<std::vector<bool> > componentsBitvecRefs(numComponents, std::vector<bool>(componentsRef.size()));
  for (std::size_t i = 0; i < componentsRef.size(); ++i)
    componentsBitvecRefs[componentsRef[i]][i] = true;

  std::vector<std::vector<bool> > componentsBitvecs(numComponents, std::vector<bool>(componentsRef.size()));
  for (std::size_t i = 0; i < components.size(); ++i)
    componentsBitvecs[components[i]][i] = true;

  for (std::size_t i = 0; i < numComponents; ++i) {
    ASSERT(std::find(componentsBitvecs.begin(), componentsBitvecs.end(), componentsBitvecRefs[i]) != componentsBitvecs.end());
    if (std::find(componentsBitvecs.begin(), componentsBitvecs.end(), componentsBitvecRefs[i]) == componentsBitvecs.end()) {
      return false;
    }

  }

  return true;
}

void test_dynamic_atom_components(SmartMol &mol, DynamicAtomComponents *dynamicAtomComponents)
{
  std::cout << write_smiles(mol) << std::endl;
  std::vector<unsigned int> componentsRef;

  componentsRef = connected_atom_components(mol);
  ASSERT(compare_components(componentsRef, dynamicAtomComponents->components()));
  if (!compare_components(componentsRef, dynamicAtomComponents->components())) {
    return;
  }


  while (num_atoms(mol)) {
    remove_atom(mol, get_atom(mol, 0));

    componentsRef = connected_atom_components(mol);
    ASSERT(compare_components(componentsRef, dynamicAtomComponents->components()));
    if (!compare_components(componentsRef, dynamicAtomComponents->components())) {
      exit(0);
      return;
    }
  }
}

void test_dynamic_atom_components(const std::string &smiles)
{
  std::cout << "test_dynamic_atom_components(" << smiles << ")" << std::endl;
  SmartMol mol;

  DynamicAtomComponents *dynamicAtomComponents = new DynamicAtomComponents(mol);
  COMPARE("Helium::atom_components", dynamicAtomComponents->name());
  mol.addAttribute(dynamicAtomComponents);

  try {
    parse_smiles(smiles, mol);
  } catch (Smiley::Exception &e) {
    std::cerr << e.what();
    ASSERT(0);
    return;
  }

  test_dynamic_atom_components(mol, dynamicAtomComponents);
}

void test_dynamic_atom_components_file(const std::string &filename)
{
  std::cout << "test_dynamic_atom_components_file(" << filename << ")" << std::endl;
  MoleculeFile file(filename);

  for (unsigned int i = 0; i < file.numMolecules(); ++i) {
    SmartMol mol;
    DynamicAtomComponents *dynamicAtomComponents = new DynamicAtomComponents(mol);
    mol.addAttribute(dynamicAtomComponents);
    file.readMolecule(mol);
    test_dynamic_atom_components(mol, dynamicAtomComponents);
  }
}

int main()
{
  test_lazy_bond_components("CC(C)C");
  test_lazy_bond_components("CC.CC");
  test_lazy_bond_components("C1CC1.CC");

  test_lazy_atom_components("CC(C)C");
  test_lazy_atom_components("CC.CC");
  test_lazy_atom_components("C1CC1.CC");

  test_dynamic_atom_components("CCC");
  test_dynamic_atom_components("CC(C)C");
  test_dynamic_atom_components("CC.CC");
  test_dynamic_atom_components("C1CC1.CC");

  test_dynamic_atom_components_file(datadir() + "100K.hel");
}

