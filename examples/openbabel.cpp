// examples/openbabel.cpp
#include <Helium/toolkits/openbabel.h>
#include <openbabel/obconversion.h>
#include <Helium/smirks.h>
#include <Helium/smiles.h>

using namespace Helium;

Smiles SMILES;

template<typename EditableMoleculeType>
void amide_formation(EditableMoleculeType &mol)
{
  // initiaize SMIRKS
  Smirks smirks;
  if (!smirks.init("[Cl:3][C:1]=[O:2].[C:4][N;H2:5]>>[O:2]=[C:1][N:5][C:4].[Cl:3]")) {
    std::cout << "Error: " << smirks.error().what() << std::endl;
    return;
  }

  // apply the SMIRKS transformation
  if (!smirks.apply(mol, RingSet<EditableMoleculeType>(mol))) {
    std::cout << "Could not apply SMIKRS, reactant SMARTS not found" << std::endl;
    return;
  }
}

void amide_formation(const std::string &amine, const std::string &acylchloride)
{
  // read SMILES using OpenBabel
  OpenBabel::OBConversion conv;
  conv.SetInFormat("smi");
  OpenBabel::OBMol mol;
  conv.ReadString(&mol, amine + "." + acylchloride);

  // apply SMIRKS
  amide_formation(mol);

  std::string product = SMILES.write(mol);
  std::cout << amine << " + " << acylchloride + " -> " + product << std::endl;
}

int main()
{
  amide_formation("CCN", "CCC(=O)Cl");
  amide_formation("CCN", "c1ccccc1CC(=O)Cl");
}
