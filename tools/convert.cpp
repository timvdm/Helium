#include "../src/molecule.h"
#include "../src/fileio.h"

#include "args.h"

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>

using namespace Helium;

/*
void write_molecule(std::ostream &os, OpenBabel::OBMol *mol)
{
  std::vector<unsigned short> indices(mol->NumAtoms());
  unsigned short numAtoms = 0;
  FOR_ATOMS_OF_MOL (a, mol) {
    indices[a->GetIndex()] = numAtoms;
    if (!a->IsHydrogen())
      ++numAtoms;
  }
  unsigned short numBonds = 0;
  FOR_BONDS_OF_MOL (b, mol)
    if (!b->GetBeginAtom()->IsHydrogen() && !b->GetEndAtom()->IsHydrogen())
      ++numBonds;

  write16<unsigned short>(os, numAtoms);
  write16<unsigned short>(os, numBonds);

  FOR_ATOMS_OF_MOL (atom, mol) {
    if (atom->IsHydrogen())
      continue;

    // always write these properies
    write8<unsigned char>(os, atom->GetAtomicNum());

    unsigned char aromaticCyclic = 0;
    if (atom->IsAromatic())
      aromaticCyclic |= 1;
    if (atom->IsInRing())
      aromaticCyclic |= 2;

    unsigned char flags = 0;
    if (aromaticCyclic)
      flags |= AromaticCyclic;
    if (atom->GetIsotope())
      flags |= Mass;
    if (atom->ExplicitHydrogenCount() + atom->ImplicitHydrogenCount())
      flags |= Hydrogens;
    if (atom->GetFormalCharge())
      flags |= Charge;

    // write flags
    write8<unsigned char>(os, flags);

    if (flags & AromaticCyclic)
      write8<unsigned char>(os, aromaticCyclic);
    if (flags & Mass)
      write8<unsigned char>(os, atom->GetIsotope());
    if (flags & Hydrogens)
      write8<unsigned char>(os, atom->ExplicitHydrogenCount() + atom->ImplicitHydrogenCount());
    if (flags & Charge)
      write8<signed char>(os, atom->GetFormalCharge());
  }

  FOR_BONDS_OF_MOL (bond, mol) {
    if (bond->GetBeginAtom()->IsHydrogen() || bond->GetEndAtom()->IsHydrogen())
      continue;

    write16<unsigned short>(os, indices[bond->GetBeginAtom()->GetIndex()]);
    write16<unsigned short>(os, indices[bond->GetEndAtom()->GetIndex()]);
    unsigned char props = bond->GetBondOrder();
    if (bond->IsAromatic())
      props |= 128;
    if (bond->IsInRing())
      props |= 64;
    write8<unsigned char>(os, props);
  }
}
*/


int main(int argc, char**argv)
{
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <in_file> <out_file>" << std::endl;
    return 0;
  }

  ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("in_file", "out_file"));
  std::string inFile = args.GetArgString("in_file");
  std::string outFile = args.GetArgString("out_file");


  std::ifstream ifs(inFile.c_str());
  std::ofstream ofs(outFile.c_str(), std::ios_base::out | std::ios_base::binary);

  unsigned int numMolecules = 0;
  write32(ofs, numMolecules);

  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv(&ifs);
  conv.SetInFormat(conv.FormatFromExt(inFile));

  while (conv.Read(&mol)) {
    ++numMolecules;
    if ((numMolecules % 1000) == 0)
      std::cout << numMolecules << std::endl;
    write_molecule(ofs, &mol);
  }

  ofs.seekp(0);
  write32(ofs, numMolecules);

}
