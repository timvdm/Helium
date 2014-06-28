/**
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
#include <Helium/molecule.h>
#include <Helium/element.h>
#include <Helium/fileio/file.h>
#include <Helium/json/json.h>

#include "args.h"
#include "progress.h"

#include <openbabel/obconversion.h>
#include <openbabel/mol.h>


using namespace Helium;

std::string normalize_smiles(const std::string &smiles)
{
  OpenBabel::OBMol obmol;
  OpenBabel::OBConversion conv;
  conv.SetInAndOutFormats("smi", "smi");
  conv.ReadString(&obmol, smiles);
  return conv.WriteString(&obmol, true);
}

void write_molecule(std::ostream &os, OpenBabel::OBMol *mol)
{
  // hydrogen atoms are not written to the file
  // create a map to map atom indices & count the number of heavy atoms
  std::vector<unsigned short> indices(mol->NumAtoms());
  unsigned short numAtoms = 0;
  FOR_ATOMS_OF_MOL (a, mol) {
    indices[a->GetIndex()] = numAtoms;
    if (!a->IsHydrogen())
      ++numAtoms;
  }
  // count the number of bonds between heavy atoms
  unsigned short numBonds = 0;
  FOR_BONDS_OF_MOL (b, mol)
    if (!b->GetBeginAtom()->IsHydrogen() && !b->GetEndAtom()->IsHydrogen())
      ++numBonds;

  // write the number of atoms & bonds
  write16<unsigned short>(os, numAtoms);
  write16<unsigned short>(os, numBonds);

  // write atoms (6 byte / atom)
  FOR_ATOMS_OF_MOL (atom, mol) {
    if (atom->IsHydrogen())
      continue;

    // write the element
    write8<unsigned char>(os, atom->GetAtomicNum());
    // write cyclic property
    write8<unsigned char>(os, atom->IsInRing());
    // write aromatic property
    write8<unsigned char>(os, atom->IsAromatic());
    // write isotope
    if (atom->GetIsotope())
      write8<unsigned char>(os, atom->GetIsotope());
    else
      write8<unsigned char>(os, Element::averageMass(atom->GetAtomicNum()));
    // write hydrogen count
    write8<unsigned char>(os, atom->ExplicitHydrogenCount() + atom->ImplicitHydrogenCount());
    // write formal charge
    write8<signed char>(os, atom->GetFormalCharge());
  }

  // write bonds (5 byte / bond
  FOR_BONDS_OF_MOL (bond, mol) {
    if (bond->GetBeginAtom()->IsHydrogen() || bond->GetEndAtom()->IsHydrogen())
      continue;

    // write source & target indices
    write16<unsigned short>(os, indices[bond->GetBeginAtom()->GetIndex()]);
    write16<unsigned short>(os, indices[bond->GetEndAtom()->GetIndex()]);
    // write bond order + aromatic & cyclic properties
    unsigned char props = bond->GetBondOrder();
    if (bond->IsAromatic())
      props |= 128;
    if (bond->IsInRing())
      props |= 64;
    write8<unsigned char>(os, props);
  }
}


int main(int argc, char**argv)
{
  if (argc < 2) {
    std::cout << "Usage: " << argv[0] << " <in_file> <out_file>" << std::endl;
    return 0;
  }

  ParseArgs args(argc, argv, ParseArgs::Args(), ParseArgs::Args("in_file", "out_file"));
  std::string inFile = args.GetArgString("in_file");
  std::string outFile = args.GetArgString("out_file");


  // open the input file
  std::ifstream ifs(inFile.c_str());

  // open the output file
  BinaryOutputFile outputFile;
  if (!outputFile.open(outFile)) {
    std::cerr << outputFile.error().what() << std::endl;
    return -1;
  }

  unsigned int numMolecules = 0;

  // open the input file using OpenBabel
  OpenBabel::OBMol mol;
  OpenBabel::OBConversion conv(&ifs);
  conv.SetInFormat(conv.FormatFromExt(inFile));

  // start converting the molecules
  std::vector<std::size_t> positions;
  while (conv.Read(&mol)) {
    ++numMolecules;
    unknown_progress("Converting molecules", numMolecules);
    positions.push_back(outputFile.stream().tellp());
    write_molecule(outputFile.stream(), &mol);
  }
  std::cout << std::endl;

  // save the stream position where the molecule positions are stored
  Json::UInt64 positionsPos = outputFile.stream().tellp();

  // write the molecule positions to the file
  outputFile.write(&positions[0], positions.size() * sizeof(std::size_t));

  // create JSON header
  Json::Value data;
  data["filetype"] = "molecules";
  data["num_molecules"] = numMolecules;
  data["molecule_indexes"] = positionsPos;

  // write JSON header
  Json::StyledWriter writer;
  outputFile.writeHeader(writer.write(data));
}
