#include <Helium/fileio/moleculefile.h>
#include <Helium/fileio/file.h>
#include <Helium/hemol.h>

#include "test.h"

using namespace Helium;

void test_binary_file()
{
  // write a file
  BinaryOutputFile out("tmp.hel");

  ASSERT(out);

  const char data[] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
  out.write(data, 10 * sizeof(char));

  std::string header = "{ foo: 42 }";
  ASSERT(out.writeHeader(header));

  out.close();

  // read the file
  BinaryInputFile in("tmp.hel");

  ASSERT(in);
  COMPARE(header, in.header());

  char value;
  for (int i = 0; i < 10; ++i) {
    ASSERT(in.read(&value, sizeof(char)));
    COMPARE(i + 1, value);
  }

}

void test_molecule_file()
{
  MoleculeFile molFile1(datadir() + "10K.hel");
  MoleculeFile molFile2(datadir() + "10K.hel");

  HeMol mol1, mol2;

  for (unsigned int i = 0; i < molFile1.numMolecules(); ++i) {
    COMPARE(molFile1.stream().tellg(), molFile2.stream().tellg());

    molFile1.readMolecule(mol1);
    molFile2.readMolecule(i, mol2);
  }

}

int main()
{
  test_binary_file();
  test_molecule_file();
}
