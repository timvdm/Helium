#include "../src/fileio/molecules.h"
#include "../src/fileio/file.h"

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

int main()
{
  test_binary_file();
}
