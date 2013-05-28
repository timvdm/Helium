#include <Helium/bitvec.h>

#include "test.h"

using namespace Helium;


int main()
{
  Word bitvec;

  hex_to_bitvec("0000000000000000", &bitvec, 1);
  COMPARE(bitvec, 0);

  hex_to_bitvec("0123456789abcdef", &bitvec, 1);
  COMPARE("0123456789abcdef", bitvec_to_hex(&bitvec, 1));

  hex_to_bitvec("01ab02cd03ef", &bitvec, 1);
  COMPARE("01ab02cd03ef0000", bitvec_to_hex(&bitvec, 1));
}
