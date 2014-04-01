#include <Helium/bitvec.h>

#include "test.h"

using namespace Helium;


int main()
{
  COMPARE(8, sizeof(Word));



  std::pair<Word*, int> bitvec = bitvec_from_hex("0000000000000000");
  COMPARE(*bitvec.first, 0);
  COMPARE(bitvec.second, 1);

  bitvec = bitvec_from_hex("0123456789abcdef");
  COMPARE("0123456789abcdef", bitvec_to_hex(bitvec.first, 1));

  bitvec = bitvec_from_hex("01ab02cd03ef");
  COMPARE("01ab02cd03ef0000", bitvec_to_hex(bitvec.first, 1));

  bitvec = bitvec_from_hex("01ab02cd03efabcdef");
  COMPARE(2, bitvec.second);

  // |          |          |          |
  //    12  40      80        14x
  Word word[3];
  bitvec_zero(word, 3);
  bitvec_set(12, word);
  bitvec_set(40, word);
  bitvec_set(80, word);
  bitvec_set(140, word);
  bitvec_set(141, word);
  bitvec_set(142, word);
  COMPARE(6, bitvec_count(word, 3));
  COMPARE(6, bitvec_count(word, 0, 192));
  COMPARE(5, bitvec_count(word, 30, 192));
  COMPARE(4, bitvec_count(word, 70, 192));
  COMPARE(3, bitvec_count(word, 90, 192));
  COMPARE(2, bitvec_count(word, 141, 192));
  COMPARE(1, bitvec_count(word, 142, 192));
  COMPARE(0, bitvec_count(word, 143, 192));

  COMPARE(6, bitvec_count(word, 0, 192));
  COMPARE(5, bitvec_count(word, 0, 142));
  COMPARE(4, bitvec_count(word, 0, 141));
  COMPARE(3, bitvec_count(word, 0, 130));
  COMPARE(2, bitvec_count(word, 0, 70));
  COMPARE(1, bitvec_count(word, 0, 30));
  COMPARE(0, bitvec_count(word, 0, 6));

  COMPARE(6, bitvec_count(word, 0, 192));
  COMPARE(4, bitvec_count(word, 30, 142));
  COMPARE(1, bitvec_count(word, 50, 100));
  COMPARE(0, bitvec_count(word, 45, 70));
  COMPARE(0, bitvec_count(word, 90, 130));

  bitvec_zero(word, 3);
  bitvec_set(4, word);
  bitvec_set(8, word);
  bitvec_set(12, word);
  bitvec_set(40, word);
  bitvec_set(80, word);
  bitvec_set(100, word);
  COMPARE(2, bitvec_count(word, 6, 30));
  COMPARE(1, bitvec_count(word, 70, 90));
  COMPARE(2, bitvec_count(word, 30, 90));

}
