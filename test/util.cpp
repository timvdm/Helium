#include <Helium/util.h>

#include "test.h"

using namespace Helium;

void test_factorial()
{
  COMPARE(1, factorial(0));
  COMPARE(1, factorial(1));
  COMPARE(2, factorial(2));
  COMPARE(6, factorial(3));
  COMPARE(24, factorial(4));
}

void test_num_combinations()
{
  COMPARE(8, num_combinations(8, 1));
  COMPARE(28, num_combinations(8, 2));
  COMPARE(56, num_combinations(8, 3));
}

void test_combinations()
{
  std::string s = "ABCDEFGH";

  // 1 element
  int count = 0;
  do {
    ++count;
  } while (next_combination(s.begin(), s.begin() + 1, s.end()));

  COMPARE(8, count);

  // 2 elements
  count = 0;
  do {
    ++count;
  } while (next_combination(s.begin(), s.begin() + 2, s.end()));

  COMPARE(28, count);

  // 3 elements
  count = 0;
  do {
    ++count;
  } while (next_combination(s.begin(), s.begin() + 3, s.end()));

  COMPARE(56, count);
}

int main()
{
  test_factorial();
  test_num_combinations();
  test_combinations();
}
