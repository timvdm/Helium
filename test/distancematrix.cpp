#include <Helium/distancematrix.h>

#include "test.h"

using namespace Helium;

int main()
{
  DistanceMatrix m1(3, 0, 99);

  COMPARE(0, m1(0, 0));
  COMPARE(99, m1(0, 1));
  COMPARE(99, m1(0, 2));
  COMPARE(99, m1(1, 0));
  COMPARE(0, m1(1, 1));
  COMPARE(99, m1(1, 2));
  COMPARE(99, m1(2, 0));
  COMPARE(99, m1(2, 1));
  COMPARE(0, m1(2, 2));
}
