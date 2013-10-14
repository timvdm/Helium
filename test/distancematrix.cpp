#include <Helium/distancematrix.h>

#include "test.h"

using namespace Helium;

int main()
{
  SymmetricDistanceMatrix m1(3, 0, 99);

  COMPARE(0, m1(0, 0));
  COMPARE(99, m1(0, 1));
  COMPARE(99, m1(0, 2));
  COMPARE(99, m1(1, 0));
  COMPARE(0, m1(1, 1));
  COMPARE(99, m1(1, 2));
  COMPARE(99, m1(2, 0));
  COMPARE(99, m1(2, 1));
  COMPARE(0, m1(2, 2));

  m1(0, 0) = 1;
  m1(1, 0) = 2;
  m1(1, 1) = 3;
  m1(2, 0) = 4;
  m1(2, 1) = 5;
  m1(2, 2) = 6;

  COMPARE(1, m1(0, 0));
  COMPARE(2, m1(0, 1));
  COMPARE(4, m1(0, 2));
  COMPARE(2, m1(1, 0));
  COMPARE(3, m1(1, 1));
  COMPARE(5, m1(1, 2));
  COMPARE(4, m1(2, 0));
  COMPARE(5, m1(2, 1));
  COMPARE(6, m1(2, 2));

  m1(0, 0) = 10;
  m1(0, 1) = 20;
  m1(0, 2) = 30;
  m1(1, 1) = 40;
  m1(1, 2) = 50;
  m1(2, 2) = 60;

  COMPARE(10, m1(0, 0));
  COMPARE(20, m1(0, 1));
  COMPARE(30, m1(0, 2));
  COMPARE(20, m1(1, 0));
  COMPARE(40, m1(1, 1));
  COMPARE(50, m1(1, 2));
  COMPARE(30, m1(2, 0));
  COMPARE(50, m1(2, 1));
  COMPARE(60, m1(2, 2));

}
