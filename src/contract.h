#ifndef HELIUM_CONTRACT_H
#define HELIUM_CONTRACT_H

#ifndef NDEBUG

#include <cassert>

#define PRE(expr) \
  assert(expr)

#define POST(expr) \
  assert(expr)

#else

#define PRE(expr)

#define POST(expr)

#endif

#endif
