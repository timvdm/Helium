#ifndef HELIUM_UTIL_H
#define HELIUM_UTIL_H

#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>

#include <iostream>

#include "util/typetraits.h"
#include "util/fileio.h"
#include "util/string.h"
#include "util/vector.h"
#include "util/functor.h"

namespace Helium {

  template<typename T>
  std::ostream& operator<<(std::ostream &os, const std::set<T> &s)
  {
    os << "set[ ";
    for (typename std::set<T>::const_iterator i = s.begin(); i != s.end(); ++i)
      os << *i << " ";
    os << "]";
    return os;
  }

  template<typename T1, typename T2>
  std::ostream& operator<<(std::ostream &os, const std::pair<T1, T2> &p)
  {
    os << "( " << p.first << " " << p.second << " )";
    return os;
  }

  inline unsigned int factorial(unsigned int n)
  {
    unsigned int result = 1;
    for (int i = 2; i <= n; ++i)
      result *= i;
    return result;
  }

  inline unsigned int num_permutations(unsigned int n)
  {
    return factorial(n);
  }

  inline unsigned int num_combinations(unsigned int n, unsigned int k)
  {
    return factorial(n) / (factorial(n - k) * factorial(k));
  }

  template <typename Iterator>
  bool next_combination(const Iterator first, Iterator k, const Iterator last)
  {
    /* Credits: Mark Nelson http://marknelson.us */
    if ((first == last) || (first == k) || (last == k))
      return false;
    Iterator i1 = first;
    Iterator i2 = last;
    ++i1;
    if (last == i1)
      return false;
    i1 = last;
    --i1;
    i1 = k;
    --i2;
    while (first != i1)
    {
      if (*--i1 < *i2)
      {
        Iterator j = k;
        while (!(*i1 < *j)) ++j;
        std::iter_swap(i1,j);
        ++i1;
        ++j;
        i2 = k;
        std::rotate(i1,j,last);
        while (last != j)
        {
          ++j;
          ++i2;
        }
        std::rotate(k,i2,last);
        return true;
      }
    }
    std::rotate(first,k,last);
    return false;
  }


}

#endif
