/*
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
#ifndef HELIUM_UTIL_H
#define HELIUM_UTIL_H

#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <algorithm>
#include <cassert>

#include <iostream>

#include <Helium/config.h>
#include <Helium/util/typetraits.h>
#include <Helium/util/fileio.h>
#include <Helium/util/string.h>
#include <Helium/util/vector.h>
#include <Helium/util/functor.h>

#define UNREACHABLE_RETURN_REF(type) \
  assert(0); \
  return *(new type);

#define ENABLE_TIMERS
#ifdef ENABLE_TIMERS
  #include <boost/timer/timer.hpp>
  #define TIMER(message) \
    std::cout << message; \
    boost::timer::auto_cpu_timer t;
#else
  #define TIMER(message)
#endif

namespace std {

  /**
   * @brief STL output stream operator for std::set.
   */
  template<typename T>
  std::ostream& operator<<(std::ostream &os, const std::set<T> &s)
  {
    os << "set[ ";
    for (typename std::set<T>::const_iterator i = s.begin(); i != s.end(); ++i)
      os << *i << " ";
    os << "]";
    return os;
  }

  /**
   * @brief STL output stream operator for std::pair.
   */
  template<typename T1, typename T2>
  std::ostream& operator<<(std::ostream &os, const std::pair<T1, T2> &p)
  {
    os << "( " << p.first << " " << p.second << " )";
    return os;
  }

}

namespace Helium {

  /**
   * @brief Compute the factorial of a mumber.
   *
   * @return The factorial of @p n.
   */
  inline unsigned int factorial(unsigned int n)
  {
    unsigned int result = 1;
    for (int i = 2; i <= n; ++i)
      result *= i;
    return result;
  }

  /**
   * @brief Check if a number is prime.
   *
   * @return True if @p n is a prime number.
   */
  inline unsigned long is_prime(unsigned long n)
  {
    for (unsigned long i = 2; i < n; ++i)
      if (!(n % i))
        return false;
    return true;
  }

  /**
   * @brief Get the largest prime less than or equal to a given number.
   *
   * @return The the largest prime less than or equal to n.
   */
  inline unsigned long previous_prime(unsigned long n)
  {
    while (!is_prime(n))
      --n;
    return n;
  }

  /**
   * @brief Compute the number of permutations of @p n obejcts.
   *
   * @return The number of permutations of @p n objects.
   */
  inline unsigned int num_permutations(unsigned int n)
  {
    return factorial(n);
  }

  /**
   * @brief Compute the number of combinations for @p n objects taken @p k at a
   * time.
   *
   * @return The number of combinations.
   */
  inline unsigned int num_combinations(unsigned int n, unsigned int k)
  {
    return factorial(n) / (factorial(n - k) * factorial(k));
  }

  /**
   * @brief Compute the next combination.
   */
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
