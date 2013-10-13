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
#ifndef HELIUM_TEST_H
#define HELIUM_TEST_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#ifdef _MSC_VER
#define FUNCTION_SIGNATURE __FUNCSIG__
#else
#define FUNCTION_SIGNATURE __PRETTY_FUNCTION__
#endif

/*inline*/ void report_error(const char* msg, const char* file, int line, const char* func_name, bool require = false)
{
    std::cout << file << ":" << line << ": " << msg << " (FAIL)" << std::endl;
    if (require)
      exit(-1);
}

template <typename T1, typename T2>
void __test_compare__(T1 a, T2 b, const char *expr, const char *file, int line, const char *func_name)
{
  if (!(a == b))
    std::cout << file << ":" << line << ": " << expr << " [" << a << " == " << b << "] (FAIL)" << std::endl;
}

#define ASSERT(exp) \
  ( (exp) ? static_cast<void>(0) : report_error(#exp, __FILE__, __LINE__, FUNCTION_SIGNATURE, false) )

#define REQUIRE(exp) \
  ( (exp) ? static_cast<void>(0) : report_error(#exp, __FILE__, __LINE__, FUNCTION_SIGNATURE, true) )

const char* __test_expr__(const char *expr) { return expr; }
#define P_EXPR(expr) __test_expr__(#expr)

#define COMPARE(a,b) \
  __test_compare__(a, b, P_EXPR( a == b ), __FILE__, __LINE__, FUNCTION_SIGNATURE)


#ifndef HEDATADIR
#define HEDATADIR "./"
#endif

std::string datadir()
{
  return HEDATADIR + std::string("/");
}

inline void compare_file(const std::string &filename, const std::string &exp)
{
  std::ifstream ifs(filename.c_str());
  ASSERT(ifs);

  std::string line;
  std::stringstream correct;
  while (std::getline(ifs, line))
    correct << line << std::endl;

  COMPARE(correct.str(), exp);
}

inline void compare_file(const std::string &filename, const std::stringstream &exp)
{
  compare_file(filename, exp.str());
}

#endif
