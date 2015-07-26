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
#ifndef HE_CONFIG_H
#define HE_CONFIG_H

// used to export symbols
#if defined(WIN32) && !defined(__MINGW32__)
  #define HE_EXPORT __declspec(dllexport)
  #define HE_IMPORT __declspec(dllimport)
  #define HE_HIDDEN
#else
  #define HE_EXPORT
  #define HE_IMPORT
  #define HE_HIDDEN
#endif

#define HEAPI HE_EXPORT

#define HEDATADIR "@HEDATADIR@"

#define HEUNUSED(var) (void)var

#cmakedefine HAVE_CPP11
#cmakedefine HAVE_OPENCL

#ifdef _MSC_VER
  // Supress warning on forcing int etc. value to bool 'true' or 'false' (performance warning)
  #pragma warning( disable : 4800 )
  // compare signed/unsigned mismatch
  #pragma warning( disable : 4018 )
#endif  // _MSC_VER

namespace Color {
#ifndef WIN32
  inline const char *white() { return "\033[0;37m"; }
  inline const char *red() { return "\033[0;31m"; }
  inline const char *green() { return "\033[0;32m"; }
  inline const char *yellow() { return "\033[0;33m"; }
  inline const char *blue() { return "\033[0;34m"; }
  inline const char *magenta() { return "\033[0;35m"; }
  inline const char *cyan() { return "\033[0;36m"; }
#else
  inline const char *white() { return ""; }
  inline const char *red() { return ""; }
  inline const char *green() { return ""; }
  inline const char *yellow() { return ""; }
  inline const char *blue() { return ""; }
  inline const char *magenta() { return ""; }
  inline const char *cyan() { return ""; }
#endif
}

#endif // HE_CONFIG_H
