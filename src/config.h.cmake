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

#cmakedefine HAVE_CPP11
#cmakedefine HAVE_OPENCL

#ifdef _MSC_VER

  // Supress warning on forcing int etc. value to bool 'true' or 'false' (performance warning)
 #pragma warning( disable : 4800 )
 // compare signed/unsigned mismatch
 #pragma warning( disable : 4018 )

#endif  // _MSC_VER

#endif // HE_CONFIG_H
