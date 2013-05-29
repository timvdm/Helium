#ifndef HELIUM_UTIL_TYPETRAITS_H
#define HELIUM_UTIL_TYPETRAITS_H

namespace Helium {

  //@cond dev

  struct NullType {};
  typedef NullType EmptyType;

  template<int Expr>
  struct StaticAssert;
  template<>
  struct StaticAssert<true> {};

  template<typename T>
  struct IsReference { enum { result = false }; };
  template<typename T>
  struct IsReference<T&> { enum { result = true }; };

  template<typename T>
  struct RemoveReference { typedef T result; };
  template<typename T>
  struct RemoveReference<T&> { typedef T result; };

  template<typename T>
  struct IsPointer { enum { result = false }; };
  template<typename T>
  struct IsPointer<T*> { enum { result = true }; };

  template<typename T>
  struct RemovePointer { typedef T result; };
  template<typename T>
  struct RemovePointer<T*> { typedef T result; };

  template<typename T>
  struct IsConst { enum { result = false }; };
  template<typename T>
  struct IsConst<const T> { enum { result = true }; };

  template<typename T>
  struct RemoveConst { typedef T result; };
  template<typename T>
  struct RemoveConst<const T> { typedef T result; };

  template<typename T1, typename T2>
  struct SameType { enum { result = false }; };
  template<typename T> 
  struct SameType<T, T> { enum { result = true }; };
  template<typename T> 
  struct SameType<const T, T> { enum { result = true }; };
  template<typename T> 
  struct SameType<T, const T> { enum { result = true }; };
  template<typename T> 
  struct SameType<T&, T> { enum { result = true }; };
  template<typename T> 
  struct SameType<T, T&> { enum { result = true }; };
  template<typename T> 
  struct SameType<T*, T> { enum { result = true }; };
  template<typename T> 
  struct SameType<T, T*> { enum { result = true }; };









  template<typename T>
  struct IsInteger { enum { result = false }; };
  /// @cond IMPL
  template<> struct IsInteger<char> { enum { result = true }; };
  template<> struct IsInteger<signed char> { enum { result = true }; };
  template<> struct IsInteger<unsigned char> { enum { result = true }; };
  template<> struct IsInteger<signed short> { enum { result = true }; };
  template<> struct IsInteger<unsigned short> { enum { result = true }; };
  template<> struct IsInteger<signed int> { enum { result = true }; };
  template<> struct IsInteger<unsigned int> { enum { result = true }; };
  template<> struct IsInteger<signed long> { enum { result = true }; };
  template<> struct IsInteger<unsigned long> { enum { result = true }; };
  /// @endcond

  /**
   * Integer to type template. This is used to turn an integer to type and is
   * mainly used to dispatch to helper functions.
   *
   * Example:
   * @code
   * template<typename T>
   * void fooHelper(T v, Int2Type<true>)
   * {
   *   // v is a pointer
   *   v->bar();
   * }
   *
   * template<typename T>
   * void fooHelper(T v, Int2Type<false>)
   * {
   *   // v is a reference
   *   v.bar();
   * }
   *
   * template<typename T>
   * void foo(T v)
   * {
   *   fooHelper(v, Int2Type<IsPointer<T>::result>());
   * }
   * @endcode
   */
  template<int N>
  struct Int2Type { enum { result = N }; };
  template<typename T>
  struct Type2Type { typedef T result; };

  template<typename IfTrueType, typename IfFalseType, int Expr>
  struct SelectHelper { typedef IfTrueType result; };
  template<typename IfTrueType, typename IfFalseType>
  struct SelectHelper<IfTrueType, IfFalseType, false> { typedef IfFalseType result; };
  /**
   * Select a type based on an value of an expression. This works like the
   * expr ? true_value : false_value statement in C++ but operates on types
   * instead.
   *
   * Example:
   * @code
   * // ? : expression with values
   * bool foo = true;
   * int bar = foo ? 0 : 42; // bar = 0;
   *
   * Select<false, int, int*>::result v; // v is of type int*
   * @endcode
   */
  template<int Expr, typename IfTrueType, typename IfFalseType>
  struct Select { typedef typename SelectHelper<IfTrueType, IfFalseType, Expr>::result result; };

  //@endcond

}

#endif
