#ifndef HELIUM_UTIL_VECTOR_H
#define HELIUM_UTIL_VECTOR_H

#include <vector>
#include <set>
#include <ostream>
#include <cassert>
#include <algorithm>

namespace std {

  //@cond dev

  template<typename T>
  std::ostream& operator<<(std::ostream &os, const std::vector<T> &v)
  {
    os << "[ ";
    for (std::size_t i = 0; i < v.size(); ++i)
      os << v[i] << " ";
    os << "]";
    return os;
  }

  template<typename T>
  bool operator<(const std::vector<T> &v1, const std::vector<T> &v2)
  {
    if (v1.size() < v2.size())
      return true;
    for (std::size_t i = 0; i < v1.size(); ++i) {
      if (v1[i] < v2[i])
        return true;
      if (v1[i] > v2[i])
        return false;
    }
    return false;
  }

  template<typename T>
  bool operator>(const std::vector<T> &v1, const std::vector<T> &v2)
  {
    if (v1.size() > v2.size())
      return true;
    for (std::size_t i = 0; i < v1.size(); ++i) {
      if (v1[i] > v2[i])
        return true;
      if (v1[i] < v2[i])
        return false;
    }
    return false;
  }

  template<typename T>
  bool operator==(const std::vector<T> &v1, const std::vector<T> &v2)
  {
    if (v1.size() != v2.size())
      return false;
    for (std::size_t i = 0; i < v1.size(); ++i)
      if (v1[i] != v2[i])
        return false;
    return true;
  }  

}

namespace Helium {

  template<typename T>
  std::size_t unique_elements(const std::vector<T> &v)
  {
    std::set<T> set;
    for (std::size_t i = 0; i < v.size(); ++i)
      set.insert(v[i]);
    return set.size();
  }

  template<typename T>
  void renumber(std::vector<T> &v)
  {
    std::set<T> values;
    for (std::size_t i = 0; i < v.size(); ++i)
      values.insert(v[i]);
    T newValue = 0;
    for (typename std::set<T>::iterator i = values.begin(); i != values.end(); ++i) {
      for (std::size_t j = 0; j < v.size(); ++j)
        if (v[j] == *i)
          v[j] = newValue;
      ++newValue;
    }
  }

  template<typename T1, typename T2>
  bool contains(const std::vector<T1> &v, const T2 &value)
  {
    return std::find(v.begin(), v.end(), value) != v.end();
  }

  template<typename T>
  std::size_t index_of(const std::vector<T> &v, const typename std::vector<T>::value_type &value)
  {
    typename std::vector<T>::const_iterator pos = std::find(v.begin(), v.end(), value);
    assert(pos != v.end());
    return pos - v.begin();
  }

  //@endcond

}

#endif
