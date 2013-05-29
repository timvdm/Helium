#ifndef HELIUM_UTIL_FUNCTOR_H
#define HELIUM_UTIL_FUNCTOR_H

#include <algorithm>

namespace Helium {

  //@cond dev

  template<typename T1, typename T2, template<typename> class Compare = std::less>
  struct compare_first
  {
    bool operator()(const std::pair<T1, T2> &left, const std::pair<T1, T2> &right) const
    {
      return compare(left.first, right.first);
    }

    Compare<T1> compare;
  };

  template<typename T1, typename T2, template<typename> class Compare = std::less>
  struct compare_second
  {
    bool operator()(const std::pair<T1, T2> &left, const std::pair<T1, T2> &right) const
    {
      return compare(left.second, right.second);
    }

    Compare<T2> compare;
  };


  template<typename T1, typename T2, template<typename> class Compare = std::equal_to>
  struct find_first
  {
    find_first(const T1 &value_) : value(value_)
    {
    }

    bool operator()(const std::pair<T1, T2> &p) const 
    {
      return compare(p.first, value);
    } 
    
    Compare<T1> compare;
    const T1 &value;
  };

  template<typename T1, typename T2, template<typename> class Compare = std::equal_to>
  struct find_second
  {
    find_second(const T2 &value_) : value(value_)
    {
    }

    bool operator()(const std::pair<T1, T2> &p) const 
    {
      return compare(p.second, value);
    } 
    
    Compare<T2> compare;
    const T2 &value;
  };


  template<typename T>
  struct ContainerSize
  {
    typedef typename T::size_type property_type;

    bool operator()(const T &c) const
    {
      return c.size();
    }
  };

  template<typename T>
  struct ContainerMinElement
  {
    typedef typename T::value_type property_type;

    bool operator()(const T &c) const
    {
      return *std::min_element(c.begin(), c.end());
    }
  };

  template<typename T>
  struct ContainerMaxElement
  {
    typedef typename T::value_type property_type;

    bool operator()(const T &c) const
    {
      return *std::max_element(c.begin(), c.end());
    }
  };

  template<typename T, template<typename> class Property = ContainerSize, template<typename> class Compare = std::less>
  struct SortContainers
  {
    bool operator()(const T &left, const T &right) const
    {
      return compare(property(left), property(right));
    }

    Property<T> property;
    Compare<typename Property<T>::property_type> compare;
  };

  //@endcond

}

#endif
