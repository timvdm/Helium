#ifndef HELIUM_TIE_H
#define HELIUM_TIE_H

namespace Helium {

  template<typename T1, typename T2>
  class tie_impl
  {
    public:
      tie_impl(T1 &first, T2 &second) : m_first(first), m_second(second)
      {
      }

      void operator=(const std::pair<T1, T2> &p)
      {
        m_first = p.first;
        m_second = p.second;
      }

    private:
      T1 &m_first;
      T2 &m_second;
  };

  template<typename T1, typename T2>
  tie_impl<T1, T2> tie(T1 &first, T2 &second)
  {
    return tie_impl<T1, T2>(first, second);
  }

}

#endif
