#ifndef HELIUM_TIMEOUT_H
#define HELIUM_TIMEOUT_H

#include <boost/timer/timer.hpp>

namespace Helium {

  struct timeout_error : std::exception {};

  class Timeout
  {
    public:
      static const boost::timer::nanosecond_type one_milisecond = 1000000L;

      /**
       * @param ms The number of miliseconds after which a timeout exception
       * will be thrown.
       */
      Timeout(unsigned int ms) : m_ns(ms * one_milisecond)
      {
      }

      void check()
      {
        boost::timer::cpu_times elapsed = timer.elapsed();
        boost::timer::nanosecond_type ns = elapsed.system + elapsed.user;

        if (ns > m_ns)
          throw timeout_error();
      }

    private:
      boost::timer::cpu_timer timer;
      unsigned long m_ns;
  };

}

#endif
