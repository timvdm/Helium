#include <boost/python.hpp>

#include "../../src/molecule.h"
#include "../../src/bitvec.h"

using namespace boost::python;

#define DATA_MEMBER_TO_FUNCTION(Klass, Return, member) \
  const Return& member(const Klass &t) \
  { \
    return t.member; \
  }

template<typename T>
std::vector<T> vector_from_list(const list &l)
{
  std::vector<T> result;

  ssize_t size = len(l);
  for (ssize_t i = 0; i < size; ++i)
    result.push_back(extract<T>(l[i]));

  return result;
}

inline object pass_through(object const& o)
{
  return o;
}

template<typename Item, typename Iterator>
struct IteratorWrapper
{
  static Item next(Helium::iterator_pair<Iterator> &iters)
  {
    if (iters.begin() == iters.end()) {
      PyErr_SetString(PyExc_StopIteration, "No more data.");
      boost::python::throw_error_already_set();
    }

    Iterator tmp(iters.begin());
    ++iters.begin();
    return *tmp;
  }

  static void wrap(const char* python_name)
  {
    class_<Helium::iterator_pair<Iterator> >(python_name, no_init)
      .def("next", next)
      .def("__iter__", pass_through)
      ;
  }

};

struct Fingerprint
{
  Fingerprint(Helium::Word *data_, int numWords_, bool owned_)
    : data(data_), numWords(numWords_), owned(owned_)
  {
  }

  Fingerprint(int numWords_) : data(new Helium::Word[numWords_]),
      numWords(numWords_), owned(true)
  {
    Helium::bitvec_zero(data, numWords);
  }

  ~Fingerprint()
  {
    if (owned) {
      delete [] data;
      data = 0;
      numWords = 0;
    }
  }

  Helium::Word *data;
  int numWords;
  bool owned;
};
