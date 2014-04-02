#include "common.h"
#include <vector>

#include "../../src/chemist/molecule.h"
#include "../../src/ring.h"

using Helium::Chemist::Molecule;

template<typename T1, typename T2>
struct PairToTuple
{
  static PyObject* convert(const std::pair<T1, T2> &p)
  {
    tuple *t = new tuple(make_tuple(p.first, p.second));
    return t->ptr();
  }
};

template<typename T>
struct VectorToList
{
  static PyObject* convert(const std::vector<T> &v)
  {
    list *l = new list();
    for (std::size_t i = 0; i < v.size(); ++i)
      l->append(T(v[i]));
    return l->ptr();
  }
};

template<typename T>
struct VectorVectorToList
{
  static PyObject* convert(const std::vector<std::vector<T> > &v)
  {
    list *l = new list();
    for (std::size_t i = 0; i < v.size(); ++i) {
      list l2;
      for (std::size_t j = 0; j < v[i].size(); ++j)
        l2.append(v[i][j]);
      l->append(v[i]);
    }
    return l->ptr();
  }
};

void export_common()
{
  to_python_converter<std::pair<unsigned int, double>, PairToTuple<unsigned int, double> >();

  to_python_converter<std::vector<bool>, VectorToList<bool> >();
  to_python_converter<std::vector<unsigned int>, VectorToList<unsigned int> >();
  to_python_converter<std::vector<std::pair<unsigned int, double> >, VectorToList<std::pair<unsigned int, double> > >();

  to_python_converter<std::vector<Molecule::atom_type>, VectorToList<Molecule::atom_type> >();
  to_python_converter<std::vector<Molecule::bond_type>, VectorToList<Molecule::bond_type> >();
  to_python_converter<std::vector<Helium::Ring<Molecule> >, VectorToList<Helium::Ring<Molecule> > >();

  to_python_converter<std::vector<std::vector<unsigned int> >, VectorVectorToList<unsigned int> >();
}
