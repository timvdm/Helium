#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/diagram.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

PyObject* generate_diagram(const Molecule &mol)
{
  std::vector<std::pair<double, double> > coords = Helium::generate_diagram(mol);

  list *l = new list();
  for (std::size_t i = 0; i < coords.size(); ++i) {
    //tuple *t = new tuple(make_tuple(coords[i].first, coords[i].second));
    //l->append(t);
    l->append(coords[i]);
  }

  return l->ptr();
}

void export_diagram()
{

  def("generate_diagram", &generate_diagram);

}
