#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/enumeratepaths.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

void export_enumerate_paths()
{

  def("enumerate_paths", &Helium::enumerate_paths<Molecule>);

}
