#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/dijkstra.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

void export_dijkstra()
{

  class_<Helium::Dijkstra<Molecule> >("Dijkstra", init<Molecule, Molecule::atom_type>())
    //.def(init<Molecule, Molecule::atom_type, bool>())
    .def("infinity", &Helium::Dijkstra<Molecule>::infinity)
    .def("distance", &Helium::Dijkstra<Molecule>::distance<Molecule::atom_type>)
    .def("path", &Helium::Dijkstra<Molecule>::path<Molecule::atom_type>)
    ;

}
