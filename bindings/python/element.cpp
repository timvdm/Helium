#include <boost/python.hpp>

#include "../../src/element.h"

using namespace boost::python;

void export_element()
{

  class_<Helium::Element>("Element", no_init)
    .def("symbol", &Helium::Element::symbol)
    .staticmethod("symbol")
    .def("element", &Helium::Element::element)
    .staticmethod("element")
    .def("averageMass", &Helium::Element::averageMass)
    .staticmethod("averageMass")
    .def("addHydrogens", &Helium::Element::addHydrogens)
    .staticmethod("addHydrogens")
    .def("valence", &Helium::Element::valence)
    .staticmethod("valence")
    ;

}
