#include <boost/python.hpp>

#include "../../src/error.h"

using namespace boost::python;

void export_error()
{

  class_<Helium::Error>("Error", no_init)
    .def("__nonzero__", &Helium::Error::operator bool) // python3: __bool__
    .def("__str__", &Helium::Error::what, return_value_policy<copy_const_reference>())
    ;

}
