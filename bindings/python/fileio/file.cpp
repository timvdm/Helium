#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/fileio/file.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

void export_file()
{

  class_<Helium::BinaryInputFile, boost::noncopyable>("BinaryInputFile")
    .def(init<const std::string&>())
    .def("open", &Helium::BinaryInputFile::open)
    .def("close", &Helium::BinaryInputFile::close)
    .def("__nonzero__", &Helium::BinaryInputFile::operator bool)
    .def("header", &Helium::BinaryInputFile::header, return_value_policy<copy_const_reference>())
    ;

}
