#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/fileio/fps.h"
#include "../common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

Fingerprint* fingerprint(const Helium::FpsFile &file, unsigned int index)
{
  return new Fingerprint(file.fingerprint(index), Helium::bitvec_num_words_for_bits(file.numBits()), false);
}

void export_fps()
{

  class_<Helium::FpsFile>("FpsFile")
    .def("load", &Helium::FpsFile::load)
    .def("type", &Helium::FpsFile::type, return_value_policy<copy_const_reference>())
    .def("numBits", &Helium::FpsFile::numBits)
    .def("numFingerprints", &Helium::FpsFile::numFingerprints)
    .def("fingerprint", &fingerprint, return_value_policy<manage_new_object>())
    ;

}
