#include <boost/python.hpp>

#include "../../src/fileio/fingerprints.h"
#include "../common.h"

using namespace boost::python;

void RowMajor_writeFingerprint(Helium::RowMajorFingerprintOutputFile &file, const Fingerprint &fp)
{
  file.writeFingerprint(fp.data);
}

void ColumnMajor_writeFingerprint(Helium::ColumnMajorFingerprintOutputFile &file, const Fingerprint &fp)
{
  file.writeFingerprint(fp.data);
}

Fingerprint* RowMajor_fingerprint(const Helium::InMemoryRowMajorFingerprintStorage &storage, unsigned int index)
{
  return new Fingerprint(storage.fingerprint(index), Helium::bitvec_num_words_for_bits(storage.numBits()), false);
}

Fingerprint* ColumnMajor_bit(const Helium::InMemoryColumnMajorFingerprintStorage &storage, unsigned int index)
{
  return new Fingerprint(storage.bit(index), Helium::bitvec_num_words_for_bits(storage.numFingerprints()), false);
}

void export_fingerprintfiles()
{

  class_<Helium::RowMajorFingerprintOutputFile, boost::noncopyable>("RowMajorFingerprintOutputFile",
      init<const std::string&, unsigned int>())
    .def("writeFingerprint", &RowMajor_writeFingerprint)
    .def("writeHeader", &Helium::RowMajorFingerprintOutputFile::writeHeader)
    .def("error", &Helium::RowMajorFingerprintOutputFile::error, return_internal_reference<>())
    ;

  class_<Helium::ColumnMajorFingerprintOutputFile, boost::noncopyable>("ColumnMajorFingerprintOutputFile",
      init<const std::string&, unsigned int, unsigned int>())
    .def("writeFingerprint", &ColumnMajor_writeFingerprint)
    .def("writeHeader", &Helium::ColumnMajorFingerprintOutputFile::writeHeader)
    .def("error", &Helium::ColumnMajorFingerprintOutputFile::error, return_internal_reference<>())
    ;

  class_<Helium::InMemoryRowMajorFingerprintStorage, boost::noncopyable>("InMemoryRowMajorFingerprintStorage")
    .def("load", &Helium::InMemoryRowMajorFingerprintStorage::load)
    .def("header", &Helium::InMemoryRowMajorFingerprintStorage::header)
    .def("numBits", &Helium::InMemoryRowMajorFingerprintStorage::numBits)
    .def("numFingerprints", &Helium::InMemoryRowMajorFingerprintStorage::numFingerprints)
    .def("fingerprint", &RowMajor_fingerprint, return_value_policy<manage_new_object>())
    .def("error", &Helium::InMemoryRowMajorFingerprintStorage::error, return_internal_reference<>())
    ;

  class_<Helium::InMemoryColumnMajorFingerprintStorage, boost::noncopyable>("InMemoryColumnMajorFingerprintStorage")
    .def("load", &Helium::InMemoryColumnMajorFingerprintStorage::load)
    .def("header", &Helium::InMemoryColumnMajorFingerprintStorage::header)
    .def("numBits", &Helium::InMemoryColumnMajorFingerprintStorage::numBits)
    .def("numFingerprints", &Helium::InMemoryColumnMajorFingerprintStorage::numFingerprints)
    .def("bit", &ColumnMajor_bit, return_value_policy<manage_new_object>())
    .def("error", &Helium::InMemoryColumnMajorFingerprintStorage::error, return_internal_reference<>())
    ;

}
