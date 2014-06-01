#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/bitvec.h"
#include "common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

void bitvec_zero(Fingerprint &fp)
{
  return Helium::bitvec_zero(fp.data, fp.numWords);
}

Fingerprint* bitvec_copy(const Fingerprint &fp)
{
  return new Fingerprint(Helium::bitvec_copy(fp.data, fp.numWords), fp.numWords, true);
}

bool bitvec_get(int index, const Fingerprint &fp)
{
  if (index >= Helium::BitsPerWord * fp.numWords)
    throw std::runtime_error("Invalid fingerprint bit index");
  return Helium::bitvec_get(index, fp.data);
}

void bitvec_set(int index, Fingerprint &fp)
{
  if (index >= Helium::BitsPerWord * fp.numWords)
    throw std::runtime_error("Invalid fingerprint bit index");
  return Helium::bitvec_set(index, fp.data);
}

void bitvec_reset(int index, Fingerprint &fp)
{
  if (index >= Helium::BitsPerWord * fp.numWords)
    throw std::runtime_error("Invalid fingerprint bit index");
  return Helium::bitvec_reset(index, fp.data);
}

bool bitvec_is_subset_superset(const Fingerprint &fp1, const Fingerprint &fp2)
{
  return Helium::bitvec_is_subset_superset(fp1.data, fp2.data, fp1.numWords);
}

int bitvec_count_1(const Fingerprint &fp)
{
  return Helium::bitvec_count(fp.data, fp.numWords);
}

int bitvec_count_2(const Fingerprint &fp, int begin, int end)
{
  if (begin >= Helium::BitsPerWord * fp.numWords)
    throw std::runtime_error("Invalid fingerprint bit index for begin");
  if (end >= Helium::BitsPerWord * fp.numWords + 1)
    throw std::runtime_error("Invalid fingerprint bit index for end");
  return Helium::bitvec_count(fp.data, begin, end);
}

int bitvec_union_count(const Fingerprint &fp1, const Fingerprint &fp2)
{
  return Helium::bitvec_union_count(fp1.data, fp2.data, fp1.numWords);
}

double bitvec_tanimoto_1(const Fingerprint &fp1, const Fingerprint &fp2)
{
  return Helium::bitvec_tanimoto(fp1.data, fp2.data, fp1.numWords);
}

double bitvec_tanimoto_2(const Fingerprint &fp1, const Fingerprint &fp2, int bitCount1, int bitCount2)
{
  return Helium::bitvec_tanimoto(fp1.data, fp2.data, bitCount1, bitCount2, fp1.numWords);
}

double bitvec_cosine_1(const Fingerprint &fp1, const Fingerprint &fp2, int bitCount1, int bitCount2)
{
  return Helium::bitvec_cosine(fp1.data, fp2.data, bitCount1, bitCount2, fp1.numWords);
}

double bitvec_cosine_2(const Fingerprint &fp1, const Fingerprint &fp2)
{
  return Helium::bitvec_cosine(fp1.data, fp2.data, fp1.numWords);
}

double bitvec_hamming_1(const Fingerprint &fp1, const Fingerprint &fp2, int bitCount1, int bitCount2)
{
  return Helium::bitvec_hamming(fp1.data, fp2.data, bitCount1, bitCount2, fp1.numWords);
}

double bitvec_hamming_2(const Fingerprint &fp1, const Fingerprint &fp2)
{
  return Helium::bitvec_hamming(fp1.data, fp2.data, fp1.numWords);
}

double bitvec_russell_rao(const Fingerprint &fp1, const Fingerprint &fp2)
{
  return Helium::bitvec_russell_rao(fp1.data, fp2.data, fp1.numWords);
}

double bitvec_forbes_1(const Fingerprint &fp1, const Fingerprint &fp2, int bitCount1, int bitCount2)
{
  return Helium::bitvec_forbes(fp1.data, fp2.data, bitCount1, bitCount2, fp1.numWords);
}

double bitvec_forbes_2(const Fingerprint &fp1, const Fingerprint &fp2)
{
  return Helium::bitvec_forbes(fp1.data, fp2.data, fp1.numWords);
}

Fingerprint* bitvec_from_binary(const std::string &binary)
{
  std::pair<Helium::Word*, int> fp = Helium::bitvec_from_binary(binary);
  return new Fingerprint(fp.first, fp.second, true);
}

std::string bitvec_to_binary(const Fingerprint &fp)
{
  return Helium::bitvec_to_binary(fp.data, fp.numWords);
}

Fingerprint* bitvec_from_hex(const std::string &hex)
{
  std::pair<Helium::Word*, int> fp = Helium::bitvec_from_hex(hex);
  return new Fingerprint(fp.first, fp.second, true);
}

std::string bitvec_to_hex(const Fingerprint &fp)
{
  return Helium::bitvec_to_hex(fp.data, fp.numWords);
}

std::string hex(const Fingerprint &fp)
{
  return Helium::bitvec_to_hex(fp.data, fp.numWords);
}

std::string bin_1(const Fingerprint &fp)
{
  return Helium::bitvec_to_binary(fp.data, fp.numWords);
}

std::string bin_2(const Fingerprint &fp, bool spaces)
{
  return Helium::bitvec_to_binary(fp.data, fp.numWords, spaces);
}

bool Fingerprint_get(const Fingerprint &fp, int index)
{
  if (index >= Helium::BitsPerWord * fp.numWords)
    throw std::runtime_error("Invalid fingerprint bit index");
  return Helium::bitvec_get(index, fp.data);
}

void Fingerprint_set(Fingerprint &fp, int index , bool value)
{
  if (index >= Helium::BitsPerWord * fp.numWords)
    throw std::runtime_error("Invalid fingerprint bit index");
  value ? Helium::bitvec_set(index, fp.data) : Helium::bitvec_reset(index, fp.data);
}

void export_bitvec()
{

  class_<Fingerprint, boost::noncopyable>("Fingerprint", init<int>())
    .def("bin", &bin_1)
    .def("bin", &bin_2)
    .def("hex", &hex)
    .def("zero", &bitvec_zero)
    .def("count", &bitvec_count_1)
    .def("count", &bitvec_count_2)
    .def("__getitem__", &Fingerprint_get)
    .def("__setitem__", &Fingerprint_set)
    .def_readonly("numWords", &Fingerprint::numWords)
    ;

  scope().attr("BitsPerWord") = Helium::BitsPerWord;

  def("bitvec_num_words_for_bits", &Helium::bitvec_num_words_for_bits);
  def("bitvec_zero", &bitvec_zero);
  def("bitvec_copy", &bitvec_copy, return_value_policy<manage_new_object>());
  def("bitvec_get", &bitvec_get);
  def("bitvec_set", &bitvec_set);
  def("bitvec_reset", &bitvec_reset);
  def("bitvec_is_subset_superset", &bitvec_is_subset_superset);
  def("bitvec_count", &bitvec_count_1);
  def("bitvec_count", &bitvec_count_2);
  def("bitvec_union_count", &bitvec_union_count);
  def("bitvec_tanimoto", &bitvec_tanimoto_1);
  def("bitvec_tanimoto", &bitvec_tanimoto_2);
  def("bitvec_cosine", &bitvec_cosine_1);
  def("bitvec_cosine", &bitvec_cosine_2);
  def("bitvec_hamming", &bitvec_hamming_1);
  def("bitvec_hamming", &bitvec_hamming_2);
  def("bitvec_russell_rao", &bitvec_russell_rao);
  def("bitvec_forbes", &bitvec_forbes_1);
  def("bitvec_forbes", &bitvec_forbes_2);
  def("bitvec_to_binary", &bin_1);
  def("bitvec_to_binary", &bin_2);
  def("bitvec_to_hex", &hex);
  def("bitvec_from_hex", &bitvec_from_hex, return_value_policy<manage_new_object>());
  def("bitvec_from_binary", &bitvec_from_binary, return_value_policy<manage_new_object>());

}
