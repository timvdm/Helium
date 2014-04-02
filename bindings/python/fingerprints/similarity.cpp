#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/fingerprints/similarity.h"
#include "../../src/fileio/fingerprints.h"
#include "../common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

std::vector<std::pair<unsigned int, double> > brute_force_similarity_search(const Fingerprint &query,
    Helium::InMemoryRowMajorFingerprintStorage &storage, double Tmin)
{
  return Helium::brute_force_similarity_search(query.data, storage, Tmin);
}

std::vector<std::pair<unsigned int, double> >
search(const Helium::SimilaritySearchIndex<Helium::InMemoryRowMajorFingerprintStorage> &storage,
    const Fingerprint &query, double threshold, unsigned int maxResults = 0)
{
  return storage.search(query.data, threshold, maxResults);
}

BOOST_PYTHON_FUNCTION_OVERLOADS(search_overloads, search, 3, 4);

void export_similarity()
{

  def("brute_force_similarity_search", &brute_force_similarity_search);

  class_<Helium::SimilaritySearchIndex<Helium::InMemoryRowMajorFingerprintStorage>, boost::noncopyable>("SimilaritySearchIndex",
      init<const Helium::InMemoryRowMajorFingerprintStorage&, int>())
    .def("search", search, search_overloads())
    ;

}
