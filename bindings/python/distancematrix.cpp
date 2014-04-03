#include <boost/python.hpp>

#include "../../src/distancematrix.h"

using namespace boost::python;

Helium::Size DistanceMatrix_call(const Helium::DistanceMatrix &m, Helium::Size i, Helium::Size j)
{
  return m(i, j);
}

std::string DistanceMatrix_str(const Helium::DistanceMatrix &m)
{
  std::stringstream ss;
  ss << m;
  return ss.str();
}

void export_distance_matrix()
{

  class_<Helium::DistanceMatrix, boost::noncopyable>("DistanceMatrix", no_init)
    .def("__str__", &DistanceMatrix_str)
    .def("dim", &Helium::DistanceMatrix::dim)
    .def("__call__", &DistanceMatrix_call)
    .def("infinity", &Helium::DistanceMatrix::infinity)
    .staticmethod("infinity")
    ;

}
