#include <boost/python.hpp>

using namespace boost::python;

void export_molecule();

BOOST_PYTHON_MODULE(helium) {

  export_molecule();

}
