#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/smarts.h"
#include "common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

bool search_1(Helium::Smarts &smarts, const Molecule &mol, Helium::NoMapping &mapping,
    const Helium::RingSet<Molecule> &rings, bool uniqueComponents = true)
{  return smarts.search(mol, mapping, rings, uniqueComponents); }

bool search_2(Helium::Smarts &smarts, const Molecule &mol, Helium::CountMapping &mapping,
    const Helium::RingSet<Molecule> &rings, bool uniqueComponents = true)
{  return smarts.search(mol, mapping, rings, uniqueComponents); }

bool search_3(Helium::Smarts &smarts, const Molecule &mol, Helium::SingleMapping &mapping,
    const Helium::RingSet<Molecule> &rings, bool uniqueComponents = true)
{  return smarts.search(mol, mapping, rings, uniqueComponents); }

bool search_4(Helium::Smarts &smarts, const Molecule &mol, Helium::MappingList &mapping,
    const Helium::RingSet<Molecule> &rings, bool uniqueComponents = true)
{  return smarts.search(mol, mapping, rings, uniqueComponents); }

BOOST_PYTHON_FUNCTION_OVERLOADS(search_overloads_1, search_1, 4, 5);
BOOST_PYTHON_FUNCTION_OVERLOADS(search_overloads_2, search_2, 4, 5);
BOOST_PYTHON_FUNCTION_OVERLOADS(search_overloads_3, search_3, 4, 5);
BOOST_PYTHON_FUNCTION_OVERLOADS(search_overloads_4, search_4, 4, 5);

DATA_MEMBER_TO_FUNCTION(Helium::SingleMapping, std::vector<unsigned int>, map);
DATA_MEMBER_TO_FUNCTION(Helium::MappingList, std::vector<std::vector<unsigned int> >, maps);

void export_smarts()
{

  class_<Helium::NoMapping>("NoMapping")
    .def_readonly("match", &Helium::NoMapping::match)
    ;

  class_<Helium::CountMapping>("CountMapping")
    .def_readonly("count", &Helium::CountMapping::count)
    ;

  class_<Helium::SingleMapping>("SingleMapping")
    .add_property("map", make_function(&map, return_value_policy<copy_const_reference>()))
    ;

  class_<Helium::MappingList>("MappingList")
    .add_property("maps", make_function(&maps, return_value_policy<copy_const_reference>()))
    ;

  class_<Helium::Smarts>("Smarts")
    .def("init", &Helium::Smarts::init)
    .def("error", &Helium::Smarts::error, return_internal_reference<>())
    .def("requiresCycles", &Helium::Smarts::requiresCycles)
    .def("requiresExplicitHydrogens", &Helium::Smarts::requiresExplicitHydrogens)
    .def("search", (bool(*)(Helium::Smarts&, const Molecule&, Helium::NoMapping&, const Helium::RingSet<Molecule>&, bool)) 0,
        search_overloads_1())
    .def("search", (bool(*)(Helium::Smarts&, const Molecule&, Helium::CountMapping&, const Helium::RingSet<Molecule>&, bool)) 0,
        search_overloads_2())
    .def("search", (bool(*)(Helium::Smarts&, const Molecule&, Helium::SingleMapping&, const Helium::RingSet<Molecule>&, bool)) 0,
        search_overloads_3())
    .def("search", (bool(*)(Helium::Smarts&, const Molecule&, Helium::MappingList&, const Helium::RingSet<Molecule>&, bool)) 0,
        search_overloads_4())
    ;

}
