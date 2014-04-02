#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/smarts.h"
#include "common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

bool findMapping_1(Helium::Smarts &smarts, const Molecule &mol, const Helium::RingSet<Molecule> &rings,
    Helium::NoMapping &mapping, bool uniqueComponents = true)
{  return smarts.findMapping(mol, rings, mapping, uniqueComponents); }

bool findMapping_2(Helium::Smarts &smarts, const Molecule &mol, const Helium::RingSet<Molecule> &rings,
    Helium::CountMapping &mapping, bool uniqueComponents = true)
{  return smarts.findMapping(mol, rings, mapping, uniqueComponents); }

bool findMapping_3(Helium::Smarts &smarts, const Molecule &mol, const Helium::RingSet<Molecule> &rings,
    Helium::SingleMapping &mapping, bool uniqueComponents = true)
{  return smarts.findMapping(mol, rings, mapping, uniqueComponents); }

bool findMapping_4(Helium::Smarts &smarts, const Molecule &mol, const Helium::RingSet<Molecule> &rings,
    Helium::MappingList &mapping, bool uniqueComponents = true)
{  return smarts.findMapping(mol, rings, mapping, uniqueComponents); }

BOOST_PYTHON_FUNCTION_OVERLOADS(findMapping_overloads_1, findMapping_1, 4, 5);
BOOST_PYTHON_FUNCTION_OVERLOADS(findMapping_overloads_2, findMapping_2, 4, 5);
BOOST_PYTHON_FUNCTION_OVERLOADS(findMapping_overloads_3, findMapping_3, 4, 5);
BOOST_PYTHON_FUNCTION_OVERLOADS(findMapping_overloads_4, findMapping_4, 4, 5);

bool findMapping_5(Helium::Smarts &smarts, const Molecule &mol,
    Helium::NoMapping &mapping, bool uniqueComponents = true)
{  return smarts.findMapping(mol, mapping, uniqueComponents); }

bool findMapping_6(Helium::Smarts &smarts, const Molecule &mol,
    Helium::CountMapping &mapping, bool uniqueComponents = true)
{  return smarts.findMapping(mol, mapping, uniqueComponents); }

bool findMapping_7(Helium::Smarts &smarts, const Molecule &mol,
    Helium::SingleMapping &mapping, bool uniqueComponents = true)
{  return smarts.findMapping(mol, mapping, uniqueComponents); }

bool findMapping_8(Helium::Smarts &smarts, const Molecule &mol,
    Helium::MappingList &mapping, bool uniqueComponents = true)
{  return smarts.findMapping(mol, mapping, uniqueComponents); }

BOOST_PYTHON_FUNCTION_OVERLOADS(findMapping_overloads_5, findMapping_5, 3, 4);
BOOST_PYTHON_FUNCTION_OVERLOADS(findMapping_overloads_6, findMapping_6, 3, 4);
BOOST_PYTHON_FUNCTION_OVERLOADS(findMapping_overloads_7, findMapping_7, 3, 4);
BOOST_PYTHON_FUNCTION_OVERLOADS(findMapping_overloads_8, findMapping_8, 3, 4);

bool find_1(Helium::Smarts &smarts, const Molecule &mol, const Helium::RingSet<Molecule> &rings,
    bool uniqueComponents = true)
{  return smarts.find(mol, rings, uniqueComponents); }

bool find_2(Helium::Smarts &smarts, const Molecule &mol, const Helium::RingSet<Molecule> &rings,
    bool uniqueComponents = true)
{  return smarts.find(mol, rings, uniqueComponents); }

bool find_3(Helium::Smarts &smarts, const Molecule &mol, const Helium::RingSet<Molecule> &rings,
    bool uniqueComponents = true)
{  return smarts.find(mol, rings, uniqueComponents); }

bool find_4(Helium::Smarts &smarts, const Molecule &mol, const Helium::RingSet<Molecule> &rings,
    bool uniqueComponents = true)
{  return smarts.find(mol, rings, uniqueComponents); }

BOOST_PYTHON_FUNCTION_OVERLOADS(find_overloads_1, find_1, 3, 4);
BOOST_PYTHON_FUNCTION_OVERLOADS(find_overloads_2, find_2, 3, 4);
BOOST_PYTHON_FUNCTION_OVERLOADS(find_overloads_3, find_3, 3, 4);
BOOST_PYTHON_FUNCTION_OVERLOADS(find_overloads_4, find_4, 3, 4);

bool find_5(Helium::Smarts &smarts, const Molecule &mol,
    bool uniqueComponents = true)
{  return smarts.find(mol, uniqueComponents); }

bool find_6(Helium::Smarts &smarts, const Molecule &mol,
    bool uniqueComponents = true)
{  return smarts.find(mol, uniqueComponents); }

bool find_7(Helium::Smarts &smarts, const Molecule &mol,
    bool uniqueComponents = true)
{  return smarts.find(mol, uniqueComponents); }

bool find_8(Helium::Smarts &smarts, const Molecule &mol,
    bool uniqueComponents = true)
{  return smarts.find(mol, uniqueComponents); }

BOOST_PYTHON_FUNCTION_OVERLOADS(find_overloads_5, find_5, 2, 3);
BOOST_PYTHON_FUNCTION_OVERLOADS(find_overloads_6, find_6, 2, 3);
BOOST_PYTHON_FUNCTION_OVERLOADS(find_overloads_7, find_7, 2, 3);
BOOST_PYTHON_FUNCTION_OVERLOADS(find_overloads_8, find_8, 2, 3);

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
    .def("findMapping", (bool(*)(Helium::Smarts&, const Molecule&, const Helium::RingSet<Molecule>&, Helium::NoMapping&, bool)) 0,
        findMapping_overloads_1())
    .def("findMapping", (bool(*)(Helium::Smarts&, const Molecule&, const Helium::RingSet<Molecule>&, Helium::CountMapping&, bool)) 0,
        findMapping_overloads_2())
    .def("findMapping", (bool(*)(Helium::Smarts&, const Molecule&, const Helium::RingSet<Molecule>&, Helium::SingleMapping&, bool)) 0,
        findMapping_overloads_3())
    .def("findMapping", (bool(*)(Helium::Smarts&, const Molecule&, const Helium::RingSet<Molecule>&, Helium::MappingList&, bool)) 0,
        findMapping_overloads_4())
    .def("findMapping", (bool(*)(Helium::Smarts&, const Molecule&, Helium::NoMapping&, bool)) 0,
        findMapping_overloads_5())
    .def("findMapping", (bool(*)(Helium::Smarts&, const Molecule&, Helium::CountMapping&, bool)) 0,
        findMapping_overloads_6())
    .def("findMapping", (bool(*)(Helium::Smarts&, const Molecule&, Helium::SingleMapping&, bool)) 0,
        findMapping_overloads_7())
    .def("findMapping", (bool(*)(Helium::Smarts&, const Molecule&, Helium::MappingList&, bool)) 0,
        findMapping_overloads_8())
    .def("find", (bool(*)(Helium::Smarts&, const Molecule&, const Helium::RingSet<Molecule>&, bool)) 0,
        find_overloads_1())
    .def("find", (bool(*)(Helium::Smarts&, const Molecule&, const Helium::RingSet<Molecule>&, bool)) 0,
        find_overloads_2())
    .def("find", (bool(*)(Helium::Smarts&, const Molecule&, const Helium::RingSet<Molecule>&, bool)) 0,
        find_overloads_3())
    .def("find", (bool(*)(Helium::Smarts&, const Molecule&, const Helium::RingSet<Molecule>&, bool)) 0,
        find_overloads_4())
    .def("find", (bool(*)(Helium::Smarts&, const Molecule&, bool)) 0,
        find_overloads_5())
    .def("find", (bool(*)(Helium::Smarts&, const Molecule&, bool)) 0,
        find_overloads_6())
    .def("find", (bool(*)(Helium::Smarts&, const Molecule&, bool)) 0,
        find_overloads_7())
    .def("find", (bool(*)(Helium::Smarts&, const Molecule&, bool)) 0,
        find_overloads_8())
    ;

}
