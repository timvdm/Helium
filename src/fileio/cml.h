#include <Helium/molecule.h>
#include <Helium/element.h>

#include <Eigen/Core>

namespace Helium {

  template<typename MoleculeType>
  void write_cml(std::ostream &os, const MoleculeType &mol, const std::vector<Eigen::Vector2d> &coords)
  {
    os << "<?xml version=\"1.0\"?>" << std::endl;
    os << "<molecule xmlns=\"http://www.xml-cml.org/schema\">" << std::endl;

    os << " <atomArray>" << std::endl;
    for (auto &atom : get_atoms(mol)) {
      os << "  <atom ";
      os << "id=\"a" << get_index(mol, atom) + 1 << "\" ";
      os << "elementType=\"" << Element::symbol(get_element(mol, atom)) << "\" ";
      os << "x2=\"" << coords[get_index(mol, atom)].x() << "\" ";
      os << "y2=\"" << coords[get_index(mol, atom)].y() << "\"";
      os << "/>" << std::endl;
    }
    os << " </atomArray>" << std::endl;

    os << " <bondArray>" << std::endl;
    for (auto &bond : get_bonds(mol)) {
      os << "  <bond ";
      os << "atomRefs2=\"a" << get_index(mol, get_source(mol, *bond)) + 1 << " a" << get_index(mol, get_target(mol, *bond)) + 1 << "\" ";
      os << "order=\"" << get_order(mol, *bond) << "\"";
      os << "/>" << std::endl;
    }
    os << " </bondArray>" << std::endl;

    os << "</molecule>" << std::endl;
  }

}
