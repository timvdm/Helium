#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/dfs.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

struct DFSVisitor
{
  typedef Helium::molecule_traits<Molecule>::atom_type atom_type;
  typedef Helium::molecule_traits<Molecule>::bond_type bond_type;

  virtual void initialize(const Molecule &mol) {}
  virtual void component(int i) {}
  virtual void atom(const Molecule &mol, atom_type prev, atom_type atom) {}
  virtual void bond(const Molecule &mol, atom_type prev, bond_type bond) {}
  virtual void backtrack(const Molecule &mol, atom_type atom) {}
  virtual void back_bond(const Molecule &mol, bond_type bond) {}
  virtual bool stop() const { return false; }
};

struct DFSVisitorWrapper : DFSVisitor, wrapper<DFSVisitor>
{
  void initialize(const Molecule &mol)
  {
    if (override f = this->get_override("initialize")) {
      f(boost::ref(mol));
      return;
    }
  }

  void component(int i) 
  {
    if (override f = this->get_override("component")) {
      f(i);
      return;
    }
  }

  void atom(const Molecule &mol, atom_type prev, atom_type atom)
  {
    if (override f = this->get_override("atom")) {
      f(boost::ref(mol), boost::ref(prev), boost::ref(atom));
      return;
    }
  }

  void bond(const Molecule &mol, atom_type prev, bond_type bond)
  {
    if (override f = this->get_override("bond")) {
      f(boost::ref(mol), boost::ref(prev), boost::ref(bond));
      return;
    }
  }

  void backtrack(const Molecule &mol, atom_type atom)
  {
    if (override f = this->get_override("backtrack")) {
      f(boost::ref(mol), boost::ref(atom));
      return;
    }
  }

  void back_bond(const Molecule &mol, bond_type bond)
  {
    if (override f = this->get_override("back_bond")) {
      f(boost::ref(mol), boost::ref(bond));
      return;
    }
  }

  bool stop() const
  {
    if (override f = this->get_override("stop"))
      return f();
    return false;
  }
};

void depth_first_search_1(const Molecule &mol, DFSVisitorWrapper &visitor)
{
  Helium::depth_first_search(mol, visitor);
}

void depth_first_search_2(const Molecule &mol, DFSVisitorWrapper &visitor,
      const std::vector<bool> &atomMask)
{
  Helium::depth_first_search_mask(mol, visitor, atomMask);
}

void depth_first_search_3(const Molecule &mol, DFSVisitorWrapper &visitor,
      const std::vector<bool> &atomMask, const std::vector<bool> &bondMask)
{
  Helium::depth_first_search_mask(mol, visitor, atomMask, bondMask);
}

void depth_first_search_4(const Molecule &mol, Molecule::atom_type atom,
      DFSVisitorWrapper &visitor)
{
  Helium::depth_first_search(mol, atom, visitor);
}

/*
  void depth_first_search_mask(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom, DFSVisitorType &visitor,
      const std::vector<bool> &atomMask)
  void depth_first_search_mask(const MoleculeType &mol,
      typename molecule_traits<MoleculeType>::atom_type atom, DFSVisitorType &visitor,
      const std::vector<bool> &atomMask, const std::vector<bool> &bondMask)
  void depth_first_search(const MoleculeType &mol, const std::vector<Index> &order, DFSVisitorType &visitor)
  void exhaustive_depth_first_search(const MoleculeType &mol, AtomType atom, DFSVisitorType &visitor)
  */



void export_dfs()
{

  class_<DFSVisitorWrapper, boost::noncopyable>("DFSVisitor")
    .def("initialize", &DFSVisitor::initialize)
    .def("component", &DFSVisitor::component)
    .def("atom", &DFSVisitor::atom)
    .def("bond", &DFSVisitor::bond)
    .def("backtrack", &DFSVisitor::backtrack)
    .def("back_bond", &DFSVisitor::back_bond)
    .def("stop", &DFSVisitor::stop)
    ;
 

  def("depth_first_search", &depth_first_search_1);
  def("depth_first_search_mask", &depth_first_search_2);
  def("depth_first_search_mask", &depth_first_search_3);
  def("depth_first_search", &depth_first_search_4);





}
