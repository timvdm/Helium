#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/bfs.h"
#include "../common.h"

using Helium::Chemist::Molecule;
using namespace boost::python;

struct BFSVisitor
{
  typedef Helium::molecule_traits<Molecule>::atom_type atom_type;
  typedef Helium::molecule_traits<Molecule>::bond_type bond_type;

  virtual void initialize(const Molecule &mol) {}
  virtual void component(int i) {}
  virtual void depth(int d) {}
  virtual void atom(const Molecule &mol, atom_type prev, atom_type atom) {}
  virtual void bond(const Molecule &mol, atom_type prev, bond_type bond) {}
  virtual void back_bond(const Molecule &mol, bond_type bond) {}
  virtual bool stop() const { return false; }
};

struct BFSVisitorWrapper : BFSVisitor, wrapper<BFSVisitor>
{
  void initialize(const Molecule &mol)
  {
    if (override f = this->get_override("initialize"))
      f(boost::ref(mol));
  }

  void component(int i)
  {
    if (override f = this->get_override("component"))
      f(i);
  }

  void depth(int d)
  {
    if (override f = this->get_override("depth"))
      f(d);
  }

  void atom(const Molecule &mol, atom_type prev, atom_type atom)
  {
    if (override f = this->get_override("atom"))
      f(boost::ref(mol), boost::ref(prev), boost::ref(atom));
  }

  void bond(const Molecule &mol, atom_type prev, bond_type bond)
  {
    if (override f = this->get_override("bond"))
      f(boost::ref(mol), boost::ref(prev), boost::ref(bond));
  }

  void back_bond(const Molecule &mol, bond_type bond)
  {
    if (override f = this->get_override("back_bond"))
      f(boost::ref(mol), boost::ref(bond));
  }

  bool stop() const
  {
    if (override f = this->get_override("stop"))
      return f();
    return false;
  }
};

struct BFSDebugVisitorWrapper : Helium::BFSDebugVisitor<Molecule>
{
  BFSDebugVisitorWrapper() : BFSDebugVisitor(ss)
  {
  }

  std::string output() const
  {
    return ss.str();
  }

  std::stringstream ss;
};

DATA_MEMBER_TO_FUNCTION(Helium::BFSAtomOrderVisitor<Molecule>, std::vector<unsigned int>, atoms)
DATA_MEMBER_TO_FUNCTION(Helium::BFSBondOrderVisitor<Molecule>, std::vector<unsigned int>, bonds)
DATA_MEMBER_TO_FUNCTION(Helium::BFSClosureRecorderVisitor<Molecule>, std::vector<unsigned int>, back_bonds)

template<typename BFSVisitorType>
void breadth_first_search_1(const Molecule &mol, BFSVisitorType &visitor)
{
  Helium::breadth_first_search(mol, visitor);
}

template<typename BFSVisitorType>
void breadth_first_search_2(const Molecule &mol, BFSVisitorType &visitor,
      const list &atomMask)
{
  if (len(atomMask) != mol.numAtoms())
    throw std::runtime_error("Invalid atom mask, size should be equal to the number of atoms");
  Helium::breadth_first_search_mask(mol, visitor, vector_from_list<bool>(atomMask));
}

template<typename BFSVisitorType>
void breadth_first_search_3(const Molecule &mol, BFSVisitorType &visitor,
      const list &atomMask, const list &bondMask)
{
  if (len(atomMask) != mol.numAtoms())
    throw std::runtime_error("Invalid atom mask, size should be equal to the number of atoms");
  if (len(bondMask) != mol.numBonds())
    throw std::runtime_error("Invalid bond mask, size should be equal to the number of bonds");
  Helium::breadth_first_search_mask(mol, visitor, vector_from_list<bool>(atomMask),
      vector_from_list<bool>(bondMask));
}

template<typename BFSVisitorType>
void breadth_first_search_4(const Molecule &mol, Molecule::atom_type atom,
      BFSVisitorType &visitor)
{
  Helium::breadth_first_search(mol, atom, visitor);
}

template<typename BFSVisitorType>
void breadth_first_search_5(const Molecule &mol, Molecule::atom_type atom,
      BFSVisitorType &visitor, const list &atomMask)
{
  if (len(atomMask) != mol.numAtoms())
    throw std::runtime_error("Invalid atom mask, size should be equal to the number of atoms");
  Helium::breadth_first_search_mask(mol, atom, visitor,
      vector_from_list<bool>(atomMask));
}

template<typename BFSVisitorType>
void breadth_first_search_6(const Molecule &mol, Molecule::atom_type atom,
    BFSVisitorType &visitor, const list &atomMask, const list &bondMask)
{
  if (len(atomMask) != mol.numAtoms())
    throw std::runtime_error("Invalid atom mask, size should be equal to the number of atoms");
  if (len(bondMask) != mol.numBonds())
    throw std::runtime_error("Invalid bond mask, size should be equal to the number of bonds");
  Helium::breadth_first_search_mask(mol, atom, visitor,
      vector_from_list<bool>(atomMask), vector_from_list<bool>(bondMask));
}

void export_bfs()
{

  class_<BFSVisitorWrapper, boost::noncopyable>("BFSVisitor")
    .def("initialize", &BFSVisitor::initialize)
    .def("component", &BFSVisitor::component)
    .def("depth", &BFSVisitor::depth)
    .def("atom", &BFSVisitor::atom)
    .def("bond", &BFSVisitor::bond)
    .def("back_bond", &BFSVisitor::back_bond)
    .def("stop", &BFSVisitor::stop)
    ;

  class_<Helium::BFSAtomOrderVisitor<Molecule>, boost::noncopyable>("BFSAtomOrderVisitor")
    .add_property("atoms", make_function(&atoms, return_value_policy<copy_const_reference>()))
    ;

  class_<Helium::BFSBondOrderVisitor<Molecule>, boost::noncopyable>("BFSBondOrderVisitor")
    .add_property("bonds", make_function(&bonds, return_value_policy<copy_const_reference>()))
    ;

  class_<Helium::BFSClosureRecorderVisitor<Molecule>, boost::noncopyable>("BFSClosureRecorderVisitor")
    .add_property("back_bonds", make_function(&back_bonds, return_value_policy<copy_const_reference>()))
    ;

  class_<BFSDebugVisitorWrapper, boost::noncopyable>("BFSDebugVisitor")
    .add_property("output", &BFSDebugVisitorWrapper::output)
    ;

  def("breadth_first_search", &breadth_first_search_1<BFSVisitorWrapper>);
  def("breadth_first_search", &breadth_first_search_1<Helium::BFSAtomOrderVisitor<Molecule> >);
  def("breadth_first_search", &breadth_first_search_1<Helium::BFSBondOrderVisitor<Molecule> >);
  def("breadth_first_search", &breadth_first_search_1<Helium::BFSClosureRecorderVisitor<Molecule> >);
  def("breadth_first_search", &breadth_first_search_1<BFSDebugVisitorWrapper>);

  def("breadth_first_search_mask", &breadth_first_search_2<BFSVisitorWrapper>);
  def("breadth_first_search_mask", &breadth_first_search_2<Helium::BFSAtomOrderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_2<Helium::BFSBondOrderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_2<Helium::BFSClosureRecorderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_2<BFSDebugVisitorWrapper>);

  def("breadth_first_search_mask", &breadth_first_search_3<BFSVisitorWrapper>);
  def("breadth_first_search_mask", &breadth_first_search_3<Helium::BFSAtomOrderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_3<Helium::BFSBondOrderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_3<Helium::BFSClosureRecorderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_3<BFSDebugVisitorWrapper>);

  def("breadth_first_search", &breadth_first_search_4<BFSVisitorWrapper>);
  def("breadth_first_search", &breadth_first_search_4<Helium::BFSAtomOrderVisitor<Molecule> >);
  def("breadth_first_search", &breadth_first_search_4<Helium::BFSBondOrderVisitor<Molecule> >);
  def("breadth_first_search", &breadth_first_search_4<Helium::BFSClosureRecorderVisitor<Molecule> >);
  def("breadth_first_search", &breadth_first_search_4<BFSDebugVisitorWrapper>);

  def("breadth_first_search_mask", &breadth_first_search_5<BFSVisitorWrapper>);
  def("breadth_first_search_mask", &breadth_first_search_5<Helium::BFSAtomOrderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_5<Helium::BFSBondOrderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_5<Helium::BFSClosureRecorderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_5<BFSDebugVisitorWrapper>);

  def("breadth_first_search_mask", &breadth_first_search_6<BFSVisitorWrapper>);
  def("breadth_first_search_mask", &breadth_first_search_6<Helium::BFSAtomOrderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_6<Helium::BFSBondOrderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_6<Helium::BFSClosureRecorderVisitor<Molecule> >);
  def("breadth_first_search_mask", &breadth_first_search_6<BFSDebugVisitorWrapper>);

}
