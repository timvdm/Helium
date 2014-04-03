#include <boost/python.hpp>

#include "../../src/chemist/molecule.h"
#include "../../src/algorithms/dfs.h"
#include "../common.h"

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
    if (override f = this->get_override("initialize"))
      f(boost::ref(mol));
  }

  void component(int i)
  {
    if (override f = this->get_override("component"))
      f(i);
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

  void backtrack(const Molecule &mol, atom_type atom)
  {
    if (override f = this->get_override("backtrack"))
      f(boost::ref(mol), boost::ref(atom));
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

struct DFSDebugVisitorWrapper : Helium::DFSDebugVisitor<Molecule>
{
  DFSDebugVisitorWrapper() : DFSDebugVisitor(ss)
  {
  }

  std::string output() const
  {
    return ss.str();
  }

  std::stringstream ss;
};

DATA_MEMBER_TO_FUNCTION(Helium::DFSAtomOrderVisitor<Molecule>, std::vector<unsigned int>, atoms);
DATA_MEMBER_TO_FUNCTION(Helium::DFSBondOrderVisitor<Molecule>, std::vector<unsigned int>, bonds);
DATA_MEMBER_TO_FUNCTION(Helium::DFSClosureRecorderVisitor<Molecule>, std::vector<unsigned int>, back_bonds);

template<typename DFSVisitorType>
void depth_first_search_1(const Molecule &mol, DFSVisitorType &visitor)
{
  Helium::depth_first_search(mol, visitor);
}

template<typename DFSVisitorType>
void depth_first_search_2(const Molecule &mol, DFSVisitorType &visitor,
      const list &atomMask)
{
  Helium::depth_first_search_mask(mol, visitor, vector_from_list<bool>(atomMask));
}

template<typename DFSVisitorType>
void depth_first_search_3(const Molecule &mol, DFSVisitorType &visitor,
      const list &atomMask, const list &bondMask)
{
  Helium::depth_first_search_mask(mol, visitor, vector_from_list<bool>(atomMask),
      vector_from_list<bool>(bondMask));
}

template<typename DFSVisitorType>
void depth_first_search_4(const Molecule &mol, Molecule::atom_type atom,
      DFSVisitorType &visitor)
{
  Helium::depth_first_search(mol, atom, visitor);
}

template<typename DFSVisitorType>
void depth_first_search_5(const Molecule &mol, Molecule::atom_type atom,
      DFSVisitorType &visitor, const list &atomMask)
{
  Helium::depth_first_search_mask(mol, atom, visitor,
      vector_from_list<bool>(atomMask));
}

template<typename DFSVisitorType>
void depth_first_search_6(const Molecule &mol, Molecule::atom_type atom,
    DFSVisitorType &visitor, const list &atomMask, const list &bondMask)
{
  Helium::depth_first_search_mask(mol, atom, visitor,
      vector_from_list<bool>(atomMask), vector_from_list<bool>(bondMask));
}

template<typename DFSVisitorType>
void depth_first_search_7(const Molecule &mol, const list &order,
    DFSVisitorType &visitor)
{
  Helium::ordered_depth_first_search(mol, vector_from_list<Helium::Index>(order),
      visitor);
}

template<typename DFSVisitorType>
void depth_first_search_8(const Molecule &mol, Molecule::atom_type atom,
    DFSVisitorType &visitor)
{
  return exhaustive_depth_first_search(mol, atom, visitor);
}

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

  class_<Helium::DFSAtomOrderVisitor<Molecule>, boost::noncopyable>("DFSAtomOrderVisitor")
    .add_property("atoms", make_function(&atoms, return_value_policy<copy_const_reference>()))
    ;

  class_<Helium::DFSBondOrderVisitor<Molecule>, boost::noncopyable>("DFSBondOrderVisitor")
    .add_property("bonds", make_function(&bonds, return_value_policy<copy_const_reference>()))
    ;

  class_<Helium::DFSClosureRecorderVisitor<Molecule>, boost::noncopyable>("DFSClosureRecorderVisitor")
    .add_property("back_bonds", make_function(&back_bonds, return_value_policy<copy_const_reference>()))
    ;

  class_<DFSDebugVisitorWrapper, boost::noncopyable>("DFSDebugVisitor")
    .add_property("output", &DFSDebugVisitorWrapper::output)
    ;

  def("depth_first_search", &depth_first_search_1<DFSVisitorWrapper>);
  def("depth_first_search", &depth_first_search_1<Helium::DFSAtomOrderVisitor<Molecule> >);
  def("depth_first_search", &depth_first_search_1<Helium::DFSBondOrderVisitor<Molecule> >);
  def("depth_first_search", &depth_first_search_1<Helium::DFSClosureRecorderVisitor<Molecule> >);
  def("depth_first_search", &depth_first_search_1<DFSDebugVisitorWrapper>);

  def("depth_first_search_mask", &depth_first_search_2<DFSVisitorWrapper>);
  def("depth_first_search_mask", &depth_first_search_2<Helium::DFSAtomOrderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_2<Helium::DFSBondOrderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_2<Helium::DFSClosureRecorderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_2<DFSDebugVisitorWrapper>);

  def("depth_first_search_mask", &depth_first_search_3<DFSVisitorWrapper>);
  def("depth_first_search_mask", &depth_first_search_3<Helium::DFSAtomOrderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_3<Helium::DFSBondOrderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_3<Helium::DFSClosureRecorderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_3<DFSDebugVisitorWrapper>);

  def("depth_first_search", &depth_first_search_4<DFSVisitorWrapper>);
  def("depth_first_search", &depth_first_search_4<Helium::DFSAtomOrderVisitor<Molecule> >);
  def("depth_first_search", &depth_first_search_4<Helium::DFSBondOrderVisitor<Molecule> >);
  def("depth_first_search", &depth_first_search_4<Helium::DFSClosureRecorderVisitor<Molecule> >);
  def("depth_first_search", &depth_first_search_4<DFSDebugVisitorWrapper>);

  def("depth_first_search_mask", &depth_first_search_5<DFSVisitorWrapper>);
  def("depth_first_search_mask", &depth_first_search_5<Helium::DFSAtomOrderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_5<Helium::DFSBondOrderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_5<Helium::DFSClosureRecorderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_5<DFSDebugVisitorWrapper>);

  def("depth_first_search_mask", &depth_first_search_6<DFSVisitorWrapper>);
  def("depth_first_search_mask", &depth_first_search_6<Helium::DFSAtomOrderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_6<Helium::DFSBondOrderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_6<Helium::DFSClosureRecorderVisitor<Molecule> >);
  def("depth_first_search_mask", &depth_first_search_6<DFSDebugVisitorWrapper>);

  def("ordered_depth_first_search", &depth_first_search_7<DFSVisitorWrapper>);
  def("ordered_depth_first_search", &depth_first_search_7<Helium::DFSAtomOrderVisitor<Molecule> >);
  def("ordered_depth_first_search", &depth_first_search_7<Helium::DFSBondOrderVisitor<Molecule> >);
  def("ordered_depth_first_search", &depth_first_search_7<Helium::DFSClosureRecorderVisitor<Molecule> >);
  def("ordered_depth_first_search", &depth_first_search_7<DFSDebugVisitorWrapper>);

  def("exhaustive_depth_first_search", &depth_first_search_8<DFSVisitorWrapper>);
  def("exhaustive_depth_first_search", &depth_first_search_8<Helium::DFSAtomOrderVisitor<Molecule> >);
  def("exhaustive_depth_first_search", &depth_first_search_8<Helium::DFSBondOrderVisitor<Molecule> >);
  def("exhaustive_depth_first_search", &depth_first_search_8<Helium::DFSClosureRecorderVisitor<Molecule> >);
  def("exhaustive_depth_first_search", &depth_first_search_8<DFSDebugVisitorWrapper>);

}
