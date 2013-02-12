#ifndef HELIUM_FILEIO_H
#define HELIUM_FILEIO_H

#include "hemol.h"
#include "util.h"

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>

#include <cassert>

namespace Helium {

  enum AtomProperty {
    AromaticCyclic = 1,
    Mass = 2,
    Hydrogens = 4,
    Charge = 8
  };

  inline void write_molecule(std::ostream &os, OpenBabel::OBMol *mol)
  {
    std::vector<unsigned short> indices(mol->NumAtoms());
    unsigned short numAtoms = 0;
    FOR_ATOMS_OF_MOL (a, mol) {
      indices[a->GetIndex()] = numAtoms;
      if (!a->IsHydrogen())
        ++numAtoms;
    }
    unsigned short numBonds = 0;
    FOR_BONDS_OF_MOL (b, mol)
      if (!b->GetBeginAtom()->IsHydrogen() && !b->GetEndAtom()->IsHydrogen())
        ++numBonds;

    write16<unsigned short>(os, numAtoms);
    write16<unsigned short>(os, numBonds);

    FOR_ATOMS_OF_MOL (atom, mol) {
      if (atom->IsHydrogen())
        continue;

      // always write these properies
      write8<unsigned char>(os, atom->GetAtomicNum());

      unsigned char aromaticCyclic = 0;
      if (atom->IsAromatic())
        aromaticCyclic |= 1;
      if (atom->IsInRing())
        aromaticCyclic |= 2;

      unsigned char flags = 0;
      if (aromaticCyclic)
        flags |= AromaticCyclic;
      if (atom->GetIsotope())
        flags |= Mass;
      if (atom->ExplicitHydrogenCount() + atom->ImplicitHydrogenCount())
        flags |= Hydrogens;
      if (atom->GetFormalCharge())
        flags |= Charge;

      // write flags
      write8<unsigned char>(os, flags);

      if (flags & AromaticCyclic)
        write8<unsigned char>(os, aromaticCyclic);
      if (flags & Mass)
        write8<unsigned char>(os, atom->GetIsotope());
      if (flags & Hydrogens)
        write8<unsigned char>(os, atom->ExplicitHydrogenCount() + atom->ImplicitHydrogenCount());
      if (flags & Charge)
        write8<signed char>(os, atom->GetFormalCharge());
    }

    FOR_BONDS_OF_MOL (bond, mol) {
      if (bond->GetBeginAtom()->IsHydrogen() || bond->GetEndAtom()->IsHydrogen())
        continue;

      write16<unsigned short>(os, indices[bond->GetBeginAtom()->GetIndex()]);
      write16<unsigned short>(os, indices[bond->GetEndAtom()->GetIndex()]);
      unsigned char props = bond->GetBondOrder();
      if (bond->IsAromatic())
        props |= 128;
      if (bond->IsInRing())
        props |= 64;
      write8<unsigned char>(os, props);
    }
  }

  template<typename MoleculeType>
  bool read_molecule(std::istream &is, MoleculeType &mol)
  {
    mol.clear();

    if (!is)
      return false;

    unsigned short numAtoms, numBonds;

    read16(is, numAtoms);
    read16(is, numBonds);

    if (!is)
      return false;

    unsigned char flags, aromaticCyclic, aromatic = 0, cyclic = 0;
    unsigned char element, mass, hydrogens;
    signed char charge;
    for (int i = 0; i < numAtoms; ++i) {
      // properties that arways there
      read8(is, element);

      // read flags
      read8(is, flags);

      if (flags & AromaticCyclic) {
        read8(is, aromaticCyclic);
        aromatic = aromaticCyclic & 1;
        cyclic = aromaticCyclic & 2;
      } else {
        aromatic = 0;
        cyclic = 0;
      }

      if (flags & Mass)
        read8(is, mass);
      else
        mass = 0;
      if (flags & Hydrogens)
        read8(is, hydrogens);
      else
        hydrogens = 0;
      if (flags & Charge)
        read8(is, charge);
      else
        charge = 0;


      HeAtom atom = mol.addAtom();
      atom.setAromatic(aromatic);
      atom.setCyclic(cyclic);
      atom.setElement(element);
      atom.setMass(mass);
      atom.setHydrogens(hydrogens);
      atom.setCharge(charge);
    }

    unsigned short source, target;
    unsigned char props;
    for (int i = 0; i < numBonds; ++i) {
      read16(is, source);
      read16(is, target);
      read8(is, props);

      HeAtom s = mol.atom(source);
      HeAtom t = mol.atom(target);

      HeBond bond = mol.addBond(s, t);
      bond.setAromatic(props & 128);
      bond.setCyclic(props & 64);
      bond.setOrder(props & 63);
    }

    return (bool)is;
  }

  inline std::string normalize_smiles(const std::string &smiles)
  {
    OpenBabel::OBMol obmol;
    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("smi", "smi");
    conv.ReadString(&obmol, smiles);
    return conv.WriteString(&obmol, true);
  }

  class MoleculeFile
  {
    public:
      MoleculeFile(const std::string &filename) : m_ifs(filename.c_str()), m_current(-1)
      {
        if (m_ifs)
          read32(m_ifs, m_numMolecules);
      }

      unsigned int numMolecules() const
      {
        return m_numMolecules;
      }

      unsigned int current() const
      {
        return m_current;
      }

      template<typename MoleculeType>
      bool read_molecule(MoleculeType &mol)
      {
        if (!m_ifs)
          return false;
        ++m_current;
        if (m_current == m_numMolecules)
          return false;
        return Helium::read_molecule(m_ifs, mol);
      }

    private:
      std::ifstream m_ifs;
      unsigned int m_numMolecules;
      unsigned int m_current;
  };




}

#endif
