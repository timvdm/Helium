/*
 * Copyright (c) 2014, Tim Vandermeersch
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the <organization> nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef HELIUM_SDF_H
#define HELIUM_SDF_H

#include <Helium/molecule.h>
#include <Helium/error.h>
#include <Helium/util.h>

#include <fstream>

namespace Helium {

  /**
   * @file .h
   * @brief .
   */
  class SDFInputFile
  {
    public:
      SDFInputFile()
      {
      }

      SDFInputFile(const std::string &filename)
      {
        open(filename);
      }

      bool open(const std::string &filename)
      {
        m_error = Error();

        m_ifs.open(filename.c_str());
        if (!m_ifs) {
          m_error = Error(make_string("Could not open ", filename));
          return false;
        }

        return true;
      }

      template<typename EditableMoleculeType>
      bool read(EditableMoleculeType &mol, std::string &title)
      {
        typedef typename molecule_traits<EditableMoleculeType>::atom_type atom_type;
        typedef typename molecule_traits<EditableMoleculeType>::bond_type bond_type;

        clear_molecule(mol);

        if (!m_ifs)
          return false;

        std::string buffer;

        //
        // read header
        //
        // title
        if (!std::getline(m_ifs, buffer)) {
          m_error = Error("Could not read title line");
          return false;
        }

        title = buffer;
        // ignore line 2-3
        std::getline(m_ifs, buffer);
        std::getline(m_ifs, buffer);

        //
        // counts line
        //
        // aaabbblllfffccc...
        // 012345678901234
        //           1
        //
        // aaa: number of atoms [generic]
        // bbb: number of bonds [generic]
        // lll: number of atom lists [query]
        // ccc: chiral flag: 0 = not chiral, 1 = chiral [generic]
        std::getline(m_ifs, buffer);
        // number of atoms
        int numAtoms;
        if (!parseNumber(buffer.substr(0, 3), numAtoms, "Could not parse number of atoms in counts line"))
          return false;
        // number of bonds
        int numBonds;
        if (!parseNumber(buffer.substr(3, 3), numBonds, "Could not parse number of bonds in counts line"))
          return false;
        bool chiralFlag;
        if (!parseNumber(buffer.substr(12, 3), chiralFlag, "Could not parse chiral flag in counts line"))
          return false;

        //
        // atom block
        //
        std::vector<int> valences;
        for (int i = 0; i < numAtoms; ++i) {
          // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
          // 012345678901234567890123456789012345678901234567890123456789012345678
          //           1         2         3         4         5         6
          //
          // x y z: atom coordinates [generic]
          // aaa: atom symbol [generic]
          // dd: mass difference [generic]
          // ccc: charge [generic]
          // sss: atom stereo parity [generic]
          // hhh: hydrogen count +1 [query]
          // bbb: stereo care box [query]
          // vvv: valence [generic]
          // HHH: H0 designator [ISIS/Desktop]
          // rrr: not used
          // iii: not used
          // mmm: atom-atom mapping number [reaction]
          // nnn: inversion/retention flag [reaction]
          // eee: exact change flag [query, reaction]
          std::getline(m_ifs, buffer);
          //std::cout << "ATOM LINE: " << buffer << std::endl;
          if (buffer.size() < 50) {
            m_error = Error(make_string("Line for atom ", i + 1, " is to short"));
            return false;
          }
          // x y z
          double x, y, z;
          if (!parseNumber(buffer.substr(0, 10), x, make_string("Could not parse x coordinate for atom ", i + 1)))
            return false;
          if (!parseNumber(buffer.substr(10, 10), y, make_string("Could not parse y coordinate for atom ", i + 1)))
            return false;
          if (!parseNumber(buffer.substr(20, 10), z, make_string("Could not parse z coordinate for atom ", i + 1)))
            return false;
          // atom symbol
          std::string symbol = buffer.substr(31, 3);
          while (symbol.front() == ' ')
            symbol = symbol.substr(1);
          while (symbol.back() == ' ')
            symbol.resize(symbol.size() - 1);
          int element = Element::element(symbol);
          // mass difference
          int mass;
          if (!parseNumber(buffer.substr(34, 2), mass, make_string("Could not parse mass difference for atom ", i + 1)))
            return false;
          mass += Element::averageMass(element);
          // charge
          int charge;
          if (!parseNumber(buffer.substr(36, 3), charge, make_string("Could not parse charge for atom ", i + 1)))
            return false;
          switch (charge) {
            case 0:
              charge = 0; // 0 -> 0
              break;
            case 1:
              charge = 3; // 0 -> 3
              break;
            case 2:
              charge = 2; // 0 -> 2
              break;
            case 3:
              charge = 1; // 0 -> 1
              break;
            case 4:
              charge = 0; // 0 -> doublet radical
              break;
            case 5:
              charge = -1; // 0 -> -1
              break;
            case 6:
              charge = -2; // 0 -> -2
              break;
            case 7:
              charge = -3; // 0 -> -3
              break;
            default:
              m_error = Error(make_string("Invalid charge for atom ", i + 1, ": ", charge));
              return false;
          }
          // atom stereo parity ignored when read
          // valence
          int valence;
          if (!parseNumber(buffer.substr(48, 3), valence, make_string("Could not parse valence for atom ", i + 1)))
            return false;
          valences.push_back(valence);

          // handle D & T
          if (symbol == "D")
            mass = 2;
          else if (symbol == "T")
            mass = 3;

          atom_type atom = add_atom(mol);
          set_element(mol, atom, element);
          set_mass(mol, atom, mass);
          set_charge(mol, atom, charge);
        }

        //
        // bond block
        //
        for (int i = 0; i < numBonds; ++i) {
          // 111222tttsssxxxrrrccc
          // 012345678901234567890
          //           1         2
          //
          // 111: source atom number [1,n]
          // 222: target atom number [1,n]
          // ttt: bond type [generic]
          // sss: bond stereo [generic]
          // xxx: not used
          // rrr: bond topology [query]
          // ccc: reaction center status [reaction]
          std::getline(m_ifs, buffer);
          if (buffer.size() < 11) {
            m_error = Error(make_string("Line for bond ", i + 1, " is to short"));
            return false;
          }
          // source atom
          int source;
          if (!parseNumber(buffer.substr(0, 3), source, make_string("Could not parse source atom number for bond ", i + 1)))
            return false;
          // target atom
          int target;
          if (!parseNumber(buffer.substr(3, 3), target, make_string("Could not parse target atom number for bond ", i + 1)))
            return false;
          // bond type
          int order;
          if (!parseNumber(buffer.substr(6, 3), order, make_string("Could not parse bond order for bond ", i + 1)))
            return false;
          // bond stereo
          int stereo;
          if (!parseNumber(buffer.substr(9, 3), stereo, make_string("Could not parse bond stereo for bond ", i + 1)))
            return false;

          if (source > num_atoms(mol)) {
            m_error = Error(make_string("Bond source index out of range for bond ", i + 1));
            return false;
          }
          if (target > num_atoms(mol)) {
            m_error = Error(make_string("Bond target index out of range for bond ", i + 1));
            return false;
          }

          bond_type bond = add_bond(mol, get_atom(mol, source - 1), get_atom(mol, target - 1));
          set_order(mol, bond, order);
        }

        // properties block
        bool firstCharge = true;
        do {
          std::getline(m_ifs, buffer);

          // check for end
          if (buffer.substr(0, 6) == "M  END")
            break;

          // charge
          if (buffer.substr(0, 6) == "M  CHG") {
            // reset all charges if this is the first charge property
            if (firstCharge) {
              FOREACH_ATOM (atom, mol)
                set_charge(mol, *atom, 0);
              firstCharge = false;
            }

            // M  CHGnn8 aaa vvv aaa vvv ...
            // 012345678901234567890123456
            //           1         2
            //
            // nn8: number of entries on line [1,8]
            // aaa: atom number
            // vvv: [-15,+15]
            int count;
            if (!parseNumber(buffer.substr(6, 3), count, "Could not parse number of entries for charge property"))
              return false;

            for (int i = 0; i < count; ++i) {
              int atom;
              if (!parseNumber(buffer.substr(10 + 8 * i, 3), atom, "Could not parse atom number for charge property"))
                return false;

              if (atom > num_atoms(mol)) {
                m_error = Error(make_string("Invalid atom number for charge property: ", atom));
                return false;
              }

              int charge;
              if (!parseNumber(buffer.substr(14 + 8 * i, 3), charge, "Could not parse charge for charge property"))
                return false;

              set_charge(mol, get_atom(mol, atom - 1), charge);
            }
          }

        } while (m_ifs);

        //
        // adjust hydrogens based on valences
        //
        FOREACH_ATOM (atom, mol) {
          if (!addHydrogens(get_element(mol, *atom)))
            continue;
          int valence = valences[get_index(mol, *atom)];
          if (valence == 15)
            continue;
          if (valence == 0)
            valence = Element::valence(get_element(mol, *atom),
                get_charge(mol, *atom), get_valence(mol, *atom));

          int numH = valence - get_valence(mol, *atom);
          if (numH > 0)
            set_hydrogens(mol, *atom, numH);
        }


        // try to read $$$$
        if (m_ifs && m_ifs.peek() == '$')
          std::getline(m_ifs, buffer);

        return true;
      }

      template<typename EditableMoleculeType>
      bool read(EditableMoleculeType &mol)
      {
        std::string title;
        return read(mol, title);
      }

      bool isGood()
      {
        // make for eof...
        m_ifs.peek();
        return m_ifs.good();
      }

      /**
       * @brief Get the last error.
       *
       * @return The last error.
       */
      const Error& error() const
      {
        return m_error;
      }


    private:
      template<typename Number>
      bool parseNumber(const std::string &str, Number &result, const std::string &error)
      {
        std::stringstream ss(str);
        if (!(ss >> result)) {
          m_error = Error(error);
          return false;
        }
        return true;
      }

      bool addHydrogens(int element) const
      {
        switch (element) {
          case 1:
          case 3:
          case 4:
          case 5:
          case 6:
          case 7:
          case 8:
          case 9:
          case 11:
          case 12:
          case 13:
          case 14:
          case 15:
          case 16:
          case 17:
          case 19:
          case 20:
          case 31:
          case 32:
          case 33:
          case 34:
          case 35:
          case 37:
          case 38:
          case 49:
          case 50:
          case 51:
          case 52:
          case 53:
          case 55:
          case 56:
          case 81:
          case 82:
          case 83:
          case 84:
          case 85:
          case 87:
          case 88:
            return true;
          default:
            return false;
        }
      }

      Error m_error; //!< The last error.
      std::ifstream m_ifs; //!< The input stream.

  };

}

#endif
