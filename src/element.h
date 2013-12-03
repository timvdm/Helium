#ifndef HELIUM_ELEMENT_H
#define HELIUM_ELEMENT_H

#include <string>

namespace Helium {

  /**
   * @brief Element properties.
   */
  class Element
  {
    public:
      /**
       * @brief Get the element symbol for an atom.
       *
       * @param element The atom element number.
       *
       * @return The element symbol.
       */
      static std::string symbol(int element)
      {
        switch (element) {
          case 1:
            return "H";
          case 2:
            return "He";
          case 3:
            return "Li";
          case 4:
            return "Be";
          case 5:
            return "B";
          case 6:
            return "C";
          case 7:
            return "N";
          case 8:
            return "O";
          case 9:
            return "F";
          case 10:
            return "Ne";
          case 11:
            return "Na";
          case 12:
            return "Mg";
          case 13:
            return "Al";
          case 14:
            return "Si";
          case 15:
            return "P";
          case 16:
            return "S";
          case 17:
            return "Cl";
          case 18:
            return "Ar";
          case 19:
            return "K";
          case 20:
            return "Ca";
          case 21:
            return "Sc";
          case 22:
            return "Ti";
          case 23:
            return "V";
          case 24:
            return "Cr";
          case 25:
            return "Mn";
          case 26:
            return "Fe";
          case 27:
            return "Co";
          case 28:
            return "Ni";
          case 29:
            return "Cu";
          case 30:
            return "Zn";
          case 31:
            return "Ga";
          case 32:
            return "Ge";
          case 33:
            return "As";
          case 34:
            return "Se";
          case 35:
            return "Br";
          case 36:
            return "Kr";
          case 37:
            return "Rb";
          case 38:
            return "Sr";
          case 39:
            return "Y";
          case 40:
            return "Zr";
          case 41:
            return "Nb";
          case 42:
            return "Mo";
          case 43:
            return "Tc";
          case 44:
            return "Ru";
          case 45:
            return "Rh";
          case 46:
            return "Pd";
          case 47:
            return "Ag";
          case 48:
            return "Cd";
          case 49:
            return "In";
          case 50:
            return "Sn";
          case 51:
            return "Sb";
          case 52:
            return "Te";
          case 53:
            return "I";
          case 54:
            return "Xe";
          case 55:
            return "Cs";
          case 56:
            return "Ba";
          case 57:
            return "La";
          case 58:
            return "Ce";
          case 59:
            return "Pr";
          case 60:
            return "Nd";
          case 61:
            return "Pm";
          case 62:
            return "Sm";
          case 63:
            return "Eu";
          case 64:
            return "Gd";
          case 65:
            return "Tb";
          case 66:
            return "Dy";
          case 67:
            return "Ho";
          case 68:
            return "Er";
          case 69:
            return "Tm";
          case 70:
            return "Yb";
          case 71:
            return "Lu";
          case 72:
            return "Hf";
          case 73:
            return "Ta";
          case 74:
            return "W";
          case 75:
            return "Re";
          case 76:
            return "Os";
          case 77:
            return "Ir";
          case 78:
            return "Pt";
          case 79:
            return "Au";
          case 80:
            return "Hg";
          case 81:
            return "Tl";
          case 82:
            return "Pb";
          case 83:
            return "Bi";
          case 84:
            return "Po";
          case 85:
            return "At";
          case 86:
            return "Rn";
          case 87:
            return "Fr";
          case 88:
            return "Ra";
          case 89:
            return "Ac";
          case 90:
            return "Th";
          case 91:
            return "Pa";
          case 92:
            return "U";
          case 93:
            return "Np";
          case 94:
            return "Pu";
          case 95:
            return "Am";
          case 96:
            return "Cm";
          case 97:
            return "Bk";
          case 98:
            return "Cf";
          case 99:
            return "Es";
          case 100:
            return "Fm";
          case 101:
            return "Md";
          case 102:
            return "No";
          case 103:
            return "Lr";
          case 104:
            return "Rf";
          case 105:
            return "Db";
          case 106:
            return "Sg";
          case 107:
            return "Bh";
          case 108:
            return "Hs";
          case 109:
            return "Mt";
          case 110:
            return "Ds";
          case 111:
            return "Rg";
          case 112:
            return "Cn";
          default:
            return "Xx";
        }
      }

      /**
       * @brief Get the valence for an atom.
       *
       * @param element The atom element number.
       *
       * @return The average mass.
       */
      static int averageMass(int element)
      {
        switch (element) {
          case 1:
            return 1;
          case 2:
            return 4;
          case 3:
            return 7;
          case 4:
            return 9;
          case 5:
            return 11;
          case 6:
            return 12;
          case 7:
            return 14;
          case 8:
            return 16;
          case 9:
            return 19;
          case 10:
            return 20;
          case 11:
            return 23;
          case 12:
            return 24;
          case 13:
            return 27;
          case 14:
            return 28;
          case 15:
            return 31;
          case 16:
            return 32;
          case 17:
            return 35;
          case 18:
            return 40;
          case 19:
            return 39;
          case 20:
            return 40;
          case 21:
            return 45;
          case 22:
            return 48;
          case 23:
            return 51;
          case 24:
            return 52;
          case 25:
            return 55;
          case 26:
            return 56;
          case 27:
            return 59;
          case 28:
            return 59;
          case 29:
            return 64;
          case 30:
            return 65;
          case 31:
            return 70;
          case 32:
            return 73;
          case 33:
            return 75;
          case 34:
            return 79;
          case 35:
            return 80;
          case 36:
            return 84;
          case 37:
            return 85;
          case 38:
            return 88;
          case 39:
            return 89;
          case 40:
            return 91;
          case 41:
            return 93;
          case 42:
            return 96;
          case 43:
            return 98;
          case 44:
            return 101;
          case 45:
            return 103;
          case 46:
            return 106;
          case 47:
            return 108;
          case 48:
            return 112;
          case 49:
            return 115;
          case 50:
            return 119;
          case 51:
            return 122;
          case 52:
            return 128;
          case 53:
            return 127;
          case 54:
            return 131;
          case 55:
            return 133;
          case 56:
            return 137;
          case 57:
            return 139;
          case 58:
            return 140;
          case 59:
            return 141;
          case 60:
            return 144;
          case 61:
            return 145;
          case 62:
            return 150;
          case 63:
            return 152;
          case 64:
            return 157;
          case 65:
            return 159;
          case 66:
            return 163;
          case 67:
            return 165;
          case 68:
            return 167;
          case 69:
            return 169;
          case 70:
            return 173;
          case 71:
            return 175;
          case 72:
            return 178;
          case 73:
            return 181;
          case 74:
            return 184;
          case 75:
            return 186;
          case 76:
            return 190;
          case 77:
            return 192;
          case 78:
            return 195;
          case 79:
            return 197;
          case 80:
            return 201;
          case 81:
            return 204;
          case 82:
            return 207;
          case 83:
            return 209;
          case 84:
            return 209;
          case 85:
            return 210;
          case 86:
            return 222;
          case 87:
            return 223;
          case 88:
            return 226;
          case 89:
            return 227;
          case 90:
            return 232;
          case 91:
            return 231;
          case 92:
            return 238;
          case 93:
            return 237;
          case 94:
            return 244;
          case 95:
            return 243;
          case 96:
            return 247;
          case 97:
            return 247;
          case 98:
            return 251;
          case 99:
            return 252;
          case 100:
            return 257;
          case 101:
            return 258;
          case 102:
            return 259;
          case 103:
            return 260;
          case 104:
            return 261;
          case 105:
            return 268;
          case 106:
            return 271;
          case 107:
            return 267;
          case 108:
            return 277;
          case 109:
            return 276;
          case 110:
            return 281;
          case 111:
            return 280;
          case 112:
            return 285;
          default:
            return 0;
        }
      }

      /**
       * @brief Check if hydrogens should be added for this atom.
       *
       * @param element The atom element number.
       *
       * @return True if hydrogens should be added.
       */
      static bool addHydrogens(int element)
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

      /**
       * @brief Get the valence for an atom
       *
       * @param element The atom element number.
       * @param charge The atom charge.
       * @param degree The current atom degree.
       *
       * @return The valence.
       */
      static int valence(int element, int charge, int degree)
      {
        if (charge < -2 || charge > 2)
          return 0;
        switch (element) {
          case 1:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 2:
            return 0;
          case 3:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 4:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                if (degree < 1)
                  return 1;
                return 0;
              case 2:
                return 0;
            }
          case 5:
            switch (charge) {
              case -2:
                if (degree < 3)
                  return 3;
                return 0;
              case -1:
                if (degree < 4)
                  return 4;
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                if (degree < 2)
                  return 2;
                return 0;
              case 2:
                if (degree < 1)
                  return 1;
                return 0;
            }
          case 6:
            switch (charge) {
              case -2:
                if (degree < 2)
                  return 2;
                return 0;
              case -1:
                if (degree < 3)
                  return 3;
                return 0;
              case 0:
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                if (degree < 3)
                  return 3;
                return 0;
              case 2:
                if (degree < 2)
                  return 2;
                return 0;
            }
          case 7:
            switch (charge) {
              case -2:
                if (degree < 1)
                  return 1;
                return 0;
              case -1:
                if (degree < 2)
                  return 2;
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                if (degree < 4)
                  return 4;
                return 0;
              case 2:
                if (degree < 3)
                  return 3;
                return 0;
            }
          case 8:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                if (degree < 1)
                  return 1;
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 2:
                if (degree < 4)
                  return 4;
                return 0;
            }
          case 9:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                if (degree < 2)
                  return 2;
                return 0;
              case 2:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
            }
          case 10:
            return 0;
          case 11:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 12:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                if (degree < 1)
                  return 1;
                return 0;
              case 2:
                return 0;
            }
          case 13:
            switch (charge) {
              case -2:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case -1:
                if (degree < 4)
                  return 4;
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                if (degree < 2)
                  return 2;
                return 0;
              case 2:
                if (degree < 1)
                  return 1;
                return 0;
            }
          case 14:
            switch (charge) {
              case -2:
                if (degree < 2)
                  return 2;
                return 0;
              case -1:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 0:
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                if (degree < 3)
                  return 3;
                return 0;
              case 2:
                if (degree < 2)
                  return 2;
                return 0;
            }
          case 15:
            switch (charge) {
              case -2:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case -1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 1:
                if (degree < 4)
                  return 4;
                return 0;
              case 2:
                if (degree < 3)
                  return 3;
                return 0;
            }
          case 16:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 2:
                if (degree < 4)
                  return 4;
                return 0;
            }
          case 17:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case 1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 2:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
            }
          case 18:
            return 0;
          case 19:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 20:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                if (degree < 1)
                  return 1;
                return 0;
              case 2:
                return 0;
            }
          case 21:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 22:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 23:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 5)
                  return 5;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 24:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 25:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 26:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 27:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 28:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 29:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 30:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 31:
            switch (charge) {
              case -2:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case -1:
                if (degree < 4)
                  return 4;
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                if (degree < 1)
                  return 1;
                return 0;
            }
          case 32:
            switch (charge) {
              case -2:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case -1:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 0:
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                if (degree < 3)
                  return 3;
                return 0;
              case 2:
                return 0;
            }
          case 33:
            switch (charge) {
              case -2:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case -1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 1:
                if (degree < 4)
                  return 4;
                return 0;
              case 2:
                if (degree < 3)
                  return 3;
                return 0;
            }
          case 34:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 2:
                if (degree < 4)
                  return 4;
                return 0;
            }
          case 35:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case 1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 2:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
            }
          case 36:
            return 0;
          case 37:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 38:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                if (degree < 1)
                  return 1;
                return 0;
              case 2:
                return 0;
            }
          case 39:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 40:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 41:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 42:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 5)
                  return 5;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 43:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 7)
                  return 7;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 44:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 45:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 46:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 47:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 48:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 49:
            switch (charge) {
              case -2:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case -1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                if (degree < 1)
                  return 1;
                return 0;
            }
          case 50:
            switch (charge) {
              case -2:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case -1:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                if (degree < 3)
                  return 3;
                return 0;
              case 2:
                return 0;
            }
          case 51:
            switch (charge) {
              case -2:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case -1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                return 0;
              case 2:
                if (degree < 3)
                  return 3;
                return 0;
            }
          case 52:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 2:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                return 0;
            }
          case 53:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case 1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 2:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
            }
          case 54:
            return 0;
          case 55:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 56:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                if (degree < 1)
                  return 1;
                return 0;
              case 2:
                return 0;
            }
          case 57:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 58:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 59:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 60:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 61:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 62:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 63:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 64:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 65:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 66:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 67:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 68:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 69:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 70:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 71:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 72:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 73:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 5)
                  return 5;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 74:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 5)
                  return 5;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 75:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                if (degree < 7)
                  return 7;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 76:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 77:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 78:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 79:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 80:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 81:
            switch (charge) {
              case -2:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case -1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 82:
            switch (charge) {
              case -2:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case -1:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                if (degree < 3)
                  return 3;
                return 0;
              case 2:
                return 0;
            }
          case 83:
            switch (charge) {
              case -2:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case -1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                return 0;
              case 2:
                if (degree < 3)
                  return 3;
                return 0;
            }
          case 84:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
              case 2:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                return 0;
            }
          case 85:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                if (degree < 7)
                  return 7;
                return 0;
              case 1:
                if (degree < 2)
                  return 2;
                if (degree < 4)
                  return 4;
                if (degree < 6)
                  return 6;
                return 0;
              case 2:
                if (degree < 3)
                  return 3;
                if (degree < 5)
                  return 5;
                return 0;
            }
          case 86:
            return 0;
          case 87:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 88:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                if (degree < 1)
                  return 1;
                return 0;
              case 2:
                return 0;
            }
          case 89:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 90:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 91:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 5)
                  return 5;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 92:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 5)
                  return 5;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 93:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 5)
                  return 5;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 94:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 5)
                  return 5;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 95:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                if (degree < 5)
                  return 5;
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 96:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 97:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 98:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 99:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 100:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 101:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 102:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 2)
                  return 2;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 103:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 3)
                  return 3;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 104:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 4)
                  return 4;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 105:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 5)
                  return 5;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 106:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 6)
                  return 6;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 107:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 7)
                  return 7;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 108:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 109:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 110:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 111:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          case 112:
            switch (charge) {
              case -2:
                return 0;
              case -1:
                return 0;
              case 0:
                if (degree < 1)
                  return 1;
                return 0;
              case 1:
                return 0;
              case 2:
                return 0;
            }
          default:
            return 0;
        }
      }
  };

}

#endif
