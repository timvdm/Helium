/*
 * Copyright (c) 2013, Tim Vandermeersch
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
#ifndef HELIUM_ELEMENT_H
#define HELIUM_ELEMENT_H

#include <Helium/config.h>
#include <string>
#include <cassert>

namespace Helium {

  /**
   * @file element.h
   * @brief Element properties.
   */

  /**
   * @class Element element.h <Helium/element.h>
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
          case 114:
            return "Fl";
          case 116:
            return "Lv";
          default:
            return "Xx";
        }
      }

      static int element(const std::string &symbol)
      {
        assert(symbol.size() >= 1 && symbol.size() <= 2);

        switch (symbol[0]) {
          case 'H':
            if (symbol.size() == 1)
              return 1;
            switch (symbol[1]) {
              case 'e':
                return 2;
              case 'f':
                return 72;
              case 'g':
                return 80;
              case 's':
                return 108;
              case 'o':
                return 67;
            }
            break;
          case 'L':
            if (symbol.size() == 1)
              break;
            switch (symbol[1]) {
              case 'i':
                return 3;
              case 'a':
                return 57;
              case 'u':
                return 71;
              case 'r':
                return 103;
              case 'v':
                return 116;
            }
            break;
          case 'B':
            if (symbol.size() == 1)
              return 5;
            switch (symbol[1]) {
              case 'e':
                return 4;
              case 'r':
                return 35;
              case 'a':
                return 56;
              case 'i':
                return 83;
              case 'h':
                return 107;
              case 'k':
                return 97;
            }
            break;
          case 'C':
            if (symbol.size() == 1)
              return 6;
            switch (symbol[1]) {
              case 'l':
                return 17;
              case 'a':
                return 20;
              case 'r':
                return 24;
              case 'o':
                return 27;
              case 'u':
                return 29;
              case 'd':
                return 48;
              case 's':
                return 55;
              case 'e':
                return 58;
              case 'm':
                return 96;
              case 'f':
                return 98;
              case 'n':
                return 112;
            }
            break;
          case 'N':
            if (symbol.size() == 1)
              return 7;
            switch (symbol[1]) {
              case 'e':
                return 10;
              case 'a':
                return 11;
              case 'i':
                return 28;
              case 'b':
                return 41;
              case 'd':
                return 60;
              case 'p':
                return 93;
              case 'o':
                return 102;
            }
            break;
          case 'O':
            if (symbol.size() == 1)
              return 8;
            if (symbol[1] == 's')
              return 76;
            break;
          case 'F':
            if (symbol.size() == 1)
              return 9;
            switch (symbol[1]) {
              case 'e':
                return 26;
              case 'r':
                return 87;
              case 'm':
                return 100;
              case 'l':
                return 114;
            }
            break;
          case 'M':
            if (symbol.size() == 1)
              break;
            switch (symbol[1]) {
              case 'n':
                return 25;
              case 'o':
                return 42;
              case 't':
                return 109;
              case 'd':
                return 101;
              case 'g':
                return 12;
            }
            break;
          case 'A':
            if (symbol.size() == 1)
              break;
            switch (symbol[1]) {
              case 'l':
                return 13;
              case 'r':
                return 18;
              case 's':
                return 33;
              case 'g':
                return 47;
              case 'u':
                return 79;
              case 't':
                return 85;
              case 'c':
                return 89;
              case 'm':
                return 95;
            }
            break;
          case 'S':
            if (symbol.size() == 1)
              return 16;
            switch (symbol[1]) {
              case 'i':
                return 14;
              case 'c':
                return 21;
              case 'e':
                return 34;
              case 'r':
                return 38;
              case 'n':
                return 50;
              case 'b':
                return 51;
              case 'g':
                return 106;
              case 'm':
                return 62;
            }
            break;
          case 'P':
            if (symbol.size() == 1)
              return 15;
            switch (symbol[1]) {
              case 'd':
                return 46;
              case 't':
                return 78;
              case 'b':
                return 82;
              case 'o':
                return 84;
              case 'r':
                return 59;
              case 'm':
                return 61;
              case 'a':
                return 91;
              case 'u':
                return 94;
            }
            break;
          case 'K':
            if (symbol.size() == 1)
              return 19;
            if (symbol[1] == 'r')
              return 36;
            break;
          case 'T':
            if (symbol.size() == 1)
              return 1;
            switch (symbol[1]) {
              case 'i':
                return 22;
              case 'c':
                return 43;
              case 'e':
                return 52;
              case 'a':
                return 73;
              case 'l':
                return 81;
              case 'b':
                return 65;
              case 'm':
                return 69;
              case 'h':
                return 90;
            }
            break;
          case 'V':
            if (symbol.size() == 1)
              return 23;
            break;
          case 'Z':
            if (symbol.size() == 1)
              break;
            switch (symbol[1]) {
              case 'n':
                return 30;
              case 'r':
                return 40;
            }
            break;
          case 'G':
            if (symbol.size() == 1)
              break;
            switch (symbol[1]) {
              case 'a':
                return 31;
              case 'e':
                return 32;
              case 'd':
                return 64;
            }
            break;
          case 'R':
            if (symbol.size() == 1)
              break;
            switch (symbol[1]) {
              case 'b':
                return 37;
              case 'u':
                return 44;
              case 'h':
                return 45;
              case 'e':
                return 75;
              case 'n':
                return 86;
              case 'a':
                return 88;
              case 'f':
                return 104;
              case 'g':
                return 111;
            }
            break;
          case 'Y':
            if (symbol.size() == 1)
              return 39;
            if (symbol[1] == 'b')
              return 70;
            break;
          case 'I':
            if (symbol.size() == 1)
              return 53;
            switch (symbol[1]) {
              case 'n':
                return 49;
              case 'r':
                return 77;
            }
            break;
          case 'X':
            if (symbol.size() == 1)
              break;
            if (symbol[1] == 'e')
              return 54;
            break;
          case 'W':
            if (symbol.size() == 1)
              return 74;
            break;
          case 'D':
            if (symbol.size() == 1)
              return 1;
            switch (symbol[1]) {
              case 'b':
                return 105;
              case 's':
                return 110;
              case 'y':
                return 66;
            }
            break;
          case 'E':
            if (symbol.size() == 1)
              break;
            switch (symbol[1]) {
              case 'u':
                return 63;
              case 'r':
                return 68;
              case 's':
                return 99;
            }
            break;
          case 'U':
            if (symbol.size() == 1)
              return 92;
            break;
        }

        return 0; // unknown
      }

      /**
       * @brief Get the average mass for an atom.
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
       * @brief Get the valence for an atom.
       *
       * This function will return the expected valence for an atom based on
       * the atom's element, charge and current valence. The values returned by
       * this function are based on the MDL and InChI valences.
       *
       * @param element The atom element number.
       * @param charge The atom charge.
       * @param currentValence The current atom valence.
       *
       * @return The valence.
       */
      static int valence(int element, int charge, int currentValence)
      {
        if (charge < -4 || charge > 6)
          return 0;

        switch (element) {
          case 1:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 2:
            return 0;
          case 3:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 4:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 1:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 5:
            switch (charge) {
              case -4:
                if (currentValence <= 1)
                  return 1;
                return 0;
              case -3:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case -2:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case -1:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 1:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 2:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 6:
            switch (charge) {
              case -3:
                if (currentValence <= 1)
                  return 1;
                return 0;
              case -2:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case -1:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 0:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 1:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 2:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 3:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 7:
            switch (charge) {
              case -2:
                if (currentValence <= 1)
                  return 1;
                return 0;
              case -1:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 1:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 2:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 3:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 4:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 8:
            switch (charge) {
              case -1:
                if (currentValence <= 1)
                  return 1;
                return 0;
              case 0:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 1:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 2:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 3:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 4:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 5:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 9:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              case 1:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 2:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 3:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 4:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 5:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 6:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 10:
            return 0;
          case 11:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 12:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 1:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 13:
            switch (charge) {
              case -4:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -3:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case -2:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case -1:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 1:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 2:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 14:
            switch (charge) {
              case -3:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -2:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case -1:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 0:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 1:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 2:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 3:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 15:
            switch (charge) {
              case -2:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 1:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 2:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 3:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 4:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 16:
            switch (charge) {
              case -1:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 1:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 2:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 3:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 4:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 5:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 17:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case 1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 2:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 3:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 4:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 5:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 6:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 18:
            return 0;
          case 19:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 20:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 1:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 21:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 22:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 23:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 5)
                  return 5;
                return 0;
              default:
                return 0;
            }
          case 24:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 25:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 26:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 27:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 28:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 29:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 2)
                  return 2;
                return 0;
              default:
                return 0;
            }
          case 30:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                return 0;
              default:
                return 0;
            }
          case 31:
            switch (charge) {
              case -4:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -3:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case -2:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case -1:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 2:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 32:
            switch (charge) {
              case -3:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -2:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case -1:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 0:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 1:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 3:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 33:
            switch (charge) {
              case -2:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 1:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 2:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 4:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 34:
            switch (charge) {
              case -1:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 1:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 2:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 3:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 5:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 35:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case 1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 2:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 3:
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 4:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 6:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 36:
            return 0;
          case 37:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 38:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 1:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 39:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 40:
            switch (charge) {
              case 0:
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 41:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              default:
                return 0;
            }
          case 42:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 43:
            switch (charge) {
              case 0:
                if (currentValence <= 7)
                  return 7;
                return 0;
              default:
                return 0;
            }
          case 44:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 45:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 46:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 47:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 48:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                return 0;
              default:
                return 0;
            }
          case 49:
            switch (charge) {
              case -4:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -3:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case -2:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case -1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 2:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 50:
            switch (charge) {
              case -3:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -2:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case -1:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 1:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 3:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 51:
            switch (charge) {
              case -2:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 2:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 4:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 52:
            switch (charge) {
              case -1:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 1:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 2:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 3:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 5:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 53:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case 1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 2:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 3:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 4:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 6:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 54:
            return 0;
          case 55:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 56:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 1:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 57:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 58:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 59:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 60:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 61:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 62:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 63:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 64:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 65:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 66:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 67:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 68:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 69:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 70:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 71:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 72:
            switch (charge) {
              case 0:
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 73:
            switch (charge) {
              case 0:
                if (currentValence <= 5)
                  return 5;
                return 0;
              default:
                return 0;
            }
          case 74:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 75:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                if (currentValence <= 7)
                  return 7;
                return 0;
              default:
                return 0;
            }
          case 76:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 77:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 78:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 79:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 80:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 2)
                  return 2;
                return 0;
              default:
                return 0;
            }
          case 81:
            switch (charge) {
              case -4:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -3:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case -2:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case -1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 0:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 82:
            switch (charge) {
              case -3:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -2:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case -1:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 1:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 3:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 83:
            switch (charge) {
              case -2:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case -1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 2:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 4:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 84:
            switch (charge) {
              case -1:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case 0:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 1:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 2:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 3:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 5:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 85:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 7)
                  return 7;
                return 0;
              case 1:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 6)
                  return 6;
                return 0;
              case 2:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 5)
                  return 5;
                return 0;
              case 3:
                if (currentValence <= 2)
                  return 2;
                if (currentValence <= 4)
                  return 4;
                return 0;
              case 4:
                if (currentValence <= 3)
                  return 3;
                return 0;
              case 6:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 86:
            return 0;
          case 87:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 88:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                return 0;
              case 1:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 89:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 90:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 91:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 5)
                  return 5;
                return 0;
              default:
                return 0;
            }
          case 92:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 93:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 94:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 95:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                if (currentValence <= 5)
                  return 5;
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 96:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 97:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 98:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 99:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 100:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 101:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 102:
            switch (charge) {
              case 0:
                if (currentValence <= 2)
                  return 2;
                return 0;
              default:
                return 0;
            }
          case 103:
            switch (charge) {
              case 0:
                if (currentValence <= 3)
                  return 3;
                return 0;
              default:
                return 0;
            }
          case 104:
            switch (charge) {
              case 0:
                if (currentValence <= 4)
                  return 4;
                return 0;
              default:
                return 0;
            }
          case 105:
            switch (charge) {
              case 0:
                if (currentValence <= 5)
                  return 5;
                return 0;
              default:
                return 0;
            }
          case 106:
            switch (charge) {
              case 0:
                if (currentValence <= 6)
                  return 6;
                return 0;
              default:
                return 0;
            }
          case 107:
            switch (charge) {
              case 0:
                if (currentValence <= 7)
                  return 7;
                return 0;
              default:
                return 0;
            }
          case 108:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 109:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 110:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 111:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          case 112:
            switch (charge) {
              case 0:
                if (currentValence <= 1)
                  return 1;
                return 0;
              default:
                return 0;
            }
          default:
            return 0;
        }
      }

  };

}

#endif
