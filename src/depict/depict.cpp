/*
 * Copyright (C) 2009-2010,2014 by Tim Vandermeersch
 * Some portions Copyright (C) 2009 by Chris Morley
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

#include <Helium/depict/depict.h>

#include <algorithm> // std::reverse
#include <iterator> // std::istream_iterator
#include <cmath>
#include <iostream>

using namespace std;

namespace Helium {

  namespace impl {

    Color elementColor(int element)
    {
      switch (element) {
        case 0:
          return Color(0.07, 0.50, 0.70);
        case 1:
          return Color(0.75, 0.75, 0.75);
        case 2:
          return Color(0.85, 1.00, 1.00);
        case 3:
          return Color(0.80, 0.50, 1.00);
        case 4:
          return Color(0.76, 1.00, 0.00);
        case 5:
          return Color(1.00, 0.71, 0.71);
        case 6:
          return Color(0.40, 0.40, 0.40);
        case 7:
          return Color(0.05, 0.05, 1.00);
        case 8:
          return Color(1.00, 0.05, 0.05);
        case 9:
          return Color(0.50, 0.70, 1.00);
        case 10:
          return Color(0.70, 0.89, 0.96);
        case 11:
          return Color(0.67, 0.36, 0.95);
        case 12:
          return Color(0.54, 1.00, 0.00);
        case 13:
          return Color(0.75, 0.65, 0.65);
        case 14:
          return Color(0.50, 0.60, 0.60);
        case 15:
          return Color(1.00, 0.50, 0.00);
        case 16:
          return Color(0.70, 0.70, 0.00);
        case 17:
          return Color(0.12, 0.94, 0.12);
        case 18:
          return Color(0.50, 0.82, 0.89);
        case 19:
          return Color(0.56, 0.25, 0.83);
        case 20:
          return Color(0.24, 1.00, 0.00);
        case 21:
          return Color(0.90, 0.90, 0.90);
        case 22:
          return Color(0.75, 0.76, 0.78);
        case 23:
          return Color(0.65, 0.65, 0.67);
        case 24:
          return Color(0.54, 0.60, 0.78);
        case 25:
          return Color(0.61, 0.48, 0.78);
        case 26:
          return Color(0.88, 0.40, 0.20);
        case 27:
          return Color(0.94, 0.56, 0.63);
        case 28:
          return Color(0.31, 0.82, 0.31);
        case 29:
          return Color(0.78, 0.50, 0.20);
        case 30:
          return Color(0.49, 0.50, 0.69);
        case 31:
          return Color(0.76, 0.56, 0.56);
        case 32:
          return Color(0.40, 0.56, 0.56);
        case 33:
          return Color(0.74, 0.50, 0.89);
        case 34:
          return Color(1.00, 0.63, 0.00);
        case 35:
          return Color(0.65, 0.16, 0.16);
        case 36:
          return Color(0.36, 0.72, 0.82);
        case 37:
          return Color(0.44, 0.18, 0.69);
        case 38:
          return Color(0.00, 1.00, 0.00);
        case 39:
          return Color(0.58, 1.00, 1.00);
        case 40:
          return Color(0.58, 0.88, 0.88);
        case 41:
          return Color(0.45, 0.76, 0.79);
        case 42:
          return Color(0.33, 0.71, 0.71);
        case 43:
          return Color(0.23, 0.62, 0.62);
        case 44:
          return Color(0.14, 0.56, 0.56);
        case 45:
          return Color(0.04, 0.49, 0.55);
        case 46:
          return Color(0.00, 0.41, 0.52);
        case 47:
          return Color(0.88, 0.88, 1.00);
        case 48:
          return Color(1.00, 0.85, 0.56);
        case 49:
          return Color(0.65, 0.46, 0.45);
        case 50:
          return Color(0.40, 0.50, 0.50);
        case 51:
          return Color(0.62, 0.39, 0.71);
        case 52:
          return Color(0.83, 0.48, 0.00);
        case 53:
          return Color(0.58, 0.00, 0.58);
        case 54:
          return Color(0.26, 0.62, 0.69);
        case 55:
          return Color(0.34, 0.09, 0.56);
        case 56:
          return Color(0.00, 0.79, 0.00);
        case 57:
          return Color(0.44, 0.83, 1.00);
        case 58:
          return Color(1.00, 1.00, 0.78);
        case 59:
          return Color(0.85, 1.00, 0.78);
        case 60:
          return Color(0.78, 1.00, 0.78);
        case 61:
          return Color(0.64, 1.00, 0.78);
        case 62:
          return Color(0.56, 1.00, 0.78);
        case 63:
          return Color(0.38, 1.00, 0.78);
        case 64:
          return Color(0.27, 1.00, 0.78);
        case 65:
          return Color(0.19, 1.00, 0.78);
        case 66:
          return Color(0.12, 1.00, 0.78);
        case 67:
          return Color(0.00, 1.00, 0.61);
        case 68:
          return Color(0.00, 0.90, 0.46);
        case 69:
          return Color(0.00, 0.83, 0.32);
        case 70:
          return Color(0.00, 0.75, 0.22);
        case 71:
          return Color(0.00, 0.67, 0.14);
        case 72:
          return Color(0.30, 0.76, 1.00);
        case 73:
          return Color(0.30, 0.65, 1.00);
        case 74:
          return Color(0.13, 0.58, 0.84);
        case 75:
          return Color(0.15, 0.49, 0.67);
        case 76:
          return Color(0.15, 0.40, 0.59);
        case 77:
          return Color(0.09, 0.33, 0.53);
        case 78:
          return Color(0.96, 0.93, 0.82);
        case 79:
          return Color(0.80, 0.82, 0.12);
        case 80:
          return Color(0.71, 0.71, 0.76);
        case 81:
          return Color(0.65, 0.33, 0.30);
        case 82:
          return Color(0.34, 0.35, 0.38);
        case 83:
          return Color(0.62, 0.31, 0.71);
        case 84:
          return Color(0.67, 0.36, 0.00);
        case 85:
          return Color(0.46, 0.31, 0.27);
        case 86:
          return Color(0.26, 0.51, 0.59);
        case 87:
          return Color(0.26, 0.00, 0.40);
        case 88:
          return Color(0.00, 0.49, 0.00);
        case 89:
          return Color(0.44, 0.67, 0.98);
        case 90:
          return Color(0.00, 0.73, 1.00);
        case 91:
          return Color(0.00, 0.63, 1.00);
        case 92:
          return Color(0.00, 0.56, 1.00);
        case 93:
          return Color(0.00, 0.50, 1.00);
        case 94:
          return Color(0.00, 0.42, 1.00);
        case 95:
          return Color(0.33, 0.36, 0.95);
        case 96:
          return Color(0.47, 0.36, 0.89);
        case 97:
          return Color(0.54, 0.31, 0.89);
        case 98:
          return Color(0.63, 0.21, 0.83);
        case 99:
          return Color(0.70, 0.12, 0.83);
        case 100:
          return Color(0.70, 0.12, 0.73);
        case 101:
          return Color(0.70, 0.05, 0.65);
        case 102:
          return Color(0.74, 0.05, 0.53);
        case 103:
          return Color(0.78, 0.00, 0.40);
        case 104:
          return Color(0.80, 0.00, 0.35);
        case 105:
          return Color(0.82, 0.00, 0.31);
        case 106:
          return Color(0.85, 0.00, 0.27);
        case 107:
          return Color(0.88, 0.00, 0.22);
        case 108:
          return Color(0.90, 0.00, 0.18);
        case 109:
          return Color(0.92, 0.00, 0.15);
        case 110:
          return Color(0.93, 0.00, 0.14);
        case 111:
          return Color(0.94, 0.00, 0.13);
        case 112:
          return Color(0.95, 0.00, 0.12);
        case 113:
          return Color(0.96, 0.00, 0.11);
        case 114:
          return Color(0.97, 0.00, 0.10);
        case 115:
          return Color(0.98, 0.00, 0.09);
        case 116:
          return Color(0.99, 0.00, 0.08);
        case 117:
          return Color(0.99, 0.00, 0.07);
        case 118:
          return Color(0.99, 0.00, 0.06);
      }

      return Color(0, 0, 0);
    }

  }


  enum {
    Left,
    Right,
    Up,
    Down
  };

  Depict::Depict(Painter *painter) : m_painter(painter), m_bondLength(40.0), m_penWidth(2.0),
          m_bondSpacing(6.0), m_bondWidth(8.0), m_fontSize(16), m_subscriptSize(13),
          m_bondColor("black"), m_options(0)
  {
  }

  Depict::~Depict()
  {
  }

  /*
  void Depict::SetDrawingTerminalCarbon(bool enabled)
  {
    d->drawTerminalC = enabled;
  }

  bool Depict::GetDrawingTerminalCarbon() const
  {
    return d->drawTerminalC;
  }
  */

  void Depict::setFontSize(int pointSize, bool subscript)
  {
    if (subscript) {
      m_subscriptSize = pointSize;
      return;
    }

    m_fontSize = pointSize;
    m_subscriptSize = (int)(0.85 * pointSize);
  }

  int Depict::fontSize(bool subscript) const
  {
    if (subscript)
      return m_subscriptSize;
    return m_fontSize;
  }


  /*
  bool Depict::AddAtomLabels(AtomLabelType type)
  {
    m_painter->SetPenColor(OBColor("red"));
    m_painter->SetFillColor(OBColor("red"));
    m_painter->setFontSize((int)(GetFontSize() * 0.8));// smaller text
    OBAtomIterator i;
    for (OBAtom *atom = d->mol->BeginAtom(i); atom; atom = d->mol->NextAtom(i)) {
      Eigen::Vector2d pos(atom->GetVector());
      std::stringstream ss;
      switch (type) {
        case AtomId:
          ss << atom->GetId();
          m_painter->drawText(pos.x(), pos.y(), ss.str());
          break;
        case AtomSymmetryClass:
          ss << GetAtomSymClass(atom);
          m_painter->drawText(pos.x(), pos.y(), ss.str());
          break;
        case AtomIndex:
          ss << atom->GetIdx();
          m_painter->drawText(pos.x(), pos.y(), ss.str());
          break;

        default:
          break;
      }
    }

    return true;
  }
  */


  void Depict::drawWobblyBond(Eigen::Vector2d begin, Eigen::Vector2d end,
      bool beginLbl, bool endLbl)
  {
    Eigen::Vector2d vb = end - begin;

    if (beginLbl)
      begin += 0.33 * vb;
    if (endLbl)
      end -= 0.33 * vb;

    vb = end - begin; // Resize the extents of the vb vector

    Eigen::Vector2d orthogonalLine(vb.y(), -vb.x());
    orthogonalLine.normalize();
    orthogonalLine *= 0.5 * m_bondWidth;

    double lines[6] = { 0.20, 0.36, 0.52, 0.68, 0.84, 1.0 };

    // This code is adapted from drawWedge():
    // What we do is just join up the opposite ends of each of the wedge strokes
    // to create a zig-zag bond

    double oldx, oldy, newx, newy;
    oldx = begin.x();
    oldy = begin.y();
    int sign = 1;
    for (int k = 0; k < 6; ++k) {
      double w = lines[k];
      newx = begin.x() + vb.x() * w + sign * orthogonalLine.x() * w;
      newy = begin.y() + vb.y() * w + sign * orthogonalLine.y() * w;
      m_painter->drawLine(oldx, oldy, newx, newy);
      oldx = newx;
      oldy = newy;
      sign = -sign;
    }
  }

  void Depict::drawWedge(Eigen::Vector2d begin, Eigen::Vector2d end,
      bool beginLbl, bool endLbl)
  {
    Eigen::Vector2d vb = end - begin;

    if (beginLbl)
      begin += 0.33 * vb;
    if (endLbl)
      end -= 0.33 * vb;

    Eigen::Vector2d orthogonalLine(vb.y(), -vb.x());
    orthogonalLine.normalize();
    orthogonalLine *= 0.5 * m_bondWidth;
    std::vector<std::pair<double,double> > points;

    points.push_back(std::pair<double,double>(begin.x(), begin.y()));
    points.push_back(std::pair<double,double>(end.x() + orthogonalLine.x(),
                                              end.y() + orthogonalLine.y()));
    points.push_back(std::pair<double,double>(end.x() - orthogonalLine.x(),
                                              end.y() - orthogonalLine.y()));
    m_painter->drawPolygon(points);
  }

  void Depict::drawHash(Eigen::Vector2d begin, Eigen::Vector2d end,
      bool beginLbl, bool endLbl)
  {
    Eigen::Vector2d vb = end - begin;

    if (beginLbl)
      begin += 0.33 * vb;
    if (endLbl)
      end -= 0.33 * vb;

    vb = end - begin; // Resize the extents of the vb vector

    Eigen::Vector2d orthogonalLine(vb.y(), -vb.x());
    orthogonalLine.normalize();
    orthogonalLine *= 0.5 * m_bondWidth;

    double lines[6] = { 0.20, 0.36, 0.52, 0.68, 0.84, 1.0 };
    double oldwidth = m_painter->penWidth();
    m_painter->setPenWidth(1);
    for (int k = 0; k < 6; ++k) {
      double w = lines[k];
      m_painter->drawLine(begin.x() + vb.x() * w + orthogonalLine.x() * w,
                        begin.y() + vb.y() * w + orthogonalLine.y() * w,
                        begin.x() + vb.x() * w - orthogonalLine.x() * w,
                        begin.y() + vb.y() * w - orthogonalLine.y() * w);
    }

    m_painter->setPenWidth(oldwidth);
  }

  void Depict::drawSimpleBond(Eigen::Vector2d begin, Eigen::Vector2d end,
      bool beginLbl, bool endLbl, int beginValence, int endValence, int order,
      bool crossed_dbl_bond)
  {
    Eigen::Vector2d vb = end - begin;
    vb.normalize();

    if (beginLbl)
      begin += 13. * vb; // Length is normally 40
    if (endLbl)
      end -= 13. * vb;

    if (order == 1 || order == 5) {
      m_painter->drawLine(begin.x(), begin.y(), end.x(), end.y());
    } else if (order == 2) {
      Eigen::Vector2d orthogonalLine(vb.y(), -vb.x());

      bool useAsymmetricDouble = m_options & Depict::AsymmetricDoubleBond;
      if (beginLbl && endLbl)
        useAsymmetricDouble = false;
      if (beginLbl && endValence == 3)
        useAsymmetricDouble = false;
      if (endLbl && beginValence == 3)
        useAsymmetricDouble = false;
      if (crossed_dbl_bond)
        useAsymmetricDouble = false; // The bond looks very strange otherwise in the case of cis

      if (!useAsymmetricDouble) {
        // style1
        //
        // -----------
        // -----------
        Eigen::Vector2d offset = orthogonalLine * 0.5 * m_bondSpacing;
        if (!crossed_dbl_bond) {
          m_painter->drawLine(begin.x() + offset.x(), begin.y() + offset.y(),
                            end.x() + offset.x(), end.y() + offset.y());
          m_painter->drawLine(begin.x() - offset.x(), begin.y() - offset.y(),
                            end.x() - offset.x(), end.y() - offset.y());
        }
        else {
          m_painter->drawLine(begin.x() + offset.x(), begin.y() + offset.y(),
                            end.x() - offset.x(), end.y() - offset.y());
          m_painter->drawLine(begin.x() - offset.x(), begin.y() - offset.y(),
                            end.x() + offset.x(), end.y() + offset.y());
        }
      } else {
        // style2
        //
        //   -------
        // -----------
        Eigen::Vector2d offset1 = orthogonalLine * /*0.5 * */ m_bondSpacing;
        Eigen::Vector2d offset2 = vb * /*0.5 * */ m_bondSpacing;
        Eigen::Vector2d offset3 = -vb * /*0.5 * */ m_bondSpacing;

        if (beginLbl)
          offset2 = Eigen::Vector2d::Zero();
        if (endLbl)
          offset3 = Eigen::Vector2d::Zero();

        m_painter->drawLine(begin.x(), begin.y(), end.x(), end.y());
        m_painter->drawLine(begin.x() + offset1.x() + offset2.x(),
                          begin.y() + offset1.y() + offset2.y(),
                          end.x() + offset1.x() + offset3.x(),
                          end.y() + offset1.y() + offset3.y());
      }
    } else if (order == 3) {
      Eigen::Vector2d orthogonalLine(vb.y(), -vb.x());
      Eigen::Vector2d offset = orthogonalLine * 0.7 * m_bondSpacing;
      m_painter->drawLine(begin.x(), begin.y(), end.x(), end.y());
      m_painter->drawLine(begin.x() + offset.x(), begin.y() + offset.y(),
                        end.x() + offset.x(), end.y() + offset.y());
      m_painter->drawLine(begin.x() - offset.x(), begin.y() - offset.y(),
                        end.x() - offset.x(), end.y() - offset.y());
    }
  }

  void Depict::drawRingBond(Eigen::Vector2d begin, Eigen::Vector2d end,
      bool beginLbl, bool endLbl, int beginValence, int endValence, int order,
      const Eigen::Vector2d &center)
  {
    if (order != 2) {
      drawSimpleBond(begin, end, beginLbl, endLbl, beginValence, endValence, order);
      return;
    }

    Eigen::Vector2d vb = (end - begin).normalized();
    Eigen::Vector2d orthogonalLine(vb.y(), -vb.x());
    Eigen::Vector2d spacing = orthogonalLine * m_bondSpacing * 1.2;
    Eigen::Vector2d offset = vb * m_bondSpacing;
    if ((begin + spacing - center).norm() > (begin - spacing - center).norm())
      spacing *= -1.0;

    Eigen::Vector2d vbb = end - begin;
    if (beginLbl)
      begin += 0.33 * vbb;
    if (endLbl)
      end -= 0.33 * vbb;
    m_painter->drawLine(begin.x(), begin.y(), end.x(), end.y());

    if (beginLbl)
      begin -= 0.10 * vbb;
    if (endLbl)
      end += 0.10 * vbb;
    m_painter->drawLine(begin.x() + spacing.x() + offset.x(), begin.y() + spacing.y() + offset.y(),
                      end.x() + spacing.x() - offset.x(), end.y() + spacing.y() - offset.y());
  }

  void Depict::drawAtomLabel(const std::string &label, int alignment, const Eigen::Vector2d &pos)
  {
   /*
    cout << "FontMetrics(" << label << "):" << endl;
    cout << "  ascent = " << metrics.ascent << endl;
    cout << "  descent = " << metrics.descent << endl;
    cout << "  width = " << metrics.width << endl;
    cout << "  height = " << metrics.height << endl;

    m_painter->SetFillColor(OBColor("white"));
    m_painter->SetPenColor(OBColor("white"));
    m_painter->drawCircle(pos.x(), pos.y(), metrics.ascent / 2);
    m_painter->SetPenColor(OBColor("black"));
    */

    // compute the total width
    double totalWidth = 0.0;
    if ((alignment == Right) || (alignment == Left) || (label.find("H") == std::string::npos)) {
      for (unsigned int i = 0; i < label.size(); ++i) {
        if (!isalpha(label[i])) {
          m_painter->setFontSize(m_subscriptSize);
          totalWidth += m_painter->fontMetrics(label.substr(i, 1)).width;
        } else {
          m_painter->setFontSize(m_fontSize);
          totalWidth += m_painter->fontMetrics(label.substr(i, 1)).width;
        }
      }
    } else {
      m_painter->setFontSize(m_fontSize);
      totalWidth = m_painter->fontMetrics(label.substr(0, label.find("H"))).width;
      double width = 0.0;
      for (unsigned int i = label.find("H"); i < label.size(); ++i) {
        if (!isalpha(label[i])) {
          m_painter->setFontSize(m_subscriptSize);
          width += m_painter->fontMetrics(label.substr(i, 1)).width;
        } else {
          m_painter->setFontSize(m_fontSize);
          width += m_painter->fontMetrics(label.substr(i, 1)).width;
        }
      }

      if (width > totalWidth)
        totalWidth = width;
    }

    m_painter->setFontSize(m_fontSize);
    FontMetrics metrics = m_painter->fontMetrics(label);


    std::string str, subscript;
    // compute the horizontal starting position
    double xOffset, yOffset, yOffsetSubscript;
    switch (alignment) {
      case Right:
        xOffset = 0.5 * m_painter->fontMetrics(label.substr(0, 1)).width -
                  m_painter->fontMetrics(label).width;
        break;
      case Left:
        xOffset = - 0.5 * m_painter->fontMetrics(label.substr(label.size()-1, 1)).width;
        break;
      case Up:
      case Down:
        if (label.find("H") != std::string::npos)
          xOffset = - 0.5 * m_painter->fontMetrics(label.substr(0, label.find("H"))).width;
        else
          xOffset = - 0.5 * totalWidth;
        break;
      default:
        xOffset = - 0.5 * totalWidth;
        break;
    }

    // compute the vertical starting position
    yOffset = 0.5 * (metrics.ascent /*- metrics.descent*/);
    yOffsetSubscript = yOffset - metrics.descent;
    double xInitial = xOffset;

    for (unsigned int i = 0; i < label.size(); ++i) {
      if (label[i] == 'H') {
        if ((alignment == Up) || (alignment == Down))
          if (!str.empty()) {
            // write the current string
            m_painter->setFontSize(m_fontSize);
            m_painter->drawText(pos.x() + xOffset, pos.y() + yOffset, str);
            if (alignment == Down) {
              yOffset += metrics.fontSize;
              yOffsetSubscript += metrics.fontSize;
            } else {
              yOffset -= metrics.fontSize;
              yOffsetSubscript -= metrics.fontSize;
            }
            xOffset = xInitial;
            str.clear();
          }
      }


      if (!isalpha(label[i])) {
        if (!str.empty()) {
          // write the current string
          m_painter->setFontSize(m_fontSize);
          FontMetrics metrics = m_painter->fontMetrics(str);
          m_painter->drawText(pos.x() + xOffset, pos.y() + yOffset, str);
          xOffset += metrics.width;
          str.clear();
        }

        subscript += label.substr(i, 1);
      } else {
        if (!subscript.empty()) {
          // write the current subscript
          m_painter->setFontSize(m_subscriptSize);
          FontMetrics metrics = m_painter->fontMetrics(subscript);
          m_painter->drawText(pos.x() + xOffset, pos.y() + yOffsetSubscript, subscript);
          xOffset += metrics.width;
          subscript.clear();
        }

        str += label.substr(i, 1);
      }
    }
    if (!str.empty()) {
      m_painter->setFontSize(m_fontSize);
      //FontMetrics metrics = m_painter->fontMetrics(str);
      m_painter->drawText(pos.x() + xOffset, pos.y() + yOffset, str);
    }
    if (!subscript.empty()) {
      m_painter->setFontSize(m_subscriptSize);
      //FontMetrics metrics = m_painter->fontMetrics(subscript);
      double yOffset = ispunct(subscript[subscript.size()-1]) || ispunct(subscript[0])
        || (subscript.size()>1 && ispunct(subscript[1]))
        ? -yOffsetSubscript : yOffsetSubscript;
      m_painter->drawText(pos.x() + xOffset, pos.y() + yOffset, subscript);
    }

  }


  /*
  void Depict::setWedgeAndHash(OBMol* mol)
  {
    // Remove any existing wedge and hash bonds
    FOR_BONDS_OF_MOL(b,mol)  {
      b->UnsetWedge();
      b->UnsetHash();
    }

    std::map<OBBond*, enum OBStereo::BondDirection> updown;
    std::map<OBBond*, OBStereo::Ref> from;
    std::map<OBBond*, OBStereo::Ref>::const_iterator from_cit;
    TetStereoToWedgeHash(*mol, updown, from);

    for(from_cit=from.begin();from_cit!=from.end();++from_cit) {
      OBBond* pbond = from_cit->first;
      if(updown[pbond]==OBStereo::UpBond)
        pbond->SetHash();
      else if(updown[pbond]==OBStereo::DownBond)
        pbond->SetWedge();
    }
  }
  */

}
