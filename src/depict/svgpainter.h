/*
 * Copyright (C) 2009 by Chris Morley
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
#ifndef HELIUM_SVGPAINTER_H
#define HELIUM_SVGPAINTER_H

#include <Helium/depict/painter.h>
#include <iostream>

namespace Helium {

  class SVGPainter : public Painter
  {
    public:
      SVGPainter(std::ostream& ofs, bool withViewBox=false,
        double width=0.0, double height=0.0, double x=0.0, double y=0.0);
      ~SVGPainter();
      //! @name Painter methods
      //@{
      void newCanvas(double width, double height);
      bool isGood() const;
      void setFontFamily(const std::string &fontFamily);
      void setFontSize(int pointSize);
      void setFillColor(const Color &color);
      void setPenColor(const Color &color);
      void setPenWidth(double width);
      double penWidth();
      void drawLine(double x1, double y1, double x2, double y2);
      void drawPolygon(const std::vector<std::pair<double,double> > &points);
      void drawCircle(double x, double y, double r);
      void drawText(double x, double y, const std::string &text);
      FontMetrics fontMetrics(const std::string &text);
      //@}

    private:
      std::string makeRGB(Color color);

      std::ostream& m_ofs;
      bool m_withViewBox;
      double m_width, m_height, m_x, m_y;
      Color m_Pencolor;
      Color m_OrigBondcolor;
      Color m_Fillcolor;
      double m_PenWidth;
      int m_fontPointSize;
      std::string m_fontFamily;
  };

}

#endif
