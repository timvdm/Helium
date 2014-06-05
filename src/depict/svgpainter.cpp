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

#include <Helium/depict/svgpainter.h>

#include <iostream>
#include <sstream>

namespace Helium {

  SVGPainter::SVGPainter(std::ostream& ofs, bool withViewBox,
    double width, double height, double x, double y)
    :  m_ofs(ofs), m_withViewBox(withViewBox), m_width(width), m_height(height),
       m_x(x), m_y(y), m_Pencolor("black"), m_Fillcolor("white"), m_PenWidth(1),
       m_fontPointSize(16)  {}

  SVGPainter::~SVGPainter()
  {
    m_ofs << "</svg>\n";
    if(m_withViewBox)
      m_ofs << "</g>\n";
  }

  void SVGPainter::newCanvas(double width, double height)
  {
    //Using withViewBox to supress xml header and xmlns attributes. May need another way.
    if(!m_withViewBox)
      m_ofs << "<?xml version=\"1.0\"?>\n";

    if(m_withViewBox)
      m_ofs << "<g transform=\"translate(" << m_x << "," << m_y << ")\">\n";

    m_ofs << "<svg ";
    if(!m_withViewBox)
      m_ofs << "xmlns=\"http://www.w3.org/2000/svg\"\n"
               "xmlns:xlink=\"http://www.w3.org/1999/xlink\" "
               "xmlns:cml=\"http://www.xml-cml.org/schema\" ";
    if(m_withViewBox)
      m_ofs << "width=\"" << m_width << "\" height=\"" << m_height << "\" "
            << "x=\"0\" y=\"0\" "
            << "viewBox=\"0 0 " << width << ' ' << height << "\"\n";
    else
      m_ofs << "width=\"" << width << "\" height=\"" << height << "\" "
            << "x=\"0\" y=\"0\" ";

    //Bond color and width are the initial m_Pencolor and m_PenWidth
    m_ofs << "font-family=\"" << m_fontFamily << "\" stroke=" << makeRGB(m_Pencolor)
          << "stroke-width=\"" << m_PenWidth << "\"  stroke-linecap=\"round\"" << ">\n";

    if(!m_withViewBox && m_Fillcolor.alpha!=0.0)//Background color for single molecule. Handled by outer svg when table.
      m_ofs << "<rect x=\"0%\" y=\"0%\" width=\"100%\" height=\"100%\" stroke-width=\"0\" fill="
            << makeRGB(m_Fillcolor) << " />\n";
    m_OrigBondcolor = m_Pencolor;
  }

  bool SVGPainter::isGood() const
  {
    return true;
  }

  void SVGPainter::setFontSize(int pointSize)
  {
    m_fontPointSize = pointSize;
  }

  void SVGPainter::setFontFamily(const std::string &fontFamily)
  {
    m_fontFamily = fontFamily;
  }

  void SVGPainter::setFillColor(const Color &color)
  {
    m_Fillcolor = color; //value when NewCanvas called used for background
  }

  void SVGPainter::setPenColor(const Color &color)
  {
    m_Pencolor = color; //value when NewCanvas called used for bonds
  }

  void SVGPainter::setPenWidth(double width)
  {
    m_PenWidth = width; //value when NewCanvas called used for bonds
  }

  double SVGPainter::penWidth()
  {
    return m_PenWidth;
  }

  void SVGPainter::drawLine(double x1, double y1, double x2, double y2)
  {
    std::streamsize oldprec = m_ofs.precision(1);
    m_ofs << std::fixed << "<line x1=\"" << x1 << "\" y1=\"" << y1 << "\" x2=\""
      << x2 << "\" y2=\"" << y2 << "\"";
    // if(m_Pencolor!=m_OrigBondcolor) // TODO: Bring this line back once Pybel is fine with this
      m_ofs << " stroke=" << makeRGB(m_Pencolor);
    m_ofs << " stroke-width=\"" << m_PenWidth << "\"";
    m_ofs << "/>\n";
    m_ofs.precision(oldprec);
  }

  void SVGPainter::drawPolygon(const std::vector<std::pair<double,double> > &points)
  {
    m_ofs << "<polygon points=\"";
      std::vector<std::pair<double,double> >::const_iterator i;
    for (i = points.begin(); i != points.end(); ++i)
      m_ofs << i->first << ' ' << i->second << ' ';
    m_ofs << "\"";
    m_ofs << " stroke-width=\"" << m_PenWidth << "\"";
    m_ofs << " fill=" << makeRGB(m_Pencolor);
    m_ofs << " stroke=" << makeRGB(m_Pencolor);
    m_ofs << "/>\n";
  }

  void SVGPainter::drawCircle(double x, double y, double r)
  {
    m_ofs << "<circle cx=\"" << x << "\" cy=\"" << y << "\" r=\"" << r << "\" "
          << " fill=\"rgb(255,255,255)\" stroke=" << makeRGB(m_Pencolor) << "stroke-width=\"" << m_PenWidth << "\" />\n";
  }

  void SVGPainter::drawText(double x, double y, const std::string &text)
  {
    m_ofs << "<text x=\"" << x << "\" y=\"" << y << "\""
      << " fill=" << makeRGB(m_Pencolor) << " stroke=" << makeRGB(m_Pencolor) << "stroke-width=\"1\" "
      << "font-size=\"" << m_fontPointSize << "\" >"
      << text << "</text>\n";
  }

  FontMetrics SVGPainter::fontMetrics(const std::string &text)
  {
    FontMetrics metrics;
    metrics.fontSize = m_fontPointSize;
    metrics.ascent   = m_fontPointSize;
    metrics.descent  = m_fontPointSize * -0.23; // Offset from baseline of bottom of text
    metrics.height   = m_fontPointSize *  1.26; // Distance between successive lines of text
    metrics.width = 0.0;
    for(std::string::size_type i = 0; i < text.size(); ++i)
      metrics.width += m_fontPointSize * (std::isalpha(text[i]) ? 0.75 : 0.5);

    return metrics;
  }

  std::string SVGPainter::makeRGB(Color color)
  {
    std::stringstream ss;
    ss << "\"rgb(" << (int)(255*color.red) << ',' << (int)(255*color.green)
       << ',' << (int)(255*color.blue) << ")\" ";
    return ss.str();
  }

}
