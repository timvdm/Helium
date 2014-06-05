/*
 * Copyright (C) 2009,2014 by Tim Vandermeersch
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
#ifndef HELIUM_PAINTER_H
#define HELIUM_PAINTER_H

#include <string>
#include <vector>
#include <sstream>

namespace Helium {

  /**
   * @class Color painter.h <Helium/depict/painter.h>
   * @brief Color class used by Depict.
   */
  struct Color
  {
    Color()
    {
      *this = Color(0.0, 0.0, 0.0);
    }

    Color(double _red, double _green, double _blue, double _alpha = 1.0) :
        red(_red), green(_green), blue(_blue), alpha(_alpha)
    {
    }

    Color(const std::string &color)
    {
      if (color[0] == '#') {
        std::stringstream ss(color.substr(1));
        unsigned c;
        ss >> std::hex >> c;
        *this = Color((c / 0x10000) / 256.0,
                      ((c % 0x10000) / 0x100 / 256.0),
                      (c % 0x100) / 256.0);
        return;
      }

      if (color == "black")
        *this = Color(0.0, 0.0, 0.0);
      else if (color == "white")
        *this = Color(1.0, 1.0, 1.0);
      else if (color == "red")
        *this = Color(1.0, 0.0, 0.0);
      else if (color == "green")
        *this = Color(0.0, 1.0, 0.0);
      else if (color == "blue")
        *this = Color(0.0, 0.0, 1.0);
      else if (color == "yellow")
        *this = Color(1.0, 1.0, 0.0);
      else if (color == "gray")
        *this = Color(0.3, 0.3, 0.3);
      else if (color == "cyan")
        *this = Color(1.0, 0.0, 1.0);
      else if (color == "purple")
        *this = Color(0.5, 0.0, 0.5);
      else if (color == "teal")
        *this = Color(0.0, 0.5, 0.5);
      else if (color == "olive")
        *this = Color(0.5, 0.5, 0.0);
      else if (color == "none")
        *this = Color(0.0, 0.0, 0.0, 0.0);
      else
        *this = Color(0.5, 0.5, 0.5);
    }

    Color(std::vector<double> vec) : red(vec[0]), green(vec[1]), blue(vec[2]), alpha(1.0){}

    bool operator !=(const Color& other)
    {
      return red!=other.red || green!=other.green || blue!=other.blue;
    }

    double red, green, blue, alpha;
  };

  /**
   * @class FontMetrics painter.h <Helium/depict/painter.h>
   * @brief Font metrics class used by Depict.
   */
  struct FontMetrics
  {
    int    fontSize;
    double ascent, descent;
    double width, height;
  };

  /**
   * @class Painter painter.h <Helium/depict/painter.h>
   * @brief Abstract painter base class used by Depict.
   */
  class Painter
  {
    public:
      /**
       * Virtual destructor to be inhereted by subclasses
       */
      virtual ~Painter() {}

      /**
       * Create a new canvas to paint on with size @p width x @p height.
       * Depict will always call NewCanvas before performing any drawing
       * operations. Painters that are capable of drawing on a previously
       * unspecified area don't need to implement this.
       */
      virtual void newCanvas(double width, double height) = 0;
      /**
       * Before Depict performes any drawing operation, this method is called
       * to check if the painter is ready to start drawing. If this method
       * returns false, drawing is aborted.
       */
      virtual bool isGood() const = 0;
      /**
       * Set the painter's font family.
       */
      virtual void setFontFamily(const std::string &fontFamily) = 0;
      /**
       * Set the painter's font point size.
       */
      virtual void setFontSize(int pointSize) = 0;
      /**
       * Set the painter's fill color.
       */
      virtual void setFillColor(const Color &color) = 0;
      /**
       * Set the painter's pen color.
       */
      virtual void setPenColor(const Color &color) = 0;
      /**
       * Set the painter's pen width.
       */
      virtual void setPenWidth(double width) = 0;
      /**
       * Get the painter's pen width.
       */
      virtual double penWidth() = 0;
      /**
       * Draw a line from @p x1, @p y1 to @p x2, @p y2. The line is drawn using
       * the current pen color and width.
       */
      virtual void drawLine(double x1, double y1, double x2, double y2) = 0;
      virtual void drawCircle(double x, double y, double r) = 0;
      /**
       * Draw a polygon by connecting consecutive points. The last point will be
       * connected to the first one. The lines are drawn using the current pen
       * color and width. The area inside the polygon is filled with the current
       * fill color.
       */
      virtual void drawPolygon(const std::vector<std::pair<double,double> > &points) = 0;
      virtual void drawText(double x, double y, const std::string &text) = 0;
      virtual FontMetrics fontMetrics(const std::string &text) = 0;
  };

}

#endif
