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

#include <Eigen/Core>

namespace Helium {

  /**
   * @file depict/painter.h
   * @brief Painter class used by Depict.
   */

  /**
   * @class Color depict/painter.h <Helium/depict/painter.h>
   * @brief Color class used by Depict.
   */
  struct Color
  {
    /**
     * @brief Default constructor.
     */
    Color()
    {
      *this = Color(0.0, 0.0, 0.0);
    }

    /**
     * @brief Constructor.
     *
     * @param _red The red component [0, 1].
     * @param _green The green component [0, 1].
     * @param _blue The blue component [0, 1].
     * @param _alpha The alpha component [0, 1].
     */
    Color(double _red, double _green, double _blue, double _alpha = 1.0) :
        red(_red), green(_green), blue(_blue), alpha(_alpha)
    {
    }

    /**
     * @brief constructor.
     *
     * Acceptable color strings: white, red, green, blue, yelloq, gray, cyan,
     * purple, teal, olive, none (black). It is also possible to specify colors
     * using hexadecimal notation using the # prefix (e.g. "#FF0000").
     *
     * @param color The color string.
     */
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

    /**
     * @brief Constructor.
     *
     * @param vec The color vector (RGB).
     */
    Color(const std::vector<double> &vec) : red(vec[0]), green(vec[1]), blue(vec[2]), alpha(1.0)
    {
    }

    /**
     * @brief Inequality comparison operator.
     *
     * @param other The other color.
     */
    bool operator !=(const Color& other)
    {
      return red!=other.red || green!=other.green || blue!=other.blue;
    }

    double red; //!< The red component.
    double green; //!< The green component.
    double blue; //!< The blue component.
    double alpha; //!< The alpha component.
  };

  /**
   * @class FontMetrics depict/painter.h <Helium/depict/painter.h>
   * @brief Font metrics class used by Depict.
   */
  struct FontMetrics
  {
    int    fontSize; //!< The font size.
    double ascent; //!< The ascent.
    double descent; //!< The descent.
    double width; //!< The width.
    double height; //!< The height.
  };

  /**
   * @class Painter depict/painter.h <Helium/depict/painter.h>
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
       * @brief Set the painter's font family.
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
       * @brief Draw a line.
       *
       * Draw a line from @p x1, @p y1 to @p x2, @p y2. The line is drawn using
       * the current pen color and width.
       */
      virtual void drawLine(double x1, double y1, double x2, double y2) = 0;

      /**
       * @brief Draw a circle.
       *
       * @param x The circle center x coordinate.
       * @param y The circle center y coordinate.
       * @param r The radius.
       */
      virtual void drawCircle(double x, double y, double r) = 0;

      /**
       * Draw a polygon by connecting consecutive points. The last point will be
       * connected to the first one. The lines are drawn using the current pen
       * color and width. The area inside the polygon is filled with the current
       * fill color.
       */
      virtual void drawPolygon(const std::vector<std::pair<double,double> > &points) = 0;

      /**
       * @brief Draw text.
       *
       * @param x The x coordinate.
       * @param y The y coordinate.
       * @param text The text.
       */
      virtual void drawText(double x, double y, const std::string &text) = 0;

      /**
       * @brief Get the font metrics for a piece of text.
       */
      virtual FontMetrics fontMetrics(const std::string &text) = 0;

      /**
       * @brief Draw a dashed line.
       */
      void drawDashedLine(double x1, double y1, double x2, double y2, double size)
      {
        Eigen::Vector2d p1(x1, y1);
        Eigen::Vector2d p2(x2, y2);
        Eigen::Vector2d p1p2 = p2 - p1;
        double norm = p1p2.norm();
        p1p2.normalize();
        Eigen::Vector2d last = p1;
        double drawn = 0.0;
        while (drawn < norm) {
          Eigen::Vector2d next = last + size * p1p2;
          if ((p1 - next).norm() >= norm)
            next = p2;
          drawLine(last.x(), last.y(), next.x(), next.y());
          last += 2 * size * p1p2;
          drawn += 2 * size;
        }
      }
  };

}

#endif
