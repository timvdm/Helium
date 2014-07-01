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
#ifndef HELIUM_DISTANCEMATRIX_H
#define HELIUM_DISTANCEMATRIX_H

#include <Helium/molecule.h>

#include <vector>
#include <limits>
#include <iostream>

namespace Helium {

  /**
   * @file distancematrix.h
   * @brief DistanceMatrix class.
   */

  /**
   * @class DistanceMatrix distancematrix.h <Helium/distancematrix.h>
   * @brief Class representing a distance matrix.
   *
   * This is a regular distance matrix which requires n * n elements of memory.
   * There is also the SymmetricDistanceMatrix which only stores non-redundant
   * entries.
   */
  class DistanceMatrix
  {
    public:
      /**
       * Constructor.
       *
       * @param n The number of atoms/vertices.
       * @param initDiagonal The initial value for diagonal entries.
       * @param initOffDiagonal The initial value for off-diagonal entries.
       */
      DistanceMatrix(Size n, Size initDiagonal = 0,
          Size initOffDiagonal = DistanceMatrix::infinity())
        : m_matrix(n * n, initOffDiagonal), m_n(n)
      {
        for (Size i = 0; i < n; ++i)
          this->operator()(i, i) = initDiagonal;
      }

      /**
       * @brief Get the dimensionality of the matrix.
       *
       * This is the number of atoms/vertices.
       *
       * @return The dimensionality.
       */
      Size dim() const
      {
        return m_n;
      }

      /**
       * @brief Get a constant reference to element (i, j).
       *
       * @param i The row index.
       * @param j The column index.
       *
       * @return A constant reference to element (i, j).
       */
      const Size& operator()(Size i, Size j) const
      {
        return m_matrix[index(i, j)];
      }

      /**
       * @brief Get a reference to element (i, j).
       *
       * @param i The row index.
       * @param j The column index.
       *
       * @return A reference to element (i, j).
       */
      Size& operator()(Size i, Size j)
      {
        return m_matrix[index(i, j)];
      }

      /**
       * @brief Get the value used for infinity.
       *
       * @return The infinity value.
       */
      static Size infinity()
      {
        return std::numeric_limits<Size>::max() - 99;
      }

    private:
      /**
       * @brief Get the row-major order index for element (i, j).
       *
       * @param i The row index.
       * @param j The column index.
       */
      std::size_t index(Size i, Size j) const
      {
        return i * m_n + j;
      }

      std::vector<Size> m_matrix; //!< The data.
      Size m_n; //!< The dimensionality.
  };

  /**
   * @brief Output operator for DistanceMatrix.
   *
   * @param os STL output stream.
   * @param D The distance matrix.
   */
  inline std::ostream& operator<<(std::ostream &os, const DistanceMatrix &D)
  {
    for (Size i = 0; i < D.dim(); ++i) {
      std::cout << "[ ";
      for (Size j = 0; j < D.dim(); ++j)
        if (D(i, j) == D.infinity())
          std::cout << "- ";
        else
          std::cout << D(i, j) << " ";
      std::cout << "]" << std::endl;
    }

    return os;
  }

  /**
   * @brief Class representing a symmetric distance matrix.
   *
   * This class only stores non-redundant entries. Updating element (i, j) also
   * updates element (j, i).
   */
  class SymmetricDistanceMatrix
  {
    public:
      /**
       * Constructor.
       *
       * @param n The number of atoms/vertices.
       * @param initDiagonal The initial value for diagonal entries.
       * @param initOffDiagonal The initial value for off-diagonal entries.
       */
      SymmetricDistanceMatrix(Size n, Size initDiagonal = 0,
          Size initOffDiagonal = DistanceMatrix::infinity())
        : m_matrix((n * n - n) / 2 + n, initOffDiagonal), m_n(n)
      {
        for (Size i = 0; i < n; ++i)
          this->operator()(i, i) = initDiagonal;
      }

      /**
       * @brief Get the dimensionality of the matrix.
       *
       * This is the number of atoms/vertices.
       *
       * @return The dimensionality.
       */
      Size dim() const
      {
        return m_n;
      }

      /**
       * @brief Get a constant reference to element (i, j).
       *
       * @param i The row index.
       * @param j The column index.
       *
       * @return A constant reference to element (i, j).
       */
      const Size& operator()(Size i, Size j) const
      {
        return m_matrix[index(i, j)];
      }

      /**
       * @brief Get a reference to element (i, j).
       *
       * @param i The row index.
       * @param j The column index.
       *
       * @return A reference to element (i, j).
       */
      Size& operator()(Size i, Size j)
      {
        return m_matrix[index(i, j)];
      }

      /**
       * @brief Get the value used for infinity.
       *
       * @return The infinity value.
       */
      static Size infinity()
      {
        return std::numeric_limits<Size>::max();
      }

    private:
      /**
       * @brief Get the row-major order index for element (i, j).
       *
       * @param i The row index.
       * @param j The column index.
       */
      std::size_t index(Size i, Size j) const
      {
        Size n = std::max(i, j);
        return n * (n + 1) / 2 + std::min(i, j);
      }

      std::vector<Size> m_matrix; //!< The data.
      Size m_n; //!< The dimensionality.
  };

  /**
   * @brief Output operator for SymmetricDistanceMatrix.
   *
   * @param os STL output stream.
   * @param D The distance matrix.
   */
  inline std::ostream& operator<<(std::ostream &os, const SymmetricDistanceMatrix &D)
  {
    for (Size i = 0; i < D.dim(); ++i) {
      std::cout << "[ ";
      for (Size j = 0; j < D.dim(); ++j)
        if (D(i, j) == D.infinity())
          std::cout << "- ";
        else
          std::cout << D(i, j) << " ";
      std::cout << "]" << std::endl;
    }

    return os;
  }

}

#endif
