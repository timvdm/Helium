#ifndef HELIUM_DISTANCEMATRIX_H
#define HELIUM_DISTANCEMATRIX_H

#include <vector>
#include <limits>

#include <Helium/molecule.h>

namespace Helium {

  class DistanceMatrix
  {
    public:
      DistanceMatrix(Size n, Size initDiagonal = 0,
          Size initOffDiagonal = std::numeric_limits<Size>::max())
        : m_matrix((n * n - n) / 2 + n, initOffDiagonal), m_n(n)
      {
        setDiagonal(initDiagonal);
      }

      void setDiagonal(Size value)
      {
        for (Size i = 0; i < m_n; ++i)
          m_matrix[index(i, i)] = value;
      }

      const Size& operator()(Size i, Size j) const
      {
        return m_matrix[index(i, j)];
      }

      Size& operator()(Size i, Size j)
      {
        return m_matrix[index(i, j)];
      }

    private:
      std::size_t index(Size i, Size j) const
      {
        Size n = std::max(i, j);
        return n * (n + 1) / 2 + std::min(i, j);
      }

      std::vector<Size> m_matrix;
      Size m_n;
  };


}

#endif
