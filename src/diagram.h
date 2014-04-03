#include <Helium/molecule.h>

#include <Helium/algorithms/cycles.h>

#include <Eigen/Core>

namespace Helium {

  typedef Eigen::Vector2d Point2D;
  typedef std::vector<Point2D> Points2D;

  namespace impl {

    double polygon_radius(int n, double side)
    {
      return side / (2.0 * std::sin(deg2rad(180.0 / n)));
    }

    std::vector<Point2D> small_ring(int n, double length)
    {
      assert(n > 2 && n < 9);
      std::vector<Point2D> points;

      double radius = polygon_radius(n, length);
      
      double angle;
      switch (n) {
        case 3:
        case 5:
        case 6:
        case 7:
          angle = 90.0;
          break;
        case 4:
          angle = 45.0;
          break;
        case 8:
          angle = 22.5;
          break;
      }
          
      for (int i = 0; i < n; ++i) {
        double x = radius * std::cos(deg2rad(angle));
        double y = radius * std::sin(deg2rad(angle));
        points.push_back(Point2D(x, y));
        angle += 360.0 / n;
      }

      assert(points.size() == n);
      return points;
    }

    Points2D chain(int n, double length)
    {
      Point2D up(length * std::cos(deg2rad(30.0)), length * std::sin(deg2rad(30.0)));
      Point2D down(length * std::cos(deg2rad(-30.0)), length * std::sin(deg2rad(-30.0)));

      Points2D points;
      points.push_back(Point2D(0.0, 0.0));
      for (int i = 1; i < n; ++i)
        if (points.size() % 2)
          points.push_back(points.back() + up);
        else
          points.push_back(points.back() + down);

      return points;
    }

  }

  template<typename MoleculeType>
  std::vector<Point2D> generate_diagram(const MoleculeType &mol)
  {
    double length = 100.0;
    RingSet<MoleculeType> rings = relevant_cycles(mol); // FIXME: param

    std::vector<Point2D> coords(num_atoms(mol), Point2D::Zero());

    for (int i = 0; i < rings.size(); ++i) {
      Points2D points = impl::small_ring(rings.ring(i).size(), length);
      for (std::size_t j = 0; j < rings.ring(i).size(); ++j)
        coords[get_index(mol, rings.ring(i).atom(j))] = points[j];
    }

    return coords;
  }

}
