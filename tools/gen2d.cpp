/**
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

#include <Helium/smiles.h>
#include <Helium/algorithms/cycles.h>
#include <Eigen/Geometry>

using namespace Helium;

double regular_polygon_radius(int n, int side)
{
  return side / (2.0 * std::sin(deg2rad(180.0 / n)));
}

std::vector<Eigen::Vector2d> generate_ring_coords(int size, int bondLength)
{
  double current = 0.0;
  double step = 360.0 / size;
  double r = regular_polygon_radius(size, bondLength);

  switch (size) {
    case 3:
      current = 90.0;
      break;
    case 4:
      current = 45.0;
      break;
    case 5:
      current = 90.0;
      break;
    case 6:
      current = 90.0;
      break;
    case 7:
      current = 90.0;
      break;
    case 8:
      current = 22.5;
      break;
    default:
      break;
  }

  std::vector<Eigen::Vector2d> coords;
  for (int i = 0; i < size; ++i) {
    double x = r * std::cos(deg2rad(current));
    double y = r * std::sin(deg2rad(current));
    current += step;
    coords.push_back(Eigen::Vector2d(x, y));
  }

  return coords;
}

void translate(std::vector<Eigen::Vector2d> &coords, const Eigen::Vector2d &offset)
{
  for (std::size_t i = 0; i < coords.size(); ++i)
    coords[i] += offset;
}

void rotate(std::vector<Eigen::Vector2d> &coords, double angle)
{
  Eigen::Rotation2Dd m(angle);
  for (std::size_t i = 0; i < coords.size(); ++i)
    coords[i] = m * coords[i];
}

void rotate_around_point(std::vector<Eigen::Vector2d> &coords, const Eigen::Vector2d &p, double angle)
{
  translate(coords, -p);
  rotate(coords, angle);
  translate(coords, p);
}

void rotate_around_point(Eigen::Vector2d &coord, const Eigen::Vector2d &p, double angle)
{
  std::vector<Eigen::Vector2d> coords(1, coord);
  rotate_around_point(coords, p, angle);
  coord = coords[0];
}

double vector_angle(const Eigen::Vector2d &a, const Eigen::Vector2d &b)
{
  return std::atan2(a.y(), a.x()) - std::atan2(b.y(), b.x()) + deg2rad(180.0);
}

double find_arc_radius(double n, double d, double s)
{
  double x = 1.0;
  double r = 1.0;
  while (r < 10 * s) {
    x = 2 * r * std::sin(n * std::asin(s / (2 * r)));
    if (d < x)
      return r;
    r += 1.0;
  }

  return 0.0;
}

double find_arc_angle(double s, double r)
{
  return 2 * std::asin(s / (2 * r));
}

class CyclicPart
{
  public:
    void addAtom(Index index, const Eigen::Vector2d &coord)
    {
      m_indexMap[index] = m_atoms.size();
      m_coords.push_back(coord);
      m_atoms.push_back(index);
    }

    void addBond(Index index)
    {
      m_bonds.push_back(index);
    }

    const std::vector<Index>& atoms() const
    {
      return m_atoms;
    }

    const std::vector<Index>& bonds() const
    {
      return m_bonds;
    }

    bool containsAtom(Index index) const
    {
      return std::find(m_atoms.begin(), m_atoms.end(), index) != m_atoms.end();
    }

    const std::vector<Eigen::Vector2d>& coords() const
    {
      return m_coords;
    }

    std::vector<Eigen::Vector2d>& coords()
    {
      return m_coords;
    }
    
    Eigen::Vector2d& coord(Index index)
    {
      return m_coords[m_indexMap[index]];
    }

    bool operator<(const CyclicPart &other) const
    {
      return m_atoms.size() < other.atoms().size();
    }
    
    bool operator>(const CyclicPart &other) const
    {
      return m_atoms.size() > other.m_atoms.size();
    }

    std::vector<Index> sharedAtoms(const CyclicPart &other) const
    {
      std::vector<Index> result;
      for (std::size_t i = 0; i < m_atoms.size(); ++i)
        if (std::find(other.atoms().begin(), other.atoms().end(), m_atoms[i]) != other.atoms().end())
          result.push_back(m_atoms[i]);
      return result;
    }
    
    std::vector<Index> sharedBonds(const CyclicPart &other) const
    {
      std::vector<Index> result;
      for (std::size_t i = 0; i < m_bonds.size(); ++i)
        if (std::find(other.bonds().begin(), other.bonds().end(), m_bonds[i]) != other.bonds().end())
          result.push_back(m_bonds[i]);
      return result;
    }
    
    template<typename MoleculeType>
    std::vector<std::pair<Index, Index> > connections(const MoleculeType &mol, const CyclicPart &other) const
    {
      std::vector<std::pair<Index, Index> > result;
      FOREACH_BOND_T (bond, mol, MoleculeType) {
        Index source = get_index(mol, get_source(mol, *bond));
        Index target = get_index(mol, get_target(mol, *bond));

        if (containsAtom(source) && other.containsAtom(target))
          result.push_back(std::make_pair(source, target));
        else if (containsAtom(target) && other.containsAtom(source))
          result.push_back(std::make_pair(target, source));
      }

      return result;
    }

    template<typename MoleculeType>
    Eigen::Vector2d attachDir(const MoleculeType &mol, typename molecule_traits<MoleculeType>::atom_type atom)
    {
      assert(std::find(m_atoms.begin(), m_atoms.end(), get_index(mol, atom)) != m_atoms.end());
      Eigen::Vector2d result;
      FOREACH_NBR (nbr, atom, mol, MoleculeType) {
        if (std::find(m_atoms.begin(), m_atoms.end(), get_index(mol, *nbr)) == m_atoms.end())
          continue;
        result += m_coords[m_indexMap[get_index(mol, *nbr)]] - m_coords[m_indexMap[get_index(mol, atom)]];
      }
      result.normalize();
      return result;
    }

    void merge(CyclicPart &other)
    {
      for (std::size_t i = 0; i < other.atoms().size(); ++i) {
        Index index = other.atoms()[i];
        if (std::find(m_atoms.begin(), m_atoms.end(), index) != m_atoms.end())
          continue;
        addAtom(index, other.coord(index));
      }
      for (std::size_t i = 0; i < other.bonds().size(); ++i) {
        Index index = other.bonds()[i];
        if (std::find(m_bonds.begin(), m_bonds.end(), index) != m_bonds.end())
          continue;
        addBond(index);
      }
    }

    /**
     * @brief Attach @p other to this part by a spiro-atom.
     */
    template<typename MoleculeType>
    void attachSpiro(const MoleculeType &mol, const typename molecule_traits<MoleculeType>::atom_type &atom,
        CyclicPart &other)
    {
      // get the attachment directions for this and the other part
      Eigen::Vector2d a = attachDir(mol, atom);
      Eigen::Vector2d b = other.attachDir(mol, atom);

      // compute the angle between the two
      double angle = std::atan2(a.y(), a.x()) - std::atan2(b.y(), b.x()) + deg2rad(180.0);
      // rotate & translate other
      rotate_around_point(other.coords(), other.coord(get_index(mol, atom)), angle);
      translate(other.coords(), coord(get_index(mol, atom)) - other.coord(get_index(mol, atom)));

      // merge other with this part
      merge(other);
    }
    
    /**
     * @brief Attach @p other to this part using a shared bond.
     */
    template<typename MoleculeType>
    void attach(const MoleculeType &mol, const typename molecule_traits<MoleculeType>::atom_type &atom1,
        const typename molecule_traits<MoleculeType>::atom_type &atom2, CyclicPart &other)
    {
      Eigen::Vector2d a = coord(get_index(mol, atom2)) - coord(get_index(mol, atom1));
      Eigen::Vector2d b = other.coord(get_index(mol, atom2)) - other.coord(get_index(mol, atom1));
      
      // compute the angle between the two
      double angle = std::atan2(a.y(), a.x()) - std::atan2(b.y(), b.x());
      // rotate & translate other
      rotate_around_point(other.coords(), other.coord(get_index(mol, atom1)), angle);
      translate(other.coords(), coord(get_index(mol, atom1)) - other.coord(get_index(mol, atom1)));

      // merge other with this part
      merge(other);
    }
    
    /**
     * @brief Connect @p other to this part by a single bond.
     */
    template<typename MoleculeType>
    void connect(const MoleculeType &mol, const typename molecule_traits<MoleculeType>::atom_type &atom1,
        const typename molecule_traits<MoleculeType>::atom_type &atom2, CyclicPart &other, double bondLength)
    {
      // get the attachment directions for this and the other part
      Eigen::Vector2d a = attachDir(mol, atom1);
      Eigen::Vector2d b = other.attachDir(mol, atom2);

      // compute the angle between the two
      double angle = std::atan2(a.y(), a.x()) - std::atan2(b.y(), b.x()) + deg2rad(180.0);
      // rotate & translate other
      rotate_around_point(other.coords(), other.coord(get_index(mol, atom2)), angle);
      translate(other.coords(), coord(get_index(mol, atom1)) - other.coord(get_index(mol, atom2)));
      translate(other.coords(), -a * bondLength);

      // merge other with this part
      merge(other);
    }

  private:
    std::map<Index, int> m_indexMap;
    std::vector<Eigen::Vector2d> m_coords;
    std::vector<Index> m_atoms;
    std::vector<Index> m_bonds;
};

class AcyclicPart
{
  public:
    void addAtom(Index index, const Eigen::Vector2d &coord)
    {
      m_indexMap[index] = m_atoms.size();
      m_coords.push_back(coord);
      m_atoms.push_back(index);
    }

    void addBond(Index index)
    {
      m_bonds.push_back(index);
    }

    const std::vector<Index>& atoms() const
    {
      return m_atoms;
    }

    const std::vector<Index>& bonds() const
    {
      return m_bonds;
    }

    bool containsAtom(Index index) const
    {
      return std::find(m_atoms.begin(), m_atoms.end(), index) != m_atoms.end();
    }

    const std::vector<Eigen::Vector2d>& coords() const
    {
      return m_coords;
    }

    std::vector<Eigen::Vector2d>& coords()
    {
      return m_coords;
    }
    
    Eigen::Vector2d& coord(Index index)
    {
      return m_coords[m_indexMap[index]];
    }

    bool operator<(const CyclicPart &other) const
    {
      return m_atoms.size() < other.atoms().size();
    }
    
    bool operator>(const CyclicPart &other) const
    {
      return m_atoms.size() > other.atoms().size();
    }

  private:
    std::map<Index, int> m_indexMap;
    std::vector<Eigen::Vector2d> m_coords;
    std::vector<Index> m_atoms;
    std::vector<Index> m_bonds;
};


struct Circle
{
  Circle(const Eigen::Vector2d &c_, double r_) : c(c_), r(r_)
  {
  }

  std::pair<Eigen::Vector2d, Eigen::Vector2d> intersection(const Circle &other)
  {
    const Eigen::Vector2d &p0 = c;
    const Eigen::Vector2d &p1 = other.c;
    // distance between p0 and p1
    double d = (p0 - p1).norm();
    // distance between p0 and line between the intersection points
    double a = (r * r - other.r * other.r + d * d) / (2 * d);
    // distance between line p0p1 and intersection points
    double h = std::sqrt(r * r - a * a);
    // intersection of line p0p1 and line between intersection points
    Eigen::Vector2d p2 = (p1 - p0) * (a / d) + p0;
    // intersection points
    double x3 = p2.x() + h * (p1.y() - p0.y()) / d;
    double y3 = p2.y() - h * (p1.x() - p0.x()) / d;
    double x4 = p2.x() - h * (p1.y() - p0.y()) / d;
    double y4 = p2.y() + h * (p1.x() - p0.x()) / d;

    return std::make_pair(Eigen::Vector2d(x3, y3), Eigen::Vector2d(x4, y4));
  }

  Eigen::Vector2d c;
  double r;
};

template<typename MoleculeType>
void place_fragments(const MoleculeType &mol, std::vector<Eigen::Vector2d> &coords)
{
  std::vector<unsigned int> components = connected_atom_components(mol);
  Size n = num_connected_components(mol);

  std::vector<double> minX(n, std::numeric_limits<double>::max());
  std::vector<double> maxX(n, std::numeric_limits<double>::min());
  std::vector<double> minY(n, std::numeric_limits<double>::max());
  std::vector<double> maxY(n, std::numeric_limits<double>::min());

  for (std::size_t i = 0; i < coords.size(); ++i) {
    unsigned int c = components[i];

    if (coords[i].x() < minX[c])
      minX[c] = coords[i].x();
    else if (coords[i].x() > maxX[c])
      maxX[c] = coords[i].x();

    if (coords[i].y() < minY[c])
      minY[c] = coords[i].y();
    else if (coords[i].y() > maxY[c])
      maxY[c] = coords[i].y();
  }

  std::cerr << "minX:" << minX << std::endl;
  std::cerr << "maxX:" << maxX << std::endl;

  double right = 0.0;
  for (Size i = 0; i < n; ++i) {
    std::cerr << "right = " << right << std::endl;
    Eigen::Vector2d offset(right + 30 - minX[i], 30 - minY[i]);
    right += 60 + maxX[i] - minX[i];
    //if (!i)
    //  right += 30;

    std::cerr << "offset: " << offset << std::endl;

    for (std::size_t j = 0; j < coords.size(); ++j)
      if (components[j] == i)
        coords[j] += offset;
  }
}

template<typename MoleculeType>
std::vector<Eigen::Vector2d> generate_coordinates(const MoleculeType &mol)
{
  typedef typename molecule_traits<MoleculeType>::bond_type bond_type;
  typedef typename molecule_traits<MoleculeType>::atom_type atom_type;
  
  double bondLength = 50.0;

  // create vector to hold the resulting coordinates
  std::vector<Eigen::Vector2d> coords(num_atoms(mol), Eigen::Vector2d::Zero());

  //
  // start by processing cyclic parts
  //
  RingSet<MoleculeType> rings = relevant_cycles(mol);

  // create a CyclicPart for each simple cycle
  std::vector<CyclicPart> cyclicParts;
  for (auto &ring : rings.rings()) {
    cyclicParts.resize(cyclicParts.size() + 1);
    CyclicPart &part = cyclicParts.back();
    
    std::vector<Eigen::Vector2d> ringCoords = generate_ring_coords(ring.size(), bondLength);
    for (int i = 0; i < ring.size(); ++i) {
      part.addAtom(get_index(mol, ring.atom(i)), ringCoords[i]);
      part.addBond(get_index(mol, ring.bond(i)));
    }
  }

  // sort the CyclicParts by size
  std::sort(cyclicParts.begin(), cyclicParts.end(), std::greater<CyclicPart>());

  // merge the CyclicParts
  //
  // cyclic parts are joined when the following situations occur:
  // - no shared atoms, connected by a single bond
  // - single shared atom (spiro)
  // - 2 shared atoms with 1 shared bond
  // - multiple shared bonds (construct arc..)
  std::vector<CyclicPart> mergedCyclicParts;
  for (auto &part : cyclicParts) {
    std::cerr << "part size: " << part.atoms().size() << std::endl;
    if (mergedCyclicParts.empty()) {
      mergedCyclicParts.push_back(part);
      continue;
    }

    for (auto &mergedPart : mergedCyclicParts) {
      // determine the shared atoms and bonds
      std::vector<Index> sharedAtoms = mergedPart.sharedAtoms(part);
      std::vector<Index> sharedBonds = mergedPart.sharedBonds(part);

      std::cerr << "mergedPart: " << mergedPart.atoms() << std::endl;
      std::cerr << "part: " << part.atoms() << std::endl;
      std::cerr << "sharedAtoms: " << sharedAtoms << std::endl;

      switch (sharedAtoms.size()) {
        case 0:
          {
            // TODO: 6-4-6 ??
            // the 2 cycles may be connected by a single bond
            std::vector<std::pair<Index, Index> > connections = mergedPart.connections(mol, part);
            if (connections.empty()) {
              // no such connection found
              mergedCyclicParts.push_back(part);
            } else {
              assert(connections.size() == 1);
              // there is such a connection
              mergedPart.connect(mol, get_atom(mol, connections[0].first), 
                  get_atom(mol, connections[0].second), part, bondLength);
            }
          }
          break;
        case 1:
          // spiro-like shared atom
          mergedPart.attachSpiro(mol, get_atom(mol, sharedAtoms[0]), part);
          break;
        case 2:
          // single shared bond
          assert(sharedBonds.size() == 1);
          mergedPart.attach(mol, get_atom(mol, sharedAtoms[0]), 
              get_atom(mol, sharedAtoms[1]), part);
          break;
        default:
          {
            // most difficult case: construct arc
            std::cerr << "shared atoms: " << sharedAtoms << std::endl;
            std::cerr << "shared bonds: " << sharedBonds << std::endl;

            // if there are no unshared atoms, no coordinates need to be determined
            if (sharedAtoms.size() == part.atoms().size())
                break;

            // determine the two outer atoms of the shared path between the new
            // cycle and the already merged cyclic part (these are the points
            // were the arc will attach to the already merged part)
            // these atoms will have degree 1 if only the sharedAtoms are considered
            // while the other atoms will have degree 2
            std::vector<int> degrees(sharedAtoms.size());
            for (std::size_t i = 0; i < sharedBonds.size(); ++i) {
              bond_type bond = get_bond(mol, sharedBonds[i]);
              Index source = get_index(mol, get_source(mol, bond));
              Index target = get_index(mol, get_target(mol, bond));

              std::vector<Index>::iterator iter = std::find(sharedAtoms.begin(), sharedAtoms.end(), source);
              if (iter != sharedAtoms.end())
                degrees[iter - sharedAtoms.begin()]++;
              iter = std::find(sharedAtoms.begin(), sharedAtoms.end(), target);
              if (iter != sharedAtoms.end())
                degrees[iter - sharedAtoms.begin()]++;
            }
            std::vector<Index> atoms;
            for (std::size_t i = 0; i < degrees.size(); ++i)
              if (degrees[i] == 1)
                atoms.push_back(sharedAtoms[i]);
            std::cerr << "atoms: " << atoms << std::endl;
            assert(atoms.size() == 2);

            // determine distance between the 2 atoms
            Eigen::Vector2d p0 = mergedPart.coord(atoms[0]);
            Eigen::Vector2d p1 = mergedPart.coord(atoms[1]);
            double d = (p0 - p1).norm();
            std::cerr << "d = " << d << std::endl;

            // find radius and angle
            int n = part.bonds().size() - sharedBonds.size();
            double r = find_arc_radius(n, d, bondLength);
            double a = find_arc_angle(bondLength, r);

            std::cerr << "r = " << r << std::endl;
            std::cerr << "a = " << a << std::endl;

            // find intersection points = arc center
            Circle c1(p0, r), c2(p1, r);
            std::pair<Eigen::Vector2d, Eigen::Vector2d> points = c1.intersection(c2);
            const Eigen::Vector2d &p2 = points.first;
            const Eigen::Vector2d &p3 = points.second;

            // determine the correct arc center
            Eigen::Vector2d p4 = p0 + bondLength * mergedPart.attachDir(mol, get_atom(mol, atoms[0]));
            Eigen::Vector2d p5 = p1 + bondLength * mergedPart.attachDir(mol, get_atom(mol, atoms[1]));
            double d1 = (p4 - p2).norm() + (p5 - p2).norm();
            double d2 = (p4 - p3).norm() + (p5 - p3).norm();
            bool moreUnshared = (part.bonds().size() - sharedBonds.size()) > sharedBonds.size();
            const Eigen::Vector2d &c = moreUnshared ? (d1 > d2 ? p2 : p3) : (d1 < d2 ? p2 : p3);

            // find a vector between the center c and p1 or p2
            // one of these vectors will be used to determine the arc points
            // by rotating it clockwise multiple times
            Eigen::Vector2d r1 = p0 - c;
            Eigen::Vector2d r2 = p1 - c;
            // copy these vectors to check which one to use
            Eigen::Vector2d v1 = r1;
            Eigen::Vector2d v2 = r2;
            // rotate v1 and v2 clockwise
            rotate_around_point(v1, Eigen::Vector2d::Zero(), a);
            rotate_around_point(v2, Eigen::Vector2d::Zero(), a);
            // determine distance of v1 and v2 to the attach point
            d1 = (v1 - mergedPart.attachDir(mol, get_atom(mol, atoms[0]))).norm();
            d2 = (v2 - mergedPart.attachDir(mol, get_atom(mol, atoms[1]))).norm();
            // choose the final v to rotate based on d1 and d2
            Eigen::Vector2d &v = d1 > d2 ? r1 : r2;

            // still need to determine if unshared atoms are in the same direction
            std::vector<Index> unsharedAtoms = part.atoms();
            while (unsharedAtoms[0] != sharedAtoms[0])
              std::rotate(unsharedAtoms.begin(), unsharedAtoms.begin() + 1, unsharedAtoms.end());
            if (std::equal(sharedAtoms.begin(), sharedAtoms.end(), unsharedAtoms.begin()))
              std::reverse(unsharedAtoms.begin(), unsharedAtoms.end());

            // construct arc
            for (std::size_t i = 0; i < unsharedAtoms.size(); ++i) {
              Index index = unsharedAtoms[i];
              if (std::find(sharedAtoms.begin(), sharedAtoms.end(), index) != sharedAtoms.end())
                continue;
              rotate_around_point(v, Eigen::Vector2d::Zero(), a);
              part.coord(index) = c + v;
              //part.coord(index) = c;
            }

            mergedPart.merge(part);



          }
          break;
      }

    }

  }





  for (auto &part : mergedCyclicParts) {
    for (int i = 0; i < part.atoms().size(); ++i) {
      Index index = part.atoms()[i];
      coords[index] = part.coord(index);
    }
  }

  place_fragments(mol, coords);



  return coords;
}

template<typename MoleculeType>
void write_svg(const MoleculeType &mol, std::vector<Eigen::Vector2d> &coords)
{
  std::cout << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\" viewBox=\"-0 0 500 500\">" << std::endl;

  std::cout << "<circle cx=\"" << coords[0].x() << "\" cy=\"" << coords[0].y() << "\" r=\"5\" fill=\"red\"/>" << std::endl;

  FOREACH_BOND_T (bond, mol, MoleculeType) {
    Index source = get_index(mol, get_source(mol, *bond));
    Index target = get_index(mol, get_target(mol, *bond));

    std::cout << "<line x1=\"" << coords[source].x() << "\" y1=\"" << coords[source].y()
              << "\" x2=\"" << coords[target].x() << "\" y2=\"" << coords[target].y()
              << "\" stroke-width=\"1\" stroke=\"black\"/>" << std::endl;
  }

  std::cout << "</svg>" << std::endl;
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " <SMILES>" << std::endl;
    return -1;
  }

  HeMol mol;
  Smiles SMILES;
  SMILES.read(argv[1], mol);
  if (SMILES.error()) {
    std::cerr << SMILES.error().what() << std::endl;
    return -1;
  }


  auto coords = generate_coordinates(mol);
  write_svg(mol, coords);


}


