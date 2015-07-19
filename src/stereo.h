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
#ifndef HELIUM_STEREOUNIT_H
#define HELIUM_STEREOUNIT_H

#include <Helium/molecule.h>
#include <Helium/ring.h>

#include <vector>
#include <limits>
#include <cassert>

namespace Helium {

  struct Stereo
  {
    typedef unsigned int Ref;

    enum Type {
      None,
      TetraPlanar = 1, //!< 4 planar ligands (i.e. square planar, cis/trans)
      TetraNonPlanar = 2, //!< 4 non-planar ligands (i.e. tetrahedral, allene-like)
      TrigonalBipyramidal = 4, //!< 5 ligands, two axial with the remaning on a plane
      Octahedral = 8, //!< 6 ligands, 3 axes
      CisTrans = 1 + 16,
      SquarePlanar = 1 + 32,
      Tetrahedral = 2 + 64,
      Allene = 2 + 128,
    };

    enum Class {
      TH1,
      TH2
    };

    enum View
    {
      ViewFrom,
      ViewTowards
    };

    enum Winding
    {
      Clockwise,
      AntiClockwise
    };

    static Ref nullRef()
    {
      return std::numeric_limits<Ref>::max();
    }

    static Ref implRef()
    {
      return std::numeric_limits<Ref>::max() - 1;
    }
  };

  /**
   *
   * Tetrahedral
   *   1 center, 4 refs: a b c d, viewing from ref a, refs b c d are clockwise
   *
   *        b  c
   *        | /
   *   a -- c1
   *         \
   *          d
   *
   * CisTrans
   *   2 centers, 4 refs: a b c d
   *
   *   a          d
   *    \        /
   *     c1 == c2
   *    /        \
   *   b          c
   *
   * SquarePlanar
   *   1 center, 4 refs: a b c d
   *
   *   a     d
   *    \   /
   *     c1
   *    /  \
   *   b    c
   *
   *
   *
   */
  class StereoStorage
  {
    public:
      StereoStorage(enum Stereo::Type type = Stereo::None) : m_type(type)
      {
        std::fill(m_centers, m_centers + 2, Stereo::nullRef());
        std::fill(m_refs, m_refs + 6, Stereo::nullRef());
      }

      ~StereoStorage()
      {
      }

      bool isValid() const
      {
        // check type
        if (m_type == Stereo::None)
          return false;

        // check center(s)
        switch (m_type) {
          case Stereo::CisTrans:
          case Stereo::Allene:
            if (m_centers[0] == Stereo::nullRef() || m_centers[1] == Stereo::nullRef())
              return false;
            if (m_centers[0] == Stereo::implRef() || m_centers[1] == Stereo::implRef())
              return false;
            break;
          default:
            if (m_centers[0] == Stereo::nullRef() || m_centers[0] == Stereo::implRef())
              return false;
        }

        // check refs
        int nullCount = std::count(m_refs, m_refs + numRefs(), Stereo::nullRef());
        if (nullCount)
          return false;

        int implCount = std::count(m_refs, m_refs + numRefs(), Stereo::implRef());
        if (implCount > 1)
          return false;

        return true;
      }

      enum Stereo::Type type() const
      {
        return m_type;
      }

      int numCenters() const
      {
        switch (m_type) {
          case Stereo::CisTrans:
          case Stereo::Allene:
            return 2;
          default:
            return 1;
        }
      }

      int numRefs() const
      {
        switch (m_type) {
          case Stereo::TrigonalBipyramidal:
            return 5;
          case Stereo::Octahedral:
            return 6;
          default:
            return 4;
        }
      }

      Stereo::Ref center(int index = 0) const
      {
        assert(index < numCenters());
        return m_centers[index];
      }

      Stereo::Ref ref(int index) const
      {
        assert(index < numRefs());
        return m_refs[index];
      }

      const Stereo::Ref* refs() const
      {
        return m_refs;
      }

      template<typename RefIter>
      static StereoStorage create(enum Stereo::Type type, Stereo::Ref center,
          RefIter beginRefs, RefIter endRefs)
      {
        StereoStorage stereo(type);
        stereo.m_centers[0] = center;
        std::copy(beginRefs, endRefs, stereo.m_refs);
        return stereo;
      }

    private:
      enum Stereo::Type m_type;
      Stereo::Ref m_centers[2]; //!< The stereo center(s)
      Stereo::Ref m_refs[6]; //!< The 4 stereo
  };

  inline std::ostream& operator<<(std::ostream &os, const StereoStorage &stereo)
  {
    os << "StereoStorage(";
    switch (stereo.type()) {
      case Stereo::Tetrahedral:
        os << "Tetrahedral, center: center: " << stereo.center() << ", refs:";
        for (std::size_t i = 0; i < 4; ++i)
          os << " " << stereo.ref(i);
        break;
      case Stereo::CisTrans:
        os << "CisTrans, center: " << stereo.center() << ", refs: ";
        break;
      case Stereo::SquarePlanar:
        os << "SquarePlanar, center: " << stereo.center() << ", refs: ";
        break;
      case Stereo::Allene:
        os << "Allene, center: " << stereo.center() << ", refs: ";
        break;
      case Stereo::TrigonalBipyramidal:
        os << "TrigonalBipyramidal, center: " << stereo.center() << ", refs: ";
        break;
      case Stereo::Octahedral:
        os << "Octahedral, center: " << stereo.center() << ", refs: ";
        break;
      default:
        os << "None";
        break;
    }
    os << ")";
    return os;
  }

  class Stereochemistry
  {
    public:
      std::size_t numStereo() const
      {
        return m_stereo.size();
      }

      std::size_t numTetrahedral() const
      {
        return m_tetrahedral.size();
      }

      bool isTetrahedral(Stereo::Ref center) const
      {
        for (auto stereo : m_tetrahedral)
          if (stereo->center() == center)
            return true;
        return false;
      }

      const StereoStorage& tetrahedral(Stereo::Ref center) const
      {
        for (auto stereo : m_tetrahedral)
          if (stereo->center() == center)
            return *stereo;
        return m_invalid;
      }

      void clear()
      {
        m_stereo.clear();
        m_tetrahedral.clear();
        m_cistrans.clear();
      }

      void add(const StereoStorage &stereo)
      {
        m_stereo.push_back(stereo);
        switch (stereo.type()) {
          case Stereo::Tetrahedral:
            m_tetrahedral.push_back(&m_stereo.back());
            break;
          default:
            break;
        }
      }


    private:
      StereoStorage m_invalid;
      std::vector<StereoStorage> m_stereo;
      std::vector<StereoStorage*> m_tetrahedral;
      std::vector<StereoStorage*> m_cistrans;
  };

  namespace impl {

    inline unsigned int stereo_num_inversions(const Stereo::Ref *refs, std::size_t size)
    {
      unsigned int result = 0;
      for (std::size_t i = 0; i < size; ++i)
        for (std::size_t j = i + 1; j < size; ++j)
          if (refs[j] < refs[i])
            ++result;
      return result;
    }

    inline bool stereo_has_even_parity(const Stereo::Ref *refs, std::size_t size)
    {
      return (stereo_num_inversions(refs, size) % 2) == 0;
    }

  } // namespace impl

  inline Stereo::Class tetrahedral_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4)
  {
    assert(std::find(strorage.refs(), storage.refs() + 4, ref1) != storage.refs() + 4);
    assert(std::find(strorage.refs(), storage.refs() + 4, ref2) != storage.refs() + 4);
    assert(std::find(strorage.refs(), storage.refs() + 4, ref3) != storage.refs() + 4);
    assert(std::find(strorage.refs(), storage.refs() + 4, ref4) != storage.refs() + 4);

    Stereo::Ref refs[4] = {ref1, ref2, ref3, ref4};
    auto pos = std::find(refs, refs + 4, storage.ref(0)) - refs;
    std::swap(refs[0], refs[pos]);

    bool invert = (pos != 0);
    bool parity1 = impl::stereo_has_even_parity(storage.refs() + 1, 3);
    bool parity2 = impl::stereo_has_even_parity(refs + 1, 3);

    if (parity1 != parity2)
      invert = !invert;
    return invert ? Stereo::TH2 : Stereo::TH1;
  }

}

#endif
