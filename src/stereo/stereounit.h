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
#include <cassert>

namespace Helium {

  class Stereo
  {
    public:
      enum Type {
        TetraPlanar, //!< 4 planar ligands (i.e. square planar, cis/trans)
        TetraNonPlanar, //!< 4 non-planar ligands (i.e. tetrahedral, allene-like)
        TrigonalBipyramidal, //!< 5 ligands, two axial with the remaning on a plane
        Octahedral, //!< 6 ligands, 3 axes
        Tetrahedral = TetraNonPlanar,
        AlleneLike = TetraNonPlanar,
        CisTrans = TetraPlanar,
        SquarePlanar = TetraPlanar
      };

      virtual ~Stereo()
      {
      }

      typedef unsigned int Ref;

      static Ref nullRef()
      {
        return std::numeric_limits<Ref>::max();
      }

      static Ref implRef()
      {
        return std::numeric_limits<Ref>::max() - 1;
      }
  };

  class StereoBase : public Stereo
  {
    public:

      StereoBase(enum Type type) : m_type(type)
      {
      }

      virtual ~StereoBase()
      {
      }

      enum Type type() const
      {
        return m_type;
      }

    private:
      enum Type m_type;
  };

  class TetraPlanarStereo : public StereoBase
  {
    public:
      TetraPlanarStereo() : StereoBase(Stereo::TetraPlanar)
      {
        m_refs[0] = Stereo::nullRef();
        m_refs[1] = Stereo::nullRef();
        m_refs[2] = Stereo::nullRef();
        m_refs[3] = Stereo::nullRef();
      }

      const std::vector<Ref>& centers() const
      {
        return m_centers;
      }

      std::vector<Ref>& centers()
      {
        return m_centers;
      }

      Ref ref(int index) const
      {
        return m_refs[index];
      }

      void setRef(int index, Ref value)
      {
        m_refs[index] = value;
      }

      void addRef(Ref value)
      {
        int i = 0;
        for (; i < 4; ++i)
          if (m_refs[i] != nullRef())
            break;
        assert(i < 4);
        m_refs[i] = value;
      }

    private:
      std::vector<Ref> m_centers; //!< The stereo center(s)
      Ref m_refs[4]; //!< The 4 stereo
  };


}

#endif
