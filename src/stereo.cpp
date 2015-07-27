/*
 * Copyright (c) 2013-2015, Tim Vandermeersch
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
#include <Helium/stereo.h>

namespace Helium {

  namespace impl {

    std::string refString(Stereo::Ref ref)
    {
      std::stringstream ss;
      if (ref == Stereo::nullRef())
        ss << "?";
      else if (ref == Stereo::implRef())
        ss << "H";
      else
        ss << ref;
      return ss.str();
    }

  } // namespace impl


  std::ostream& operator<<(std::ostream &os, const StereoStorage &stereo)
  {
    os << "StereoStorage(";
    switch (stereo.type()) {
      case Stereo::Tetrahedral:
        os << "Tetrahedral, center: " << stereo.center() << ", refs:";
        for (std::size_t i = 0; i < 4; ++i)
          os << " " << impl::refString(stereo.ref(i));
        break;
      case Stereo::CisTrans:
        os << "CisTrans, center: " << stereo.center() << ", refs:";
        for (std::size_t i = 0; i < 4; ++i)
          os << " " << impl::refString(stereo.ref(i));
        break;
      case Stereo::SquarePlanar:
        os << "SquarePlanar, center: " << stereo.center() << ", refs:";
        for (std::size_t i = 0; i < 4; ++i)
          os << " " << impl::refString(stereo.ref(i));
        break;
      case Stereo::Allene:
        os << "Allene, center: " << stereo.center() << ", refs:";
        for (std::size_t i = 0; i < 4; ++i)
          os << " " << impl::refString(stereo.ref(i));
        break;
      case Stereo::TrigonalBipyramidal:
        os << "TrigonalBipyramidal, center: " << stereo.center() << ", refs:";
        for (std::size_t i = 0; i < 5; ++i)
          os << " " << impl::refString(stereo.ref(i));
        break;
      case Stereo::Octahedral:
        os << "Octahedral, center: " << stereo.center() << ", refs:";
        for (std::size_t i = 0; i < 6; ++i)
          os << " " << impl::refString(stereo.ref(i));
        break;
      default:
        os << "None";
        break;
    }
    os << ")";
    return os;
  }

  namespace impl {

    unsigned int stereo_num_inversions(const Stereo::Ref *refs, std::size_t size)
    {
      unsigned int result = 0;
      for (std::size_t i = 0; i < size; ++i)
        for (std::size_t j = i + 1; j < size; ++j)
          if (refs[j] < refs[i])
            ++result;
      return result;
    }

    bool stereo_has_even_parity(const Stereo::Ref *refs, std::size_t size)
    {
      return (stereo_num_inversions(refs, size) % 2) == 0;
    }

  } // namespace impl

  Stereo::Class tetrahedral_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4)
  {
    assert(std::find(storage.refs(), storage.refs() + 4, ref1) != storage.refs() + 4);
    assert(std::find(storage.refs(), storage.refs() + 4, ref2) != storage.refs() + 4);
    assert(std::find(storage.refs(), storage.refs() + 4, ref3) != storage.refs() + 4);
    assert(std::find(storage.refs(), storage.refs() + 4, ref4) != storage.refs() + 4);

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


  namespace impl {

    // abcd = bcda = cdab = dabc
    bool stereo_compare_directed_cycles(const Stereo::Ref *refs, Stereo::Ref ref1, Stereo::Ref ref2,
        Stereo::Ref ref3, Stereo::Ref ref4)
    {
      Stereo::Ref ptr[4] = {ref1, ref2, ref3, ref4};

      // reduce problem to abcd = abcd
      while (ptr[0] != refs[0])
        std::rotate(ptr, ptr + 1, ptr + 4);

      if (std::equal(refs + 1, refs + 4, ptr + 1))
        return true;

      return false;
    }

    // abcd = bcda = cdab = dabc = dcba = cbad = badc = adcb
    bool stereo_compare_undirected_cycles(const Stereo::Ref *refs, Stereo::Ref ref1, Stereo::Ref ref2,
        Stereo::Ref ref3, Stereo::Ref ref4)
    {
      Stereo::Ref ptr[4] = {ref1, ref2, ref3, ref4};

      // reduce problem to abcd or adcb
      while (ptr[0] != refs[0])
        std::rotate(ptr, ptr + 1, ptr + 4);

      // handle abcd = abcd
      if (std::equal(refs + 1, refs + 4, ptr + 1))
        return true;

      // handle abcd = adcb
      std::reverse(ptr + 1, ptr + 4);
      if (std::equal(refs + 1, refs + 4, ptr + 1))
        return true;

      return false;
    }

  }

  Stereo::Class squareplanar_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4)
  {
    assert(std::find(storage.refs(), storage.refs() + 4, ref1) != storage.refs() + 4);
    assert(std::find(storage.refs(), storage.refs() + 4, ref2) != storage.refs() + 4);
    assert(std::find(storage.refs(), storage.refs() + 4, ref3) != storage.refs() + 4);
    assert(std::find(storage.refs(), storage.refs() + 4, ref4) != storage.refs() + 4);

    // assume U-shape and check
    if (impl::stereo_compare_undirected_cycles(storage.refs(), ref1, ref2, ref3, ref4))
      return Stereo::SP1;

    // assume 4-shape and check
    if (impl::stereo_compare_undirected_cycles(storage.refs(), ref1, ref3, ref2, ref4))
      return Stereo::SP2;

    // assume Z-shape and check
    assert(impl::stereo_compare_undirected_cycles(storage.refs(), ref2, ref1, ref3, ref4));

    return Stereo::SP3;
  }

  Stereo::Class trigonalbipyramidal_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4, Stereo::Ref ref5)
  {
    assert(std::find(storage.refs(), storage.refs() + 5, ref1) != storage.refs() + 5);
    assert(std::find(storage.refs(), storage.refs() + 5, ref2) != storage.refs() + 5);
    assert(std::find(storage.refs(), storage.refs() + 5, ref3) != storage.refs() + 5);
    assert(std::find(storage.refs(), storage.refs() + 5, ref4) != storage.refs() + 5);
    assert(std::find(storage.refs(), storage.refs() + 5, ref5) != storage.refs() + 5);

    Stereo::Ref refs[5] = {ref1, ref2, ref3, ref4, ref5};

    // determine axis to use
    std::size_t a = std::find(refs, refs + 5, storage.ref(0)) - refs;
    std::size_t e = std::find(refs, refs + 5, storage.ref(4)) - refs;

    // determine plane refs
    std::vector<Stereo::Ref> planeRefs;
    for (std::size_t i = 0; i < 5; ++i) {
      if (i == a || i == e)
        continue;
      planeRefs.push_back(refs[i]);
    }

    if (a > e) {
      std::swap(a, e);
      //std::swap(planeRefs[1], planeRefs[2]);
      std::reverse(planeRefs.begin(), planeRefs.end());
    }

    // reduce to bcd or bdc
    while (storage.ref(1) != planeRefs[0])
      std::rotate(planeRefs.begin(), planeRefs.begin() + 1, planeRefs.end());

    Stereo::Class classes[2];
    if (a == 0 && e == 4) {
      classes[0] = Stereo::TB1;
      classes[1] = Stereo::TB2;
    } else if (a == 0 && e == 3) {
      classes[0] = Stereo::TB3;
      classes[1] = Stereo::TB4;
    } else if (a == 0 && e == 2) {
      classes[0] = Stereo::TB5;
      classes[1] = Stereo::TB6;
    } else if (a == 0 && e == 1) {
      classes[0] = Stereo::TB7;
      classes[1] = Stereo::TB8;
    } else if (a == 1 && e == 4) {
      classes[0] = Stereo::TB9;
      classes[1] = Stereo::TB11;
    } else if (a == 1 && e == 3) {
      classes[0] = Stereo::TB10;
      classes[1] = Stereo::TB12;
    } else if (a == 1 && e == 2) {
      classes[0] = Stereo::TB13;
      classes[1] = Stereo::TB14;
    } else if (a == 2 && e == 4) {
      classes[0] = Stereo::TB15;
      classes[1] = Stereo::TB20;
    } else if (a == 2 && e == 3) {
      classes[0] = Stereo::TB16;
      classes[1] = Stereo::TB19;
    } else if (a == 3 && e == 4) {
      classes[0] = Stereo::TB17;
      classes[1] = Stereo::TB18;
    }

    if (planeRefs[1] == storage.ref(2))
      return classes[0];

    assert(planeRefs[2] == storage.ref(2));
    return classes[1];
  }

  namespace impl {

    int octahedral_permutations[30][24*6] = {
      // OH1
      {0, 1, 2, 3, 4, 5, 0, 2, 3, 4, 1, 5, 0, 3, 4, 1, 2, 5, 0, 4, 1, 2, 3, 5, 1, 0, 4, 5, 2, 3, 1, 2, 0, 4, 5, 3, 1, 4, 5, 2, 0, 3, 1, 5, 2, 0, 4, 3, 2, 0, 1, 5, 3, 4, 2, 1, 5, 3, 0, 4, 2, 3, 0, 1, 5, 4, 2, 5, 3, 0, 1, 4, 3, 0, 2, 5, 4, 1, 3, 2, 5, 4, 0, 1, 3, 4, 0, 2, 5, 1, 3, 5, 4, 0, 2, 1, 4, 0, 3, 5, 1, 2, 4, 1, 0, 3, 5, 2, 4, 3, 5, 1, 0, 2, 4, 5, 1, 0, 3, 2, 5, 1, 4, 3, 2, 0, 5, 2, 1, 4, 3, 0, 5, 3, 2, 1, 4, 0, 5, 4, 3, 2, 1, 0},
      // OH2
      {0, 1, 4, 3, 2, 5, 0, 2, 1, 4, 3, 5, 0, 3, 2, 1, 4, 5, 0, 4, 3, 2, 1, 5, 1, 0, 2, 5, 4, 3, 1, 2, 5, 4, 0, 3, 1, 4, 0, 2, 5, 3, 1, 5, 4, 0, 2, 3, 2, 0, 3, 5, 1, 4, 2, 1, 0, 3, 5, 4, 2, 3, 5, 1, 0, 4, 2, 5, 1, 0, 3, 4, 3, 0, 4, 5, 2, 1, 3, 2, 0, 4, 5, 1, 3, 4, 5, 2, 0, 1, 3, 5, 2, 0, 4, 1, 4, 0, 1, 5, 3, 2, 4, 1, 5, 3, 0, 2, 4, 3, 0, 1, 5, 2, 4, 5, 3, 0, 1, 2, 5, 1, 2, 3, 4, 0, 5, 2, 3, 4, 1, 0, 5, 3, 4, 1, 2, 0, 5, 4, 1, 2, 3, 0},
      // OH3
      {0, 1, 2, 3, 5, 4, 0, 2, 3, 4, 5, 1, 0, 3, 4, 1, 5, 2, 0, 4, 1, 2, 5, 3, 1, 0, 4, 5, 3, 2, 1, 2, 0, 4, 3, 5, 1, 4, 5, 2, 3, 0, 1, 5, 2, 0, 3, 4, 2, 0, 1, 5, 4, 3, 2, 1, 5, 3, 4, 0, 2, 3, 0, 1, 4, 5, 2, 5, 3, 0, 4, 1, 3, 0, 2, 5, 1, 4, 3, 2, 5, 4, 1, 0, 3, 4, 0, 2, 1, 5, 3, 5, 4, 0, 1, 2, 4, 0, 3, 5, 2, 1, 4, 1, 0, 3, 2, 5, 4, 3, 5, 1, 2, 0, 4, 5, 1, 0, 2, 3, 5, 1, 4, 3, 0, 2, 5, 2, 1, 4, 0, 3, 5, 3, 2, 1, 0, 4, 5, 4, 3, 2, 0, 1},
      // OH4
      {0, 1, 2, 4, 3, 5, 0, 2, 3, 1, 4, 5, 0, 3, 4, 2, 1, 5, 0, 4, 1, 3, 2, 5, 1, 0, 4, 2, 5, 3, 1, 2, 0, 5, 4, 3, 1, 4, 5, 0, 2, 3, 1, 5, 2, 4, 0, 3, 2, 0, 1, 3, 5, 4, 2, 1, 5, 0, 3, 4, 2, 3, 0, 5, 1, 4, 2, 5, 3, 1, 0, 4, 3, 0, 2, 4, 5, 1, 3, 2, 5, 0, 4, 1, 3, 4, 0, 5, 2, 1, 3, 5, 4, 2, 0, 1, 4, 0, 3, 1, 5, 2, 4, 1, 0, 5, 3, 2, 4, 3, 5, 0, 1, 2, 4, 5, 1, 3, 0, 2, 5, 1, 4, 2, 3, 0, 5, 2, 1, 3, 4, 0, 5, 3, 2, 4, 1, 0, 5, 4, 3, 1, 2, 0},
      // OH5
      {0, 1, 2, 4, 5, 3, 0, 2, 3, 1, 5, 4, 0, 3, 4, 2, 5, 1, 0, 4, 1, 3, 5, 2, 1, 0, 4, 2, 3, 5, 1, 2, 0, 5, 3, 4, 1, 4, 5, 0, 3, 2, 1, 5, 2, 4, 3, 0, 2, 0, 1, 3, 4, 5, 2, 1, 5, 0, 4, 3, 2, 3, 0, 5, 4, 1, 2, 5, 3, 1, 4, 0, 3, 0, 2, 4, 1, 5, 3, 2, 5, 0, 1, 4, 3, 4, 0, 5, 1, 2, 3, 5, 4, 2, 1, 0, 4, 0, 3, 1, 2, 5, 4, 1, 0, 5, 2, 3, 4, 3, 5, 0, 2, 1, 4, 5, 1, 3, 2, 0, 5, 1, 4, 2, 0, 3, 5, 2, 1, 3, 0, 4, 5, 3, 2, 4, 0, 1, 5, 4, 3, 1, 0, 2},
      // OH6
      {0, 1, 2, 5, 3, 4, 0, 2, 3, 5, 4, 1, 0, 3, 4, 5, 1, 2, 0, 4, 1, 5, 2, 3, 1, 0, 4, 3, 5, 2, 1, 2, 0, 3, 4, 5, 1, 4, 5, 3, 2, 0, 1, 5, 2, 3, 0, 4, 2, 0, 1, 4, 5, 3, 2, 1, 5, 4, 3, 0, 2, 3, 0, 4, 1, 5, 2, 5, 3, 4, 0, 1, 3, 0, 2, 1, 5, 4, 3, 2, 5, 1, 4, 0, 3, 4, 0, 1, 2, 5, 3, 5, 4, 1, 0, 2, 4, 0, 3, 2, 5, 1, 4, 1, 0, 2, 3, 5, 4, 3, 5, 2, 1, 0, 4, 5, 1, 2, 0, 3, 5, 1, 4, 0, 3, 2, 5, 2, 1, 0, 4, 3, 5, 3, 2, 0, 1, 4, 5, 4, 3, 0, 2, 1},
      // OH7
      {0, 1, 2, 5, 4, 3, 0, 2, 3, 5, 1, 4, 0, 3, 4, 5, 2, 1, 0, 4, 1, 5, 3, 2, 1, 0, 4, 3, 2, 5, 1, 2, 0, 3, 5, 4, 1, 4, 5, 3, 0, 2, 1, 5, 2, 3, 4, 0, 2, 0, 1, 4, 3, 5, 2, 1, 5, 4, 0, 3, 2, 3, 0, 4, 5, 1, 2, 5, 3, 4, 1, 0, 3, 0, 2, 1, 4, 5, 3, 2, 5, 1, 0, 4, 3, 4, 0, 1, 5, 2, 3, 5, 4, 1, 2, 0, 4, 0, 3, 2, 1, 5, 4, 1, 0, 2, 5, 3, 4, 3, 5, 2, 0, 1, 4, 5, 1, 2, 3, 0, 5, 1, 4, 0, 2, 3, 5, 2, 1, 0, 3, 4, 5, 3, 2, 0, 4, 1, 5, 4, 3, 0, 1, 2},
      // OH8
      {0, 1, 3, 2, 4, 5, 0, 2, 4, 3, 1, 5, 0, 3, 1, 4, 2, 5, 0, 4, 2, 1, 3, 5, 1, 0, 5, 4, 2, 3, 1, 2, 4, 0, 5, 3, 1, 4, 2, 5, 0, 3, 1, 5, 0, 2, 4, 3, 2, 0, 5, 1, 3, 4, 2, 1, 3, 5, 0, 4, 2, 3, 1, 0, 5, 4, 2, 5, 0, 3, 1, 4, 3, 0, 5, 2, 4, 1, 3, 2, 4, 5, 0, 1, 3, 4, 2, 0, 5, 1, 3, 5, 0, 4, 2, 1, 4, 0, 5, 3, 1, 2, 4, 1, 3, 0, 5, 2, 4, 3, 1, 5, 0, 2, 4, 5, 0, 1, 3, 2, 5, 1, 3, 4, 2, 0, 5, 2, 4, 1, 3, 0, 5, 3, 1, 2, 4, 0, 5, 4, 2, 3, 1, 0},
      // OH9
      {0, 1, 3, 2, 5, 4, 0, 2, 4, 3, 5, 1, 0, 3, 1, 4, 5, 2, 0, 4, 2, 1, 5, 3, 1, 0, 5, 4, 3, 2, 1, 2, 4, 0, 3, 5, 1, 4, 2, 5, 3, 0, 1, 5, 0, 2, 3, 4, 2, 0, 5, 1, 4, 3, 2, 1, 3, 5, 4, 0, 2, 3, 1, 0, 4, 5, 2, 5, 0, 3, 4, 1, 3, 0, 5, 2, 1, 4, 3, 2, 4, 5, 1, 0, 3, 4, 2, 0, 1, 5, 3, 5, 0, 4, 1, 2, 4, 0, 5, 3, 2, 1, 4, 1, 3, 0, 2, 5, 4, 3, 1, 5, 2, 0, 4, 5, 0, 1, 2, 3, 5, 1, 3, 4, 0, 2, 5, 2, 4, 1, 0, 3, 5, 3, 1, 2, 0, 4, 5, 4, 2, 3, 0, 1},
      // OH10
      {0, 1, 3, 4, 2, 5, 0, 2, 4, 1, 3, 5, 0, 3, 1, 2, 4, 5, 0, 4, 2, 3, 1, 5, 1, 0, 5, 2, 4, 3, 1, 2, 4, 5, 0, 3, 1, 4, 2, 0, 5, 3, 1, 5, 0, 4, 2, 3, 2, 0, 5, 3, 1, 4, 2, 1, 3, 0, 5, 4, 2, 3, 1, 5, 0, 4, 2, 5, 0, 1, 3, 4, 3, 0, 5, 4, 2, 1, 3, 2, 4, 0, 5, 1, 3, 4, 2, 5, 0, 1, 3, 5, 0, 2, 4, 1, 4, 0, 5, 1, 3, 2, 4, 1, 3, 5, 0, 2, 4, 3, 1, 0, 5, 2, 4, 5, 0, 3, 1, 2, 5, 1, 3, 2, 4, 0, 5, 2, 4, 3, 1, 0, 5, 3, 1, 4, 2, 0, 5, 4, 2, 1, 3, 0},
      // OH11
      {0, 1, 3, 4, 5, 2, 0, 2, 4, 1, 5, 3, 0, 3, 1, 2, 5, 4, 0, 4, 2, 3, 5, 1, 1, 0, 5, 2, 3, 4, 1, 2, 4, 5, 3, 0, 1, 4, 2, 0, 3, 5, 1, 5, 0, 4, 3, 2, 2, 0, 5, 3, 4, 1, 2, 1, 3, 0, 4, 5, 2, 3, 1, 5, 4, 0, 2, 5, 0, 1, 4, 3, 3, 0, 5, 4, 1, 2, 3, 2, 4, 0, 1, 5, 3, 4, 2, 5, 1, 0, 3, 5, 0, 2, 1, 4, 4, 0, 5, 1, 2, 3, 4, 1, 3, 5, 2, 0, 4, 3, 1, 0, 2, 5, 4, 5, 0, 3, 2, 1, 5, 1, 3, 2, 0, 4, 5, 2, 4, 3, 0, 1, 5, 3, 1, 4, 0, 2, 5, 4, 2, 1, 0, 3},
      // OH12
      {0, 1, 3, 5, 2, 4, 0, 2, 4, 5, 3, 1, 0, 3, 1, 5, 4, 2, 0, 4, 2, 5, 1, 3, 1, 0, 5, 3, 4, 2, 1, 2, 4, 3, 0, 5, 1, 4, 2, 3, 5, 0, 1, 5, 0, 3, 2, 4, 2, 0, 5, 4, 1, 3, 2, 1, 3, 4, 5, 0, 2, 3, 1, 4, 0, 5, 2, 5, 0, 4, 3, 1, 3, 0, 5, 1, 2, 4, 3, 2, 4, 1, 5, 0, 3, 4, 2, 1, 0, 5, 3, 5, 0, 1, 4, 2, 4, 0, 5, 2, 3, 1, 4, 1, 3, 2, 0, 5, 4, 3, 1, 2, 5, 0, 4, 5, 0, 2, 1, 3, 5, 1, 3, 0, 4, 2, 5, 2, 4, 0, 1, 3, 5, 3, 1, 0, 2, 4, 5, 4, 2, 0, 3, 1},
      // OH13
      {0, 1, 3, 5, 4, 2, 0, 2, 4, 5, 1, 3, 0, 3, 1, 5, 2, 4, 0, 4, 2, 5, 3, 1, 1, 0, 5, 3, 2, 4, 1, 2, 4, 3, 5, 0, 1, 4, 2, 3, 0, 5, 1, 5, 0, 3, 4, 2, 2, 0, 5, 4, 3, 1, 2, 1, 3, 4, 0, 5, 2, 3, 1, 4, 5, 0, 2, 5, 0, 4, 1, 3, 3, 0, 5, 1, 4, 2, 3, 2, 4, 1, 0, 5, 3, 4, 2, 1, 5, 0, 3, 5, 0, 1, 2, 4, 4, 0, 5, 2, 1, 3, 4, 1, 3, 2, 5, 0, 4, 3, 1, 2, 0, 5, 4, 5, 0, 2, 3, 1, 5, 1, 3, 0, 2, 4, 5, 2, 4, 0, 3, 1, 5, 3, 1, 0, 4, 2, 5, 4, 2, 0, 1, 3},
      // OH14
      {0, 1, 4, 2, 3, 5, 0, 2, 1, 3, 4, 5, 0, 3, 2, 4, 1, 5, 0, 4, 3, 1, 2, 5, 1, 0, 2, 4, 5, 3, 1, 2, 5, 0, 4, 3, 1, 4, 0, 5, 2, 3, 1, 5, 4, 2, 0, 3, 2, 0, 3, 1, 5, 4, 2, 1, 0, 5, 3, 4, 2, 3, 5, 0, 1, 4, 2, 5, 1, 3, 0, 4, 3, 0, 4, 2, 5, 1, 3, 2, 0, 5, 4, 1, 3, 4, 5, 0, 2, 1, 3, 5, 2, 4, 0, 1, 4, 0, 1, 3, 5, 2, 4, 1, 5, 0, 3, 2, 4, 3, 0, 5, 1, 2, 4, 5, 3, 1, 0, 2, 5, 1, 2, 4, 3, 0, 5, 2, 3, 1, 4, 0, 5, 3, 4, 2, 1, 0, 5, 4, 1, 3, 2, 0},
      // OH15
      {0, 1, 4, 2, 5, 3, 0, 2, 1, 3, 5, 4, 0, 3, 2, 4, 5, 1, 0, 4, 3, 1, 5, 2, 1, 0, 2, 4, 3, 5, 1, 2, 5, 0, 3, 4, 1, 4, 0, 5, 3, 2, 1, 5, 4, 2, 3, 0, 2, 0, 3, 1, 4, 5, 2, 1, 0, 5, 4, 3, 2, 3, 5, 0, 4, 1, 2, 5, 1, 3, 4, 0, 3, 0, 4, 2, 1, 5, 3, 2, 0, 5, 1, 4, 3, 4, 5, 0, 1, 2, 3, 5, 2, 4, 1, 0, 4, 0, 1, 3, 2, 5, 4, 1, 5, 0, 2, 3, 4, 3, 0, 5, 2, 1, 4, 5, 3, 1, 2, 0, 5, 1, 2, 4, 0, 3, 5, 2, 3, 1, 0, 4, 5, 3, 4, 2, 0, 1, 5, 4, 1, 3, 0, 2},
      // OH16
      {0, 1, 4, 3, 5, 2, 0, 2, 1, 4, 5, 3, 0, 3, 2, 1, 5, 4, 0, 4, 3, 2, 5, 1, 1, 0, 2, 5, 3, 4, 1, 2, 5, 4, 3, 0, 1, 4, 0, 2, 3, 5, 1, 5, 4, 0, 3, 2, 2, 0, 3, 5, 4, 1, 2, 1, 0, 3, 4, 5, 2, 3, 5, 1, 4, 0, 2, 5, 1, 0, 4, 3, 3, 0, 4, 5, 1, 2, 3, 2, 0, 4, 1, 5, 3, 4, 5, 2, 1, 0, 3, 5, 2, 0, 1, 4, 4, 0, 1, 5, 2, 3, 4, 1, 5, 3, 2, 0, 4, 3, 0, 1, 2, 5, 4, 5, 3, 0, 2, 1, 5, 1, 2, 3, 0, 4, 5, 2, 3, 4, 0, 1, 5, 3, 4, 1, 0, 2, 5, 4, 1, 2, 0, 3},
      // OH17
      {0, 1, 4, 5, 2, 3, 0, 2, 1, 5, 3, 4, 0, 3, 2, 5, 4, 1, 0, 4, 3, 5, 1, 2, 1, 0, 2, 3, 4, 5, 1, 2, 5, 3, 0, 4, 1, 4, 0, 3, 5, 2, 1, 5, 4, 3, 2, 0, 2, 0, 3, 4, 1, 5, 2, 1, 0, 4, 5, 3, 2, 3, 5, 4, 0, 1, 2, 5, 1, 4, 3, 0, 3, 0, 4, 1, 2, 5, 3, 2, 0, 1, 5, 4, 3, 4, 5, 1, 0, 2, 3, 5, 2, 1, 4, 0, 4, 0, 1, 2, 3, 5, 4, 1, 5, 2, 0, 3, 4, 3, 0, 2, 5, 1, 4, 5, 3, 2, 1, 0, 5, 1, 2, 0, 4, 3, 5, 2, 3, 0, 1, 4, 5, 3, 4, 0, 2, 1, 5, 4, 1, 0, 3, 2},
      // OH18
      {0, 1, 4, 5, 3, 2, 0, 2, 1, 5, 4, 3, 0, 3, 2, 5, 1, 4, 0, 4, 3, 5, 2, 1, 1, 0, 2, 3, 5, 4, 1, 2, 5, 3, 4, 0, 1, 4, 0, 3, 2, 5, 1, 5, 4, 3, 0, 2, 2, 0, 3, 4, 5, 1, 2, 1, 0, 4, 3, 5, 2, 3, 5, 4, 1, 0, 2, 5, 1, 4, 0, 3, 3, 0, 4, 1, 5, 2, 3, 2, 0, 1, 4, 5, 3, 4, 5, 1, 2, 0, 3, 5, 2, 1, 0, 4, 4, 0, 1, 2, 5, 3, 4, 1, 5, 2, 3, 0, 4, 3, 0, 2, 1, 5, 4, 5, 3, 2, 0, 1, 5, 1, 2, 0, 3, 4, 5, 2, 3, 0, 4, 1, 5, 3, 4, 0, 1, 2, 5, 4, 1, 0, 2, 3},
      // OH19
      {0, 1, 5, 2, 3, 4, 0, 2, 5, 3, 4, 1, 0, 3, 5, 4, 1, 2, 0, 4, 5, 1, 2, 3, 1, 0, 3, 4, 5, 2, 1, 2, 3, 0, 4, 5, 1, 4, 3, 5, 2, 0, 1, 5, 3, 2, 0, 4, 2, 0, 4, 1, 5, 3, 2, 1, 4, 5, 3, 0, 2, 3, 4, 0, 1, 5, 2, 5, 4, 3, 0, 1, 3, 0, 1, 2, 5, 4, 3, 2, 1, 5, 4, 0, 3, 4, 1, 0, 2, 5, 3, 5, 1, 4, 0, 2, 4, 0, 2, 3, 5, 1, 4, 1, 2, 0, 3, 5, 4, 3, 2, 5, 1, 0, 4, 5, 2, 1, 0, 3, 5, 1, 0, 4, 3, 2, 5, 2, 0, 1, 4, 3, 5, 3, 0, 2, 1, 4, 5, 4, 0, 3, 2, 1},
      // OH20
      {0, 1, 5, 2, 4, 3, 0, 2, 5, 3, 1, 4, 0, 3, 5, 4, 2, 1, 0, 4, 5, 1, 3, 2, 1, 0, 3, 4, 2, 5, 1, 2, 3, 0, 5, 4, 1, 4, 3, 5, 0, 2, 1, 5, 3, 2, 4, 0, 2, 0, 4, 1, 3, 5, 2, 1, 4, 5, 0, 3, 2, 3, 4, 0, 5, 1, 2, 5, 4, 3, 1, 0, 3, 0, 1, 2, 4, 5, 3, 2, 1, 5, 0, 4, 3, 4, 1, 0, 5, 2, 3, 5, 1, 4, 2, 0, 4, 0, 2, 3, 1, 5, 4, 1, 2, 0, 5, 3, 4, 3, 2, 5, 0, 1, 4, 5, 2, 1, 3, 0, 5, 1, 0, 4, 2, 3, 5, 2, 0, 1, 3, 4, 5, 3, 0, 2, 4, 1, 5, 4, 0, 3, 1, 2},
      // OH21
      {0, 1, 5, 3, 2, 4, 0, 2, 5, 4, 3, 1, 0, 3, 5, 1, 4, 2, 0, 4, 5, 2, 1, 3, 1, 0, 3, 5, 4, 2, 1, 2, 3, 4, 0, 5, 1, 4, 3, 2, 5, 0, 1, 5, 3, 0, 2, 4, 2, 0, 4, 5, 1, 3, 2, 1, 4, 3, 5, 0, 2, 3, 4, 1, 0, 5, 2, 5, 4, 0, 3, 1, 3, 0, 1, 5, 2, 4, 3, 2, 1, 4, 5, 0, 3, 4, 1, 2, 0, 5, 3, 5, 1, 0, 4, 2, 4, 0, 2, 5, 3, 1, 4, 1, 2, 3, 0, 5, 4, 3, 2, 1, 5, 0, 4, 5, 2, 0, 1, 3, 5, 1, 0, 3, 4, 2, 5, 2, 0, 4, 1, 3, 5, 3, 0, 1, 2, 4, 5, 4, 0, 2, 3, 1},
      // OH22
      {0, 1, 5, 3, 4, 2, 0, 2, 5, 4, 1, 3, 0, 3, 5, 1, 2, 4, 0, 4, 5, 2, 3, 1, 1, 0, 3, 5, 2, 4, 1, 2, 3, 4, 5, 0, 1, 4, 3, 2, 0, 5, 1, 5, 3, 0, 4, 2, 2, 0, 4, 5, 3, 1, 2, 1, 4, 3, 0, 5, 2, 3, 4, 1, 5, 0, 2, 5, 4, 0, 1, 3, 3, 0, 1, 5, 4, 2, 3, 2, 1, 4, 0, 5, 3, 4, 1, 2, 5, 0, 3, 5, 1, 0, 2, 4, 4, 0, 2, 5, 1, 3, 4, 1, 2, 3, 5, 0, 4, 3, 2, 1, 0, 5, 4, 5, 2, 0, 3, 1, 5, 1, 0, 3, 2, 4, 5, 2, 0, 4, 3, 1, 5, 3, 0, 1, 4, 2, 5, 4, 0, 2, 1, 3},
      // OH23
      {0, 1, 5, 4, 2, 3, 0, 2, 5, 1, 3, 4, 0, 3, 5, 2, 4, 1, 0, 4, 5, 3, 1, 2, 1, 0, 3, 2, 4, 5, 1, 2, 3, 5, 0, 4, 1, 4, 3, 0, 5, 2, 1, 5, 3, 4, 2, 0, 2, 0, 4, 3, 1, 5, 2, 1, 4, 0, 5, 3, 2, 3, 4, 5, 0, 1, 2, 5, 4, 1, 3, 0, 3, 0, 1, 4, 2, 5, 3, 2, 1, 0, 5, 4, 3, 4, 1, 5, 0, 2, 3, 5, 1, 2, 4, 0, 4, 0, 2, 1, 3, 5, 4, 1, 2, 5, 0, 3, 4, 3, 2, 0, 5, 1, 4, 5, 2, 3, 1, 0, 5, 1, 0, 2, 4, 3, 5, 2, 0, 3, 1, 4, 5, 3, 0, 4, 2, 1, 5, 4, 0, 1, 3, 2},
      // OH24
      {0, 1, 5, 4, 3, 2, 0, 2, 5, 1, 4, 3, 0, 3, 5, 2, 1, 4, 0, 4, 5, 3, 2, 1, 1, 0, 3, 2, 5, 4, 1, 2, 3, 5, 4, 0, 1, 4, 3, 0, 2, 5, 1, 5, 3, 4, 0, 2, 2, 0, 4, 3, 5, 1, 2, 1, 4, 0, 3, 5, 2, 3, 4, 5, 1, 0, 2, 5, 4, 1, 0, 3, 3, 0, 1, 4, 5, 2, 3, 2, 1, 0, 4, 5, 3, 4, 1, 5, 2, 0, 3, 5, 1, 2, 0, 4, 4, 0, 2, 1, 5, 3, 4, 1, 2, 5, 3, 0, 4, 3, 2, 0, 1, 5, 4, 5, 2, 3, 0, 1, 5, 1, 0, 2, 3, 4, 5, 2, 0, 3, 4, 1, 5, 3, 0, 4, 1, 2, 5, 4, 0, 1, 2, 3},
      // OH25
      {0, 5, 1, 2, 3, 4, 0, 5, 2, 3, 4, 1, 0, 5, 3, 4, 1, 2, 0, 5, 4, 1, 2, 3, 1, 3, 0, 4, 5, 2, 1, 3, 2, 0, 4, 5, 1, 3, 4, 5, 2, 0, 1, 3, 5, 2, 0, 4, 2, 4, 0, 1, 5, 3, 2, 4, 1, 5, 3, 0, 2, 4, 3, 0, 1, 5, 2, 4, 5, 3, 0, 1, 3, 1, 0, 2, 5, 4, 3, 1, 2, 5, 4, 0, 3, 1, 4, 0, 2, 5, 3, 1, 5, 4, 0, 2, 4, 2, 0, 3, 5, 1, 4, 2, 1, 0, 3, 5, 4, 2, 3, 5, 1, 0, 4, 2, 5, 1, 0, 3, 5, 0, 1, 4, 3, 2, 5, 0, 2, 1, 4, 3, 5, 0, 3, 2, 1, 4, 5, 0, 4, 3, 2, 1},
      // OH26
      {0, 5, 1, 2, 4, 3, 0, 5, 2, 3, 1, 4, 0, 5, 3, 4, 2, 1, 0, 5, 4, 1, 3, 2, 1, 3, 0, 4, 2, 5, 1, 3, 2, 0, 5, 4, 1, 3, 4, 5, 0, 2, 1, 3, 5, 2, 4, 0, 2, 4, 0, 1, 3, 5, 2, 4, 1, 5, 0, 3, 2, 4, 3, 0, 5, 1, 2, 4, 5, 3, 1, 0, 3, 1, 0, 2, 4, 5, 3, 1, 2, 5, 0, 4, 3, 1, 4, 0, 5, 2, 3, 1, 5, 4, 2, 0, 4, 2, 0, 3, 1, 5, 4, 2, 1, 0, 5, 3, 4, 2, 3, 5, 0, 1, 4, 2, 5, 1, 3, 0, 5, 0, 1, 4, 2, 3, 5, 0, 2, 1, 3, 4, 5, 0, 3, 2, 4, 1, 5, 0, 4, 3, 1, 2},
      // OH27
      {0, 5, 1, 3, 2, 4, 0, 5, 2, 4, 3, 1, 0, 5, 3, 1, 4, 2, 0, 5, 4, 2, 1, 3, 1, 3, 0, 5, 4, 2, 1, 3, 2, 4, 0, 5, 1, 3, 4, 2, 5, 0, 1, 3, 5, 0, 2, 4, 2, 4, 0, 5, 1, 3, 2, 4, 1, 3, 5, 0, 2, 4, 3, 1, 0, 5, 2, 4, 5, 0, 3, 1, 3, 1, 0, 5, 2, 4, 3, 1, 2, 4, 5, 0, 3, 1, 4, 2, 0, 5, 3, 1, 5, 0, 4, 2, 4, 2, 0, 5, 3, 1, 4, 2, 1, 3, 0, 5, 4, 2, 3, 1, 5, 0, 4, 2, 5, 0, 1, 3, 5, 0, 1, 3, 4, 2, 5, 0, 2, 4, 1, 3, 5, 0, 3, 1, 2, 4, 5, 0, 4, 2, 3, 1},
      // OH28
      {0, 5, 1, 3, 4, 2, 0, 5, 2, 4, 1, 3, 0, 5, 3, 1, 2, 4, 0, 5, 4, 2, 3, 1, 1, 3, 0, 5, 2, 4, 1, 3, 2, 4, 5, 0, 1, 3, 4, 2, 0, 5, 1, 3, 5, 0, 4, 2, 2, 4, 0, 5, 3, 1, 2, 4, 1, 3, 0, 5, 2, 4, 3, 1, 5, 0, 2, 4, 5, 0, 1, 3, 3, 1, 0, 5, 4, 2, 3, 1, 2, 4, 0, 5, 3, 1, 4, 2, 5, 0, 3, 1, 5, 0, 2, 4, 4, 2, 0, 5, 1, 3, 4, 2, 1, 3, 5, 0, 4, 2, 3, 1, 0, 5, 4, 2, 5, 0, 3, 1, 5, 0, 1, 3, 2, 4, 5, 0, 2, 4, 3, 1, 5, 0, 3, 1, 4, 2, 5, 0, 4, 2, 1, 3},
      // OH29
      {0, 5, 1, 4, 2, 3, 0, 5, 2, 1, 3, 4, 0, 5, 3, 2, 4, 1, 0, 5, 4, 3, 1, 2, 1, 3, 0, 2, 4, 5, 1, 3, 2, 5, 0, 4, 1, 3, 4, 0, 5, 2, 1, 3, 5, 4, 2, 0, 2, 4, 0, 3, 1, 5, 2, 4, 1, 0, 5, 3, 2, 4, 3, 5, 0, 1, 2, 4, 5, 1, 3, 0, 3, 1, 0, 4, 2, 5, 3, 1, 2, 0, 5, 4, 3, 1, 4, 5, 0, 2, 3, 1, 5, 2, 4, 0, 4, 2, 0, 1, 3, 5, 4, 2, 1, 5, 0, 3, 4, 2, 3, 0, 5, 1, 4, 2, 5, 3, 1, 0, 5, 0, 1, 2, 4, 3, 5, 0, 2, 3, 1, 4, 5, 0, 3, 4, 2, 1, 5, 0, 4, 1, 3, 2},
      // OH30
      {0, 5, 1, 4, 3, 2, 0, 5, 2, 1, 4, 3, 0, 5, 3, 2, 1, 4, 0, 5, 4, 3, 2, 1, 1, 3, 0, 2, 5, 4, 1, 3, 2, 5, 4, 0, 1, 3, 4, 0, 2, 5, 1, 3, 5, 4, 0, 2, 2, 4, 0, 3, 5, 1, 2, 4, 1, 0, 3, 5, 2, 4, 3, 5, 1, 0, 2, 4, 5, 1, 0, 3, 3, 1, 0, 4, 5, 2, 3, 1, 2, 0, 4, 5, 3, 1, 4, 5, 2, 0, 3, 1, 5, 2, 0, 4, 4, 2, 0, 1, 5, 3, 4, 2, 1, 5, 3, 0, 4, 2, 3, 0, 1, 5, 4, 2, 5, 3, 0, 1, 5, 0, 1, 2, 3, 4, 5, 0, 2, 3, 4, 1, 5, 0, 3, 4, 1, 2, 5, 0, 4, 1, 2, 3}
    };

  }

  Stereo::Class octahedral_class(const StereoStorage &storage, Stereo::Ref ref1,
      Stereo::Ref ref2, Stereo::Ref ref3, Stereo::Ref ref4, Stereo::Ref ref5, Stereo::Ref ref6)
  {
    assert(std::find(storage.refs(), storage.refs() + 6, ref1) != storage.refs() + 6);
    assert(std::find(storage.refs(), storage.refs() + 6, ref2) != storage.refs() + 6);
    assert(std::find(storage.refs(), storage.refs() + 6, ref3) != storage.refs() + 6);
    assert(std::find(storage.refs(), storage.refs() + 6, ref4) != storage.refs() + 6);
    assert(std::find(storage.refs(), storage.refs() + 6, ref5) != storage.refs() + 6);
    assert(std::find(storage.refs(), storage.refs() + 6, ref6) != storage.refs() + 6);

    Stereo::Ref refs[6] = {ref1, ref2, ref3, ref4, ref5, ref6};

    for (int oh = 0; oh < 30; ++oh) {
      int *perms = impl::octahedral_permutations[oh];
      for (int i = 0; i < 24; ++i) {
        int off = i * 6;
        int a = perms[off + 0];
        int b = perms[off + 1];
        int c = perms[off + 2];
        int d = perms[off + 3];
        int e = perms[off + 4];
        int f = perms[off + 5];

        Stereo::Ref refs2[6] = {storage.ref(a), storage.ref(b), storage.ref(c),
          storage.ref(d), storage.ref(e), storage.ref(f)};

        if (std::equal(refs2, refs2 + 6, refs))
          return static_cast<Stereo::Class>(Stereo::OH1 + oh);
      }
    }

    return Stereo::OH1;
  }


}
