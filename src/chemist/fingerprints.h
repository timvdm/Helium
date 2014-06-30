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
#ifndef HELIUM_CHEMIST_FINGERPRINTS_H
#define HELIUM_CHEMIST_FINGERPRINTS_H

#include <Helium/bitvec.h>
#include <Helium/fingerprints/similarity.h>
#include <Helium/fingerprints/fingerprints.h>

namespace Helium {

  namespace Chemist {

    using Helium::bitvec_num_words_for_bits;
    using Helium::bitvec_zero;
    using Helium::bitvec_copy;
    using Helium::bitvec_get;
    using Helium::bitvec_set;
    using Helium::bitvec_reset;
    using Helium::bitvec_is_subset_superset;
    using Helium::bitvec_count;
    using Helium::bitvec_union_count;
    using Helium::bitvec_tanimoto;
    using Helium::bitvec_cosine;
    using Helium::bitvec_hamming;
    using Helium::bitvec_russell_rao;
    using Helium::bitvec_forbes;
    using Helium::bitvec_write_size;
    using Helium::bitvec_read_size;
    using Helium::bitvec_write;
    using Helium::bitvec_read;
    using Helium::bitvec_to_binary;
    using Helium::bitvec_from_binary;
    using Helium::bitvec_print;
    using Helium::bitvec_to_hex;
    using Helium::bitvec_from_hex;

    using Helium::SimilaritySearchIndex;

  }

}

#endif
