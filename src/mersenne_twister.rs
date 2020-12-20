/*
   Rust port of the original mt19937ar Mersenne Twister.

   Porting, additional comments and documentation by Ronald Volgers.

   The original copyright notice and BSD 3-clause license is reproduced below.

   ---

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

use std::cmp::max;
use std::fmt;

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908b0df;
const UPPER_MASK: u32 = 0x80000000;
const LOWER_MASK: u32 = 0x7fffffff;

/// `self.index` is set to this value after seeding. It is also encountered during normal
/// operation, when a block of random numbers has been consumed and `genrand_int32` needs
/// to generate a new one before returning the next random number.
const INITIAL_INDEX: usize = N;

/// `self.index` is set to this value at creation time. It indicates the random number
/// generator has not been seeded. If this value is encountered when trying to generate
/// a random number, the implementation will panic.
const UNSEEDED_INDEX: usize = N + 1;

/// Implementation of [Mersenne Twister] (`mt19937ar`).
///
/// This is a direct Rust port of the C implementation known as [`mt19937ar.c`].
///
/// It tries to stay as close as possible to the original in terms of function naming
/// and overall interface. This makes it easier to port libraries that use the Mersenne
/// Twister internally.
///
/// The code remains functionally the same but has been cleaned up and re-commented
/// to reflect my own understanding of the original code and algorithm. Some comments from
/// the original have been omitted. Errors are mine.
///
/// [Mersenne Twister]: http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
/// [`mt19937ar.c`]: http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/CODES/mt19937ar.c
#[derive(Clone)]
pub struct MersenneTwister {
    index: usize,
    state: [u32; N],
}

impl fmt::Debug for MersenneTwister {
    fn fmt(&self, fmt: &mut fmt::Formatter) -> fmt::Result {
        fmt.debug_struct("MersenneTwister")
            .field("index", &self.index)
            .field("state", &&self.state[..])
            .finish()
    }
}

impl MersenneTwister {
    /// Creates a new, unseeded random number generator.
    ///
    /// Seed it using [`init_by_array`], [`init_genrand`], or the deprecated and buggy
    /// [`sgenrand`].
    ///
    /// Trying to generate random numbers without seeding will cause a panic.
    /// The upstream implementation automatically seeds the generator using a hardcoded
    /// seed instead. To use that seed, you can manually call `init_genrand(5489)`.
    ///
    /// [`init_by_array`]: #method.init_by_array
    /// [`init_genrand`]: #method.init_genrand
    /// [`sgenrand`]: #method.sgenrand
    pub fn new() -> MersenneTwister {
        MersenneTwister {
            index: UNSEEDED_INDEX,
            state: [0; N],
        }
    }

    /// This is the core random number generation function that all the others are
    /// based on.
    ///
    /// It generates a random `u32` value.
    pub fn genrand_int32(&mut self) -> u32 {

        let mt = &mut self.state;

        if self.index >= N {

            // This deviates from the original implementation, which will simply call
            // init_genrand(5489) in this situation instead.
            assert!(self.index != UNSEEDED_INDEX,
                "random number generator must be seeded before use");

            #[cfg(feature = "prefetch")]
            #[cfg(target_feature = "sse")]
            unsafe {
                // 64 is min possible cache line size
                for i in (0..N).step_by(64) {
                    core::arch::x86_64::_mm_prefetch((mt as *const _ as *const i8).add(i), core::arch::x86_64::_MM_HINT_T2);
                }
            }

            // Produce an entire block of outputs at once.
            // The calculation for each entry is the same, but the loop is split into
            // three separate parts according to which indexing calculations will wrap,
            // to avoid needing modulo operations or conditionals in the loop itself.

            // Process initial items, where kk+M and kk+1 both don't wrap
            for kk in 0..(N - M) {
                let y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M] ^ (y >> 1) ^ ((y & 1) * MATRIX_A);
            }

            // Process next items, where kk+M wraps but kk+1 doesn't
            for kk in (N - M)..(N - 1) {
                let y = (mt[kk] & UPPER_MASK) | (mt[kk + 1] & LOWER_MASK);
                mt[kk] = mt[kk + M - N] ^ (y >> 1) ^ ((y & 1) * MATRIX_A);
            }

            // Process final item, where both kk+M and kk+1 wrap
            let kk = N - 1;
            let y = (mt[kk] & UPPER_MASK) | (mt[kk + 1 - N] & LOWER_MASK);
            mt[kk] = mt[kk + M - N] ^ (y >> 1) ^ ((y & 1) * MATRIX_A);

            self.index = 0;
        }

        let y = mt[self.index];
        self.index += 1;

        // This is the "tempering" operation to improve the distribution of bits in
        // the output. This transformation is *not* intended to prevent recovery of
        // the state from the output; in fact this transformation is fully invertible
        // so it does nothing to prevent that.
        let mut y = y;
        y ^= y >> 11;
        y ^= (y << 7) & 0x9d2c5680;
        y ^= (y << 15) & 0xefc60000;
        y ^= y >> 18;
        y
    }

    /// Reset the random number generator with a seed of arbitary length.
    ///
    /// The seed is specified as a slice of `u32` values. It cannot be specified as an
    /// iterator because if the seed is shorter than 624 values each value may
    /// be accessed more than once.
    pub fn init_by_array(&mut self, init_key: &[u32]) {

        // Reset self.state to some decently varied constant data.
        // This also sets self.index = N, which ensures that genrand_int32 will run
        // its random number generation loop when it is called for the first time
        // after seeding.
        self.init_genrand(19650218);

        let mt = &mut self.state;

        // Cycles through 1..N, which has length N-1, not N.
        // The zero index is "skipped over" and handled specially.
        let mut i: usize = 1;

        // Cycles through the valid indices of init_key.
        let mut j: usize = 0;

        for _ in 0..max(N, init_key.len()) {
            let prev = mt[i - 1] ^ (mt[i - 1] >> 30);
            mt[i] = (mt[i] ^ prev.wrapping_mul(1664525))
                .wrapping_add(init_key[j])
                .wrapping_add(j as u32);

            i += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }

            j += 1;
            if j >= init_key.len() {
                j = 0;
            }
        }

        for _ in 0..(N - 1) {
            let prev = mt[i - 1] ^ (mt[i - 1] >> 30);
            mt[i] = (mt[i] ^ prev.wrapping_mul(1566083941))
                .wrapping_sub(i as u32);

            i += 1;
            if i >= N {
                mt[0] = mt[N - 1];
                i = 1;
            }
        }

        mt[0] = 0x80000000;
    }

    /// Reset the random number generator using a single 32 bit seed.
    ///
    /// This is generally not recommended since it means all future output
    /// of the random number generator can be reproduced by knowing, guessing,
    /// or accidentally picking the same 32 bit seed value.
    ///
    /// If you want to seed the random number generator with more than 32 bits
    /// of data, see [`init_by_array`].
    ///
    /// [`init_by_array`]: #method.init_by_array
    pub fn init_genrand(&mut self, s: u32) {

        let mt = &mut self.state;

        mt[0] = s;

        for i in 1..N {
            let prev = mt[i - 1] ^ (mt[i - 1] >> 30);
            mt[i] = prev.wrapping_mul(1812433253)
                .wrapping_add(i as u32);
        }

        self.index = INITIAL_INDEX;  // == N
    }

    /// Generates a random `value: f64` such that `0. <= value && value < 1.`,
    /// using 53 bits of randomness (the maximum possible for a `f64` value).
    pub fn genrand_res53(&mut self) -> f64 {

        let a = self.genrand_int32() >> 5;  // 32 - 5 = 27 bits
        let b = self.genrand_int32() >> 6;  // 32 - 6 = 26 bits

        // This is essentially doing ((a << 26) + b) / 2**53, but combining
        // a and b is done in floating point because on 32 bits systems the
        // CPU often has native support for f64 but not u64.

        // Combining a and b in a u64 instead of a f64 is about 25% faster in
        // my testing, but changing anything about the generation of floating
        // point numbers is extremely risky since even the smallest change
        // can cause a number to be rounded up instead of down or vice versa
        // in downstream calculations. So, let's just do it the way the original
        // implementation does.

        // Likewise, the multiplication-by-reciprocal instead of division is
        // how it's done in the original code.

        // 2**26 == 67108864
        // 2**53 == 9007199254740992
        (a as f64 * 67108864.0 + b as f64) * (1.0 / 9007199254740992.0)
    }

    /// Generates a random `value: i32` such that `0 <= value && value <= 0x7fffffff`.
    ///
    /// If you want a `u32`, see [`genrand_int32`].
    ///
    /// [`genrand_int32`]: #method.genrand_int32
    pub fn genrand_int31(&mut self) -> i32 {
        (self.genrand_int32() >> 1) as i32
    }

    /// Generates a random `value: f64` such that `0. <= value && value <= 1.`,
    /// using 32 bits worth of randomness.
    ///
    /// If you want `value < 1` instead of `value <= 1`, see [`genrand_real2`]
    /// (or [`genrand_res53`] if you also want the maximum amount of randomness).
    ///
    /// [`genrand_real2`]: #method.genrand_real2
    /// [`genrand_res53`]: #method.genrand_res53
    pub fn genrand_real1(&mut self) -> f64 {
        // 4294967295 == 2**32 - 1
        self.genrand_int32() as f64 / 4294967295.0
    }

    /// Generates a random `value: f64` such that `0. <= value && value < 1.`,
    /// using 32 bits of randomness.
    ///
    /// If you want the maximum amount of randomness, see [`genrand_res53`].
    ///
    /// [`genrand_res53`]: #method.genrand_res53
    pub fn genrand_real2(&mut self) -> f64 {
        // 4294967296 == 2**32
        self.genrand_int32() as f64 / 4294967296.0
    }

    /// Generates a random `value: f64` such that `0. < value && value < 1.`,
    /// using 32 bits of randomness.
    ///
    /// If you want `0 <= value` instead of `0 < value`, see [`genrand_real2`]
    /// (or [`genrand_res53`] if you also want the maximum amount of randomness).
    ///
    /// [`genrand_real2`]: #method.genrand_real2
    /// [`genrand_res53`]: #method.genrand_res53
    pub fn genrand_real3(&mut self) -> f64 {
        // 4294967296 == 2**32
        (self.genrand_int32() as f64 + 0.5) / 4294967296.0
    }

    /// Reset the random number generator using a single 32 bit seed.
    ///
    /// This old seed function can cause very poor quality random numbers and
    /// is mainly included for historical purposes. It is only present in older
    /// versions of the Mersenne Twister code base, not the `mt19937ar` version
    /// this module is based on.
    ///
    /// Most software uses [`init_by_array`] or [`init_genrand`] instead, which were
    /// introduced in the `mt19937ar` version.
    ///
    /// [`init_by_array`]: #method.init_by_array
    /// [`init_genrand`]: #method.init_getrand
    #[deprecated = "Can lead to very poor quality random numbers. \
                    Use `init_by_array` or `init_genrand` instead."]
    pub fn sgenrand(&mut self, seed: u32) {

        let mt = &mut self.state;

        mt[0] = seed;

        for i in 1..N {
            mt[i] = mt[i - 1].wrapping_mul(69069);
        }

    }

}