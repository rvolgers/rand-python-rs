mod mersenne_twister;

use std::slice;

pub use mersenne_twister::MersenneTwister;

type BigInt = u64; // could use `num` module to make this generic so you can use bigints etc


pub struct PythonRandom {
    mersenne_twister: MersenneTwister,
}


impl PythonRandom {
    pub fn new(mt: MersenneTwister) -> PythonRandom {
        PythonRandom {
            mersenne_twister: mt,
        }
    }

    pub fn random(&mut self) -> f64 {
        self.mersenne_twister.genrand_res53()
    }

    pub fn getrandbits(&mut self, k: u32) -> BigInt {
        assert!(0 < k && k <= BigInt::from(0u8).leading_zeros());

        if k <= 32 {
            return BigInt::from(self.mersenne_twister.genrand_int32() >> (32 - k));
        }

        let mut tmp = BigInt::from(0u8);
        let mut k = k;
        let mut shift = 0;
        while k > 0 {
            tmp |= BigInt::from(self.mersenne_twister.genrand_int32() >> 32u32.saturating_sub(k)) << shift;
            k = k.saturating_sub(32);
            shift += 32;
        }

        tmp
    }

    pub fn randbelow(&mut self, n: BigInt) -> BigInt {
        let n_bits = BigInt::from(0u8).leading_zeros() - n.leading_zeros();
        let mut r = self.getrandbits(n_bits);
        while r >= n {
            r = self.getrandbits(n_bits);
        }
        r
    }

    pub fn seed_u32(&mut self, s: u32) {
        self.mersenne_twister.init_by_array(slice::from_ref(&s));
    }

    pub fn expovariate(&mut self, lambda: f64) -> f64 {
        -(1.0 - self.mersenne_twister.genrand_res53()).ln() / lambda
    }

    pub fn shuffle<T>(&mut self, x: &mut [T]) {
        for i in (1..x.len()).rev() {
            let j = self.randbelow(BigInt::from(i as u64) + 1) as usize;
            x.swap(i, j);
        }
    }

    pub fn randint(&mut self, start: BigInt, stop: BigInt) -> BigInt {
        self.randrange(start, stop + 1)
    }

    pub fn randrange(&mut self, start: BigInt, stop: BigInt) -> BigInt {
        start + self.randbelow(stop - start)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sanity_checks() {
        // Known-good values generated using the following python program:
        //
        // import random
        // random.seed(63245986)
        // print(random.getstate())
        // print(random.random())
        // print(random.randrange(0, 100000))
        // print(random.expovariate(1.0 / 15000.0))
        //
        // tmp = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        // random.shuffle(tmp)
        //
        // print(tmp)

        let mt = MersenneTwister::new();

        let mut rand = PythonRandom::new(mt);

        rand.seed_u32(63245986);

        // println!("{:?}", &rand);

        assert_eq!(rand.random(), 0.5213761361171212);

        assert_eq!(rand.randbelow(100000u64), 58671);

        assert_eq!(rand.expovariate(1.0 / 15000.0), 13775.46713470634);

        let mut list: [u64; 10] = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];

        rand.shuffle(&mut list[..]);

        assert_eq!(&list, &[10, 3, 6, 1, 8, 5, 7, 4, 2, 9]);
    }
}
