//! MT19937 and the distributions used by the reference GNU C++ implementation.
//!
//! Keeping the reference generator makes debugging stochastic equivalence much
//! easier. Search results may still diverge after a floating point tie.

#[derive(Clone)]
pub struct Rng {
    state: [u32; 624],
    index: usize,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        let mut state = [0; 624];
        state[0] = seed as u32;
        for i in 1..624 {
            state[i] = 1812433253_u32
                .wrapping_mul(state[i - 1] ^ (state[i - 1] >> 30))
                .wrapping_add(i as u32);
        }
        Self { state, index: 624 }
    }

    pub fn next_u32(&mut self) -> u32 {
        if self.index == 624 {
            for i in 0..624 {
                let x = (self.state[i] & 0x8000_0000) | (self.state[(i + 1) % 624] & 0x7fff_ffff);
                self.state[i] = self.state[(i + 397) % 624]
                    ^ (x >> 1)
                    ^ if x & 1 != 0 { 0x9908_b0df } else { 0 };
            }
            self.index = 0;
        }
        let mut x = self.state[self.index];
        self.index += 1;
        x ^= x >> 11;
        x ^= (x << 7) & 0x9d2c_5680;
        x ^= (x << 15) & 0xefc6_0000;
        x ^ (x >> 18)
    }

    /// Uniform in [0, 1), matching libstdc++'s generate_canonical<double, 53>.
    pub fn random01(&mut self) -> f64 {
        let low = self.next_u32() as f64;
        let high = self.next_u32() as f64;
        ((low + high * 4294967296.0) / 18446744073709551616.0)
            .min(f64::from_bits(1.0_f64.to_bits() - 1))
    }

    pub fn uniform(&mut self, min: f64, max: f64) -> f64 {
        self.random01() * (max - min) + min
    }

    /// Uniform integer in [0, bound), using libstdc++'s unbiased reduction.
    pub fn usize(&mut self, bound: usize) -> usize {
        assert!(
            bound > 0 && bound <= u32::MAX as usize,
            "invalid random integer bound"
        );
        let bound = bound as u32;
        let mut product = self.next_u32() as u64 * bound as u64;
        if (product as u32) < bound {
            let threshold = bound.wrapping_neg() % bound;
            while (product as u32) < threshold {
                product = self.next_u32() as u64 * bound as u64;
            }
        }
        (product >> 32) as usize
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_mt19937_reference_vector() {
        let mut rng = Rng::new(5489);
        let expected = [
            3499211612, 581869302, 3890346734, 3586334585, 545404204, 4161255391, 3922919429,
            949333985, 2715962298, 1323567403,
        ];
        for x in expected {
            assert_eq!(rng.next_u32(), x);
        }
    }

    #[test]
    fn matches_gnu_cpp_mixed_distribution_vectors() {
        // std::mt19937(0), with interleaved real and integer distributions.
        let expected = [
            (0x3fe2f89545f18fe1, [0, 1, 4, 535, 1170127712]),
            (0x3fdb1d290276395f, [0, 1, 2, 273, 638950699]),
            (0x3fad097bb1c89679, [0, 0, 2, 298, 1700216562]),
            (0x3fe0ecb50af9fd49, [0, 1, 2, 577, 1795465483]),
            (0x3fd597e612048bdb, [0, 1, 0, 229, 1788037496]),
        ];
        let mut rng = Rng::new(0);
        for (bits, integers) in expected {
            assert_eq!(rng.random01().to_bits(), bits);
            for (bound, integer) in [1, 2, 7, 624, 2147483647].into_iter().zip(integers) {
                assert_eq!(rng.usize(bound), integer);
            }
        }
    }

    #[test]
    fn ranges_and_reseed_are_deterministic() {
        let mut a = Rng::new(19);
        let mut b = Rng::new(19);
        for _ in 0..2000 {
            let x = a.random01();
            assert!((0.0..1.0).contains(&x));
            assert_eq!(x, b.random01());
            for bound in [1, 2, 7, 624, 2_147_483_649] {
                let x = a.usize(bound);
                assert!(x < bound);
                assert_eq!(x, b.usize(bound));
            }
        }
    }
}
