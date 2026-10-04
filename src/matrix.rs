//! Dense row-major matrices with allocation-free product updates.
use crate::complex::Complex;
use std::ops::{Index, IndexMut};

#[derive(Clone, Debug, PartialEq)]
pub struct Matrix {
    pub n: usize,
    pub data: Vec<Complex>,
}
impl Matrix {
    pub fn zero(n: usize) -> Self {
        Self {
            n,
            data: vec![Complex::ZERO; n * n],
        }
    }
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zero(n);
        for i in 0..n {
            m[(i, i)] = Complex::ONE;
        }
        m
    }
    pub fn from_data(n: usize, data: Vec<Complex>) -> Result<Self, String> {
        if data.len() != n.checked_mul(n).ok_or("Matrix dimension overflow")? {
            return Err("Matrix data has the wrong length".into());
        }
        Ok(Self { n, data })
    }
    #[inline]
    pub fn mul(&self, rhs: &Self) -> Self {
        let mut out = Self::zero(self.n);
        self.mul_into(rhs, &mut out);
        out
    }
    /// The inner loop walks contiguous rows. Exact zero skips make elementary
    /// quantum gates inexpensive without discarding small floating point terms.
    pub fn mul_into(&self, rhs: &Self, out: &mut Self) {
        assert_eq!(self.n, rhs.n, "Matrix dimension mismatch");
        assert_eq!(self.n, out.n, "Output matrix dimension mismatch");
        let length = self
            .n
            .checked_mul(self.n)
            .expect("Matrix dimension overflow");
        assert_eq!(self.data.len(), length, "Invalid left matrix storage");
        assert_eq!(rhs.data.len(), length, "Invalid right matrix storage");
        assert_eq!(out.data.len(), length, "Invalid output matrix storage");
        if length == 0 {
            return;
        }
        #[cfg(target_arch = "x86_64")]
        if self.n >= 4 && std::arch::is_x86_feature_detected!("avx") {
            // SAFETY: runtime feature detection guards AVX instructions; all
            // slices have equal n*n lengths and Complex is repr(C), two f64s.
            unsafe {
                mul_avx(self.n, &self.data, &rhs.data, &mut out.data);
            }
            return;
        }
        self.mul_into_scalar(rhs, out);
    }
    fn mul_into_scalar(&self, rhs: &Self, out: &mut Self) {
        let n = self.n;
        out.data.fill(Complex::ZERO);
        for (left, row) in self.data.chunks_exact(n).zip(out.data.chunks_exact_mut(n)) {
            for (k, &a) in left.iter().enumerate() {
                if a == Complex::ZERO {
                    continue;
                }
                let right = &rhs.data[k * n..(k + 1) * n];
                // Real factors include H, CX, and all intermediate permutation
                // matrices, and need half as many arithmetic operations.
                if a.im == 0.0 {
                    for (dst, b) in row.iter_mut().zip(right) {
                        dst.re += a.re * b.re;
                        dst.im += a.re * b.im;
                    }
                } else if a.re == 0.0 {
                    for (dst, b) in row.iter_mut().zip(right) {
                        dst.re -= a.im * b.im;
                        dst.im += a.im * b.re;
                    }
                } else {
                    for (dst, b) in row.iter_mut().zip(right) {
                        dst.re += a.re * b.re - a.im * b.im;
                        dst.im += a.re * b.im + a.im * b.re;
                    }
                }
            }
        }
    }
    pub fn adjoint(&self) -> Self {
        let mut out = Self::zero(self.n);
        for i in 0..self.n {
            for j in 0..self.n {
                out[(j, i)] = self[(i, j)].conj();
            }
        }
        out
    }
    pub fn conjugate(&self) -> Self {
        Self {
            n: self.n,
            data: self.data.iter().map(|z| z.conj()).collect(),
        }
    }
    pub fn transpose(&self) -> Self {
        let mut out = Self::zero(self.n);
        for i in 0..self.n {
            for j in 0..self.n {
                out[(j, i)] = self[(i, j)];
            }
        }
        out
    }
    /// Tensor an identity to the left, preserving little-endian qubit indices.
    pub fn kron_identity(&self, n_qubits: usize) -> Self {
        let n = 1usize
            .checked_shl(n_qubits as u32)
            .expect("Too many qubits");
        assert!(self.n.is_power_of_two() && self.n <= n);
        let mut out = Self::zero(n);
        for block in 0..n / self.n {
            for i in 0..self.n {
                let row = (block * self.n + i) * n + block * self.n;
                out.data[row..row + self.n]
                    .copy_from_slice(&self.data[i * self.n..(i + 1) * self.n]);
            }
        }
        out
    }
    /// `order[q]` is the new position of old qubit q, matching Utils::changeQubits.
    pub fn permuted(&self, order: &[usize]) -> Self {
        assert_eq!(self.n, 1usize << order.len());
        let map = bit_permutation(order);
        let mut out = Self::zero(self.n);
        for i in 0..self.n {
            for j in 0..self.n {
                out[(map[i], map[j])] = self[(i, j)];
            }
        }
        out
    }
    pub fn norm_diff(&self, other: &Self) -> f64 {
        assert_eq!(self.n, other.n);
        self.data
            .iter()
            .zip(&other.data)
            .map(|(&a, &b)| (a - b).norm_sqr())
            .sum::<f64>()
            .sqrt()
    }
    pub fn max_abs_diff(&self, other: &Self) -> f64 {
        assert_eq!(self.n, other.n);
        self.data
            .iter()
            .zip(&other.data)
            .map(|(&a, &b)| (a - b).abs())
            .fold(0.0, f64::max)
    }
    pub fn approximately_equal(&self, other: &Self, epsilon: f64) -> bool {
        self.n == other.n
            && self
                .data
                .iter()
                .zip(&other.data)
                .all(|(&a, &b)| (a - b).abs() <= epsilon)
    }
    pub fn trace_conjugate_product(&self, other: &Self) -> Complex {
        assert_eq!(self.n, other.n);
        self.data
            .iter()
            .zip(&other.data)
            .fold(Complex::ZERO, |sum, (&a, &b)| sum + a.conj() * b)
    }
}
impl Index<(usize, usize)> for Matrix {
    type Output = Complex;
    #[inline]
    fn index(&self, (i, j): (usize, usize)) -> &Complex {
        &self.data[i * self.n + j]
    }
}
impl IndexMut<(usize, usize)> for Matrix {
    #[inline]
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut Complex {
        &mut self.data[i * self.n + j]
    }
}

pub fn bit_permutation(order: &[usize]) -> Vec<usize> {
    let mut seen = vec![false; order.len()];
    for &q in order {
        assert!(q < order.len() && !seen[q], "Not a qubit permutation");
        seen[q] = true;
    }
    (0..1usize << order.len())
        .map(|i| {
            order
                .iter()
                .enumerate()
                .fold(0, |mapped, (q, &to)| mapped | (((i >> q) & 1) << to))
        })
        .collect()
}
/// Lexicographic permutations, including one empty permutation.
pub fn permutations(n: usize) -> Vec<Vec<usize>> {
    let mut p: Vec<usize> = (0..n).collect();
    let mut result = vec![p.clone()];
    while next_permutation(&mut p) {
        result.push(p.clone());
    }
    result
}
pub fn next_permutation(p: &mut [usize]) -> bool {
    let Some(i) = (1..p.len()).rev().find(|&i| p[i - 1] < p[i]) else {
        return false;
    };
    let j = (i..p.len()).rev().find(|&j| p[i - 1] < p[j]).unwrap();
    p.swap(i - 1, j);
    p[i..].reverse();
    true
}

// A runtime-dispatched SIMD kernel keeps binaries portable while accelerating
// complex rows on AVX machines. No FMA is used: products and additions retain
// the same rounding boundaries as scalar complex arithmetic.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx")]
unsafe fn mul_avx(n: usize, a: &[Complex], b: &[Complex], out: &mut [Complex]) {
    use std::arch::x86_64::*;
    out.fill(Complex::ZERO);
    let bp = b.as_ptr().cast::<f64>();
    let op = out.as_mut_ptr().cast::<f64>();
    for i in 0..n {
        let dst = op.add(2 * i * n);
        for k in 0..n {
            let coefficient = a[i * n + k];
            if coefficient == Complex::ZERO {
                continue;
            }
            let source = bp.add(2 * k * n);
            let re = _mm256_set1_pd(coefficient.re);
            let im = _mm256_set1_pd(coefficient.im);
            let mut j = 0;
            if coefficient.im == 0.0 {
                while j + 8 <= 2 * n {
                    let b0 = _mm256_loadu_pd(source.add(j));
                    let b1 = _mm256_loadu_pd(source.add(j + 4));
                    let c0 = _mm256_loadu_pd(dst.add(j));
                    let c1 = _mm256_loadu_pd(dst.add(j + 4));
                    _mm256_storeu_pd(dst.add(j), _mm256_add_pd(c0, _mm256_mul_pd(re, b0)));
                    _mm256_storeu_pd(dst.add(j + 4), _mm256_add_pd(c1, _mm256_mul_pd(re, b1)));
                    j += 8;
                }
            } else {
                while j + 8 <= 2 * n {
                    let b0 = _mm256_loadu_pd(source.add(j));
                    let b1 = _mm256_loadu_pd(source.add(j + 4));
                    let p0 = _mm256_addsub_pd(
                        _mm256_mul_pd(re, b0),
                        _mm256_mul_pd(im, _mm256_permute_pd(b0, 5)),
                    );
                    let p1 = _mm256_addsub_pd(
                        _mm256_mul_pd(re, b1),
                        _mm256_mul_pd(im, _mm256_permute_pd(b1, 5)),
                    );
                    _mm256_storeu_pd(dst.add(j), _mm256_add_pd(_mm256_loadu_pd(dst.add(j)), p0));
                    _mm256_storeu_pd(
                        dst.add(j + 4),
                        _mm256_add_pd(_mm256_loadu_pd(dst.add(j + 4)), p1),
                    );
                    j += 8;
                }
            }
            // General matrix dimensions are supported too, although circuit
            // dimensions are powers of two and never take this remainder.
            for col in j / 2..n {
                out[i * n + col] += coefficient * b[k * n + col];
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn dense(n: usize, seed: usize) -> Matrix {
        Matrix {
            n,
            data: (0..n * n)
                .map(|x| {
                    Complex::new(
                        ((x * 17 + seed) % 29) as f64 / 29.0,
                        ((x * 13 + seed) % 23) as f64 / 23.0,
                    )
                })
                .collect(),
        }
    }
    fn naive(a: &Matrix, b: &Matrix) -> Matrix {
        let mut out = Matrix::zero(a.n);
        for i in 0..a.n {
            for j in 0..a.n {
                for k in 0..a.n {
                    out[(i, j)] += a[(i, k)] * b[(k, j)];
                }
            }
        }
        out
    }
    #[test]
    fn product_matches_independent_triple_loop() {
        for n in [1, 2, 4, 8, 16, 32] {
            let a = dense(n, 7);
            let b = dense(n, 11);
            assert!(a.mul(&b).max_abs_diff(&naive(&a, &b)) < 1e-12);
        }
    }
    #[test]
    fn dispatched_kernel_matches_scalar() {
        for n in [3, 4, 5, 8, 16, 32, 64] {
            let a = dense(n, 12);
            let b = dense(n, 19);
            let mut scalar = Matrix::zero(n);
            a.mul_into_scalar(&b, &mut scalar);
            assert!(a.mul(&b).max_abs_diff(&scalar) < 1e-12);
        }
    }
    #[test]
    fn malformed_public_storage_panics_before_simd() {
        let valid = Matrix::identity(4);
        for length in [0, 15, 17] {
            let invalid = Matrix {
                n: 4,
                data: vec![Complex::ZERO; length],
            };
            assert!(std::panic::catch_unwind(|| invalid.mul(&valid)).is_err());
            assert!(std::panic::catch_unwind(|| valid.mul(&invalid)).is_err());
            assert!(std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                let mut output = invalid.clone();
                valid.mul_into(&valid, &mut output);
            }))
            .is_err());
        }
        let invalid = Matrix {
            n: usize::MAX,
            data: Vec::new(),
        };
        assert!(std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            let mut output = invalid.clone();
            invalid.mul_into(&invalid, &mut output);
        }))
        .is_err());
        assert_eq!(Matrix::zero(0).mul(&Matrix::zero(0)), Matrix::zero(0));
    }
    #[test]
    fn sparse_and_real_products() {
        for n in [2, 4, 8, 16] {
            let mut a = Matrix::identity(n);
            a[(0, 0)] = Complex::I;
            a[(1, 1)] = Complex::new(0.5, 0.0);
            let b = dense(n, 9);
            assert!(a.mul(&b).max_abs_diff(&naive(&a, &b)) < 1e-14);
        }
    }
    #[test]
    fn tensor_embedding_and_permutation() {
        let x = Matrix {
            n: 2,
            data: vec![Complex::ZERO, Complex::ONE, Complex::ONE, Complex::ZERO],
        };
        let x0 = x.kron_identity(3);
        let x2 = x0.permuted(&[2, 1, 0]);
        for col in 0..8 {
            assert_eq!(x0[(col ^ 1, col)], Complex::ONE);
            assert_eq!(x2[(col ^ 4, col)], Complex::ONE);
        }
        assert_eq!(x2.permuted(&[2, 1, 0]), x0);
    }
    #[test]
    fn every_three_qubit_permutation_is_invertible() {
        let a = dense(8, 4);
        for p in permutations(3) {
            let mut inv = vec![0; 3];
            for (i, &q) in p.iter().enumerate() {
                inv[q] = i;
            }
            assert_eq!(a.permuted(&p).permuted(&inv), a);
        }
    }
    #[test]
    fn adjoint_reverses_products() {
        let a = dense(4, 1);
        let b = dense(4, 2);
        assert!(
            a.mul(&b)
                .adjoint()
                .max_abs_diff(&b.adjoint().mul(&a.adjoint()))
                < 1e-14
        );
    }
    #[test]
    fn permutation_counts() {
        for (n, count) in [(0, 1), (1, 1), (2, 2), (3, 6), (4, 24)] {
            assert_eq!(permutations(n).len(), count);
        }
    }
}
