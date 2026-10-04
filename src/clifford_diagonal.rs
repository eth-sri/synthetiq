//! Diagonalization by native Clifford basis changes inferred from a target's
//! commuting Pauli operators. Only the basis is determined here; the diagonal
//! operator remains a stochastic phase-synthesis problem.
use crate::{
    affine::inverse_word, circuit::Circuit, gates::GateLibrary, matrix::Matrix,
    phase::NativeClifford,
};

#[derive(Clone, Copy, Debug)]
struct Pauli {
    x: usize,
    z: usize,
}
impl Pauli {
    fn xor(self, other: Self) -> Self {
        Self {
            x: self.x ^ other.x,
            z: self.z ^ other.z,
        }
    }
    fn symplectic(self, other: Self) -> bool {
        ((self.x & other.z).count_ones() + (self.z & other.x).count_ones()) % 2 == 1
    }
    fn bits(self, q: usize) -> usize {
        self.x | (self.z << q)
    }
}
fn insert(basis: &mut [usize], mut vector: usize) -> bool {
    for (bit, row) in basis.iter_mut().enumerate().rev() {
        if vector & (1 << bit) == 0 {
            continue;
        }
        if *row == 0 {
            *row = vector;
            return true;
        }
        vector ^= *row;
    }
    false
}
fn commutes(matrix: &Matrix, p: Pauli) -> bool {
    let sign = |bits: usize| {
        if (bits & p.z).count_ones().is_multiple_of(2) {
            1.0
        } else {
            -1.0
        }
    };
    (0..matrix.n).all(|i| {
        (0..matrix.n).all(|j| {
            (matrix[(i, j ^ p.x)] * sign(j) - matrix[(i ^ p.x, j)] * sign(i ^ p.x)).norm_sqr()
                < 1e-18
        })
    })
}
fn independent_commuting_paulis(matrix: &Matrix, q: usize) -> Option<Vec<Pauli>> {
    let mut basis = vec![0; 2 * q];
    let mut span = Vec::new();
    for bits in 1..1usize << (2 * q) {
        let p = Pauli {
            x: bits & (matrix.n - 1),
            z: bits >> q,
        };
        if commutes(matrix, p) && insert(&mut basis, p.bits(q)) {
            span.push(p);
        }
    }
    let mut isotropic = Vec::new();
    // Symplectic Gram-Schmidt: a radical vector can be retained directly;
    // from each anticommuting pair retain one and orthogonalize the remainder.
    while let Some(v) = span.pop() {
        if let Some(index) = span.iter().position(|&w| v.symplectic(w)) {
            let w = span.swap_remove(index);
            for u in &mut span {
                let a = u.symplectic(w);
                let b = u.symplectic(v);
                if a {
                    *u = u.xor(v);
                }
                if b {
                    *u = u.xor(w);
                }
            }
        }
        isotropic.push(v);
    }
    if isotropic.len() < q {
        None
    } else {
        isotropic.truncate(q);
        Some(isotropic)
    }
}

/// Returns W such that W U W† is diagonal, using only verified native Cliffords.
pub(crate) fn basis_word(
    matrix: &Matrix,
    lib: &GateLibrary,
    native: &NativeClifford,
) -> Option<Vec<usize>> {
    let q = lib.n_qubits;
    if q > 6 {
        return None;
    }
    let paulis = independent_commuting_paulis(matrix, q)?;
    let word = diagonalizing_word(paulis, native, false)?;
    let forward = Circuit::new(word.clone(), q, lib);
    let backward = Circuit::new(inverse_word(&word, lib)?, q, lib);
    let diagonal = forward.matrix().mul(matrix).mul(backward.matrix());
    (0..matrix.n)
        .all(|i| (0..matrix.n).all(|j| i == j || diagonal[(i, j)].norm_sqr() < 1e-18))
        .then_some(word)
}

/// A basis change for ordinary gate-by-gate annealing, with no T budget.
/// In time order the original operator is prefix -> reduced -> suffix.
pub(crate) struct Reduction {
    pub matrix: Matrix,
    pub prefix: Vec<usize>,
    pub suffix: Vec<usize>,
}

pub(crate) fn reduce(matrix: &Matrix, lib: &GateLibrary) -> Option<Reduction> {
    if lib.n_qubits > 4
        || (0..matrix.n).all(|i| (0..matrix.n).all(|j| i == j || matrix[(i, j)].norm_sqr() < 1e-18))
        || !matrix
            .mul(&matrix.adjoint())
            .approximately_equal(&Matrix::identity(matrix.n), 1e-9)
    {
        return None;
    }
    let native = NativeClifford::new(lib)?;
    let prefix = basis_word(matrix, lib, &native)?;
    let suffix = inverse_word(&prefix, lib)?;
    let forward = Circuit::new(prefix.clone(), lib.n_qubits, lib);
    let backward = Circuit::new(suffix.clone(), lib.n_qubits, lib);
    let reduced = forward.matrix().mul(matrix).mul(backward.matrix());
    Some(Reduction {
        matrix: reduced,
        prefix,
        suffix,
    })
}

fn diagonalizing_word(
    mut paulis: Vec<Pauli>,
    native: &NativeClifford,
    ordered: bool,
) -> Option<Vec<usize>> {
    let q = native.native.len();

    let mut word = Vec::new();
    let h = |bit: usize, ps: &mut [Pauli], word: &mut Vec<usize>| {
        for p in ps {
            let difference = ((p.x ^ p.z) >> bit) & 1;
            p.x ^= difference << bit;
            p.z ^= difference << bit;
        }
        word.push(native.native[bit][0]);
    };
    let sdg = |bit: usize, ps: &mut [Pauli], word: &mut Vec<usize>| {
        for p in ps {
            p.z ^= p.x & (1 << bit);
        }
        word.push(native.native[bit][2]);
    };
    let cx = |a: usize, b: usize, ps: &mut [Pauli], word: &mut Vec<usize>| {
        for p in ps {
            p.x ^= ((p.x >> a) & 1) << b;
            p.z ^= ((p.z >> b) & 1) << a;
        }
        word.push(native.cx[a][b].unwrap());
    };
    for k in 0..paulis.len() {
        let pivot = (k..q).find(|&b| (paulis[k].x | paulis[k].z) & (1 << b) != 0)?;
        if pivot != k {
            cx(k, pivot, &mut paulis, &mut word);
            cx(pivot, k, &mut paulis, &mut word);
            cx(k, pivot, &mut paulis, &mut word);
        }
        for bit in k..q {
            if paulis[k].x & (1 << bit) != 0 {
                if paulis[k].z & (1 << bit) != 0 {
                    sdg(bit, &mut paulis, &mut word);
                }
                h(bit, &mut paulis, &mut word);
            }
        }
        for bit in k + 1..q {
            if paulis[k].z & (1 << bit) != 0 {
                cx(bit, k, &mut paulis, &mut word);
            }
        }
        if ordered {
            for bit in 0..k {
                if paulis[k].z & (1 << bit) != 0 {
                    cx(bit, k, &mut paulis, &mut word);
                }
            }
        }
        debug_assert_eq!(paulis[k].x, 0);
        debug_assert_eq!(paulis[k].z, 1 << k);
        if !ordered {
            for row in k + 1..paulis.len() {
                if paulis[row].z & (1 << k) != 0 {
                    paulis[row] = paulis[row].xor(paulis[k]);
                }
            }
        }
    }
    Some(word)
}

pub(crate) fn diagonalize_pauli_axes(
    axes: &[(usize, usize)],
    native: &NativeClifford,
) -> Option<Vec<usize>> {
    let q = native.native.len();
    if axes.len() > q {
        return None;
    }
    let paulis: Vec<_> = axes.iter().map(|&(x, z)| Pauli { x, z }).collect();
    let mut basis = vec![0; 2 * q];
    for (i, &p) in paulis.iter().enumerate() {
        if p.x >= 1 << q
            || p.z >= 1 << q
            || !insert(&mut basis, p.bits(q))
            || paulis[..i].iter().any(|&other| p.symplectic(other))
        {
            return None;
        }
    }
    diagonalizing_word(paulis, native, true)
}

fn recognize_pauli(matrix: &Matrix, q: usize) -> Option<Pauli> {
    let n = matrix.n;
    let x = (0..n).find(|&i| (matrix[(i, 0)].norm_sqr() - 1.0).abs() < 1e-9)?;
    let base = matrix[(x, 0)];
    let mut z = 0;
    for bit in 0..q {
        let j = 1 << bit;
        let relative = matrix[(x ^ j, j)] / base;
        if (relative + crate::complex::Complex::ONE).norm_sqr() < 1e-18 {
            z |= j;
        } else if (relative - crate::complex::Complex::ONE).norm_sqr() >= 1e-18 {
            return None;
        }
    }
    for i in 0..n {
        for j in 0..n {
            let expected = if i == j ^ x {
                base * if (j & z).count_ones().is_multiple_of(2) {
                    1.0
                } else {
                    -1.0
                }
            } else {
                crate::complex::Complex::ZERO
            };
            if (matrix[(i, j)] - expected).norm_sqr() >= 1e-18 {
                return None;
            }
        }
    }
    Some(Pauli { x, z })
}

/// Independent Clifford changes on the input and output reduce any recognized
/// semi-Clifford operator to a monomial one. This does not discover its phases.
pub(crate) fn two_sided_basis(
    matrix: &Matrix,
    lib: &GateLibrary,
    native: &NativeClifford,
) -> Option<(Vec<usize>, Vec<usize>, Matrix)> {
    let q = lib.n_qubits;
    if q > 4 {
        return None;
    }
    let n = matrix.n;
    let adjoint = matrix.adjoint();
    if !matrix
        .mul(&adjoint)
        .approximately_equal(&Matrix::identity(n), 1e-9)
    {
        return None;
    }
    let mut basis = vec![0; 2 * q];
    let mut span = Vec::new();
    for bits in 1..1usize << (2 * q) {
        let p = Pauli {
            x: bits & (n - 1),
            z: bits >> q,
        };
        let mut product = Matrix::zero(n);
        for i in 0..n {
            for j in 0..n {
                product[(i, j)] = matrix[(i, j ^ p.x)]
                    * if (j & p.z).count_ones().is_multiple_of(2) {
                        1.0
                    } else {
                        -1.0
                    };
            }
        }
        if let Some(image) = recognize_pauli(&product.mul(&adjoint), q) {
            if insert(&mut basis, p.bits(q)) {
                span.push((p, image));
            }
        }
    }
    let mut isotropic = Vec::new();
    while let Some(v) = span.pop() {
        if let Some(index) = span.iter().position(|w| v.0.symplectic(w.0)) {
            let w = span.swap_remove(index);
            for u in &mut span {
                let a = u.0.symplectic(w.0);
                let b = u.0.symplectic(v.0);
                if a {
                    *u = (u.0.xor(v.0), u.1.xor(v.1));
                }
                if b {
                    *u = (u.0.xor(w.0), u.1.xor(w.1));
                }
            }
        }
        isotropic.push(v);
    }
    if isotropic.len() < q {
        return None;
    }
    isotropic.truncate(q);
    let input = diagonalizing_word(isotropic.iter().map(|p| p.0).collect(), native, false)?;
    let output = diagonalizing_word(isotropic.iter().map(|p| p.1).collect(), native, false)?;
    let a = Circuit::new(input.clone(), q, lib);
    let b = Circuit::new(output.clone(), q, lib);
    let transformed = b.matrix().mul(matrix).mul(&a.matrix().adjoint());
    for j in 0..n {
        if (0..n)
            .filter(|&i| (transformed[(i, j)].norm_sqr() - 1.0).abs() < 1e-9)
            .count()
            != 1
        {
            return None;
        }
    }
    Some((input, output, transformed))
}

/// Candidate single-T basis changes exposed by axes (P +/- Q)/sqrt(2)
/// commuting with the target. Their phase functions are still synthesized by SA.
pub(crate) fn near_pauli_bases(
    matrix: &Matrix,
    lib: &GateLibrary,
    native: &NativeClifford,
) -> Vec<(Vec<usize>, Matrix)> {
    use crate::complex::Complex;
    use std::collections::{HashMap, HashSet};
    let q = lib.n_qubits;
    if q > 4 {
        return Vec::new();
    }
    let rank = independent_commuting_paulis(matrix, q).map_or(0, |p| p.len());
    if rank == q {
        return Vec::new();
    }
    let n = matrix.n;
    let mut commutators = HashMap::<Vec<(i64, i64)>, Vec<Pauli>>::new();
    let mut axes = HashSet::new();
    for bits in 1..1usize << (2 * q) {
        let p = Pauli {
            x: bits & (n - 1),
            z: bits >> q,
        };
        let factor = [Complex::ONE, Complex::I, -Complex::ONE, -Complex::I]
            [(p.x & p.z).count_ones() as usize % 4];
        let sign = |x: usize| {
            if (x & p.z).count_ones().is_multiple_of(2) {
                1.0
            } else {
                -1.0
            }
        };
        let mut values = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                values.push(
                    (matrix[(i, j ^ p.x)] * sign(j) - matrix[(i ^ p.x, j)] * sign(i ^ p.x))
                        * factor,
                );
            }
        }
        let Some(first) = values.iter().find(|v| v.norm_sqr() > 1e-18) else {
            continue;
        };
        let sign = if first.re.abs() > 1e-9 {
            first.re.signum()
        } else {
            first.im.signum()
        };
        let key: Vec<_> = values
            .iter()
            .map(|v| {
                (
                    (v.re * sign * 1e8).round() as i64,
                    (v.im * sign * 1e8).round() as i64,
                )
            })
            .collect();
        let group = commutators.entry(key).or_default();
        for &other in group.iter() {
            if p.symplectic(other) {
                axes.insert(p.xor(other).bits(q));
            }
        }
        group.push(p);
    }
    // Sort before consuming randomness elsewhere; HashSet traversal must not
    // make seeded synthesis depend on the process's hash randomization.
    let mut axes: Vec<_> = axes.into_iter().collect();
    axes.sort_unstable();
    let mut candidates = Vec::new();
    for bits in axes {
        let p = Pauli {
            x: bits & (n - 1),
            z: bits >> q,
        };
        let pivot = (p.x | p.z).trailing_zeros() as usize;
        let mut basis = Vec::new();
        for bit in 0..q {
            if p.x & (1 << bit) != 0 {
                if p.z & (1 << bit) != 0 {
                    basis.push(native.native[bit][2]);
                }
                basis.push(native.native[bit][0]);
            }
        }
        for bit in 0..q {
            if bit != pivot && (p.x | p.z) & (1 << bit) != 0 {
                basis.push(native.cx[bit][pivot].unwrap());
            }
        }
        let inverse = inverse_word(&basis, lib).unwrap();
        for t in [3, 4] {
            let mut word = basis.clone();
            word.push(native.native[pivot][t]);
            word.extend_from_slice(&inverse);
            let rotation = Circuit::new(word.clone(), q, lib);
            let transformed = rotation
                .matrix()
                .adjoint()
                .mul(matrix)
                .mul(rotation.matrix());
            let new_rank = independent_commuting_paulis(&transformed, q).map_or(0, |p| p.len());
            if new_rank > rank {
                candidates.push((new_rank, word, transformed));
            }
        }
    }
    candidates.sort_by_key(|(rank, word, _)| (std::cmp::Reverse(*rank), word.len()));
    candidates.truncate(16);
    candidates
        .into_iter()
        .map(|(_, word, matrix)| (word, matrix))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{complex::Complex, rng::Rng};
    #[test]
    fn budget_free_basis_restores_arbitrary_angles_without_extra_gates() {
        let mut rng = Rng::new(826513);
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let ids: Vec<_> = lib
                .all
                .iter()
                .copied()
                .filter(|&g| matches!(lib.gates[g].name.as_str(), "h" | "s" | "sdg" | "cx"))
                .collect();
            for _ in 0..12 {
                let basis = Circuit::new(
                    (0..12).map(|_| ids[rng.usize(ids.len())]).collect(),
                    q,
                    &lib,
                );
                let mut diagonal = Matrix::zero(1 << q);
                for i in 0..diagonal.n {
                    let angle = rng.random01() * 6.0;
                    diagonal[(i, i)] = Complex::new(angle.cos(), angle.sin());
                }
                assert!(reduce(&diagonal, &lib).is_none());
                let original = basis.matrix().adjoint().mul(&diagonal).mul(basis.matrix());
                if let Some(r) = reduce(&original, &lib) {
                    assert!(r.prefix.iter().chain(&r.suffix).all(|g| ids.contains(g)));
                    let a = Circuit::new(r.prefix, q, &lib);
                    let b = Circuit::new(r.suffix, q, &lib);
                    assert!(
                        b.matrix()
                            .mul(&r.matrix)
                            .mul(a.matrix())
                            .max_abs_diff(&original)
                            < 1e-12
                    );
                } else {
                    assert!((0..original.n)
                        .all(|i| (0..original.n)
                            .all(|j| i == j || original[(i, j)].norm_sqr() < 1e-18)));
                }
            }
        }
    }

    #[test]
    fn unrelated_clifford_bases_restore_random_diagonal_operators() {
        let mut rng = Rng::new(250311);
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let native = NativeClifford::new(&lib).unwrap();
            let ids: Vec<_> = lib
                .all
                .iter()
                .copied()
                .filter(|&id| matches!(lib.gates[id].name.as_str(), "h" | "s" | "sdg" | "cx"))
                .collect();
            for _ in 0..12 {
                let left = Circuit::new(
                    (0..20).map(|_| ids[rng.usize(ids.len())]).collect(),
                    q,
                    &lib,
                );
                let right = Circuit::new(
                    (0..20).map(|_| ids[rng.usize(ids.len())]).collect(),
                    q,
                    &lib,
                );
                let mut diagonal = Matrix::zero(1 << q);
                for x in 0..diagonal.n {
                    let angle = rng.random01() * 6.0;
                    diagonal[(x, x)] = Complex::new(angle.cos(), angle.sin());
                }
                let original = left.matrix().mul(&diagonal).mul(right.matrix());
                let (a, b, reduced) = two_sided_basis(&original, &lib, &native).unwrap();
                let a = Circuit::new(a, q, &lib);
                let b = Circuit::new(b, q, &lib);
                assert!(
                    b.matrix()
                        .adjoint()
                        .mul(&reduced)
                        .mul(a.matrix())
                        .max_abs_diff(&original)
                        < 1e-12
                );
                assert_eq!(
                    a.count(&["t".into(), "tdg".into()], &lib)
                        + b.count(&["t".into(), "tdg".into()], &lib),
                    0
                );
            }
        }
    }
    #[test]
    fn entangled_clifford_bases_restore_random_diagonal_operators() {
        let mut rng = Rng::new(33809);
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let native = NativeClifford::new(&lib).unwrap();
            let cliffords: Vec<_> = lib
                .all
                .iter()
                .copied()
                .filter(|&id| matches!(lib.gates[id].name.as_str(), "h" | "s" | "sdg" | "cx"))
                .collect();
            for _ in 0..20 {
                let c = Circuit::new(
                    (0..24)
                        .map(|_| cliffords[rng.usize(cliffords.len())])
                        .collect(),
                    q,
                    &lib,
                );
                let mut d = Matrix::zero(1 << q);
                for x in 0..d.n {
                    let angle = std::f64::consts::FRAC_PI_4 * rng.usize(8) as f64;
                    d[(x, x)] = Complex::new(angle.cos(), angle.sin());
                }
                let target = c.matrix().mul(&d).mul(&c.matrix().adjoint());
                let word = basis_word(&target, &lib, &native).unwrap();
                let w = Circuit::new(word.clone(), q, &lib);
                let inverse = Circuit::new(inverse_word(&word, &lib).unwrap(), q, &lib);
                let normalized = w.matrix().mul(&target).mul(inverse.matrix());
                for i in 0..d.n {
                    for j in 0..d.n {
                        if i != j {
                            assert!(normalized[(i, j)].norm_sqr() < 1e-18);
                        }
                    }
                }
                assert!(
                    inverse
                        .matrix()
                        .mul(&normalized)
                        .mul(w.matrix())
                        .max_abs_diff(&target)
                        < 1e-12
                );
                assert_eq!(w.count(&["t".into(), "tdg".into()], &lib), 0);
            }
        }
    }
    #[test]
    fn non_pauli_eigenbasis_and_inadequate_commutants_are_rejected() {
        let lib = GateLibrary::load(1, "data/gates/CliffordT", "").unwrap();
        let native = NativeClifford::new(&lib).unwrap();
        let h = &lib.gates[lib.find("h", &[0]).unwrap()].matrix;
        assert!(basis_word(h, &lib, &native).is_none());
    }
}
