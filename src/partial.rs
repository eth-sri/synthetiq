//! Partial specifications and the original phase-insensitive distance functions.
use crate::{complex::Complex, matrix::Matrix};
use std::{fs, path::Path};

#[derive(Clone, Debug)]
pub struct PartialMatrix {
    pub name: String,
    pub n_qubits: usize,
    pub matrix: Matrix,
    pub cover: Vec<bool>,
    pub squared_norm: f64,
    // Only covered entries are visited in the annealing hot path.
    entries: Vec<(usize, Complex)>,
    normalization: f64,
}

impl PartialMatrix {
    pub fn new(matrix: Matrix, cover: Vec<bool>, name: impl Into<String>) -> Result<Self, String> {
        if !matrix.n.is_power_of_two()
            || matrix.n.checked_mul(matrix.n) != Some(matrix.data.len())
            || cover.len() != matrix.data.len()
        {
            return Err(
                "specification must be square, power-of-two sized, with a matching cover".into(),
            );
        }
        let mut matrix = matrix;
        let mut entries = Vec::new();
        let mut squared_norm = 0.0;
        for (i, value) in matrix.data.iter_mut().enumerate() {
            if !value.re.is_finite() || !value.im.is_finite() {
                return Err("specification contains non-finite values".into());
            }
            if cover[i] {
                squared_norm += value.norm_sqr();
                entries.push((i, value.conj()));
            } else {
                *value = Complex::new(0.0, 0.0);
            }
        }
        let normalization = (entries.len() as f64).sqrt().sqrt();
        Ok(Self {
            n_qubits: matrix.n.trailing_zeros() as usize,
            name: name.into(),
            matrix,
            cover,
            squared_norm,
            entries,
            normalization,
        })
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self, String> {
        let path = path.as_ref();
        let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        Self::parse(&text).map_err(|e| format!("{}: {e}", path.display()))
    }

    pub fn parse(text: &str) -> Result<Self, String> {
        let mut fields = Vec::new();
        let mut rest = text.trim_start();
        while !rest.is_empty() {
            let end = if rest.starts_with('(') {
                rest.find(')').ok_or("unterminated complex value")? + 1
            } else {
                rest.find(char::is_whitespace).unwrap_or(rest.len())
            };
            fields.push(&rest[..end]);
            rest = rest[end..].trim_start();
        }
        let mut tokens = fields.into_iter();
        let name = tokens.next().ok_or("missing specification name")?;
        let nq: u32 = tokens
            .next()
            .ok_or("missing qubit count")?
            .parse()
            .map_err(|_| "invalid qubit count")?;
        let n = 1usize.checked_shl(nq).ok_or("qubit count is too large")?;
        let len = n.checked_mul(n).ok_or("matrix size overflow")?;
        // Check the file size before allocating for an untrusted dimension.
        if len > text.len() {
            return Err("matrix dimensions exceed file contents".into());
        }
        let mut data = Vec::with_capacity(len);
        for _ in 0..len {
            let s = tokens.next().ok_or("missing matrix entry")?;
            let value = s.parse::<Complex>()?;
            data.push(value);
        }
        let mut cover = Vec::with_capacity(len);
        for _ in 0..len {
            cover.push(match tokens.next() {
                Some("0") => false,
                Some("1") => true,
                _ => return Err("missing or invalid cover entry".into()),
            });
        }
        if tokens.next().is_some() {
            return Err("unexpected data after specification".into());
        }
        Self::new(Matrix { n, data }, cover, name)
    }

    pub fn write(&self, path: impl AsRef<Path>) -> Result<(), String> {
        use std::fmt::Write;
        let mut s = format!("{}\n{}\n", self.name, self.n_qubits);
        for row in self.matrix.data.chunks(self.matrix.n) {
            for z in row {
                write!(&mut s, "({:.17e},{:.17e}) ", z.re, z.im).unwrap();
            }
            s.push('\n');
        }
        for row in self.cover.chunks(self.matrix.n) {
            for &v in row {
                s.push_str(if v { "1 " } else { "0 " });
            }
            s.push('\n');
        }
        fs::write(path, s).map_err(|e| e.to_string())
    }

    pub fn n_constraints(&self) -> usize {
        self.entries.len()
    }

    pub fn cost(&self, candidate: &Matrix, simple: bool) -> f64 {
        assert_eq!(self.matrix.n, candidate.n);
        let mut circ_size = 0.0;
        let mut overlap = Complex::new(0.0, 0.0);
        for &(i, conjugate) in &self.entries {
            let value = candidate.data[i];
            circ_size += value.norm_sqr();
            overlap += conjugate * value;
        }
        if simple {
            (1.0 - overlap.abs() / (self.entries.len() as f64).sqrt())
                .max(0.0)
                .sqrt()
                * std::f64::consts::SQRT_2
        } else {
            (self.squared_norm + circ_size - 2.0 * overlap.abs())
                .max(0.0)
                .sqrt()
                / self.normalization
        }
    }

    pub fn exact_cost(&self, candidate: &Matrix, epsilon: f64) -> f64 {
        if self.cost(candidate, false) / std::f64::consts::SQRT_2 <= epsilon {
            0.0
        } else {
            1.0
        }
    }

    pub fn add_ancillas(&self, count: usize) -> Result<Self, String> {
        if count == 0 {
            return Ok(self.clone());
        }
        let n = self
            .matrix
            .n
            .checked_shl(count.try_into().map_err(|_| "too many ancillas")?)
            .ok_or("too many ancillas")?;
        let mut matrix = Matrix::zero(n);
        let mut cover = vec![false; n * n];
        let old = self.matrix.n;
        for r in 0..old {
            for c in 0..old {
                matrix.data[r * n + c] = self.matrix.data[r * old + c];
                cover[r * n + c] = self.cover[r * old + c];
            }
        }
        // Preserve MatrixGenerator::addAncilla's bottomLeftCorner(old,n-old),
        // including its asymmetric rectangle when adding multiple ancillas.
        for r in n - old..n {
            for c in 0..n - old {
                cover[r * n + c] = true;
            }
        }
        Self::new(matrix, cover, "matrix")
    }

    pub fn add_dirty(&self, count: usize) -> Result<Self, String> {
        if count == 0 {
            return Ok(self.clone());
        }
        let blocks = 1usize
            .checked_shl(count.try_into().map_err(|_| "too many dirty qubits")?)
            .ok_or("too many dirty qubits")?;
        let old = self.matrix.n;
        let n = old.checked_mul(blocks).ok_or("matrix size overflow")?;
        let mut matrix = Matrix::zero(n);
        let mut cover = vec![false; n * n];
        for b in 0..blocks {
            for r in 0..old {
                for c in 0..old {
                    let i = (b * old + r) * n + b * old + c;
                    matrix.data[i] = self.matrix.data[r * old + c];
                    cover[i] = self.cover[r * old + c];
                }
            }
        }
        Self::new(matrix, cover, "matrix")
    }
}

/// Map a new basis index to the old basis index, matching PartialMatrix.cpp.
pub fn basis_permutation(n: usize, permutation: &[usize]) -> Vec<usize> {
    (0..n)
        .map(|i| {
            permutation
                .iter()
                .enumerate()
                .fold(0, |j, (bit, &to)| j | (((i >> bit) & 1) << to))
        })
        .collect()
}

pub fn next_permutation(values: &mut [usize]) -> bool {
    let Some(i) = (1..values.len()).rfind(|&i| values[i - 1] < values[i]) else {
        return false;
    };
    let j = (i..values.len())
        .rfind(|&j| values[j] > values[i - 1])
        .unwrap();
    values.swap(i - 1, j);
    values[i..].reverse();
    true
}

#[derive(Clone, Debug)]
pub struct Target {
    pub original: PartialMatrix,
    pub variants: Vec<PartialMatrix>,
    pub permutations: Vec<Vec<usize>>,
    pub inverses: Vec<bool>,
}

impl Target {
    pub fn new(original: PartialMatrix, independent: bool, inverse: bool) -> Self {
        let mut out = Self {
            original,
            variants: Vec::new(),
            permutations: Vec::new(),
            inverses: Vec::new(),
        };
        let mut p: Vec<_> = (0..out.original.n_qubits).collect();
        let n = out.original.matrix.n;
        loop {
            let change = basis_permutation(n, &p);
            let mut matrix = Matrix::zero(n);
            let mut cover = vec![false; n * n];
            for r in 0..n {
                for c in 0..n {
                    matrix.data[r * n + c] = out.original.matrix.data[change[r] * n + change[c]];
                    cover[r * n + c] = out.original.cover[change[r] * n + change[c]];
                }
            }
            if independent || out.variants.is_empty() {
                out.add(
                    PartialMatrix::new(matrix.clone(), cover.clone(), "matrix").unwrap(),
                    &p,
                    false,
                );
            }
            if inverse {
                // C++ deliberately uses conjugate(), with a transposed cover.
                for r in 0..n {
                    for c in 0..r {
                        cover.swap(r * n + c, c * n + r);
                    }
                }
                for z in &mut matrix.data {
                    *z = z.conj();
                }
                out.add(
                    PartialMatrix::new(matrix, cover, "matrix").unwrap(),
                    &p,
                    true,
                );
            }
            if (!independent && !inverse) || !next_permutation(&mut p) {
                break;
            }
        }
        out
    }

    fn add(&mut self, m: PartialMatrix, permutation: &[usize], inverse: bool) {
        if self.variants.iter().any(|v| {
            v.cover == m.cover
                && v.matrix
                    .data
                    .iter()
                    .zip(&m.matrix.data)
                    .map(|(&a, &b)| (a - b).norm_sqr())
                    .sum::<f64>()
                    < 1e-12
        }) {
            return;
        }
        self.variants.push(m);
        self.permutations.push(permutation.to_vec());
        self.inverses.push(inverse);
    }

    pub fn closest(&self, candidate: &Matrix, simple: bool) -> usize {
        let mut best = 0;
        let mut cost = f64::INFINITY;
        for (i, target) in self.variants.iter().enumerate() {
            let value = target.cost(candidate, simple);
            if value < cost {
                cost = value;
                best = i;
            }
        }
        best
    }
    pub fn cost(&self, candidate: &Matrix, simple: bool) -> f64 {
        self.variants
            .iter()
            .map(|m| m.cost(candidate, simple))
            .fold(f64::INFINITY, f64::min)
    }
    pub fn exact_cost(&self, candidate: &Matrix, epsilon: f64) -> f64 {
        if self
            .variants
            .iter()
            .any(|m| m.exact_cost(candidate, epsilon) == 0.0)
        {
            0.0
        } else {
            1.0
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn identity() -> PartialMatrix {
        PartialMatrix::new(Matrix::identity(2), vec![true; 4], "id").unwrap()
    }
    #[test]
    fn phase_invariance() {
        let p = identity();
        let mut m = Matrix::identity(2);
        for z in &mut m.data {
            *z *= Complex::new(0.0, 1.0);
        }
        assert_eq!(p.exact_cost(&m, 1e-6), 0.0);
        assert!(p.cost(&m, false) < 1e-7);
    }
    #[test]
    fn masks_ignore_unspecified_entries() {
        let p = PartialMatrix::new(
            Matrix::identity(2),
            vec![true, false, false, false],
            "partial",
        )
        .unwrap();
        let mut m = Matrix::identity(2);
        m.data[3] = Complex::new(100.0, 50.0);
        assert_eq!(p.cost(&m, false), 0.0);
    }
    #[test]
    fn original_ancilla_and_dirty_covers() {
        let p = identity();
        assert_eq!(p.add_ancillas(1).unwrap().n_constraints(), 8);
        assert_eq!(p.add_ancillas(2).unwrap().n_constraints(), 16);
        assert_eq!(p.add_dirty(2).unwrap().n_constraints(), 16);
    }
    #[test]
    fn symmetric_variants_deduplicate() {
        let p = PartialMatrix::new(Matrix::identity(8), vec![true; 64], "id").unwrap();
        assert_eq!(Target::new(p, true, true).variants.len(), 1);
    }
    #[test]
    fn rejects_malformed_public_matrix_storage() {
        for len in [0, 3, 5, 8] {
            let matrix = Matrix {
                n: 2,
                data: vec![Complex::ZERO; len],
            };
            assert!(PartialMatrix::new(matrix, vec![true; len], "invalid").is_err());
        }
    }
    #[test]
    fn parse_validation() {
        assert!(PartialMatrix::parse("id 1 (1,0) 0 0 1 1 1 1 1").is_ok());
        assert!(PartialMatrix::parse("id 1 (1, 0) (0,0) (0) 1 1 1 1 1").is_ok());
        assert!(PartialMatrix::parse("id 1 1").is_err());
        assert!(PartialMatrix::parse("id 900 1").is_err());
        assert!(PartialMatrix::parse("id 1 NaN 0 0 1 1 1 1 1").is_err());
    }
    #[test]
    fn permutation_cycle_convention() {
        assert_eq!(
            basis_permutation(8, &[1, 2, 0]),
            vec![0, 2, 4, 6, 1, 3, 5, 7]
        );
    }
}
