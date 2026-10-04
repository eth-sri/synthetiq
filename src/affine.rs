//! Reversible affine basis changes for monomial specifications.
//! Gaussian elimination emits only native CX and X=H S S H. For up to three
//! qubits, an exhaustive affine normalization changes the coordinates of the
//! synthesis problem; the nonlinear part is still synthesized by annealing.
use crate::{circuit::Circuit, gates::GateLibrary, matrix::Matrix};

pub struct AffineReduction {
    /// left * original * right
    pub matrix: Matrix,
    pub left: Vec<usize>,
    pub right: Vec<usize>,
}
fn linear(columns: &[usize], input: usize) -> usize {
    columns.iter().enumerate().fold(0, |y, (bit, &v)| {
        y ^ if input & (1 << bit) != 0 { v } else { 0 }
    })
}
fn inverse_linear(columns: &[usize]) -> Option<Vec<usize>> {
    let n = 1 << columns.len();
    let mut inverse = vec![usize::MAX; n];
    for input in 0..n {
        let output = linear(columns, input);
        if inverse[output] != usize::MAX {
            return None;
        }
        inverse[output] = input;
    }
    Some(inverse)
}
fn affine_columns(map: &[usize], q: usize) -> Option<Vec<usize>> {
    let columns: Vec<_> = (0..q).map(|bit| map[1 << bit] ^ map[0]).collect();
    (0..map.len())
        .all(|x| map[x] == (map[0] ^ linear(&columns, x)))
        .then_some(columns)
}
fn synthesize(map: &[usize], native: &[[usize; 5]], cx: &[Vec<Option<usize>>]) -> Vec<usize> {
    let q = native.len();
    let offset = map[0];
    let columns = affine_columns(map, q).expect("affine map");
    let mut rows: Vec<usize> = (0..q)
        .map(|r| (0..q).fold(0, |bits, c| bits | (((columns[c] >> r) & 1) << c)))
        .collect();
    let shortest = crate::linear_clifford::shortest(&rows, cx);
    let mut elimination = Vec::new();
    let apply = |control: usize, to: usize, rows: &mut [usize], gates: &mut Vec<usize>| {
        rows[to] ^= rows[control];
        gates.push(cx[control][to].unwrap());
    };
    for column in 0..q {
        if rows[column] & (1 << column) == 0 {
            let pivot = (column + 1..q)
                .find(|&r| rows[r] & (1 << column) != 0)
                .expect("invertible map");
            apply(column, pivot, &mut rows, &mut elimination);
            apply(pivot, column, &mut rows, &mut elimination);
            apply(column, pivot, &mut rows, &mut elimination);
        }
        for r in 0..q {
            if r != column && rows[r] & (1 << column) != 0 {
                apply(column, r, &mut rows, &mut elimination);
            }
        }
    }
    elimination.reverse();
    if let Some(word) = shortest {
        elimination = word;
    }
    for (bit, ids) in native.iter().enumerate() {
        if offset & (1 << bit) != 0 {
            elimination.extend([ids[0], ids[1], ids[1], ids[0]]);
        }
    }
    elimination
}
fn bases(q: usize) -> Vec<Vec<usize>> {
    fn visit(q: usize, columns: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        if columns.len() == q {
            out.push(columns.clone());
            return;
        }
        let span: Vec<_> = (0..1usize << columns.len())
            .map(|x| linear(columns, x))
            .collect();
        for v in 1..1usize << q {
            if span.contains(&v) {
                continue;
            }
            columns.push(v);
            visit(q, columns, out);
            columns.pop();
        }
    }
    let mut result = Vec::new();
    visit(q, &mut Vec::new(), &mut result);
    result
}

pub fn inverse_word(word: &[usize], lib: &GateLibrary) -> Option<Vec<usize>> {
    word.iter()
        .rev()
        .map(|&id| {
            let gate = &lib.gates[id];
            let name = match gate.name.as_str() {
                "h" | "cx" | "id" => gate.name.as_str(),
                "s" => "sdg",
                "sdg" => "s",
                "t" => "tdg",
                "tdg" => "t",
                _ => return None,
            };
            lib.find(name, &gate.qubits)
        })
        .collect()
}

/// A single nonlinear output direction can be isolated by affine coordinates.
/// This includes affine changes of a controlled-X with any Boolean controls.
fn shear_coordinates(permutation: &[usize], q: usize) -> Option<(Vec<usize>, Vec<usize>)> {
    let n = permutation.len();
    let mut anf = permutation.to_vec();
    for bit in 0..q {
        for mask in 0..n {
            if mask & (1 << bit) != 0 {
                anf[mask] ^= anf[mask ^ (1 << bit)];
            }
        }
    }
    let direction = (1usize..n)
        .filter(|&m| m.count_ones() >= 2)
        .map(|m| anf[m])
        .find(|&v| v != 0)?;
    if (1usize..n).any(|m| m.count_ones() >= 2 && anf[m] != 0 && anf[m] != direction) {
        return None;
    }
    let input_direction =
        (1..n).find(|&v| (0..n).all(|x| permutation[x ^ v] ^ permutation[x] == direction))?;
    let mut columns = vec![input_direction];
    for bit in 0..q {
        let unit = 1 << bit;
        if !(0..1usize << columns.len()).any(|x| linear(&columns, x) == unit) {
            columns.push(unit);
        }
    }
    let right: Vec<_> = (0..n).map(|x| linear(&columns, x)).collect();
    let origin = permutation[0];
    let image: Vec<_> = (0..q)
        .map(|bit| permutation[right[1 << bit]] ^ origin)
        .collect();
    let inverse = inverse_linear(&image)?;
    let left = (0..n).map(|x| inverse[x ^ origin]).collect();
    Some((left, right))
}

pub fn reduce(
    matrix: &Matrix,
    lib: &GateLibrary,
    native: &[[usize; 5]],
    cx: &[Vec<Option<usize>>],
) -> Option<AffineReduction> {
    reduce_impl(matrix, lib, native, cx, false)
}
/// Complete four-qubit affine normalization is a more expensive optional
/// preprocessing step for the general permutation-search portfolio.
pub fn reduce_full(
    matrix: &Matrix,
    lib: &GateLibrary,
    native: &[[usize; 5]],
    cx: &[Vec<Option<usize>>],
) -> Option<AffineReduction> {
    reduce_impl(matrix, lib, native, cx, true)
}
fn reduce_impl(
    matrix: &Matrix,
    lib: &GateLibrary,
    native: &[[usize; 5]],
    cx: &[Vec<Option<usize>>],
    four_qubit: bool,
) -> Option<AffineReduction> {
    let q = lib.n_qubits;
    let n = matrix.n;
    let mut permutation = Vec::with_capacity(n);
    let mut seen = vec![false; n];
    for input in 0..n {
        let output = (0..n).find(|&r| (matrix[(r, input)].norm_sqr() - 1.0).abs() < 1e-9)?;
        if seen[output] || (0..n).any(|r| r != output && matrix[(r, input)].norm_sqr() > 1e-18) {
            return None;
        }
        seen[output] = true;
        permutation.push(output);
    }
    let identity: Vec<_> = (0..n).collect();
    let shear = if q >= 4 {
        shear_coordinates(&permutation, q)
    } else {
        None
    };
    let (left_map, right_map) = if let Some(columns) = affine_columns(&permutation, q) {
        let inverse = inverse_linear(&columns)?;
        (
            (0..n).map(|x| inverse[x ^ permutation[0]]).collect(),
            identity,
        )
    } else if let Some(shear) = shear {
        shear
    } else if q <= 3 || (four_qubit && q == 4) {
        let mut best = vec![usize::MAX; n];
        let mut best_left = Vec::new();
        let mut best_right = Vec::new();
        let mut best_cost = usize::MAX;
        for columns in bases(q) {
            let lin: Vec<_> = (0..n).map(|x| linear(&columns, x)).collect();
            for offset in 0..n {
                let right: Vec<_> = lin.iter().map(|&x| x ^ offset).collect();
                let origin = permutation[right[0]];
                let image: Vec<_> = (0..q)
                    .map(|bit| permutation[right[1 << bit]] ^ origin)
                    .collect();
                let Some(inverse) = inverse_linear(&image) else {
                    continue;
                };
                let left: Vec<_> = (0..n).map(|x| inverse[x ^ origin]).collect();
                let candidate: Vec<_> = (0..n).map(|x| left[permutation[right[x]]]).collect();
                if candidate > best {
                    continue;
                }
                let cost =
                    synthesize(&left, native, cx).len() + synthesize(&right, native, cx).len();
                if candidate < best || cost < best_cost {
                    best = candidate;
                    best_left = left;
                    best_right = right;
                    best_cost = cost;
                }
            }
        }
        if best_left.is_empty() {
            return None;
        }
        (best_left, best_right)
    } else {
        return None;
    };
    let left = synthesize(&left_map, native, cx);
    let right = synthesize(&right_map, native, cx);
    let l = Circuit::new(left.clone(), q, lib);
    let r = Circuit::new(right.clone(), q, lib);
    let reduced = l.matrix().mul(matrix).mul(r.matrix());
    Some(AffineReduction {
        matrix: reduced,
        left,
        right,
    })
}

/// Fit the uniquely determined affine action visible in unit-magnitude columns
/// of a partly dense operator. Removing it exposes conjugated controlled gates.
pub(crate) fn reduce_sparse_columns(
    matrix: &Matrix,
    lib: &GateLibrary,
    native: &[[usize; 5]],
    cx: &[Vec<Option<usize>>],
) -> Option<AffineReduction> {
    let q = lib.n_qubits;
    let n = matrix.n;
    if q > 6 {
        return None;
    }
    let mut rows = vec![0usize; q + 1];
    let mut values = vec![0usize; q + 1];
    let mut constrained = 0;
    for x in 0..n {
        let Some(y) = (0..n).find(|&r| (matrix[(r, x)].norm_sqr() - 1.0).abs() < 1e-9) else {
            continue;
        };
        if (0..n).any(|r| r != y && matrix[(r, x)].norm_sqr() > 1e-18) {
            continue;
        }
        constrained += 1;
        let mut bits = x | (1 << q);
        let mut value = y;
        for pivot in (0..=q).rev() {
            if bits & (1 << pivot) == 0 {
                continue;
            }
            if rows[pivot] == 0 {
                rows[pivot] = bits;
                values[pivot] = value;
                bits = 0;
                value = 0;
                break;
            }
            bits ^= rows[pivot];
            value ^= values[pivot];
        }
        if bits == 0 && value != 0 {
            return None;
        }
    }
    if constrained == n || rows.contains(&0) {
        return None;
    }
    let mut coefficients = vec![0usize; q + 1];
    for pivot in 0..=q {
        coefficients[pivot] = values[pivot];
        for bit in 0..pivot {
            if rows[pivot] & (1 << bit) != 0 {
                coefficients[pivot] ^= coefficients[bit];
            }
        }
    }
    let inverse = inverse_linear(&coefficients[..q])?;
    let left_map: Vec<_> = (0..n).map(|x| inverse[x ^ coefficients[q]]).collect();
    let left = synthesize(&left_map, native, cx);
    let circuit = Circuit::new(left.clone(), q, lib);
    Some(AffineReduction {
        matrix: circuit.matrix().mul(matrix),
        left,
        right: Vec::new(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{complex::Complex, rng::Rng};
    fn metadata(lib: &GateLibrary) -> (Vec<[usize; 5]>, Vec<Vec<Option<usize>>>) {
        let native = (0..lib.n_qubits)
            .map(|q| ["h", "s", "sdg", "t", "tdg"].map(|name| lib.find(name, &[q]).unwrap()))
            .collect();
        let cx = (0..lib.n_qubits)
            .map(|a| {
                (0..lib.n_qubits)
                    .map(|b| {
                        if a != b {
                            lib.find("cx", &[a, b])
                        } else {
                            None
                        }
                    })
                    .collect()
            })
            .collect();
        (native, cx)
    }
    #[test]
    fn complete_four_qubit_normalization_matches_affine_wrapped_u1() {
        let lib = GateLibrary::load(4, "data/gates/CliffordT", "").unwrap();
        let (native, cx) = metadata(&lib);
        let original = crate::partial::PartialMatrix::read("data/input/64/comparison/U1.txt")
            .unwrap()
            .matrix;
        let expected = reduce_full(&original, &lib, &native, &cx).unwrap();
        let mut rng = Rng::new(46225);
        for _ in 0..8 {
            let mut words = [Vec::new(), Vec::new()];
            for word in &mut words {
                for _ in 0..12 {
                    let a = rng.usize(4);
                    let mut b = rng.usize(3);
                    if b >= a {
                        b += 1;
                    }
                    word.push(cx[a][b].unwrap());
                    if rng.random01() < 0.2 {
                        let ids = native[b];
                        word.extend([ids[0], ids[1], ids[1], ids[0]]);
                    }
                }
            }
            let left = Circuit::new(words[0].clone(), 4, &lib);
            let right = Circuit::new(words[1].clone(), 4, &lib);
            let matrix = left.matrix().mul(&original).mul(right.matrix());
            let result = reduce_full(&matrix, &lib, &native, &cx).unwrap();
            assert!(result.matrix.max_abs_diff(&expected.matrix) < 1e-12);
            let l = Circuit::new(inverse_word(&result.left, &lib).unwrap(), 4, &lib);
            let r = Circuit::new(inverse_word(&result.right, &lib).unwrap(), 4, &lib);
            assert!(
                l.matrix()
                    .mul(&result.matrix)
                    .mul(r.matrix())
                    .max_abs_diff(&matrix)
                    < 1e-12
            );
        }
    }
    #[test]
    fn affine_circuits_match_all_two_qubit_affine_maps_and_random_larger_maps() {
        let mut rng = Rng::new(376);
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let (native, cx) = metadata(&lib);
            let n = 1 << q;
            let all = bases(q);
            let sample: Vec<_> = if q <= 2 {
                all
            } else {
                (0..30).map(|_| all[rng.usize(all.len())].clone()).collect()
            };
            for columns in sample {
                for offset in 0..n {
                    let map: Vec<_> = (0..n).map(|x| offset ^ linear(&columns, x)).collect();
                    let word = synthesize(&map, &native, &cx);
                    let c = Circuit::new(word.clone(), q, &lib);
                    let mut expected = Matrix::zero(n);
                    for x in 0..n {
                        expected[(map[x], x)] = Complex::ONE;
                    }
                    assert!(expected.max_abs_diff(c.matrix()) < 1e-12);
                    let inverse = Circuit::new(inverse_word(&word, &lib).unwrap(), q, &lib);
                    assert!(
                        inverse
                            .matrix()
                            .mul(c.matrix())
                            .max_abs_diff(&Matrix::identity(n))
                            < 1e-12
                    );
                    assert_eq!(c.count(&["t".into(), "tdg".into()], &lib), 0);
                }
            }
        }
    }
    #[test]
    fn larger_affine_wrapped_shears_are_normalized_and_restored() {
        let mut rng = Rng::new(18924);
        for q in [4, 5] {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let (native, cx) = metadata(&lib);
            for _ in 0..16 {
                let n = 1 << q;
                let mut matrix = Matrix::zero(n);
                for x in 0usize..n {
                    let flip = ((x >> 1) & 1) & ((x >> 2) & 1);
                    matrix[(x ^ flip, x)] = Complex::ONE;
                }
                let mut left = Vec::new();
                let mut right = Vec::new();
                for word in [&mut left, &mut right] {
                    for _ in 0..16 {
                        let a = rng.usize(q);
                        let mut b = rng.usize(q - 1);
                        if b >= a {
                            b += 1;
                        }
                        word.push(cx[a][b].unwrap());
                        if rng.usize(4) == 0 {
                            word.extend([native[a][0], native[a][1], native[a][1], native[a][0]]);
                        }
                    }
                }
                let l = Circuit::new(left, q, &lib);
                let r = Circuit::new(right, q, &lib);
                let original = l.matrix().mul(&matrix).mul(r.matrix());
                let reduced = reduce(&original, &lib, &native, &cx).unwrap();
                for x in 0..n {
                    for y in 0..n {
                        if x != y && x ^ y != 1 {
                            assert!(reduced.matrix[(y, x)].norm_sqr() < 1e-20);
                        }
                    }
                }
                let li = Circuit::new(inverse_word(&reduced.left, &lib).unwrap(), q, &lib);
                let ri = Circuit::new(inverse_word(&reduced.right, &lib).unwrap(), q, &lib);
                assert!(
                    li.matrix()
                        .mul(&reduced.matrix)
                        .mul(ri.matrix())
                        .max_abs_diff(&original)
                        < 1e-12
                );
            }
        }
    }

    #[test]
    fn normalization_preserves_the_original_operator_and_is_affine_invariant() {
        let lib = GateLibrary::load(3, "data/gates/CliffordT", "").unwrap();
        let (native, cx) = metadata(&lib);
        let all = bases(3);
        let mut rng = Rng::new(1539);
        for _ in 0..8 {
            let mut permutation: Vec<_> = (0..8).collect();
            for i in (1..8).rev() {
                permutation.swap(i, rng.usize(i + 1));
            }
            let mut matrix = Matrix::zero(8);
            for i in 0..8 {
                matrix[(permutation[i], i)] = Complex::ONE;
            }
            let reduced = reduce(&matrix, &lib, &native, &cx).unwrap();
            let l = Circuit::new(inverse_word(&reduced.left, &lib).unwrap(), 3, &lib);
            let r = Circuit::new(inverse_word(&reduced.right, &lib).unwrap(), 3, &lib);
            assert!(
                l.matrix()
                    .mul(&reduced.matrix)
                    .mul(r.matrix())
                    .max_abs_diff(&matrix)
                    < 1e-12
            );
            for _ in 0..3 {
                let a = &all[rng.usize(all.len())];
                let b = &all[rng.usize(all.len())];
                let x = rng.usize(8);
                let y = rng.usize(8);
                let am: Vec<_> = (0..8).map(|v| x ^ linear(a, v)).collect();
                let bm: Vec<_> = (0..8).map(|v| y ^ linear(b, v)).collect();
                let ac = Circuit::new(synthesize(&am, &native, &cx), 3, &lib);
                let bc = Circuit::new(synthesize(&bm, &native, &cx), 3, &lib);
                let transformed = ac.matrix().mul(&matrix).mul(bc.matrix());
                let again = reduce(&transformed, &lib, &native, &cx).unwrap();
                assert!(again.matrix.max_abs_diff(&reduced.matrix) < 1e-12);
            }
        }
    }
}
