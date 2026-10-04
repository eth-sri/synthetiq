//! T-layer scheduling across native Clifford gates. Clifford conjugation exposes
//! exact Pauli rotations; randomized topological packing preserves their order
//! whenever rotations do not commute.
use crate::{
    circuit::Circuit, complex::Complex, gates::GateLibrary, matrix::Matrix, phase::NativeClifford,
    rng::Rng, search::QualityMetric,
};
use std::collections::HashMap;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
struct Rotation {
    x: usize,
    z: usize,
    phase: u8,
}
impl Rotation {
    fn commutes(self, other: Self) -> bool {
        ((self.x & other.z).count_ones() + (self.z & other.x).count_ones()).is_multiple_of(2)
    }
    fn matrix(self, q: usize) -> Matrix {
        let mut m = Matrix::zero(1 << q);
        let factor = [Complex::ONE, Complex::I, -Complex::ONE, -Complex::I]
            [(self.x & self.z).count_ones() as usize % 4];
        for j in 0..m.n {
            m[(j ^ self.x, j)] = factor
                * if (j & self.z).count_ones().is_multiple_of(2) {
                    1.
                } else {
                    -1.
                };
        }
        m
    }
}
fn identify(matrix: &Matrix, q: usize, phase: u8) -> Option<Rotation> {
    let n = matrix.n;
    let x = (0..n).find(|&i| (matrix[(i, 0)].norm_sqr() - 1.).abs() < 1e-9)?;
    let base = matrix[(x, 0)];
    let mut z = 0;
    for bit in 0..q {
        let j = 1 << bit;
        let ratio = matrix[(x ^ j, j)] / base;
        if (ratio + Complex::ONE).norm_sqr() < 1e-18 {
            z |= j;
        } else if (ratio - Complex::ONE).norm_sqr() > 1e-18 {
            return None;
        }
    }
    let mut result = Rotation { x, z, phase };
    let expected = result.matrix(q);
    let sign = base / expected[(x, 0)];
    if (sign + Complex::ONE).norm_sqr() < 1e-18 {
        result.phase = (8 - phase) % 8;
    } else if (sign - Complex::ONE).norm_sqr() > 1e-18 {
        return None;
    }
    if matrix
        .data
        .iter()
        .zip(&expected.data)
        .any(|(&a, &b)| (a - b * sign).norm_sqr() > 1e-18)
    {
        return None;
    }
    Some(result)
}
fn split(circuit: &Circuit, lib: &GateLibrary) -> Option<(Vec<Rotation>, Vec<usize>)> {
    let q = lib.n_qubits;
    let n = 1 << q;
    let mut clifford = Matrix::identity(n);
    let mut word = Vec::new();
    let mut rotations = Vec::<Rotation>::new();
    for &id in &circuit.gates {
        let g = &lib.gates[id];
        match g.name.as_str() {
            "id" => (),
            "h" | "s" | "sdg" | "cx" => {
                clifford = g.matrix.mul(&clifford);
                word.push(id);
            }
            "t" | "tdg" => {
                let mut zc = clifford.clone();
                let bit = g.qubits[0];
                for row in 0..n {
                    if row & (1 << bit) != 0 {
                        for col in 0..n {
                            zc[(row, col)] = -zc[(row, col)];
                        }
                    }
                }
                let axis = clifford.adjoint().mul(&zc);
                let phase = if g.name == "t" { 1 } else { 7 };
                let rotation = identify(&axis, q, phase)?;
                let mut previous = None;
                for (index, &old) in rotations.iter().enumerate().rev() {
                    if old.x == rotation.x && old.z == rotation.z {
                        previous = Some(index);
                        break;
                    }
                    if !old.commutes(rotation) {
                        break;
                    }
                }
                if let Some(index) = previous {
                    rotations[index].phase = (rotations[index].phase + rotation.phase) % 8;
                    if rotations[index].phase == 0 {
                        rotations.remove(index);
                    }
                } else {
                    rotations.push(rotation);
                }
            }
            _ => return None,
        }
    }
    Some((rotations, word))
}
fn independent(basis: &mut [usize], mut vector: usize) -> bool {
    for bit in (0..basis.len()).rev() {
        if vector & (1 << bit) == 0 {
            continue;
        }
        if basis[bit] == 0 {
            basis[bit] = vector;
            return true;
        }
        vector ^= basis[bit];
    }
    false
}
fn emit(group: &[Rotation], lib: &GateLibrary, native: &NativeClifford) -> Option<Vec<usize>> {
    let q = lib.n_qubits;
    let axes: Vec<_> = group.iter().map(|r| (r.x, r.z)).collect();
    let mut word = crate::clifford_diagonal::diagonalize_pauli_axes(&axes, native)?;
    let basis = Circuit::new(word.clone(), q, lib);
    let adjoint = basis.matrix().adjoint();
    for (bit, &rotation) in group.iter().enumerate() {
        let diagonal = basis.matrix().mul(&rotation.matrix(q)).mul(&adjoint);
        let phase = if diagonal[(0, 0)].re > 0. {
            rotation.phase
        } else {
            (8 - rotation.phase) % 8
        };
        let ids = native.native[bit];
        match phase {
            0 => (),
            1 => word.push(ids[3]),
            2 => word.push(ids[1]),
            3 => word.extend([ids[1], ids[3]]),
            4 => word.extend([ids[1], ids[1]]),
            5 => word.extend([ids[2], ids[4]]),
            6 => word.push(ids[2]),
            7 => word.push(ids[4]),
            _ => unreachable!(),
        }
    }
    word.extend(crate::affine::inverse_word(&basis.gates, lib)?);
    Some(word)
}
fn phase_equal(a: &Matrix, b: &Matrix) -> bool {
    let z = a.trace_conjugate_product(b);
    if z.abs() < 1e-12 {
        return false;
    }
    let phase = z / z.abs();
    a.data
        .iter()
        .zip(&b.data)
        .all(|(&x, &y)| (x * phase - y).norm_sqr() < 1e-18)
}
fn key(c: &Circuit, lib: &GateLibrary) -> [usize; 3] {
    let names = ["t".into(), "tdg".into()];
    [c.depth(&names, lib), c.count(&names, lib), c.non_identity()]
}

fn metric_key(c: &Circuit, lib: &GateLibrary, metric: QualityMetric) -> [f64; 4] {
    let [d, t, g] = key(c, lib).map(|x| x as f64);
    match metric {
        QualityMetric::TCount => [t, d, g, c.cost()],
        QualityMetric::TDepth => [d, t, g, c.cost()],
        QualityMetric::GateCount => [g, t, d, c.cost()],
        QualityMetric::WeightedCost => [c.cost(), t, d, g],
    }
}

pub fn run(circuit: &mut Circuit, lib: &GateLibrary, rng: &mut Rng, goal: usize) -> bool {
    run_with_metric(circuit, lib, rng, QualityMetric::TDepth, Some(goal))
}

pub fn run_with_metric(
    circuit: &mut Circuit,
    lib: &GateLibrary,
    rng: &mut Rng,
    metric: QualityMetric,
    depth_goal: Option<usize>,
) -> bool {
    let q = lib.n_qubits;
    if q > 4 {
        return false;
    }
    let Some(native) = NativeClifford::new(lib) else {
        return false;
    };
    let mut original = circuit.clone();
    original.expand(lib);
    let Some((rotations, clifford)) = split(&original, lib) else {
        return false;
    };
    if rotations.is_empty() && key(&original, lib)[1] == 0 {
        return false;
    }
    let count = rotations.len();
    let mut successors = vec![Vec::new(); count];
    let mut incoming = vec![0; count];
    for i in 0..count {
        for j in i + 1..count {
            if !rotations[i].commutes(rotations[j]) {
                successors[i].push(j);
                incoming[j] += 1;
            }
        }
    }
    let mut best_key = metric_key(&original, lib, metric);
    let mut best = None;
    let mut cache = HashMap::<Vec<Rotation>, Vec<usize>>::new();
    for _ in 0..128 {
        let mut degree = incoming.clone();
        let mut done = vec![false; count];
        let mut remaining = count;
        let mut output = Vec::new();
        while remaining > 0 {
            let mut available: Vec<_> =
                (0..count).filter(|&i| !done[i] && degree[i] == 0).collect();
            for i in (1..available.len()).rev() {
                let j = rng.usize(i + 1);
                available.swap(i, j);
            }
            let mut basis = vec![0; 2 * q];
            let mut group = Vec::new();
            let mut indices = Vec::new();
            for i in available {
                let r = rotations[i];
                if group.len() < q && independent(&mut basis, r.x | (r.z << q)) {
                    group.push(r);
                    indices.push(i);
                }
            }
            if group.is_empty() {
                return false;
            }
            let word = if let Some(word) = cache.get(&group) {
                word.clone()
            } else {
                let Some(word) = emit(&group, lib, &native) else {
                    return false;
                };
                cache.insert(group, word.clone());
                word
            };
            output.extend(word);
            for i in indices {
                done[i] = true;
                remaining -= 1;
                for &j in &successors[i] {
                    degree[j] -= 1;
                }
            }
        }
        output.extend_from_slice(&clifford);
        let candidate = Circuit::new(output, q, lib);
        let candidate_key = metric_key(&candidate, lib, metric);
        if candidate_key < best_key && phase_equal(original.matrix(), candidate.matrix()) {
            best_key = candidate_key;
            best = Some(candidate);
        }
        if matches!(metric, QualityMetric::TDepth)
            && depth_goal.is_some_and(|goal| best_key[0] <= goal as f64)
            && best.is_some()
        {
            break;
        }
    }
    if let Some(best) = best {
        *circuit = best;
        true
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn budget_free_folding_preserves_operator_and_selected_resource_objective() {
        let mut source = Rng::new(568304);
        for q in 1..=3 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            for _ in 0..8 {
                let original = crate::search::random_circuit(&lib, &mut source, 16, 0.1);
                for metric in [
                    QualityMetric::WeightedCost,
                    QualityMetric::TCount,
                    QualityMetric::TDepth,
                    QualityMetric::GateCount,
                ] {
                    let mut candidate = original.clone();
                    run_with_metric(&mut candidate, &lib, &mut source, metric, None);
                    assert!(phase_equal(original.matrix(), candidate.matrix()));
                    assert!(
                        metric_key(&candidate, &lib, metric) <= metric_key(&original, &lib, metric)
                    );
                    assert!(candidate
                        .gates
                        .iter()
                        .all(|&g| lib.gates[g].decomposition.is_empty()));
                }
            }
        }
    }

    #[test]
    fn independent_pauli_layers_match_their_exact_rotation_matrices() {
        let mut rng = Rng::new(11346);
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let native = NativeClifford::new(&lib).unwrap();
            let ids: Vec<_> = lib
                .all
                .iter()
                .copied()
                .filter(|&id| matches!(lib.gates[id].name.as_str(), "h" | "s" | "sdg" | "cx"))
                .collect();
            for _ in 0..20 {
                let c = Circuit::new(
                    (0..20).map(|_| ids[rng.usize(ids.len())]).collect(),
                    q,
                    &lib,
                );
                let mut group = Vec::new();
                let mut expected = Matrix::identity(1 << q);
                for bit in 0..q {
                    let z = Rotation {
                        x: 0,
                        z: 1 << bit,
                        phase: 1,
                    }
                    .matrix(q);
                    let axis = c.matrix().adjoint().mul(&z).mul(c.matrix());
                    let phase = 1 + rng.usize(7) as u8;
                    let r = identify(&axis, q, phase).unwrap();
                    let p = r.matrix(q);
                    let angle = std::f64::consts::FRAC_PI_4 * r.phase as f64;
                    let root = Complex::new(angle.cos(), angle.sin());
                    let mut rotation = Matrix::identity(1 << q);
                    for i in 0..rotation.data.len() {
                        rotation.data[i] = (rotation.data[i] * (Complex::ONE + root)
                            + p.data[i] * (Complex::ONE - root))
                            * 0.5;
                    }
                    expected = rotation.mul(&expected);
                    group.push(r);
                }
                let circuit = Circuit::new(emit(&group, &lib, &native).unwrap(), q, &lib);
                assert!(phase_equal(&expected, circuit.matrix()));
                assert!(key(&circuit, &lib)[0] <= 1);
            }
        }
    }
    #[test]
    fn random_clifford_t_rewrites_preserve_the_full_operator_and_t_count() {
        let mut source = Rng::new(32630);
        let mut rng = Rng::new(52417);
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            for _ in 0..20 {
                let original = crate::search::random_circuit(&lib, &mut source, 30, 0.1);
                let mut candidate = original.clone();
                run(&mut candidate, &lib, &mut rng, 0);
                assert!(phase_equal(original.matrix(), candidate.matrix()));
                assert!(key(&candidate, &lib) <= key(&original, &lib));
                assert!(key(&candidate, &lib)[1] <= key(&original, &lib)[1]);
            }
        }
    }
    #[test]
    fn commuting_equal_pauli_rotations_cancel_across_clifford_gates() {
        let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
        let native = NativeClifford::new(&lib).unwrap();
        let a = native.native[0];
        let b = native.native[1];
        let word = vec![a[3], b[0], a[4], b[0]];
        let mut circuit = Circuit::new(word, 2, &lib);
        let original = circuit.matrix().clone();
        assert!(run(&mut circuit, &lib, &mut Rng::new(31), 0));
        assert_eq!(key(&circuit, &lib)[0], 0);
        assert_eq!(key(&circuit, &lib)[1], 0);
        assert!(phase_equal(&original, circuit.matrix()));
    }
}
