//! Exact parity-phase folding inside native Hadamard-free blocks. This is a
//! postprocessing rewrite; synthesis and coefficient discovery remain annealed.
use crate::{
    circuit::Circuit, gates::GateLibrary, phase::NativeClifford, phase_layers::PhaseLayout,
    quality::DepthModel, rng::Rng, search::QualityMetric,
};
use std::collections::HashSet;

fn subspaces(q: usize) -> Vec<(u64, usize)> {
    let n = 1usize << q;
    let mut spaces = vec![(1u64, 0usize)];
    let mut seen = HashSet::from([1u64]);
    let mut i = 0;
    while i < spaces.len() {
        let (space, dimension) = spaces[i];
        i += 1;
        for v in 1..n {
            if space & (1 << v) != 0 {
                continue;
            }
            let mut extended = space;
            for x in 0..n {
                if space & (1 << x) != 0 {
                    extended |= 1 << (x ^ v);
                }
            }
            if seen.insert(extended) {
                spaces.push((extended, dimension + 1));
            }
        }
    }
    spaces
}
fn layer_bound(coefficients: &[u8], spaces: &[(u64, usize)]) -> usize {
    let odd = coefficients
        .iter()
        .enumerate()
        .filter(|&(_, v)| v % 2 == 1)
        .fold(0u64, |bits, (m, _)| bits | (1 << (m + 1)));
    spaces
        .iter()
        .filter(|(_, dimension)| *dimension > 0)
        .map(|&(space, dimension)| ((space & odd).count_ones() as usize).div_ceil(dimension))
        .max()
        .unwrap_or(0)
}
fn phase_equal(a: &crate::matrix::Matrix, b: &crate::matrix::Matrix) -> bool {
    let z = b.trace_conjugate_product(a);
    if z.abs() < 1e-12 {
        return false;
    }
    let phase = z / z.abs();
    a.data
        .iter()
        .zip(&b.data)
        .all(|(&x, &y)| (x - phase * y).norm_sqr() < 1e-18)
}
fn key(c: &Circuit, lib: &GateLibrary, metric: QualityMetric) -> [f64; 4] {
    let names = ["t".into(), "tdg".into()];
    let t = c.count(&names, lib) as f64;
    let d = c.depth(&names, lib) as f64;
    let g = c.non_identity() as f64;
    let cost = c.cost();
    match metric {
        QualityMetric::TCount => [t, d, g, cost],
        QualityMetric::TDepth => [d, t, g, cost],
        QualityMetric::GateCount => [g, t, d, cost],
        QualityMetric::WeightedCost => [cost, t, d, g],
    }
}

pub fn run(circuit: &mut Circuit, lib: &GateLibrary, rng: &mut Rng, metric: QualityMetric) -> bool {
    let q = lib.n_qubits;
    if q > 4 {
        return false;
    }
    let Some(native) = NativeClifford::new(lib) else {
        return false;
    };
    let layout = PhaseLayout {
        native: &native.native,
        cx: &native.cx,
    };
    let spaces = subspaces(q);
    let depth = DepthModel::new(lib, false);
    let mut original = circuit.clone();
    original.expand(lib);
    let mut output = Vec::new();
    let mut begin = 0;
    let eligible = |id: usize| {
        matches!(
            lib.gates[id].name.as_str(),
            "id" | "cx" | "s" | "sdg" | "t" | "tdg"
        )
    };
    while begin < original.gates.len() {
        if !eligible(original.gates[begin]) {
            output.push(original.gates[begin]);
            begin += 1;
            continue;
        }
        let mut end = begin;
        let mut rows: Vec<_> = (0..q).map(|bit| 1 << bit).collect();
        let mut coefficients = vec![0u8; (1 << q) - 1];
        let mut tcount = 0;
        while end < original.gates.len() && eligible(original.gates[end]) {
            let g = &lib.gates[original.gates[end]];
            match g.name.as_str() {
                "cx" => rows[g.qubits[1]] ^= rows[g.qubits[0]],
                "id" => (),
                _ => {
                    let phase = match g.name.as_str() {
                        "s" => 2,
                        "sdg" => 6,
                        "t" => 1,
                        "tdg" => 7,
                        _ => unreachable!(),
                    };
                    tcount += usize::from(phase % 2 == 1);
                    let index = rows[g.qubits[0]] - 1;
                    coefficients[index] = (coefficients[index] + phase) % 8;
                }
            }
            end += 1;
        }
        let block = &original.gates[begin..end];
        let bound = layer_bound(&coefficients, &spaces);
        let count = coefficients.iter().filter(|&&c| c % 2 == 1).count();
        if count == tcount && depth.depth(block) <= bound {
            output.extend_from_slice(block);
        } else {
            let mut folded = layout.optimize(&coefficients, (&[], &[]), lib, rng, bound);
            folded.extend(layout.basis_word(rows));
            output.extend(folded);
        }
        begin = end;
    }
    let candidate = Circuit::new(output, q, lib);
    if key(&candidate, lib, metric) < key(&original, lib, metric)
        && phase_equal(original.matrix(), candidate.matrix())
    {
        *circuit = candidate;
        true
    } else {
        false
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn parity_subspaces_have_the_expected_small_dimensions() {
        for (q, count) in [(1, 2), (2, 5), (3, 16), (4, 67)] {
            let spaces = subspaces(q);
            assert_eq!(spaces.len(), count);
            for (space, rank) in spaces {
                assert_eq!(space.count_ones(), 1 << rank);
                assert_ne!(space & 1, 0);
            }
        }
    }
    #[test]
    fn relative_phase_toffoli_is_rewritten_to_two_t_layers() {
        let lib = GateLibrary::load(3, "data/gates/CliffordT", "").unwrap();
        let n = NativeClifford::new(&lib).unwrap();
        let t = n.native[2];
        let word = vec![
            t[0],
            t[3],
            n.cx[1][2].unwrap(),
            t[4],
            n.cx[0][2].unwrap(),
            t[3],
            n.cx[1][2].unwrap(),
            t[4],
            t[0],
        ];
        for seed in 0..16 {
            let mut c = Circuit::new(word.clone(), 3, &lib);
            let original = c.matrix().clone();
            assert!(run(
                &mut c,
                &lib,
                &mut Rng::new(seed),
                QualityMetric::TDepth
            ));
            assert_eq!(c.count(&["t".into(), "tdg".into()], &lib), 4);
            assert_eq!(c.depth(&["t".into(), "tdg".into()], &lib), 2);
            assert!(phase_equal(&original, c.matrix()));
        }
    }
    #[test]
    fn random_native_block_rewrites_preserve_the_complete_operator() {
        let mut rng = Rng::new(23907);
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            for _ in 0..24 {
                let original = crate::search::random_circuit(&lib, &mut rng, 30, 0.1);
                for metric in [QualityMetric::TCount, QualityMetric::TDepth] {
                    let mut candidate = original.clone();
                    run(&mut candidate, &lib, &mut rng, metric);
                    assert!(phase_equal(original.matrix(), candidate.matrix()));
                    assert!(key(&candidate, &lib, metric) <= key(&original, &lib, metric));
                }
            }
        }
    }
}
