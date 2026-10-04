//! Random benchmark circuits and full, isometric, or arbitrary specifications.
use crate::{
    circuit::Circuit, gates::GateLibrary, partial::PartialMatrix, resynthesis::Resynth, rng::Rng,
    search::random_gate,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[repr(u8)]
pub enum SpecificationKind {
    Full = 0,
    Isometry = 1,
    Random = 2,
}

impl TryFrom<usize> for SpecificationKind {
    type Error = String;
    fn try_from(value: usize) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::Full),
            1 => Ok(Self::Isometry),
            2 => Ok(Self::Random),
            _ => Err("specification type must be 0 (full), 1 (isometry), or 2 (random)".into()),
        }
    }
}

/// Generate until local simplification retains at least the requested size.
/// As in the reference, return the original circuit, not its simplified copy.
pub fn generate_circuit(
    lib: &GateLibrary,
    rng: &mut Rng,
    n_gates: usize,
    extra_gates: usize,
    identity_probability: f64,
    ensure_non_identity: bool,
) -> Circuit {
    assert!(
        !ensure_non_identity
            || (identity_probability < 1.0
                && lib
                    .all
                    .iter()
                    .any(|&id| lib.gates[id].name != lib.gates[lib.identity].name)),
        "non-identity generation requires a non-identity gate and probability below one"
    );
    loop {
        let gates = (0..n_gates + extra_gates)
            .map(|_| loop {
                let gate = random_gate(lib, rng, identity_probability, 1.0, 0.5);
                if !ensure_non_identity || lib.gates[gate].name != lib.gates[lib.identity].name {
                    break gate;
                }
            })
            .collect();
        let original = Circuit::new(gates, lib.n_qubits, lib);
        let mut simplified = original.clone();
        Resynth::default().run(&mut simplified, lib);
        if simplified.non_identity() >= n_gates {
            return original;
        }
    }
}

pub fn generate_matrix(
    circuit: &Circuit,
    kind: SpecificationKind,
    rng: &mut Rng,
) -> Result<PartialMatrix, String> {
    let matrix = circuit.matrix();
    let n = matrix.n;
    if n < 2 && kind != SpecificationKind::Full {
        return Err("a proper partial specification requires at least one qubit".into());
    }
    let mut cover = vec![kind == SpecificationKind::Full; n * n];
    match kind {
        SpecificationKind::Full => {}
        SpecificationKind::Isometry => loop {
            let probability = rng.random01();
            let selected: Vec<_> = (0..n).filter(|_| rng.random01() > probability).collect();
            if selected.is_empty() || selected.len() == n {
                continue;
            }
            for col in selected {
                for row in 0..n {
                    cover[row * n + col] = true;
                }
            }
            break;
        },
        SpecificationKind::Random => {
            // The C++ generator consumes a threshold before its retry loop.
            rng.random01();
            loop {
                let probability = rng.random01();
                for entry in &mut cover {
                    *entry = probability < rng.random01();
                }
                let count = cover.iter().filter(|&&x| x).count();
                if count > 0 && count < n * n {
                    break;
                }
            }
        }
    }
    PartialMatrix::new(matrix.clone(), cover, "matrix")
}

#[cfg(test)]
mod tests {
    use super::*;
    fn circuit(n: usize) -> Circuit {
        let lib = GateLibrary::load(n, "data/gates/CliffordT", "").unwrap();
        let h = lib.find("h", &[0]).unwrap();
        Circuit::new(vec![h], n, &lib)
    }
    #[test]
    fn full_specs_keep_every_entry_and_consume_no_randomness() {
        let circuit = circuit(2);
        let mut rng = Rng::new(7);
        let target = generate_matrix(&circuit, SpecificationKind::Full, &mut rng).unwrap();
        assert_eq!(target.n_constraints(), 16);
        assert_eq!(target.exact_cost(circuit.matrix(), 1e-6), 0.0);
        assert_eq!(rng.next_u32(), Rng::new(7).next_u32());
    }
    #[test]
    fn isometries_select_whole_columns_and_never_empty_or_full() {
        for n_qubits in 1..=4 {
            let circuit = circuit(n_qubits);
            let n = 1 << n_qubits;
            for seed in 0..32 {
                let target =
                    generate_matrix(&circuit, SpecificationKind::Isometry, &mut Rng::new(seed))
                        .unwrap();
                assert!(target.n_constraints() > 0 && target.n_constraints() < n * n);
                for row in 1..n {
                    assert_eq!(&target.cover[..n], &target.cover[row * n..row * n + n]);
                }
                assert_eq!(target.exact_cost(circuit.matrix(), 1e-6), 0.0);
            }
        }
    }
    #[test]
    fn random_specs_are_proper_subsets_and_keep_covered_values() {
        for n_qubits in 1..=4 {
            let circuit = circuit(n_qubits);
            let n = 1 << n_qubits;
            for seed in 0..32 {
                let target =
                    generate_matrix(&circuit, SpecificationKind::Random, &mut Rng::new(seed))
                        .unwrap();
                assert!(target.n_constraints() > 0 && target.n_constraints() < n * n);
                for i in 0..n * n {
                    if target.cover[i] {
                        assert_eq!(target.matrix.data[i], circuit.matrix().data[i]);
                    } else {
                        assert_eq!(target.matrix.data[i].norm_sqr(), 0.0);
                    }
                }
                assert_eq!(target.exact_cost(circuit.matrix(), 1e-6), 0.0);
            }
        }
    }
    #[test]
    fn generated_circuits_meet_simplified_length_requirement() {
        let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
        for seed in 0..12 {
            let mut circuit = generate_circuit(&lib, &mut Rng::new(seed), 3, 1, 0.0, true);
            assert_eq!(circuit.gates.len(), 4);
            assert_eq!(circuit.non_identity(), 4);
            Resynth::default().run(&mut circuit, &lib);
            assert!(circuit.non_identity() >= 3);
        }
    }
}
