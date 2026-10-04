//! Deterministic local rewriting and commuting-gate ordering from C++.
use crate::{circuit::Circuit, gates::GateLibrary};
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Resynth {
    pub max_gate_mult: usize,
    pub watch_depth: bool,
    pub depth_gates: Vec<String>,
    correct_order: bool,
}

pub type Resynthesize = Resynth;

impl Default for Resynth {
    fn default() -> Self {
        Self::new(12, false, vec!["t".into(), "tdg".into()])
    }
}

impl Resynth {
    pub fn new(max_gate_mult: usize, watch_depth: bool, depth_gates: Vec<String>) -> Self {
        Self {
            max_gate_mult,
            watch_depth,
            depth_gates,
            correct_order: false,
        }
    }

    pub fn run(&self, circuit: &mut Circuit, lib: &GateLibrary) {
        // The reference stops its product cache while rewriting. Work directly
        // on gate IDs, then rebuild once, instead of recomputing after swaps.
        let mut gates = circuit.gates.clone();
        let mut commutation = HashMap::new();
        self.pass(&mut gates, lib, false, &mut commutation);
        gates.reverse();
        self.pass(&mut gates, lib, true, &mut commutation);
        gates.reverse();
        *circuit = Circuit::new(gates, circuit.n_qubits, lib);
    }

    /// Remove empty slots and repeat correct-order rewriting to a fixed point.
    /// Each pass is accepted only if it preserves the complete circuit matrix.
    pub fn run_optimized(&self, circuit: &mut Circuit, lib: &GateLibrary) {
        let mut optimizer = self.clone();
        optimizer.correct_order = true;
        for _ in 0..4 {
            let before = circuit.clone();
            let ids = circuit
                .gates
                .iter()
                .copied()
                .filter(|&g| g != lib.identity)
                .collect();
            *circuit = Circuit::new(ids, circuit.n_qubits, lib);
            optimizer.run(circuit, lib);
            if !before.matrix().approximately_equal(circuit.matrix(), 1e-7) {
                *circuit = before;
                break;
            }
            if circuit.non_identity() >= before.non_identity()
                && circuit.cost() >= before.cost() - 1e-10
            {
                break;
            }
        }
    }

    fn pass(
        &self,
        gates: &mut [usize],
        lib: &GateLibrary,
        second: bool,
        cache: &mut HashMap<(usize, usize), bool>,
    ) {
        let mut index = 0;
        while index + 1 < gates.len() {
            let change = self.change(gates, index, lib, second, cache);
            if change < 0 {
                index = index.saturating_sub((-change) as usize);
            } else {
                index += change as usize;
            }
        }
    }

    fn change(
        &self,
        gates: &mut [usize],
        index: usize,
        lib: &GateLibrary,
        second: bool,
        cache: &mut HashMap<(usize, usize), bool>,
    ) -> isize {
        let id_name = &lib.gates[lib.identity].name;
        let first = gates[index];
        let second_id = gates[index + 1];
        if lib.gates[second_id].name == *id_name {
            return 1;
        }
        if !second {
            let mut product = lib.gates[second_id].matrix.clone();
            let mut cost = lib.gates[second_id].cost;
            for extra in 0..self.max_gate_mult.saturating_sub(1) {
                if extra > index {
                    break;
                }
                let earlier = &lib.gates[gates[index - extra]];
                if earlier.name == *id_name {
                    continue;
                }
                // Preserve the reference multiplication order, including its
                // nonstandard backwards accumulation for noncommuting gates.
                product = if self.correct_order {
                    product.mul(&earlier.matrix)
                } else {
                    earlier.matrix.mul(&product)
                };
                cost += earlier.cost;
                for &replacement in &lib.all {
                    let gate = &lib.gates[replacement];
                    if cost > gate.cost && product.approximately_equal(&gate.matrix, 1e-6) {
                        gates[index - extra] = replacement;
                        gates[index - extra + 1..index + 2].fill(lib.identity);
                        return -1 - extra as isize;
                    }
                }
            }
        }
        if lib.gates[first].name == *id_name {
            gates.swap(index, index + 1);
            return -1;
        }
        if commutes(first, second_id, lib, cache)
            && self.priority_change(gates, index, lib, second, cache)
        {
            gates.swap(index, index + 1);
            return -1;
        }
        1
    }

    fn priority_change(
        &self,
        gates: &mut [usize],
        index: usize,
        lib: &GateLibrary,
        second: bool,
        cache: &mut HashMap<(usize, usize), bool>,
    ) -> bool {
        if self.watch_depth && second {
            let original = gate_depth(gates, lib, &self.depth_gates);
            gates.swap(index, index + 1);
            let candidate = gate_depth(gates, lib, &self.depth_gates);
            gates.swap(index, index + 1);
            if candidate > original {
                return false;
            }
            if candidate < original {
                return true;
            }
        }
        let commuting1 = before_commutes(gates, index, lib, cache) as isize;
        let commuting2 = before_commutes(gates, index + 1, lib, cache) as isize;
        if commuting1 < commuting2 - 1 {
            return true;
        }
        if commuting1 > commuting2 - 1 {
            return false;
        }
        let first = &lib.gates[gates[index]];
        let second = &lib.gates[gates[index + 1]];
        (first.qubits.len(), &first.name, &first.qubits)
            > (second.qubits.len(), &second.name, &second.qubits)
    }
}

fn commutes(
    a: usize,
    b: usize,
    lib: &GateLibrary,
    cache: &mut HashMap<(usize, usize), bool>,
) -> bool {
    let key = (a.min(b), a.max(b));
    *cache.entry(key).or_insert_with(|| {
        let a = &lib.gates[a].matrix;
        let b = &lib.gates[b].matrix;
        // C++ uses exact matrix equality here, not Utils::matricesEqual.
        a.mul(b).data == b.mul(a).data
    })
}

fn before_commutes(
    gates: &[usize],
    index: usize,
    lib: &GateLibrary,
    cache: &mut HashMap<(usize, usize), bool>,
) -> usize {
    let mut count = 0;
    for i in (0..index).rev() {
        if !commutes(gates[index], gates[i], lib, cache) {
            break;
        }
        count += 1;
    }
    count
}

fn gate_depth(gates: &[usize], lib: &GateLibrary, names: &[String]) -> usize {
    let mut depths = vec![0; lib.n_qubits];
    for &id in gates {
        let gate = &lib.gates[id];
        let max = gate.qubits.iter().map(|&q| depths[q]).max().unwrap_or(0);
        let increment = usize::from(names.contains(&gate.name));
        for &q in &gate.qubits {
            depths[q] = max + increment;
        }
    }
    depths.into_iter().max().unwrap_or(0)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn lib() -> GateLibrary {
        GateLibrary::load(2, "data/gates/CliffordT", "").unwrap()
    }
    #[test]
    fn inverse_pairs_cancel_and_repeated_t_gates_combine() {
        let lib = lib();
        for pair in [["h", "h"], ["s", "sdg"], ["t", "tdg"]] {
            let gates = pair.map(|name| lib.find(name, &[0]).unwrap()).to_vec();
            let mut circuit = Circuit::new(gates, 2, &lib);
            let before = circuit.matrix().clone();
            Resynth::default().run(&mut circuit, &lib);
            assert_eq!(circuit.non_identity(), 0);
            assert!(circuit.matrix().approximately_equal(&before, 1e-6));
        }
        let t = lib.find("t", &[0]).unwrap();
        let mut circuit = Circuit::new(vec![t, t], 2, &lib);
        Resynth::default().run(&mut circuit, &lib);
        assert_eq!(circuit.non_identity(), 1);
        assert_eq!(circuit.count(&["s".into()], &lib), 1);
    }
    #[test]
    fn cancellation_crosses_commuting_gates_and_identity_slots() {
        let lib = lib();
        let h = lib.find("h", &[0]).unwrap();
        let t = lib.find("t", &[1]).unwrap();
        let mut circuit = Circuit::new(vec![lib.identity, h, t, lib.identity, h], 2, &lib);
        let before = circuit.matrix().clone();
        Resynth::new(12, true, vec!["t".into(), "tdg".into()]).run(&mut circuit, &lib);
        assert_eq!(circuit.non_identity(), 1);
        assert!(circuit.matrix().approximately_equal(&before, 1e-6));
    }
    #[test]
    fn empty_and_single_gate_rewrites_are_well_defined() {
        let lib = lib();
        for gates in [
            vec![],
            vec![lib.find("h", &[0]).unwrap()],
            vec![lib.identity],
        ] {
            let mut circuit = Circuit::new(gates.clone(), 2, &lib);
            Resynth::default().run(&mut circuit, &lib);
            assert_eq!(circuit.gates, gates);
        }
    }
    #[test]
    fn depth_helper_matches_circuit_depth() {
        let lib = lib();
        let gates = vec![
            lib.find("t", &[0]).unwrap(),
            lib.find("cx", &[0, 1]).unwrap(),
            lib.find("tdg", &[1]).unwrap(),
        ];
        let names = vec!["t".into(), "tdg".into()];
        let circuit = Circuit::new(gates.clone(), 2, &lib);
        assert_eq!(
            gate_depth(&gates, &lib, &names),
            circuit.depth(&names, &lib)
        );
        assert_eq!(gate_depth(&gates, &lib, &names), 2);
    }
    #[test]
    fn optimized_rewrite_uses_chronological_product_and_preserves_unitary() {
        let mut lib = lib();
        let h = lib.find("h", &[0]).unwrap();
        let s = lib.find("s", &[0]).unwrap();
        let combined = crate::gates::Gate {
            name: "sh".into(),
            matrix: lib.gates[s].matrix.mul(&lib.gates[h].matrix),
            qubits: vec![0],
            cost: 0.001,
            decomposition: vec![],
        };
        let id = lib.gates.len();
        lib.gates.push(combined);
        lib.all.push(id);
        let mut c = Circuit::new(vec![lib.identity, h, lib.identity, s], 2, &lib);
        let before = c.matrix().clone();
        Resynth::default().run_optimized(&mut c, &lib);
        assert!(before.approximately_equal(c.matrix(), 1e-12));
        assert_eq!(c.non_identity(), 1);
    }
}
