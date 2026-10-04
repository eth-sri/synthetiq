//! Fixed-width circuits and incremental binary matrix-product caching.
use crate::gates::{is_qasm_header, parse_gate_line, qasm_qubit_count, GateLibrary};
use crate::matrix::Matrix;

#[derive(Clone, Debug)]
pub struct Circuit {
    pub gates: Vec<usize>,
    pub n_qubits: usize,
    // Pair products are stored in reverse temporal order, like the original
    // BinaryMatrixComputer, including its grouping for odd circuit lengths.
    tree: Vec<Vec<Matrix>>,
    active: Vec<Vec<usize>>,
    cost: f64,
    non_identity: usize,
}
/// Reusable storage for the O(log L) cache nodes touched by one mutation.
pub struct MutationBackup {
    nodes: Vec<Matrix>,
    active: Vec<usize>,
    position: usize,
    gate: usize,
    cost: f64,
    non_identity: usize,
}
impl Circuit {
    pub fn new(gates: Vec<usize>, n_qubits: usize, lib: &GateLibrary) -> Self {
        assert_eq!(n_qubits, lib.n_qubits);
        let cost = gates.iter().map(|&id| lib.gates[id].cost).sum();
        let non_identity = gates
            .iter()
            .filter(|&&id| lib.gates[id].name != "id")
            .count();
        let mut result = Self {
            gates,
            n_qubits,
            tree: Vec::new(),
            active: Vec::new(),
            cost,
            non_identity,
        };
        result.build_tree(lib);
        result
    }
    pub fn from_qasm(text: &str, lib: &GateLibrary) -> Result<Self, String> {
        let n_qubits = qasm_qubit_count(text)?;
        if n_qubits != lib.n_qubits {
            return Err(format!(
                "QASM uses {n_qubits} qubits, but gate library uses {}",
                lib.n_qubits
            ));
        }
        let mut gates = Vec::new();
        for (line_number, line) in text.lines().enumerate() {
            if is_qasm_header(line) {
                continue;
            }
            if let Some((name, qubits)) = parse_gate_line(line)? {
                let id = lib.find(&name, &qubits).ok_or_else(|| {
                    format!("Unknown gate on QASM line {}: {line}", line_number + 1)
                })?;
                gates.push(id);
            }
        }
        Ok(Self::new(gates, n_qubits, lib))
    }
    pub fn matrix(&self) -> &Matrix {
        &self.tree.last().unwrap()[0]
    }
    pub fn cost(&self) -> f64 {
        self.cost
    }
    pub fn non_identity(&self) -> usize {
        self.non_identity
    }
    pub fn replace(&mut self, position: usize, id: usize, lib: &GateLibrary) {
        let old = self.gates[position];
        if old == id {
            return;
        }
        self.cost += lib.gates[id].cost - lib.gates[old].cost;
        self.non_identity = self.non_identity + usize::from(lib.gates[id].name != "id")
            - usize::from(lib.gates[old].name != "id");
        self.gates[position] = id;
        let mut node = self.tree[0].len() - 1 - position / 2;
        self.update_pair(position / 2, lib);
        for depth in 1..self.tree.len() {
            node /= 2;
            self.update_node(depth, node);
        }
    }
    pub fn mutation_backup(&self) -> MutationBackup {
        MutationBackup {
            nodes: self
                .tree
                .iter()
                .map(|level| Matrix::zero(level[0].n))
                .collect(),
            active: vec![0; self.tree.len()],
            position: 0,
            gate: 0,
            cost: 0.0,
            non_identity: 0,
        }
    }
    pub fn replace_recorded(
        &mut self,
        position: usize,
        id: usize,
        lib: &GateLibrary,
        backup: &mut MutationBackup,
    ) {
        backup.position = position;
        backup.gate = self.gates[position];
        backup.cost = self.cost;
        backup.non_identity = self.non_identity;
        let mut node = self.tree[0].len() - 1 - position / 2;
        for depth in 0..self.tree.len() {
            backup.nodes[depth]
                .data
                .copy_from_slice(&self.tree[depth][node].data);
            backup.active[depth] = self.active[depth][node];
            node /= 2;
        }
        self.replace(position, id, lib);
    }
    pub fn restore_recorded(&mut self, backup: &MutationBackup) {
        self.gates[backup.position] = backup.gate;
        self.cost = backup.cost;
        self.non_identity = backup.non_identity;
        let mut node = self.tree[0].len() - 1 - backup.position / 2;
        for depth in 0..self.tree.len() {
            self.tree[depth][node]
                .data
                .copy_from_slice(&backup.nodes[depth].data);
            self.active[depth][node] = backup.active[depth];
            node /= 2;
        }
    }
    pub fn rebuild(&mut self, lib: &GateLibrary) {
        self.cost = self.gates.iter().map(|&id| lib.gates[id].cost).sum();
        self.non_identity = self
            .gates
            .iter()
            .filter(|&&id| lib.gates[id].name != "id")
            .count();
        self.build_tree(lib);
    }
    pub fn count(&self, names: &[String], lib: &GateLibrary) -> usize {
        self.gates
            .iter()
            .filter(|&&id| names.contains(&lib.gates[id].name))
            .count()
    }
    /// Non-counted gates still synchronize the depths of their acting qubits.
    pub fn depth(&self, names: &[String], lib: &GateLibrary) -> usize {
        let mut depths = vec![0; self.n_qubits];
        for &id in &self.gates {
            let gate = &lib.gates[id];
            let depth = gate.qubits.iter().map(|&q| depths[q]).max().unwrap_or(0)
                + usize::from(names.contains(&gate.name));
            for &q in &gate.qubits {
                depths[q] = depth;
            }
        }
        depths.into_iter().max().unwrap_or(0)
    }
    pub fn qasm(&self, lib: &GateLibrary) -> String {
        use std::fmt::Write;
        let mut result = format!(
            "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg qubits[{}];\n",
            self.n_qubits
        );
        for &id in &self.gates {
            let gate = &lib.gates[id];
            if gate.name == "id" {
                continue;
            }
            result.push_str(&gate.name);
            for (i, q) in gate.qubits.iter().enumerate() {
                write!(result, " qubits[{q}]").unwrap();
                if i + 1 < gate.qubits.len() {
                    result.push(',');
                }
            }
            result.push_str(";\n");
        }
        result
    }
    pub fn inverse(&mut self, lib: &GateLibrary) -> Result<(), String> {
        let inverse: Result<Vec<_>, _> = self
            .gates
            .iter()
            .rev()
            .map(|&id| {
                lib.inverse_gate(id).ok_or_else(|| {
                    format!(
                        "No conjugate gate available for {} {:?}",
                        lib.gates[id].name, lib.gates[id].qubits
                    )
                })
            })
            .collect();
        self.gates = inverse?;
        self.rebuild(lib);
        Ok(())
    }
    pub fn permute(&mut self, order: &[usize], lib: &GateLibrary) -> Result<(), String> {
        if order.len() != self.n_qubits {
            return Err("Qubit permutation has wrong length".into());
        }
        let mut sorted = order.to_vec();
        sorted.sort_unstable();
        if sorted != (0..self.n_qubits).collect::<Vec<_>>() {
            return Err("Invalid qubit permutation".into());
        }
        let gates: Result<Vec<_>, _> = self
            .gates
            .iter()
            .map(|&id| {
                if lib.gates[id].name == "id" {
                    return Ok(id);
                }
                let gate = &lib.gates[id];
                let qubits: Vec<_> = gate.qubits.iter().map(|&q| order[q]).collect();
                lib.find_all(&gate.name, &qubits).ok_or_else(|| {
                    format!("No permuted gate available for {} {:?}", gate.name, qubits)
                })
            })
            .collect();
        self.gates = gates?;
        self.rebuild(lib);
        Ok(())
    }
    pub fn expand(&mut self, lib: &GateLibrary) {
        self.gates = self
            .gates
            .iter()
            .flat_map(|&id| {
                let gate = &lib.gates[id];
                if gate.decomposition.is_empty() {
                    vec![id]
                } else {
                    gate.decomposition.clone()
                }
            })
            .collect();
        self.rebuild(lib);
    }
    pub fn rotate(&mut self, lib: &GateLibrary) {
        self.gates.reverse();
        self.rebuild(lib);
    }
    fn build_tree(&mut self, lib: &GateLibrary) {
        let n = 1 << self.n_qubits;
        self.tree.clear();
        self.active.clear();
        let mut size = self.gates.len().div_ceil(2).max(1);
        loop {
            self.tree
                .push((0..size).map(|_| Matrix::identity(n)).collect());
            self.active.push(vec![0; size]);
            if size == 1 {
                break;
            }
            size = size.div_ceil(2);
        }
        for pair in 0..self.gates.len().div_ceil(2) {
            self.update_pair(pair, lib);
        }
        for depth in 1..self.tree.len() {
            for node in 0..self.tree[depth].len() {
                self.update_node(depth, node);
            }
        }
    }
    fn update_pair(&mut self, pair: usize, lib: &GateLibrary) {
        let position = pair * 2;
        let node = self.tree[0].len() - 1 - pair;
        let first = self.gates[position];
        let first_active = usize::from(first != lib.identity);
        let out = &mut self.tree[0][node];
        if position + 1 == self.gates.len() {
            out.data.copy_from_slice(&lib.gates[first].matrix.data);
            self.active[0][node] = first_active;
        } else {
            let second = self.gates[position + 1];
            let second_active = usize::from(second != lib.identity);
            multiply_or_copy(
                &lib.gates[second].matrix,
                &lib.gates[first].matrix,
                second_active,
                first_active,
                out,
            );
            self.active[0][node] = first_active + second_active;
        }
    }
    fn update_node(&mut self, depth: usize, node: usize) {
        let (before, after) = self.tree.split_at_mut(depth);
        let previous = &before[depth - 1];
        let out = &mut after[0][node];
        let left = 2 * node;
        let left_active = self.active[depth - 1][left];
        if left + 1 == previous.len() {
            out.data.copy_from_slice(&previous[left].data);
            self.active[depth][node] = left_active;
        } else {
            let right_active = self.active[depth - 1][left + 1];
            multiply_or_copy(
                &previous[left],
                &previous[left + 1],
                left_active,
                right_active,
                out,
            );
            self.active[depth][node] = left_active + right_active;
        }
    }
}
fn multiply_or_copy(a: &Matrix, b: &Matrix, a_active: usize, b_active: usize, out: &mut Matrix) {
    if a_active == 0 {
        out.data.copy_from_slice(&b.data);
    } else if b_active == 0 {
        out.data.copy_from_slice(&a.data);
    } else {
        a.mul_into(b, out);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::complex::Complex;
    use crate::matrix::permutations;
    fn lib(n: usize) -> GateLibrary {
        GateLibrary::load(n, "data/gates/basic_gates", "").unwrap()
    }
    fn linear(c: &Circuit, lib: &GateLibrary) -> Matrix {
        c.gates
            .iter()
            .fold(Matrix::identity(1 << c.n_qubits), |m, &id| {
                lib.gates[id].matrix.mul(&m)
            })
    }
    fn gate(lib: &GateLibrary, name: &str, qubits: &[usize]) -> usize {
        lib.find(name, qubits).unwrap()
    }
    #[test]
    fn empty_and_single_circuits() {
        let lib = lib(2);
        let empty = Circuit::new(Vec::new(), 2, &lib);
        assert_eq!(empty.matrix(), &Matrix::identity(4));
        assert_eq!(empty.cost(), 0.0);
        assert_eq!(empty.non_identity(), 0);
        for &id in &lib.all {
            let c = Circuit::new(vec![id], 2, &lib);
            assert_eq!(c.matrix(), &lib.gates[id].matrix);
        }
    }
    #[test]
    fn cached_products_and_thousands_of_mutations() {
        for n in 1..=4 {
            let lib = lib(n);
            for len in [1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33] {
                let ids = (0..len).map(|i| (i * 13 + 7) % lib.gates.len()).collect();
                let mut c = Circuit::new(ids, n, &lib);
                assert!(c.matrix().max_abs_diff(&linear(&c, &lib)) < 1e-12);
                for step in 0..48 {
                    let pos = (step * 17 + 3) % len;
                    let id = (step * 19 + 11) % lib.gates.len();
                    let old = c.gates[pos];
                    c.replace(pos, id, &lib);
                    assert!(
                        c.matrix().max_abs_diff(&linear(&c, &lib)) < 1e-12,
                        "n={n}, len={len}, step={step}"
                    );
                    assert!(
                        (c.cost() - c.gates.iter().map(|&id| lib.gates[id].cost).sum::<f64>())
                            .abs()
                            < 1e-10
                    );
                    assert_eq!(
                        c.non_identity(),
                        c.gates
                            .iter()
                            .filter(|&&id| lib.gates[id].name != "id")
                            .count()
                    );
                    if step % 4 == 0 {
                        c.replace(pos, old, &lib);
                        assert!(c.matrix().max_abs_diff(&linear(&c, &lib)) < 1e-12);
                    }
                }
            }
        }
    }
    #[test]
    fn bell_and_ghz_state_vectors() {
        for n in 2..=5 {
            let lib = lib(n);
            let mut ids = vec![gate(&lib, "h", &[0])];
            for q in 1..n {
                ids.push(gate(&lib, "cx", &[q - 1, q]));
            }
            let c = Circuit::new(ids, n, &lib);
            for i in 0..1 << n {
                let expected = if i == 0 || i == (1 << n) - 1 {
                    Complex::from(0.5f64.sqrt())
                } else {
                    Complex::ZERO
                };
                assert!((c.matrix()[(i, 0)] - expected).abs() < 1e-15);
            }
        }
    }
    #[test]
    fn qasm_roundtrip_and_identity_elision() {
        let lib = lib(3);
        let text="OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg q[3];\nh q[2];\ncx q[2],q[0];\nt q[1];\nid q[0];\n";
        let c = Circuit::from_qasm(text, &lib).unwrap();
        let parsed = Circuit::from_qasm(&c.qasm(&lib), &lib).unwrap();
        assert_eq!(parsed.gates.len(), 3);
        assert_eq!(parsed.matrix(), c.matrix());
        assert_eq!(c.cost(), parsed.cost());
        assert!(Circuit::from_qasm(&text.replace("qreg q[3]", "qreg q[2]"), &lib).is_err());
    }
    #[test]
    fn inverse_and_all_qubit_permutations() {
        let lib = lib(3);
        let ids = vec![
            gate(&lib, "h", &[0]),
            gate(&lib, "cx", &[0, 2]),
            gate(&lib, "t", &[1]),
            gate(&lib, "sdg", &[2]),
            0,
        ];
        let c = Circuit::new(ids, 3, &lib);
        let mut inverse = c.clone();
        inverse.inverse(&lib).unwrap();
        assert!(inverse.matrix().max_abs_diff(&c.matrix().adjoint()) < 1e-14);
        for p in permutations(3) {
            let mut permuted = c.clone();
            permuted.permute(&p, &lib).unwrap();
            assert!(permuted.matrix().max_abs_diff(&c.matrix().permuted(&p)) < 1e-14);
        }
    }
    #[test]
    fn weighted_depth_synchronizes_noncounted_gates() {
        let lib = lib(3);
        let ids = vec![
            gate(&lib, "t", &[0]),
            gate(&lib, "t", &[0]),
            gate(&lib, "cx", &[0, 1]),
            gate(&lib, "tdg", &[1]),
            gate(&lib, "t", &[2]),
        ];
        let c = Circuit::new(ids, 3, &lib);
        let names = vec!["t".into(), "tdg".into()];
        assert_eq!(c.count(&names, &lib), 4);
        assert_eq!(c.depth(&names, &lib), 3);
        assert_eq!(c.depth(&[], &lib), 0);
    }
    #[test]
    fn composite_expansion_preserves_matrix_and_cost() {
        let lib =
            GateLibrary::load(3, "data/gates/basic_gates", "data/gates/composite_ccx").unwrap();
        let mut c = Circuit::new(vec![0, lib.composite[0], lib.composite[3]], 3, &lib);
        let old = c.matrix().clone();
        let cost = c.cost();
        c.expand(&lib);
        assert!(old.max_abs_diff(c.matrix()) < 1e-13);
        assert!((cost - c.cost()).abs() < 1e-13);
        assert!(c
            .gates
            .iter()
            .all(|&id| lib.gates[id].decomposition.is_empty()));
    }
    #[test]
    fn rejected_transactions_restore_all_cached_state() {
        for nq in 1..=4 {
            let lib = lib(nq);
            for length in [1, 2, 3, 8, 17, 33] {
                let mut c = Circuit::new(vec![lib.identity; length], nq, &lib);
                let mut backup = c.mutation_backup();
                for step in 0..100 {
                    let position = (step * 17 + 3) % length;
                    let id = lib.all[(step * 7 + 1) % lib.all.len()];
                    let before = c.clone();
                    c.replace_recorded(position, id, &lib, &mut backup);
                    if step % 3 != 0 {
                        c.restore_recorded(&backup);
                        assert_eq!(c.matrix(), before.matrix());
                        assert_eq!(c.gates, before.gates);
                        assert_eq!(c.cost(), before.cost());
                        assert_eq!(c.non_identity(), before.non_identity());
                    }
                    assert!(c.matrix().max_abs_diff(&linear(&c, &lib)) < 1e-12);
                }
            }
        }
    }
}
