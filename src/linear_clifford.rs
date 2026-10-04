//! Small Clifford rewrites: shortest native CX words for linear reversible maps.
//! The table is a reusable finite rewrite library, never a target synthesis path.
use crate::{circuit::Circuit, gates::GateLibrary, phase::NativeClifford};
use std::sync::OnceLock;

struct Table {
    parent: Vec<u16>,
    step: Vec<u8>,
    identity: usize,
}
impl Table {
    fn new(q: usize) -> Self {
        if q == 4 {
            let bytes = include_bytes!(concat!(env!("OUT_DIR"), "/linear4-parent.bin"));
            let parent = bytes
                .chunks_exact(2)
                .map(|pair| u16::from_le_bytes([pair[0], pair[1]]))
                .collect();
            let step = include_bytes!(concat!(env!("OUT_DIR"), "/linear4-step.bin")).to_vec();
            let identity = (0..q).fold(0, |state, bit| state | (1 << (bit * q + bit)));
            Self {
                parent,
                step,
                identity,
            }
        } else {
            Self::generate(q)
        }
    }
    fn generate(q: usize) -> Self {
        let identity = (0..q).fold(0, |state, bit| state | (1 << (bit * q + bit)));
        let mut parent = vec![u16::MAX; 1 << (q * q)];
        let mut step = vec![0; parent.len()];
        let mut queue = Vec::new();
        parent[identity] = identity as u16;
        queue.push(identity);
        let mask = (1 << q) - 1;
        let mut head = 0;
        while head < queue.len() {
            let state = queue[head];
            head += 1;
            for control in 0..q {
                let source = (state >> (control * q)) & mask;
                for target in 0..q {
                    if control == target {
                        continue;
                    }
                    let next = state ^ (source << (target * q));
                    if parent[next] == u16::MAX {
                        parent[next] = state as u16;
                        step[next] = (control * q + target) as u8;
                        queue.push(next);
                    }
                }
            }
        }
        Self {
            parent,
            step,
            identity,
        }
    }
    fn word(&self, mut state: usize, q: usize, cx: &[Vec<Option<usize>>]) -> Option<Vec<usize>> {
        if self.parent[state] == u16::MAX {
            return None;
        }
        let mut result = Vec::new();
        while state != self.identity {
            let code = self.step[state] as usize;
            result.push(cx[code / q][code % q]?);
            state = self.parent[state] as usize;
        }
        result.reverse();
        Some(result)
    }
}

pub(crate) fn shortest(rows: &[usize], cx: &[Vec<Option<usize>>]) -> Option<Vec<usize>> {
    static TABLES: [OnceLock<Table>; 4] = [const { OnceLock::new() }; 4];
    let q = rows.len();
    if !(1..=4).contains(&q) || rows.iter().any(|&row| row >= 1 << q) {
        return None;
    }
    let state = rows
        .iter()
        .enumerate()
        .fold(0, |key, (r, &row)| key | (row << (r * q)));
    TABLES[q - 1]
        .get_or_init(|| Table::new(q))
        .word(state, q, cx)
}

struct SingleTable {
    words: Vec<Vec<usize>>,
    next: Vec<[usize; 3]>,
}
impl SingleTable {
    fn new() -> Self {
        // Signed Pauli images of X and Z determine a Clifford up to global phase.
        let transform = |state: [i8; 2], gate: usize| {
            let images = match gate {
                0 => [3, -2, 1],
                1 => [2, -1, 3],
                _ => [-2, 1, 3],
            };
            state.map(|axis| axis.signum() * images[axis.unsigned_abs() as usize - 1])
        };
        let mut states = vec![[1, 3]];
        let mut words = vec![Vec::new()];
        let mut head = 0;
        while head < states.len() {
            for gate in 0..3 {
                let next = transform(states[head], gate);
                if !states.contains(&next) {
                    states.push(next);
                    let mut word = words[head].clone();
                    word.push(gate);
                    words.push(word);
                }
            }
            head += 1;
        }
        assert_eq!(states.len(), 24);
        let next = states
            .iter()
            .map(|&state| {
                std::array::from_fn(|gate| {
                    states
                        .iter()
                        .position(|&other| other == transform(state, gate))
                        .unwrap()
                })
            })
            .collect();
        Self { words, next }
    }
}

fn single_words(input: &[usize], lib: &GateLibrary, native: &NativeClifford) -> Vec<usize> {
    static TABLE: OnceLock<SingleTable> = OnceLock::new();
    let table = TABLE.get_or_init(SingleTable::new);
    let mut pending = vec![Vec::new(); lib.n_qubits];
    let mut states = vec![0; lib.n_qubits];
    let mut output = Vec::with_capacity(input.len());
    let flush =
        |bit: usize, pending: &mut [Vec<usize>], states: &mut [usize], output: &mut Vec<usize>| {
            let word: Vec<_> = table.words[states[bit]]
                .iter()
                .map(|&gate| native.native[bit][gate])
                .collect();
            if word.len() <= pending[bit].len()
                && word.iter().map(|&id| lib.gates[id].cost).sum::<f64>()
                    <= pending[bit]
                        .iter()
                        .map(|&id| lib.gates[id].cost)
                        .sum::<f64>()
                        + 1e-12
            {
                output.extend(word);
                pending[bit].clear();
            } else {
                output.append(&mut pending[bit]);
            }
            states[bit] = 0;
        };
    for &id in input {
        if id == lib.identity {
            continue;
        }
        let gate = &lib.gates[id];
        if gate.qubits.len() == 1 {
            let bit = gate.qubits[0];
            if let Some(code) = native.native[bit][..3].iter().position(|&g| g == id) {
                states[bit] = table.next[states[bit]][code];
                pending[bit].push(id);
                continue;
            }
        }
        for &bit in &gate.qubits {
            flush(bit, &mut pending, &mut states, &mut output);
        }
        output.push(id);
    }
    for bit in 0..lib.n_qubits {
        flush(bit, &mut pending, &mut states, &mut output);
    }
    output
}

fn same_phase(a: &crate::matrix::Matrix, b: &crate::matrix::Matrix) -> bool {
    let z = a.trace_conjugate_product(b);
    if z.abs() < 1e-12 {
        return a.approximately_equal(b, 1e-10);
    }
    let phase = z / z.abs();
    a.data
        .iter()
        .zip(&b.data)
        .all(|(&x, &y)| (y - phase * x).norm_sqr() < 1e-20)
}

/// Compress contiguous CX blocks, accepting only nonincreasing native resource
/// counts and cost. Full-circuit depth accounts for changed dependencies.
pub(crate) fn fold(circuit: &mut Circuit, lib: &GateLibrary) -> bool {
    if lib.n_qubits > 6 {
        return false;
    }
    let Some(native) = NativeClifford::new(lib) else {
        return false;
    };
    let compact = single_words(&circuit.gates, lib, &native);
    let input = &compact;
    let mut output = Vec::with_capacity(input.len());
    let mut begin = 0;
    while begin < input.len() {
        if lib.n_qubits > 4 || lib.gates[input[begin]].name != "cx" {
            output.push(input[begin]);
            begin += 1;
            continue;
        }
        let mut rows: Vec<_> = (0..lib.n_qubits).map(|q| 1 << q).collect();
        let mut end = begin;
        while end < input.len() && lib.gates[input[end]].name == "cx" {
            let bits = &lib.gates[input[end]].qubits;
            rows[bits[1]] ^= rows[bits[0]];
            end += 1;
        }
        let replacement = shortest(&rows, &native.cx).expect("invertible native CX block");
        if replacement.len() < end - begin {
            output.extend(replacement);
        } else {
            output.extend_from_slice(&input[begin..end]);
        }
        begin = end;
    }
    if output.len() >= circuit.gates.len() {
        return false;
    }
    let candidate = Circuit::new(output, lib.n_qubits, lib);
    let names = ["t".into(), "tdg".into()];
    if candidate.cost() <= circuit.cost() + 1e-12
        && candidate.depth(&names, lib) <= circuit.depth(&names, lib)
        && same_phase(candidate.matrix(), circuit.matrix())
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
    fn compiled_table_matches_independent_runtime_breadth_first_search() {
        let compiled = Table::new(4);
        let runtime = Table::generate(4);
        assert_eq!(compiled.parent, runtime.parent);
        assert_eq!(compiled.step, runtime.step);
        assert_eq!(compiled.identity, runtime.identity);
    }

    #[test]
    fn single_clifford_compaction_preserves_mixed_circuits_up_to_global_phase() {
        let mut rng = crate::rng::Rng::new(77134);
        let mut improved = 0;
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let native = NativeClifford::new(&lib).unwrap();
            for _ in 0..100 {
                let mut word = Vec::new();
                for _ in 0..30 {
                    let bit = rng.usize(q);
                    word.push(native.native[bit][rng.usize(5)]);
                    if q > 1 && rng.random01() < 0.1 {
                        word.push(native.cx[0][1].unwrap());
                    }
                }
                let bit = rng.usize(q);
                for _ in 0..3 {
                    word.extend([native.native[bit][0], native.native[bit][1]]);
                }
                let before = Circuit::new(word, q, &lib);
                let mut after = before.clone();
                improved += usize::from(fold(&mut after, &lib));
                let names = ["t".into(), "tdg".into()];
                assert!(same_phase(before.matrix(), after.matrix()));
                assert_eq!(after.count(&names, &lib), before.count(&names, &lib));
                assert!(after.depth(&names, &lib) <= before.depth(&names, &lib));
                assert!(after.non_identity() <= before.non_identity());
                assert!(after.cost() <= before.cost() + 1e-12);
            }
        }
        assert!(improved > 300);
    }

    #[test]
    fn tables_cover_every_invertible_map_and_reconstruct_it() {
        for (q, count) in [(1, 1), (2, 6), (3, 168), (4, 20160)] {
            let table = Table::new(q);
            assert_eq!(
                table.parent.iter().filter(|&&p| p != u16::MAX).count(),
                count
            );
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let native = NativeClifford::new(&lib).unwrap();
            for state in 0..table.parent.len() {
                let Some(word) = table.word(state, q, &native.cx) else {
                    continue;
                };
                for input in 0usize..(1 << q) {
                    let mut output = input;
                    for &id in &word {
                        let bits = &lib.gates[id].qubits;
                        output ^= ((output >> bits[0]) & 1) << bits[1];
                    }
                    let expected = (0..q).fold(0, |value, r| {
                        let row = (state >> (r * q)) & ((1 << q) - 1);
                        value | (((row & input).count_ones() as usize % 2) << r)
                    });
                    assert_eq!(output, expected);
                }
            }
        }
    }
    #[test]
    fn block_compression_preserves_all_native_resources() {
        let mut rng = crate::rng::Rng::new(632911);
        for q in 2..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let native = NativeClifford::new(&lib).unwrap();
            for _ in 0..100 {
                let mut word = Vec::new();
                for _ in 0..3 {
                    word.push(native.native[rng.usize(q)][3]);
                    for _ in 0..15 {
                        let a = rng.usize(q);
                        let mut b = rng.usize(q - 1);
                        if b >= a {
                            b += 1;
                        }
                        word.push(native.cx[a][b].unwrap());
                    }
                }
                let before = Circuit::new(word, q, &lib);
                let mut after = before.clone();
                fold(&mut after, &lib);
                let names = ["t".into(), "tdg".into()];
                assert!(after.matrix().approximately_equal(before.matrix(), 1e-10));
                assert_eq!(after.count(&names, &lib), before.count(&names, &lib));
                assert!(after.depth(&names, &lib) <= before.depth(&names, &lib));
                assert!(after.non_identity() <= before.non_identity());
                assert!(after.cost() <= before.cost() + 1e-12);
            }
        }
    }
}
