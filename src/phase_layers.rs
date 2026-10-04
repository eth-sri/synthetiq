//! Native Clifford basis changes schedule independent parity phases in parallel.
//! Coefficients come from the stochastic phase search; randomized packing only
//! chooses an equivalent gate layout and never changes those coefficients.
use crate::{gates::GateLibrary, quality::DepthModel, rng::Rng};

pub struct PhaseLayout<'a> {
    pub native: &'a [[usize; 5]],
    pub cx: &'a [Vec<Option<usize>>],
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
fn shuffle<T>(values: &mut [T], rng: &mut Rng) {
    for i in (1..values.len()).rev() {
        values.swap(i, rng.usize(i + 1));
    }
}
impl PhaseLayout<'_> {
    /// Row masks define y_i = parity(row_i & x). Return a native CX word for y.
    pub(crate) fn basis_word(&self, mut rows: Vec<usize>) -> Vec<usize> {
        if let Some(word) = crate::linear_clifford::shortest(&rows, self.cx) {
            return word;
        }
        let n = rows.len();
        let mut gates = Vec::new();
        let apply = |a: usize, b: usize, rows: &mut [usize], gates: &mut Vec<usize>| {
            rows[b] ^= rows[a];
            gates.push(self.cx[a][b].unwrap());
        };
        for column in 0..n {
            if rows[column] & (1 << column) == 0 {
                let pivot = (column + 1..n)
                    .find(|&r| rows[r] & (1 << column) != 0)
                    .unwrap();
                apply(column, pivot, &mut rows, &mut gates);
                apply(pivot, column, &mut rows, &mut gates);
                apply(column, pivot, &mut rows, &mut gates);
            }
            for r in 0..n {
                if r != column && rows[r] & (1 << column) != 0 {
                    apply(column, r, &mut rows, &mut gates);
                }
            }
        }
        gates.reverse();
        gates
    }
    fn compress_cx_blocks(&self, gates: &[usize], lib: &GateLibrary) -> Vec<usize> {
        let mut result = Vec::new();
        let mut begin = 0;
        while begin < gates.len() {
            if lib.gates[gates[begin]].name != "cx" {
                result.push(gates[begin]);
                begin += 1;
                continue;
            }
            let mut end = begin;
            let mut rows: Vec<_> = (0..self.native.len()).map(|bit| 1 << bit).collect();
            while end < gates.len() && lib.gates[gates[end]].name == "cx" {
                let qubits = &lib.gates[gates[end]].qubits;
                rows[qubits[1]] ^= rows[qubits[0]];
                end += 1;
            }
            let replacement = self.basis_word(rows);
            if replacement.len() < end - begin {
                result.extend(replacement);
            } else {
                result.extend_from_slice(&gates[begin..end]);
            }
            begin = end;
        }
        result
    }
    fn append_phase(&self, gates: &mut Vec<usize>, wire: usize, phase: u8) {
        let ids = self.native[wire];
        match phase {
            0 => (),
            1 => gates.push(ids[3]),
            2 => gates.push(ids[1]),
            3 => gates.extend([ids[1], ids[3]]),
            4 => gates.extend([ids[1], ids[1]]),
            5 => gates.extend([ids[2], ids[4]]),
            6 => gates.push(ids[2]),
            7 => gates.push(ids[4]),
            _ => unreachable!(),
        }
    }
    fn append_groups(
        &self,
        gates: &mut Vec<usize>,
        masks: &[usize],
        coefficients: &[u8],
        rng: &mut Rng,
    ) {
        let n = self.native.len();
        let mut groups: Vec<(Vec<usize>, Vec<usize>)> = Vec::new();
        for &mask in masks {
            let group_index = groups.iter().position(|(_, basis)| {
                let mut copy = basis.clone();
                insert(&mut copy, mask)
            });
            if let Some(index) = group_index {
                groups[index].0.push(mask);
                assert!(insert(&mut groups[index].1, mask));
            } else {
                let mut basis = vec![0; n];
                assert!(insert(&mut basis, mask));
                groups.push((vec![mask], basis));
            }
        }
        for (group, mut basis) in groups {
            let mut extension: Vec<_> = (0..n).map(|bit| 1 << bit).collect();
            shuffle(&mut extension, rng);
            let mut complete = group.clone();
            for vector in extension {
                if insert(&mut basis, vector) {
                    complete.push(vector);
                }
            }
            let mut wires: Vec<_> = (0..n).collect();
            shuffle(&mut wires, rng);
            let mut rows = vec![0; n];
            for (i, &mask) in complete.iter().enumerate() {
                rows[wires[i]] = mask;
            }
            let word = self.basis_word(rows);
            gates.extend(&word);
            for (i, &mask) in group.iter().enumerate() {
                self.append_phase(gates, wires[i], coefficients[mask - 1]);
            }
            gates.extend(word.into_iter().rev());
        }
    }
    pub fn optimize(
        &self,
        coefficients: &[u8],
        wrappers: (&[usize], &[usize]),
        lib: &GateLibrary,
        rng: &mut Rng,
        depth_limit: usize,
    ) -> Vec<usize> {
        let (prefix, suffix) = wrappers;
        let model = DepthModel::new(lib, false);
        let mut odd: Vec<_> = (1..=coefficients.len())
            .filter(|&m| coefficients[m - 1] % 2 == 1)
            .collect();
        let mut even: Vec<_> = (1..=coefficients.len())
            .filter(|&m| coefficients[m - 1] != 0 && coefficients[m - 1].is_multiple_of(2))
            .collect();
        let mut best = Vec::new();
        let mut best_key = (usize::MAX, usize::MAX);
        for attempt in 0..128 {
            shuffle(&mut odd, rng);
            shuffle(&mut even, rng);
            let mut gates = prefix.to_vec();
            if attempt % 2 == 0 {
                self.append_groups(&mut gates, &even, coefficients, rng);
            }
            self.append_groups(&mut gates, &odd, coefficients, rng);
            if attempt % 2 != 0 {
                self.append_groups(&mut gates, &even, coefficients, rng);
            }
            gates.extend(suffix);
            let key = (model.depth(&gates), gates.len());
            if key < best_key {
                best_key = key;
                best = gates;
            }
            if best_key.0 <= depth_limit {
                break;
            }
        }
        let compact = self.compress_cx_blocks(&best, lib);
        if compact.len() < best.len() && model.depth(&compact) <= best_key.0 {
            compact
        } else {
            best
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{circuit::Circuit, complex::Complex, matrix::Matrix};
    #[test]
    fn randomized_parallel_layouts_preserve_every_phase_and_t_count() {
        let mut rng = Rng::new(37891);
        for n in 1..=4 {
            let lib = GateLibrary::load(n, "data/gates/CliffordT", "").unwrap();
            let native: Vec<_> = (0..n)
                .map(|q| ["h", "s", "sdg", "t", "tdg"].map(|name| lib.find(name, &[q]).unwrap()))
                .collect();
            let cx: Vec<_> = (0..n)
                .map(|a| (0..n).map(|b| lib.find("cx", &[a, b])).collect())
                .collect();
            let layout = PhaseLayout {
                native: &native,
                cx: &cx,
            };
            for _ in 0..32 {
                let coeff: Vec<_> = (1..1 << n).map(|_| rng.usize(8) as u8).collect();
                let gates = layout.optimize(&coeff, (&[], &[]), &lib, &mut rng, n);
                let c = Circuit::new(gates, n, &lib);
                let mut expected = Matrix::zero(1 << n);
                for x in 0usize..1 << n {
                    let phase: usize = (1usize..1 << n)
                        .filter(|&m| (x & m).count_ones() % 2 == 1)
                        .map(|m| coeff[m - 1] as usize)
                        .sum();
                    let angle = std::f64::consts::FRAC_PI_4 * (phase % 8) as f64;
                    expected[(x, x)] = Complex::new(angle.cos(), angle.sin());
                }
                assert!(c.matrix().max_abs_diff(&expected) < 1e-12);
                assert_eq!(
                    c.count(&["t".into(), "tdg".into()], &lib),
                    coeff.iter().filter(|&&c| c % 2 == 1).count()
                );
            }
        }
    }
    #[test]
    fn all_seven_three_qubit_parities_fit_three_t_layers() {
        let lib = GateLibrary::load(3, "data/gates/CliffordT", "").unwrap();
        let native: Vec<_> = (0..3)
            .map(|q| ["h", "s", "sdg", "t", "tdg"].map(|name| lib.find(name, &[q]).unwrap()))
            .collect();
        let cx: Vec<_> = (0..3)
            .map(|a| (0..3).map(|b| lib.find("cx", &[a, b])).collect())
            .collect();
        let layout = PhaseLayout {
            native: &native,
            cx: &cx,
        };
        for seed in 0..32 {
            let gates = layout.optimize(&[1; 7], (&[], &[]), &lib, &mut Rng::new(seed), 3);
            assert_eq!(DepthModel::new(&lib, false).depth(&gates), 3);
        }
    }
}
