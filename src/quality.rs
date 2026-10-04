//! Exact T-depth proposals using max-plus prefix/suffix contractions.
use crate::gates::GateLibrary;

/// Row i lists (input qubit j, added T depth) for the longest causal path j -> i.
pub struct DepthModel {
    n: usize,
    transitions: Vec<Vec<Vec<(usize, usize)>>>,
}
impl DepthModel {
    pub fn new(lib: &GateLibrary, expand: bool) -> Self {
        let n = lib.n_qubits;
        let transitions = lib
            .gates
            .iter()
            .enumerate()
            .map(|(id, gate)| {
                let sequence: &[usize] = if expand && !gate.decomposition.is_empty() {
                    &gate.decomposition
                } else {
                    std::slice::from_ref(&id)
                };
                let mut paths = vec![vec![None; n]; n];
                for (q, row) in paths.iter_mut().enumerate() {
                    row[q] = Some(0);
                }
                for &id in sequence {
                    let g = &lib.gates[id];
                    let weight = usize::from(matches!(g.name.as_str(), "t" | "tdg"));
                    let row: Vec<_> = (0..n)
                        .map(|source| {
                            g.qubits
                                .iter()
                                .filter_map(|&q| paths[q][source])
                                .max()
                                .map(|x| x + weight)
                        })
                        .collect();
                    for &q in &g.qubits {
                        paths[q].clone_from(&row);
                    }
                }
                paths
                    .into_iter()
                    .map(|row| {
                        row.into_iter()
                            .enumerate()
                            .filter_map(|(q, w)| w.map(|v| (q, v)))
                            .collect()
                    })
                    .collect()
            })
            .collect();
        Self { n, transitions }
    }
    fn push(&self, gate: usize, prefix: &[usize], out: &mut [usize]) {
        for (dst, row) in out.iter_mut().zip(&self.transitions[gate]) {
            *dst = row.iter().map(|&(q, w)| prefix[q] + w).max().unwrap_or(0);
        }
    }
    fn pull(&self, gate: usize, suffix: &[usize], out: &mut [usize]) {
        out.fill(0);
        for (r, row) in self.transitions[gate].iter().enumerate() {
            for &(q, w) in row {
                out[q] = out[q].max(suffix[r] + w);
            }
        }
    }
    fn combine(&self, gate: usize, prefix: &[usize], suffix: &[usize]) -> usize {
        self.transitions[gate]
            .iter()
            .enumerate()
            .flat_map(|(r, row)| row.iter().map(move |&(q, w)| prefix[q] + w + suffix[r]))
            .max()
            .unwrap_or(0)
    }
    pub fn depth(&self, gates: &[usize]) -> usize {
        let mut prefix = vec![0; self.n];
        let mut scratch = prefix.clone();
        for &id in gates {
            self.push(id, &prefix, &mut scratch);
            std::mem::swap(&mut prefix, &mut scratch);
        }
        prefix.into_iter().max().unwrap_or(0)
    }
    pub fn workspace(&self, slots: usize) -> DepthWorkspace {
        DepthWorkspace {
            n: self.n,
            fixed: vec![0; (slots + 1) * self.n],
            moving: vec![0; self.n],
            scratch: vec![0; self.n],
            forward: true,
        }
    }
}

pub struct DepthWorkspace {
    n: usize,
    fixed: Vec<usize>,
    moving: Vec<usize>,
    scratch: Vec<usize>,
    forward: bool,
}
impl DepthWorkspace {
    pub fn reset(&mut self, model: &DepthModel, gates: &[usize], forward: bool) {
        self.forward = forward;
        self.fixed.fill(0);
        self.moving.fill(0);
        if forward {
            for p in (0..gates.len()).rev() {
                let (before, after) = self.fixed.split_at_mut((p + 1) * self.n);
                model.pull(gates[p], &after[..self.n], &mut before[p * self.n..]);
            }
        } else {
            for (p, &gate) in gates.iter().enumerate() {
                let (before, after) = self.fixed.split_at_mut((p + 1) * self.n);
                model.push(gate, &before[p * self.n..], &mut after[..self.n]);
            }
        }
    }
    pub fn candidate(&self, model: &DepthModel, position: usize, gate: usize) -> usize {
        if self.forward {
            model.combine(
                gate,
                &self.moving,
                &self.fixed[(position + 1) * self.n..(position + 2) * self.n],
            )
        } else {
            model.combine(
                gate,
                &self.fixed[position * self.n..(position + 1) * self.n],
                &self.moving,
            )
        }
    }
    pub fn swapped(
        &mut self,
        model: &DepthModel,
        position: usize,
        current: usize,
        neighbor: usize,
    ) -> usize {
        if self.forward {
            model.push(neighbor, &self.moving, &mut self.scratch);
            model.combine(
                current,
                &self.scratch,
                &self.fixed[(position + 2) * self.n..(position + 3) * self.n],
            )
        } else {
            model.push(
                current,
                &self.fixed[(position - 1) * self.n..position * self.n],
                &mut self.scratch,
            );
            model.combine(neighbor, &self.scratch, &self.moving)
        }
    }
    pub fn advance(&mut self, model: &DepthModel, gate: usize) {
        if self.forward {
            model.push(gate, &self.moving, &mut self.scratch);
        } else {
            model.pull(gate, &self.moving, &mut self.scratch);
        }
        std::mem::swap(&mut self.moving, &mut self.scratch);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{circuit::Circuit, rng::Rng, search::random_circuit};

    #[test]
    fn candidate_and_swap_depths_match_complete_circuit_scheduling() {
        for nq in 1..=4 {
            let mut lib =
                GateLibrary::load(nq, "data/gates/CliffordT", "data/gates/composite_rccx").unwrap();
            lib.add_clifford_macros();
            for expand in [false, true] {
                let model = DepthModel::new(&lib, expand);
                let mut rng = Rng::new(573);
                let mut reference = Circuit::new(vec![lib.identity], nq, &lib);
                let mut reference_depth = |ids: &[usize]| {
                    reference.gates = ids
                        .iter()
                        .flat_map(|&id| {
                            if expand && !lib.gates[id].decomposition.is_empty() {
                                lib.gates[id].decomposition.clone()
                            } else {
                                vec![id]
                            }
                        })
                        .collect();
                    // Circuit::depth is independent of matrix/product caches.
                    reference.depth(&["t".into(), "tdg".into()], &lib)
                };
                for _ in 0..8 {
                    for forward in [true, false] {
                        let mut gates = random_circuit(&lib, &mut rng, 11, 0.1).gates;
                        assert_eq!(model.depth(&gates), reference_depth(&gates));
                        let mut work = model.workspace(gates.len());
                        work.reset(&model, &gates, forward);
                        for offset in 0..gates.len() {
                            let pos = if forward {
                                offset
                            } else {
                                gates.len() - 1 - offset
                            };
                            for &candidate in &lib.all {
                                gates[pos] = candidate;
                                assert_eq!(
                                    work.candidate(&model, pos, candidate),
                                    reference_depth(&gates)
                                );
                            }
                            gates[pos] = lib.all[rng.usize(lib.all.len())];
                            if offset + 1 < gates.len() && rng.random01() < 0.5 {
                                let next = if forward { pos + 1 } else { pos - 1 };
                                let depth = work.swapped(&model, pos, gates[pos], gates[next]);
                                gates.swap(pos, next);
                                assert_eq!(depth, reference_depth(&gates));
                            }
                            work.advance(&model, gates[pos]);
                        }
                    }
                }
            }
        }
    }
}
