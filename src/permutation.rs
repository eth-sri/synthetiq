//! Gibbs annealing of monomial circuit skeletons followed by stochastic phase
//! correction. Every skeleton move expands to the supplied native gate set.
use crate::{
    circuit::Circuit,
    complex::Complex,
    gates::GateLibrary,
    matrix::Matrix,
    partial::{PartialMatrix, Target},
    phase::{NativeClifford, PhaseKernel},
    rng::Rng,
    search::{SearchOptions, SearchResult},
};

use std::time::Instant;

#[derive(Clone, Copy)]
struct CorrectionBudget {
    used: usize,
    limit: usize,
    deadline: Option<Instant>,
}
fn expired(deadline: Option<Instant>) -> bool {
    deadline.is_some_and(|end| Instant::now() >= end)
}

struct Monomial {
    map: Vec<usize>,
    inverse: Vec<usize>,
    phases: Vec<Complex>,
    word: Vec<usize>,
    tcount: usize,
}
fn monomial(matrix: &Matrix, word: Vec<usize>, lib: &GateLibrary) -> Option<Monomial> {
    let n = matrix.n;
    let mut map = Vec::new();
    let mut phases = Vec::new();
    let mut inverse = vec![usize::MAX; n];
    for input in 0..n {
        let output = (0..n).find(|&r| (matrix[(r, input)].norm_sqr() - 1.0).abs() < 1e-9)?;
        if inverse[output] != usize::MAX
            || (0..n).any(|r| r != output && matrix[(r, input)].norm_sqr() > 1e-18)
        {
            return None;
        }
        inverse[output] = input;
        map.push(output);
        phases.push(matrix[(output, input)]);
    }
    let tcount = word
        .iter()
        .filter(|&&id| matches!(lib.gates[id].name.as_str(), "t" | "tdg"))
        .count();
    Some(Monomial {
        map,
        inverse,
        phases,
        word,
        tcount,
    })
}
fn same_gate(a: &Monomial, b: &Monomial) -> bool {
    if a.map != b.map {
        return false;
    }
    let phase = a.phases[0] / b.phases[0];
    a.phases
        .iter()
        .zip(&b.phases)
        .all(|(&x, &y)| (x - phase * y).norm_sqr() < 1e-18)
}

pub struct PermutationKernel {
    moves: Vec<Monomial>,
    target_map: Vec<usize>,
    target_phases: Vec<Complex>,
    t_limit: usize,
    budget_floor: Option<usize>,
    qubits: usize,
    native: NativeClifford,
}
impl PermutationKernel {
    pub fn new(target: &Target, lib: &GateLibrary, t_limit: usize) -> Option<Self> {
        if lib.n_qubits > 6
            || target.original.cover.iter().any(|&b| !b)
            || target.inverses.iter().any(|&b| b)
        {
            return None;
        }
        let native = NativeClifford::new(lib)?;
        let specification = monomial(&target.original.matrix, Vec::new(), lib)?;
        let identity = Matrix::identity(1 << lib.n_qubits);
        let mut moves = vec![monomial(&identity, Vec::new(), lib)?];
        for &id in &lib.all {
            let gate = &lib.gates[id];
            let word = if gate.decomposition.is_empty() {
                vec![id]
            } else {
                gate.decomposition.clone()
            };
            let Some(g) = monomial(&gate.matrix, word, lib) else {
                continue;
            };
            if g.tcount > t_limit
                || g.map.iter().enumerate().all(|(x, &y)| x == y)
                || moves.iter().any(|m| same_gate(m, &g))
            {
                continue;
            }
            moves.push(g);
        }
        for ids in &native.native {
            let word = vec![ids[0], ids[1], ids[1], ids[0]];
            let circuit = Circuit::new(word.clone(), lib.n_qubits, lib);
            let g = monomial(circuit.matrix(), word, lib)?;
            if !moves.iter().any(|m| same_gate(m, &g)) {
                moves.push(g);
            }
        }
        // A standard relative-phase controlled-controlled-X is another finite
        // native macro. Include both control orders because their phases differ.
        if t_limit >= 4 {
            for a in 0..lib.n_qubits {
                for b in 0..lib.n_qubits {
                    for to in 0..lib.n_qubits {
                        if a == b || a == to || b == to {
                            continue;
                        }
                        let ids = native.native[to];
                        let word = vec![
                            ids[0],
                            ids[3],
                            native.cx[b][to]?,
                            ids[4],
                            native.cx[a][to]?,
                            ids[3],
                            native.cx[b][to]?,
                            ids[4],
                            ids[0],
                        ];
                        let circuit = Circuit::new(word.clone(), lib.n_qubits, lib);
                        let g = monomial(circuit.matrix(), word, lib)?;
                        debug_assert!((0..g.map.len())
                            .all(|x| g.map[x] == (x ^ ((((x >> a) & 1) & ((x >> b) & 1)) << to))));
                        if !moves.iter().any(|m| same_gate(m, &g)) {
                            moves.push(g);
                        }
                    }
                }
            }
        }
        if !moves.iter().any(|m| m.tcount > 0) {
            return None;
        }
        Some(Self {
            moves,
            target_map: specification.map,
            target_phases: specification.phases,
            t_limit,
            budget_floor: None,
            qubits: lib.n_qubits,
            native,
        })
    }
    pub fn vary_t_budget(mut self) -> Self {
        self.budget_floor = Some(self.t_limit / 2);
        self
    }
    fn word(&self, gates: &[usize]) -> Vec<usize> {
        gates
            .iter()
            .flat_map(|&id| self.moves[id].word.iter().copied())
            .collect()
    }
    fn correct_phases(
        &self,
        gates: &[usize],
        budget: CorrectionBudget,
        target: &Target,
        lib: &GateLibrary,
        rng: &mut Rng,
        options: &SearchOptions,
    ) -> Option<SearchResult> {
        let n = self.target_map.len();
        let mut map: Vec<_> = (0..n).collect();
        let mut phases = vec![Complex::ONE; n];
        for &id in gates {
            let g = &self.moves[id];
            for x in 0..n {
                phases[x] *= g.phases[map[x]];
                map[x] = g.map[map[x]];
            }
        }
        if map != self.target_map {
            return None;
        }
        let remaining = budget.limit.checked_sub(budget.used)?;
        let differences: Vec<_> = (0..n).map(|x| self.target_phases[x] / phases[x]).collect();
        let mut cut_map: Vec<_> = (0..n).collect();
        let mut seen = Vec::new();
        for cut in 0..=gates.len() {
            if expired(budget.deadline) {
                return None;
            }
            if cut > 0 {
                for x in &mut cut_map {
                    *x = self.moves[gates[cut - 1]].map[*x];
                }
            }
            if seen.contains(&cut_map) {
                continue;
            }
            seen.push(cut_map.clone());
            let mut diagonal = vec![Complex::ZERO; n];
            for x in 0..n {
                diagonal[cut_map[x]] = differences[x];
            }
            let Some(kernel) = PhaseKernel::from_diagonal(&diagonal, &self.native, remaining)
            else {
                continue;
            };
            let mut residual = Matrix::zero(n);
            for x in 0..n {
                residual[(x, x)] = diagonal[x];
            }
            let partial =
                PartialMatrix::new(residual, vec![true; n * n], "phase correction").ok()?;
            let phase_target = Target::new(partial, false, false);
            for _ in 0..4 {
                if expired(budget.deadline) {
                    return None;
                }
                let result = kernel.run(&phase_target, lib, rng, options);
                if !result.found {
                    continue;
                }
                let mut word = self.word(&gates[..cut]);
                word.extend(result.circuit.gates);
                word.extend(self.word(&gates[cut..]));
                let circuit = Circuit::new(word, self.qubits, lib);
                if target
                    .original
                    .exact_cost(circuit.matrix(), options.epsilon)
                    != 0.0
                {
                    continue;
                }
                let eq = target.original.cost(circuit.matrix(), options.simple_cost);
                return Some(SearchResult {
                    circuit,
                    best_eq: eq,
                    best_energy: eq.powf(2.0 * options.cost_power),
                    found: true,
                    steps: result.steps,
                    accepted: result.accepted,
                });
            }
        }
        None
    }
    pub fn run(
        &self,
        target: &Target,
        lib: &GateLibrary,
        rng: &mut Rng,
        options: &SearchOptions,
    ) -> SearchResult {
        self.run_until(target, lib, rng, options, None)
    }
    pub(crate) fn run_until(
        &self,
        target: &Target,
        lib: &GateLibrary,
        rng: &mut Rng,
        options: &SearchOptions,
        deadline: Option<Instant>,
    ) -> SearchResult {
        let limit = self.budget_floor.map_or(self.t_limit, |floor| {
            floor + rng.usize(self.t_limit - floor + 1)
        });
        self.run_limit(target, lib, rng, options, (limit, deadline))
    }
    fn run_limit(
        &self,
        target: &Target,
        lib: &GateLibrary,
        rng: &mut Rng,
        options: &SearchOptions,
        (limit, deadline): (usize, Option<Instant>),
    ) -> SearchResult {
        let n = self.target_map.len();
        let slots = 2 * self.qubits + limit / 4;
        let mut tcount = 0;
        let mut gates: Vec<_> = (0..slots)
            .map(|_| {
                let id = if rng.random01() < options.pid {
                    0
                } else {
                    1 + rng.usize(self.moves.len() - 1)
                };
                if tcount + self.moves[id].tcount > limit {
                    0
                } else {
                    tcount += self.moves[id].tcount;
                    id
                }
            })
            .collect();
        let mut best = gates.clone();
        let mut best_mismatch = usize::MAX;
        let mut prefix: Vec<_> = (0..n).collect();
        let mut suffix: Vec<_> = (0..n).collect();
        let mut scratch = vec![0; n];
        let mut energies = vec![0.0; self.moves.len()];
        let mut weights = energies.clone();
        let mut steps = 0;
        let mut accepted = 0;
        let mut last_improvement = 0;
        for &id in &gates[1..] {
            for value in &mut suffix {
                *value = self.moves[id].map[*value];
            }
        }
        let sweeps = (options.iterations_factor * self.qubits as f64).ceil() as usize;
        'sweeps: for sweep in 0..sweeps {
            let forward = sweep % 2 == 0;
            let temperature = options.start_temp_base / (n as f64).sqrt()
                * options.sweep_cooling.powi(sweep as i32);
            for offset in 0..slots {
                if expired(deadline) || options.is_cancelled() {
                    break 'sweeps;
                }
                let pos = if forward { offset } else { slots - 1 - offset };
                steps += 1;
                let old = gates[pos];
                let others = tcount - self.moves[old].tcount;
                for (id, g) in self.moves.iter().enumerate() {
                    if others + g.tcount > limit {
                        energies[id] = f64::INFINITY;
                        continue;
                    }
                    let mismatch = (0..n)
                        .filter(|&x| suffix[g.map[prefix[x]]] != self.target_map[x])
                        .count();
                    if mismatch == 0 {
                        gates[pos] = id;
                        if let Some(mut result) = self.correct_phases(
                            &gates,
                            CorrectionBudget {
                                used: others + g.tcount,
                                limit,
                                deadline,
                            },
                            target,
                            lib,
                            rng,
                            options,
                        ) {
                            result.steps += steps;
                            result.accepted += accepted;
                            return result;
                        }
                        gates[pos] = old;
                    }
                    energies[id] = mismatch as f64 / n as f64;
                }
                let min = energies.iter().copied().fold(f64::INFINITY, f64::min);
                let mut sum = 0.0;
                for id in 0..self.moves.len() {
                    let prior = if id == 0 {
                        options.pid
                    } else {
                        (1.0 - options.pid) / (self.moves.len() - 1) as f64
                    };
                    weights[id] = (-(energies[id] - min) / temperature.max(1e-14)).exp()
                        * prior.max(1e-12).powf(options.gate_prior);
                    sum += weights[id];
                }
                let mut draw = rng.random01() * sum;
                let mut chosen = self.moves.len() - 1;
                for (id, &w) in weights.iter().enumerate() {
                    draw -= w;
                    if draw <= 0.0 {
                        chosen = id;
                        break;
                    }
                }
                gates[pos] = chosen;
                tcount = others + self.moves[chosen].tcount;
                accepted += usize::from(chosen != old);
                let mismatch = (energies[chosen] * n as f64).round() as usize;
                if mismatch < best_mismatch {
                    best_mismatch = mismatch;
                    best.clone_from(&gates);
                    last_improvement = sweep;
                }
                if offset + 1 < slots {
                    if forward {
                        for x in 0..n {
                            prefix[x] = self.moves[chosen].map[prefix[x]];
                            scratch[x] = suffix[self.moves[gates[pos + 1]].inverse[x]];
                        }
                        std::mem::swap(&mut scratch, &mut suffix);
                    } else {
                        for x in 0..n {
                            prefix[x] = self.moves[gates[pos - 1]].inverse[prefix[x]];
                            scratch[x] = suffix[self.moves[chosen].map[x]];
                        }
                        std::mem::swap(&mut scratch, &mut suffix);
                    }
                }
            }
            if options.stall_sweeps > 0 && sweep >= last_improvement + options.stall_sweeps {
                break;
            }
        }
        let circuit = Circuit::new(self.word(&best), self.qubits, lib);
        let eq = target.original.cost(circuit.matrix(), options.simple_cost);
        let found = target
            .original
            .exact_cost(circuit.matrix(), options.epsilon)
            == 0.0;
        SearchResult {
            circuit,
            best_eq: eq,
            best_energy: eq.powf(2.0 * options.cost_power),
            found,
            steps,
            accepted,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn expired_deadline_prevents_permutation_and_phase_proposals() {
        let lib =
            GateLibrary::load(4, "data/gates/CliffordT", "data/gates/composite_rccx").unwrap();
        let target = Target::new(
            PartialMatrix::read("data/input/64/comparison/U1.txt").unwrap(),
            false,
            false,
        );
        let kernel = PermutationKernel::new(&target, &lib, 11).unwrap();
        let options = SearchOptions {
            iterations_factor: 1e9,
            ..SearchOptions::tuned()
        };
        let result = kernel.run_until(
            &target,
            &lib,
            &mut Rng::new(31),
            &options,
            Some(Instant::now()),
        );
        assert_eq!(result.steps, 0);
        assert_eq!(result.accepted, 0);
    }

    #[test]
    fn all_monomial_moves_have_exact_native_decompositions() {
        for q in 3..=4 {
            let lib =
                GateLibrary::load(q, "data/gates/CliffordT", "data/gates/composite_rccx").unwrap();
            let target = Target::new(
                PartialMatrix::new(Matrix::identity(1 << q), vec![true; 1 << (2 * q)], "id")
                    .unwrap(),
                false,
                false,
            );
            let kernel = PermutationKernel::new(&target, &lib, 15).unwrap();
            for g in kernel.moves {
                let c = Circuit::new(g.word, q, &lib);
                let mut expected = Matrix::zero(1 << q);
                for x in 0..expected.n {
                    expected[(g.map[x], x)] = g.phases[x];
                    assert_eq!(g.inverse[g.map[x]], x);
                }
                assert!(expected.max_abs_diff(c.matrix()) < 1e-12);
                assert_eq!(g.tcount, c.count(&["t".into(), "tdg".into()], &lib));
            }
        }
    }
    #[test]
    fn direct_two_macro_phase_correction_diagnostic() {
        let lib =
            GateLibrary::load(4, "data/gates/CliffordT", "data/gates/composite_rccx").unwrap();
        let target = Target::new(
            PartialMatrix::read("data/input/64/comparison/U1.txt").unwrap(),
            false,
            false,
        );
        let kernel = PermutationKernel::new(&target, &lib, 11).unwrap();
        let mut matched = 0;
        let mut corrected = 0;
        let mut rng = Rng::new(723);
        let options = SearchOptions::tuned();
        for a in 0..kernel.moves.len() {
            for b in 0..kernel.moves.len() {
                if !(0..16)
                    .all(|x| kernel.moves[b].map[kernel.moves[a].map[x]] == kernel.target_map[x])
                {
                    continue;
                }
                matched += 1;
                let count = kernel.moves[a].tcount + kernel.moves[b].tcount;
                if kernel
                    .correct_phases(
                        &[a, b],
                        CorrectionBudget {
                            used: count,
                            limit: 11,
                            deadline: None,
                        },
                        &target,
                        &lib,
                        &mut rng,
                        &options,
                    )
                    .is_some()
                {
                    corrected += 1;
                }
            }
        }
        eprintln!("two-macro matches {matched}, corrected {corrected}");
        assert!(matched > 0);
    }

    #[test]
    fn skeleton_and_phase_search_reach_the_development_u1_t_goal() {
        let lib =
            GateLibrary::load(4, "data/gates/CliffordT", "data/gates/composite_rccx").unwrap();
        let target = Target::new(
            PartialMatrix::read("data/input/64/comparison/U1.txt").unwrap(),
            false,
            false,
        );
        let kernel = PermutationKernel::new(&target, &lib, 11).unwrap();
        let options = SearchOptions::tuned();
        for seed in 0..4 {
            let mut rng = Rng::new(seed);
            let mut solved = false;
            for _ in 0..100 {
                let result = kernel.run(&target, &lib, &mut rng, &options);
                if result.found {
                    assert_eq!(
                        target.original.exact_cost(result.circuit.matrix(), 1e-6),
                        0.0
                    );
                    assert!(result.circuit.count(&["t".into(), "tdg".into()], &lib) <= 11);
                    solved = true;
                    break;
                }
            }
            assert!(solved, "seed {seed}");
        }
    }
}
