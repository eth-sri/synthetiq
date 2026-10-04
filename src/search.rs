//! Reference simulated annealing, mutation probabilities, and adaptive lengths.
use crate::{circuit::Circuit, gates::GateLibrary, partial::Target, rng::Rng};

#[derive(Clone, Copy, Debug, Default)]
pub enum QualityMetric {
    TCount,
    TDepth,
    GateCount,
    #[default]
    WeightedCost,
}

#[derive(Clone, Debug)]
pub struct SearchOptions {
    pub pid: f64,
    pub pcomp: f64,
    pub proba_name: f64,
    pub enable_permutations: bool,
    pub simple_cost: bool,
    pub iterations_factor: f64,
    pub start_temp_base: f64,
    pub n_norm: f64,
    pub epsilon: f64,
    pub optimized: bool,
    pub cache_rollback: bool,
    pub cost_power: f64,
    pub sweep_cooling: f64,
    pub stall_sweeps: usize,
    pub gate_prior: f64,
    pub sweep_permutations: bool,
    pub improved_post: bool,
    pub quality_weight: f64,
    pub swap_probability: f64,
    pub quality_slack: Option<f64>,
    pub quality_goal: Option<f64>,
    pub quality_metric: QualityMetric,
    pub quality_expand: bool,
    pub masked_sweep: bool,
    pub masked_proposals: usize,
    pub cancellation: Option<std::sync::Arc<std::sync::atomic::AtomicBool>>,
}

impl Default for SearchOptions {
    fn default() -> Self {
        Self {
            pid: 0.3,
            pcomp: 0.2,
            proba_name: 0.5,
            enable_permutations: true,
            simple_cost: false,
            iterations_factor: 40.0,
            start_temp_base: 0.1,
            n_norm: 80.0,
            epsilon: 1e-6,
            optimized: false,
            cache_rollback: true,
            cost_power: 0.5,
            sweep_cooling: 0.95,
            stall_sweeps: 24,
            gate_prior: 0.0,
            sweep_permutations: true,
            improved_post: false,
            quality_weight: 0.0,
            swap_probability: 0.0,
            quality_slack: None,
            quality_goal: None,
            quality_metric: QualityMetric::WeightedCost,
            quality_expand: false,
            masked_sweep: false,
            masked_proposals: 1,
            cancellation: None,
        }
    }
}

impl SearchOptions {
    #[inline]
    pub fn is_cancelled(&self) -> bool {
        self.cancellation
            .as_ref()
            .is_some_and(|flag| flag.load(std::sync::atomic::Ordering::Relaxed))
    }
    /// Global settings selected on the predeclared tuning split.
    pub fn tuned() -> Self {
        Self {
            optimized: true,
            cost_power: 1.0,
            gate_prior: 1.0,
            stall_sweeps: 12,
            sweep_permutations: false,
            ..Self::default()
        }
    }
}

#[derive(Clone, Debug)]
pub struct GateScheme {
    pub min_start_gates: usize,
    pub max_start_gates: usize,
    pub min_gates: usize,
    pub max_gates: usize,
    pub current_best_gates: usize,
    pub start_best_gates: usize,
    pub beta: f64,
    pub min_factor: f64,
    pub max_factor: f64,
}

impl Default for GateScheme {
    fn default() -> Self {
        Self::new(30, 120, 0.05, 2.5, 3.5)
    }
}

impl GateScheme {
    pub fn new(min: usize, max: usize, beta: f64, min_factor: f64, max_factor: f64) -> Self {
        Self {
            min_start_gates: min,
            max_start_gates: max,
            min_gates: min,
            max_gates: max,
            current_best_gates: max,
            start_best_gates: max,
            beta,
            min_factor,
            max_factor,
        }
    }

    pub fn set_min_start_gates(&mut self, gates: usize) {
        self.min_start_gates = gates;
        self.reset();
    }

    pub fn set_max_start_gates(&mut self, gates: usize) {
        self.max_start_gates = gates;
        self.reset();
    }

    pub fn set_start_best_gates(&mut self, gates: usize) {
        self.start_best_gates = gates;
    }

    pub fn reset(&mut self) {
        self.min_gates = self.min_start_gates;
        self.max_gates = self.max_start_gates;
        self.current_best_gates = self.start_best_gates;
    }

    pub fn update(&mut self, found: usize) {
        self.current_best_gates = self.current_best_gates.min(found);
        fn adjust(current: usize, desired: f64, beta: f64) -> usize {
            if desired > current as f64 {
                current + ((beta * (desired - current as f64) as usize as f64) as usize).max(1)
            } else {
                desired as usize
            }
        }
        self.max_gates = adjust(
            self.max_gates,
            self.max_factor * self.current_best_gates as f64,
            self.beta,
        );
        self.min_gates = adjust(
            self.min_gates,
            self.min_factor * self.current_best_gates as f64,
            self.beta,
        );
    }

    /// The reference intentionally samples (min, max], rather than [min, max).
    pub fn get_start_gates(&self, rng: &mut Rng) -> usize {
        if self.max_gates == 0 {
            1
        } else if self.max_gates <= self.min_gates {
            // The C++ distribution is undefined for equal bounds. A fixed
            // length is useful and is the natural extension of this option.
            self.max_gates
        } else {
            rng.usize(self.max_gates - self.min_gates) + self.min_gates + 1
        }
    }
}

pub struct SearchResult {
    pub circuit: Circuit,
    pub best_eq: f64,
    pub best_energy: f64,
    pub found: bool,
    pub steps: usize,
    pub accepted: usize,
}

/// Select a gate with exactly the reference two-stage mixture.
pub fn random_gate(
    lib: &GateLibrary,
    rng: &mut Rng,
    pid: f64,
    pcomp: f64,
    proba_name: f64,
) -> usize {
    if rng.random01() < pid {
        return lib.identity;
    }
    if rng.random01() < proba_name {
        let probability = lib.basic_by_name.len() as f64
            / (lib.basic_by_name.len() as f64 + pcomp * lib.composite_by_name.len() as f64);
        let groups = if rng.random01() < probability {
            &lib.basic_by_name
        } else {
            &lib.composite_by_name
        };
        if groups.is_empty() {
            return lib.identity;
        }
        let group = &groups[rng.usize(groups.len())];
        group[rng.usize(group.len())]
    } else {
        let probability =
            lib.basic.len() as f64 / (lib.basic.len() as f64 + pcomp * lib.composite.len() as f64);
        let gates = if rng.random01() < probability {
            &lib.basic
        } else {
            &lib.composite
        };
        if gates.is_empty() {
            return lib.identity;
        }
        gates[rng.usize(gates.len())]
    }
}

pub fn random_circuit(lib: &GateLibrary, rng: &mut Rng, n_gates: usize, pid: f64) -> Circuit {
    // Initial circuits use pcomp=1 independently of the mutation setting.
    let gates = (0..n_gates)
        .map(|_| random_gate(lib, rng, pid, 1.0, 0.5))
        .collect();
    Circuit::new(gates, lib.n_qubits, lib)
}

pub fn temperature(start: f64, accepted: usize, n_gates: usize, normalizer: f64) -> f64 {
    start * (-(accepted as f64) / (n_gates as f64 * normalizer)).exp()
}

pub fn accept_mutation(uniform: f64, candidate: f64, current: f64, temperature: f64) -> bool {
    candidate <= current || uniform <= (-(candidate - current) / temperature).exp()
}

pub fn search(
    target: &Target,
    lib: &GateLibrary,
    rng: &mut Rng,
    n_gates: usize,
    options: &SearchOptions,
) -> Result<SearchResult, String> {
    let initial = random_circuit(lib, rng, n_gates, options.pid);
    search_from(target, lib, rng, initial, options)
}

pub fn search_from(
    target: &Target,
    lib: &GateLibrary,
    rng: &mut Rng,
    candidate: Circuit,
    options: &SearchOptions,
) -> Result<SearchResult, String> {
    search_from_observed(target, lib, rng, candidate, options, |_, _, _| {})
}

/// Observe the actual annealing schedule without changing random draws or moves.
/// The no-observer entry point specializes this callback away in release builds.
pub fn search_from_observed(
    target: &Target,
    lib: &GateLibrary,
    rng: &mut Rng,
    mut candidate: Circuit,
    options: &SearchOptions,
    mut observe: impl FnMut(usize, usize, f64),
) -> Result<SearchResult, String> {
    let equality = |circuit: &Circuit| {
        if options.enable_permutations {
            target.cost(circuit.matrix(), options.simple_cost)
        } else {
            target.original.cost(circuit.matrix(), options.simple_cost)
        }
    };
    let initial_cost = equality(&candidate);
    let mut result = SearchResult {
        circuit: candidate.clone(),
        best_eq: initial_cost,
        best_energy: initial_cost,
        found: false,
        steps: 0,
        accepted: 0,
    };
    let mut current_energy = initial_cost;
    let n_gates = candidate.gates.len();
    let iterations =
        (options.iterations_factor * lib.n_qubits as f64 * n_gates as f64).ceil() as usize;
    let start_temp = options.start_temp_base / ((1_usize << lib.n_qubits) as f64).sqrt();
    let mut improved_match = false;
    let mut undo = options.cache_rollback.then(|| candidate.mutation_backup());
    for _ in 0..iterations {
        if options.is_cancelled() {
            break;
        }
        result.steps += 1;
        let temp = temperature(start_temp, result.accepted, n_gates, options.n_norm);
        observe(result.steps - 1, result.accepted, temp);
        let position = rng.usize(n_gates);
        let old_id = candidate.gates[position];
        let new_id = random_gate(lib, rng, options.pid, options.pcomp, options.proba_name);
        let old_gate = &lib.gates[old_id];
        let new_gate = &lib.gates[new_id];
        if old_gate.name == new_gate.name && old_gate.qubits == new_gate.qubits {
            continue;
        }
        if let Some(undo) = &mut undo {
            candidate.replace_recorded(position, new_id, lib, undo);
        } else {
            candidate.replace(position, new_id, lib);
        }
        let candidate_energy = equality(&candidate);
        // This draw is consumed even for downhill moves, as in C++.
        let uniform = rng.random01();
        if accept_mutation(uniform, candidate_energy, current_energy, temp) {
            result.accepted += 1;
            if candidate_energy < result.best_energy {
                result.best_energy = candidate_energy;
                result.best_eq = candidate_energy;
                result.circuit = candidate.clone();
                if target.exact_cost(result.circuit.matrix(), options.epsilon) < 1e-3 {
                    improved_match = true;
                    break;
                }
            }
            current_energy = candidate_energy;
        } else {
            if let Some(undo) = &undo {
                candidate.restore_recorded(undo);
            } else {
                candidate.replace(position, old_id, lib);
            }
        }
    }
    if improved_match {
        let mut best_index = 0;
        let mut best_exact = f64::INFINITY;
        for (index, variant) in target.variants.iter().enumerate() {
            let exact = variant.exact_cost(result.circuit.matrix(), options.epsilon);
            if exact < best_exact {
                best_index = index;
                best_exact = exact;
            }
        }
        if target.inverses[best_index] {
            result.circuit.inverse(lib)?;
        }
        result
            .circuit
            .permute(&target.permutations[best_index], lib)?;
    }
    result.found = target.exact_cost(result.circuit.matrix(), options.epsilon) < 1e-3;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cancelled_annealers_stop_before_mutation_proposals() {
        let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
        let flag = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(true));
        let options = SearchOptions {
            cancellation: Some(flag.clone()),
            ..SearchOptions::tuned()
        };
        let target = Target::new(
            crate::partial::PartialMatrix::new(
                crate::matrix::Matrix::identity(4),
                vec![true; 16],
                "cancel",
            )
            .unwrap(),
            false,
            false,
        );
        let kernel = crate::anneal::SweepKernel::new(&lib, &options);
        let result = kernel
            .run(&target, &lib, &mut Rng::new(17), 20, &options)
            .unwrap();
        assert_eq!(result.steps, 0);
        let initial = random_circuit(&lib, &mut Rng::new(18), 20, 0.3);
        let result = search_from(&target, &lib, &mut Rng::new(19), initial, &options).unwrap();
        assert_eq!(result.steps, 0);
        let partial = Target::new(
            crate::partial::PartialMatrix::new(
                crate::matrix::Matrix::identity(4),
                (0..16).map(|i| i % 3 != 0).collect(),
                "cancel mask",
            )
            .unwrap(),
            false,
            false,
        );
        let result = kernel
            .run(&partial, &lib, &mut Rng::new(20), 20, &options)
            .unwrap();
        assert_eq!(result.steps, 0);
        flag.store(false, std::sync::atomic::Ordering::Relaxed);
        assert!(!options.is_cancelled());
    }

    #[test]
    fn temperature_depends_on_accepted_moves() {
        assert_eq!(temperature(0.1, 0, 20, 80.0), 0.1);
        assert!((temperature(0.1, 1600, 20, 80.0) - 0.1 / std::f64::consts::E).abs() < 1e-16);
    }

    #[test]
    fn metropolis_boundary_cases() {
        assert!(accept_mutation(0.999, 0.1, 0.2, 0.001));
        assert!(accept_mutation(1.0, 0.2, 0.2, 0.001));
        assert!(!accept_mutation(0.5, 0.3, 0.2, 0.01));
        assert!(accept_mutation(0.00001, 0.3, 0.2, 0.01));
        assert!(!accept_mutation(0.1, 0.3, 0.2, 0.0));
    }

    #[test]
    fn adaptive_lengths_match_reference_rounding() {
        let mut scheme = GateScheme::default();
        scheme.update(20);
        assert_eq!((scheme.min_gates, scheme.max_gates), (31, 70));
        scheme.update(20);
        assert_eq!((scheme.min_gates, scheme.max_gates), (32, 70));
        scheme.update(4);
        assert_eq!((scheme.min_gates, scheme.max_gates), (10, 14));
        scheme.update(6);
        assert_eq!((scheme.min_gates, scheme.max_gates), (10, 14));
        scheme.reset();
        assert_eq!((scheme.min_gates, scheme.max_gates), (30, 120));
    }

    #[test]
    fn adaptive_lengths_use_open_lower_bound() {
        let scheme = GateScheme::new(3, 4, 0.05, 2.5, 3.5);
        let mut rng = Rng::new(0);
        for _ in 0..100 {
            assert_eq!(scheme.get_start_gates(&mut rng), 4);
        }
        assert_eq!(
            GateScheme::new(0, 0, 0.05, 2.5, 3.5).get_start_gates(&mut rng),
            1
        );
    }
}

#[cfg(test)]
mod integration_tests {
    use super::*;
    use crate::{matrix::Matrix, partial::PartialMatrix};
    fn library(n: usize) -> GateLibrary {
        GateLibrary::load(n, "data/gates/CliffordT", "").unwrap()
    }

    #[test]
    fn recovers_every_one_qubit_gate_for_multiple_seeds() {
        let lib = library(1);
        for &id in &lib.all {
            let target = Target::new(
                PartialMatrix::new(lib.gates[id].matrix.clone(), vec![true; 4], "gate").unwrap(),
                true,
                false,
            );
            for seed in 0..8 {
                let options = SearchOptions {
                    iterations_factor: 300.0,
                    ..SearchOptions::default()
                };
                let result = search(&target, &lib, &mut Rng::new(seed), 1, &options).unwrap();
                assert!(result.found, "gate {} seed {seed}", lib.gates[id].name);
                assert_eq!(
                    target.original.exact_cost(result.circuit.matrix(), 1e-6),
                    0.0
                );
                assert!(result.best_eq < 1e-6);
            }
        }
    }

    #[test]
    fn identical_seeds_replay_the_same_search() {
        let lib = library(2);
        let expected = Circuit::new(
            vec![
                lib.find("h", &[0]).unwrap(),
                lib.find("cx", &[0, 1]).unwrap(),
            ],
            2,
            &lib,
        );
        let target = Target::new(
            PartialMatrix::new(expected.matrix().clone(), vec![true; 16], "bell").unwrap(),
            true,
            false,
        );
        let a = search(
            &target,
            &lib,
            &mut Rng::new(7),
            5,
            &SearchOptions::default(),
        )
        .unwrap();
        let b = search(
            &target,
            &lib,
            &mut Rng::new(7),
            5,
            &SearchOptions::default(),
        )
        .unwrap();
        assert_eq!(a.circuit.gates, b.circuit.gates);
        assert_eq!(
            (a.best_eq, a.steps, a.accepted),
            (b.best_eq, b.steps, b.accepted)
        );
        assert!(a.best_eq <= target.cost(&Matrix::identity(4), false));
    }

    #[test]
    fn disabled_permutations_uses_original_for_search_cost() {
        let lib = library(2);
        let h0 = lib.find("h", &[0]).unwrap();
        let h1 = lib.find("h", &[1]).unwrap();
        let target = Target::new(
            PartialMatrix::new(lib.gates[h0].matrix.clone(), vec![true; 16], "h0").unwrap(),
            true,
            false,
        );
        let options = SearchOptions {
            iterations_factor: 0.0,
            ..SearchOptions::default()
        };
        let initial = Circuit::new(vec![h1], 2, &lib);
        let a = search_from(&target, &lib, &mut Rng::new(0), initial.clone(), &options).unwrap();
        let b = search_from(
            &target,
            &lib,
            &mut Rng::new(0),
            initial,
            &SearchOptions {
                enable_permutations: false,
                ..options
            },
        )
        .unwrap();
        assert!(a.best_eq < 1e-6);
        assert!(b.best_eq > 0.5);
    }

    #[test]
    fn empty_search_returns_identity_without_random_draws() {
        let lib = library(1);
        let target = Target::new(
            PartialMatrix::new(Matrix::identity(2), vec![true; 4], "id").unwrap(),
            false,
            false,
        );
        let mut rng = Rng::new(42);
        let result = search(&target, &lib, &mut rng, 0, &SearchOptions::default()).unwrap();
        assert!(result.found);
        assert_eq!(result.steps, 0);
        assert_eq!(rng.next_u32(), Rng::new(42).next_u32());
    }

    #[test]
    fn identity_probability_one_never_produces_other_gates() {
        let lib = library(2);
        let circuit = random_circuit(&lib, &mut Rng::new(42), 50, 1.0);
        assert_eq!(circuit.non_identity(), 0);
        assert_eq!(circuit.cost(), 0.0);
        assert!(circuit
            .matrix()
            .approximately_equal(&Matrix::identity(4), 0.0));
    }
}
