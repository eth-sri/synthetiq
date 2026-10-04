//! Coordinate Metropolis sweeps for arbitrary masks. Prefix and suffix products
//! turn each proposal into one sparse and one dense product, while evaluating
//! the complete original masked objective, including its varying norm.
use super::*;

impl SweepKernel {
    pub(super) fn run_masked(
        &self,
        target: &Target,
        lib: &GateLibrary,
        rng: &mut Rng,
        n_gates: usize,
        options: &SearchOptions,
    ) -> Result<SearchResult, String> {
        let equality = |matrix: &Matrix| {
            if options.enable_permutations && options.sweep_permutations {
                target.cost(matrix, options.simple_cost)
            } else {
                target.original.cost(matrix, options.simple_cost)
            }
        };
        let mut initial = random_circuit(lib, rng, n_gates, options.pid);
        if n_gates == 0 {
            return crate::search::search_from(target, lib, rng, initial, options);
        }
        let cap = options
            .quality_slack
            .zip(options.quality_goal)
            .map(|(s, g)| s + g);
        let mut quality = self.circuit_quality(&initial.gates);
        if let Some(cap) = cap {
            let mut positions: Vec<_> = (0..n_gates)
                .filter(|&p| self.qualities[initial.gates[p]] > 0.0)
                .collect();
            while quality > cap + 1e-12 && !positions.is_empty() {
                let p = positions.swap_remove(rng.usize(positions.len()));
                initial.gates[p] = lib.identity;
                quality = self.circuit_quality(&initial.gates);
            }
            initial.rebuild(lib);
        }
        let mut gates = initial.gates.clone();
        let mut best_gates = gates.clone();
        let mut current_eq = equality(initial.matrix());
        let energy = |eq: f64, quality: f64| {
            eq.powf(2.0 * options.cost_power) + options.quality_weight * quality
        };
        let mut current_energy = energy(current_eq, quality);
        let mut best_energy = current_energy;
        let mut best_eq = current_eq;
        let n = initial.matrix().n;
        let mut prefix = Matrix::identity(n);
        let mut suffix = Matrix::identity(n);
        let mut scratch = Matrix::zero(n);
        let mut product = Matrix::zero(n);
        let mut steps = 0;
        let mut accepted = 0;
        let mut last_improvement = 0;
        let mut depth_work = self.depth.as_ref().map(|d| d.workspace(n_gates));
        let sweeps = (options.iterations_factor * lib.n_qubits as f64).ceil() as usize;
        for sweep in 0..sweeps {
            if options.is_cancelled() {
                break;
            }
            let forward = sweep % 2 == 0;
            if let (Some(model), Some(work)) = (&self.depth, &mut depth_work) {
                work.reset(model, &gates, forward);
            }
            if sweep % 8 == 0 {
                prefix = Matrix::identity(n);
                suffix = Matrix::identity(n);
                for &id in &gates[1..] {
                    self.sparse[id].left(&suffix, &mut scratch);
                    std::mem::swap(&mut suffix, &mut scratch);
                }
            }
            let temp = options.start_temp_base / (n as f64).sqrt()
                * options.sweep_cooling.powi(sweep as i32);
            for offset in 0..n_gates {
                let pos = if forward {
                    offset
                } else {
                    n_gates - 1 - offset
                };
                for _ in 0..options.masked_proposals {
                    steps += 1;
                    let id = crate::search::random_gate(
                        lib,
                        rng,
                        options.pid,
                        options.pcomp,
                        options.proba_name,
                    );
                    if id == gates[pos] {
                        continue;
                    }
                    let candidate_quality =
                        if let (Some(model), Some(work)) = (&self.depth, &depth_work) {
                            work.candidate(model, pos, id) as f64
                        } else {
                            quality - self.qualities[gates[pos]] + self.qualities[id]
                        };
                    if cap.is_some_and(|cap| candidate_quality > cap + 1e-12) {
                        continue;
                    }
                    self.sparse[id].left(&prefix, &mut scratch);
                    suffix.mul_into(&scratch, &mut product);
                    let eq = equality(&product);
                    let candidate_energy = energy(eq, candidate_quality);
                    if crate::search::accept_mutation(
                        rng.random01(),
                        candidate_energy,
                        current_energy,
                        temp,
                    ) {
                        gates[pos] = id;
                        current_eq = eq;
                        current_energy = candidate_energy;
                        quality = candidate_quality;
                        accepted += 1;
                        if current_energy + 1e-12 < best_energy {
                            best_energy = current_energy;
                            best_eq = current_eq;
                            best_gates.clone_from(&gates);
                            last_improvement = sweep;
                        }
                        if eq <= options.epsilon * std::f64::consts::SQRT_2 + 1e-7 {
                            let circuit = Circuit::new(gates.clone(), lib.n_qubits, lib);
                            if let Some(circuit) = correct(circuit, target, lib, options.epsilon)? {
                                return Ok(SearchResult {
                                    circuit,
                                    best_eq: eq,
                                    best_energy: current_energy,
                                    found: true,
                                    steps,
                                    accepted,
                                });
                            }
                        }
                    }
                }
                if offset + 1 < n_gates {
                    let chosen = gates[pos];
                    if let (Some(model), Some(work)) = (&self.depth, &mut depth_work) {
                        work.advance(model, chosen);
                    }
                    if forward {
                        self.sparse[chosen].left(&prefix, &mut scratch);
                        std::mem::swap(&mut prefix, &mut scratch);
                        self.adjoints[gates[pos + 1]].right(&suffix, &mut scratch);
                        std::mem::swap(&mut suffix, &mut scratch);
                    } else {
                        self.adjoints[gates[pos - 1]].left(&prefix, &mut scratch);
                        std::mem::swap(&mut prefix, &mut scratch);
                        self.sparse[chosen].right(&suffix, &mut scratch);
                        std::mem::swap(&mut suffix, &mut scratch);
                    }
                }
            }
            if options.stall_sweeps > 0 && sweep >= last_improvement + options.stall_sweeps {
                break;
            }
        }
        let circuit = Circuit::new(best_gates, lib.n_qubits, lib);
        let circuit = correct(circuit.clone(), target, lib, options.epsilon)?.unwrap_or(circuit);
        let found = target
            .original
            .exact_cost(circuit.matrix(), options.epsilon)
            == 0.0;
        Ok(SearchResult {
            circuit,
            best_eq,
            best_energy,
            found,
            steps,
            accepted,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn coordinate_products_match_complete_circuits_in_both_directions() {
        let options = SearchOptions::tuned();
        let mut rng = Rng::new(8971);
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let kernel = SweepKernel::new(&lib, &options);
            let mut gates = random_circuit(&lib, &mut rng, 9, 0.2).gates;
            let n = 1 << q;
            let mut prefix = Matrix::identity(n);
            let mut suffix = Matrix::identity(n);
            let mut scratch = Matrix::zero(n);
            for &id in &gates[1..] {
                kernel.sparse[id].left(&suffix, &mut scratch);
                std::mem::swap(&mut suffix, &mut scratch);
            }
            for sweep in 0..12 {
                let forward = sweep % 2 == 0;
                for offset in 0..gates.len() {
                    let pos = if forward {
                        offset
                    } else {
                        gates.len() - 1 - offset
                    };
                    let id = lib.all[rng.usize(lib.all.len())];
                    gates[pos] = id;
                    let candidate = suffix.mul(&lib.gates[id].matrix).mul(&prefix);
                    let full = Circuit::new(gates.clone(), q, &lib);
                    assert!(candidate.max_abs_diff(full.matrix()) < 1e-12);
                    if offset + 1 == gates.len() {
                        continue;
                    }
                    if forward {
                        kernel.sparse[id].left(&prefix, &mut scratch);
                        std::mem::swap(&mut prefix, &mut scratch);
                        kernel.adjoints[gates[pos + 1]].right(&suffix, &mut scratch);
                        std::mem::swap(&mut suffix, &mut scratch);
                    } else {
                        kernel.adjoints[gates[pos - 1]].left(&prefix, &mut scratch);
                        std::mem::swap(&mut prefix, &mut scratch);
                        kernel.sparse[id].right(&suffix, &mut scratch);
                        std::mem::swap(&mut suffix, &mut scratch);
                    }
                }
            }
        }
    }
    #[test]
    fn arbitrary_mask_solutions_and_capped_failures_match_the_original_objective() {
        let mut rng = Rng::new(22319);
        let mut solved = 0;
        for q in 1..=3 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let options = SearchOptions {
                masked_sweep: true,
                masked_proposals: 8,
                quality_goal: Some(1.0),
                quality_slack: Some(0.0),
                quality_metric: QualityMetric::TCount,
                ..SearchOptions::tuned()
            };
            let kernel = SweepKernel::new(&lib, &options);
            for _ in 0..12 {
                let source = random_circuit(&lib, &mut rng, 4, 0.2);
                let mut cover: Vec<_> = (0..1usize << (2 * q)).map(|_| rng.usize(3) != 0).collect();
                cover[0] = true;
                let target = Target::new(
                    PartialMatrix::new(source.matrix().clone(), cover, "mask").unwrap(),
                    false,
                    false,
                );
                let r = kernel
                    .run_masked(&target, &lib, &mut rng, 8, &options)
                    .unwrap();
                let eq = target
                    .original
                    .cost(r.circuit.matrix(), options.simple_cost);
                assert!((eq - r.best_eq).abs() < 1e-6);
                assert!((eq.powf(2.0 * options.cost_power) - r.best_energy).abs() < 1e-9);
                assert!(r.circuit.count(&["t".into(), "tdg".into()], &lib) <= 1);
                if r.found {
                    solved += 1;
                    assert_eq!(
                        target
                            .original
                            .exact_cost(r.circuit.matrix(), options.epsilon),
                        0.0
                    );
                }
            }
        }
        assert!(solved >= 12, "only {solved} solutions");
    }
}
