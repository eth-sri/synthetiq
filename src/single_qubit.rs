//! Annealing over native Clifford+T syllables for full single-qubit targets.
//! A fixed T budget avoids spending proposals on canceling native gate pairs.
use crate::{
    circuit::Circuit,
    complex::Complex,
    gates::GateLibrary,
    matrix::Matrix,
    partial::Target,
    phase::NativeClifford,
    rng::Rng,
    search::{SearchOptions, SearchResult},
};

type Mat = [Complex; 4];
const ID: Mat = [Complex::ONE, Complex::ZERO, Complex::ZERO, Complex::ONE];
fn mul(a: Mat, b: Mat) -> Mat {
    [
        a[0] * b[0] + a[1] * b[2],
        a[0] * b[1] + a[1] * b[3],
        a[2] * b[0] + a[3] * b[2],
        a[2] * b[1] + a[3] * b[3],
    ]
}
fn overlap(a: Mat, b: Mat) -> Complex {
    (0..4).fold(Complex::ZERO, |z, k| z + a[k].conj() * b[k])
}
struct Syllable {
    matrix: Mat,
    word: Vec<usize>,
}
fn syllable(word: Vec<usize>, lib: &GateLibrary) -> Syllable {
    let circuit = Circuit::new(word.clone(), 1, lib);
    Syllable {
        matrix: circuit.matrix().data.clone().try_into().unwrap(),
        word,
    }
}
pub struct SingleQubitKernel {
    cliffords: Vec<Syllable>,
    syllables: Vec<Syllable>,
    target: Mat,
    limit: usize,
}
impl SingleQubitKernel {
    pub fn new(target: &Target, lib: &GateLibrary, limit: usize) -> Option<Self> {
        if lib.n_qubits != 1
            || limit > 128
            || target.original.cover.iter().any(|&b| !b)
            || target.inverses.iter().any(|&b| b)
        {
            return None;
        }
        let native = NativeClifford::new(lib)?;
        let ids = native.native[0];
        let matrix = &target.original.matrix;
        if !matrix
            .mul(&matrix.adjoint())
            .approximately_equal(&Matrix::identity(2), 1e-9)
        {
            return None;
        }
        let mut cliffords = vec![syllable(Vec::new(), lib)];
        let mut next = 0;
        while next < cliffords.len() {
            for id in [ids[0], ids[1], ids[2]] {
                let mut word = cliffords[next].word.clone();
                word.push(id);
                let candidate = syllable(word, lib);
                if !cliffords
                    .iter()
                    .any(|c| (overlap(c.matrix, candidate.matrix).abs() - 2.0).abs() < 1e-10)
                {
                    cliffords.push(candidate);
                }
            }
            next += 1;
            if cliffords.len() > 24 {
                return None;
            }
        }
        if cliffords.len() != 24 {
            return None;
        }
        let syllables = vec![
            syllable(vec![ids[3], ids[0]], lib),
            syllable(vec![ids[3], ids[0], ids[1]], lib),
            syllable(vec![ids[3]], lib),
        ];
        Some(Self {
            cliffords,
            syllables,
            target: matrix.data.clone().try_into().ok()?,
            limit,
        })
    }
    fn choices(&self, pos: usize, slots: usize) -> &[Syllable] {
        if pos == 0 {
            &self.cliffords
        } else if pos + 1 == slots {
            &self.syllables
        } else {
            &self.syllables[..2]
        }
    }
    fn circuit(&self, gates: &[usize], lib: &GateLibrary) -> Circuit {
        let word = gates
            .iter()
            .enumerate()
            .flat_map(|(pos, &id)| self.choices(pos, gates.len())[id].word.iter().copied())
            .collect();
        Circuit::new(word, 1, lib)
    }
    fn best_clifford(&self, body: Mat) -> (usize, f64) {
        let t = self.target;
        let adjoint = [t[0].conj(), t[2].conj(), t[1].conj(), t[3].conj()];
        let e = mul(adjoint, body);
        let mut best = 0.;
        let mut chosen = 0;
        for (id, c) in self.cliffords.iter().enumerate() {
            let c = c.matrix;
            let norm = (e[0] * c[0] + e[1] * c[2] + e[2] * c[1] + e[3] * c[3]).norm_sqr();
            if norm > best {
                best = norm;
                chosen = id;
            }
        }
        (chosen, (2. - best.sqrt()).max(0.))
    }
    pub fn run(
        &self,
        target: &Target,
        lib: &GateLibrary,
        rng: &mut Rng,
        options: &SearchOptions,
    ) -> SearchResult {
        let count = if self.limit == 0 || rng.random01() < 0.8 {
            self.limit
        } else {
            rng.usize(self.limit)
        };
        let slots = count + 1;
        let mut gates: Vec<_> = (0..slots)
            .map(|p| rng.usize(self.choices(p, slots).len()))
            .collect();
        let mut prefix = vec![ID; slots + 1];
        let mut suffix = prefix.clone();
        let mut best = gates.clone();
        let mut best_energy = f64::INFINITY;
        let mut last_improvement = 0;
        let mut steps = 0;
        let mut accepted = 0;
        if count == 0 {
            best[0] = self.best_clifford(ID).0;
        }
        let sweeps = options.iterations_factor.ceil() as usize;
        for sweep in 0..sweeps {
            if options.is_cancelled() {
                break;
            }
            let forward = sweep % 2 == 0;
            // The terminal Clifford is minimized jointly with each proposed
            // syllable. Products contain only the remaining non-Clifford word.
            for pos in 1..slots {
                prefix[pos + 1] = mul(self.choices(pos, slots)[gates[pos]].matrix, prefix[pos]);
            }
            for pos in (1..slots).rev() {
                suffix[pos] = mul(suffix[pos + 1], self.choices(pos, slots)[gates[pos]].matrix);
            }
            let temperature = options.start_temp_base / std::f64::consts::SQRT_2
                * options.sweep_cooling.powi(sweep as i32);
            for offset in 0..count {
                let pos = if forward {
                    offset + 1
                } else {
                    slots - 1 - offset
                };
                let choices = self.choices(pos, slots);
                let old = gates[pos];
                steps += 1;
                let mut energies = [f64::INFINITY; 3];
                let mut cliffords = [0; 3];
                for (id, choice) in choices.iter().enumerate() {
                    let candidate = mul(suffix[pos + 1], mul(choice.matrix, prefix[pos]));
                    let (clifford, eq2) = self.best_clifford(candidate);
                    cliffords[id] = clifford;
                    energies[id] = eq2.powf(options.cost_power);
                    if eq2 <= 2.0 * options.epsilon * options.epsilon + 1e-12 {
                        gates[pos] = id;
                        gates[0] = clifford;
                        let circuit = self.circuit(&gates, lib);
                        if target
                            .original
                            .exact_cost(circuit.matrix(), options.epsilon)
                            == 0.0
                        {
                            let eq = target.original.cost(circuit.matrix(), options.simple_cost);
                            return SearchResult {
                                circuit,
                                best_eq: eq,
                                best_energy: energies[id],
                                found: true,
                                steps,
                                accepted,
                            };
                        }
                        gates[pos] = old;
                    }
                }
                let min = energies[..choices.len()]
                    .iter()
                    .copied()
                    .fold(f64::INFINITY, f64::min);
                let mut weights = [0.; 3];
                let mut sum = 0.;
                for id in 0..choices.len() {
                    weights[id] = (-(energies[id] - min) / temperature.max(1e-14)).exp();
                    sum += weights[id];
                }
                let mut draw = rng.random01() * sum;
                let mut chosen = choices.len() - 1;
                for (id, &weight) in weights[..choices.len()].iter().enumerate() {
                    draw -= weight;
                    if draw <= 0.0 {
                        chosen = id;
                        break;
                    }
                }
                gates[pos] = chosen;
                gates[0] = cliffords[chosen];
                accepted += usize::from(chosen != old);
                if energies[chosen] + 1e-12 < best_energy {
                    best_energy = energies[chosen];
                    best.clone_from(&gates);
                    last_improvement = sweep;
                }
                if forward {
                    prefix[pos + 1] = mul(choices[chosen].matrix, prefix[pos]);
                } else {
                    suffix[pos] = mul(suffix[pos + 1], choices[chosen].matrix);
                }
            }
            if count == 0
                || (options.stall_sweeps > 0 && sweep >= last_improvement + options.stall_sweeps)
            {
                break;
            }
        }
        let circuit = self.circuit(&best, lib);
        let eq = target.original.cost(circuit.matrix(), options.simple_cost);
        let found = target
            .original
            .exact_cost(circuit.matrix(), options.epsilon)
            == 0.0;
        SearchResult {
            circuit,
            best_eq: eq,
            best_energy: eq.powf(2. * options.cost_power),
            found,
            steps,
            accepted,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partial::PartialMatrix;
    #[test]
    fn every_syllable_preserves_the_native_t_budget() {
        let lib = GateLibrary::load(1, "data/gates/CliffordT", "").unwrap();
        let target = Target::new(
            PartialMatrix::new(Matrix::identity(2), vec![true; 4], "id").unwrap(),
            false,
            false,
        );
        let kernel = SingleQubitKernel::new(&target, &lib, 16).unwrap();
        let names = ["t".into(), "tdg".into()];
        let mut rng = Rng::new(2282);
        for count in 0..=32 {
            let gates: Vec<_> = (0..=count)
                .map(|p| rng.usize(kernel.choices(p, count + 1).len()))
                .collect();
            let circuit = kernel.circuit(&gates, &lib);
            assert_eq!(circuit.count(&names, &lib), count);
            let mut matrix = ID;
            for (p, &id) in gates.iter().enumerate() {
                matrix = mul(kernel.choices(p, count + 1)[id].matrix, matrix);
            }
            assert!(matrix
                .iter()
                .zip(&circuit.matrix().data)
                .all(|(&a, &b)| (a - b).abs() < 1e-12));
        }
    }
    #[test]
    fn annealed_syllables_reproduce_short_random_targets() {
        let lib = GateLibrary::load(1, "data/gates/CliffordT", "").unwrap();
        let mut source = Rng::new(18631);
        let mut rng = Rng::new(75138);
        for length in 1..=15 {
            let circuit = crate::search::random_circuit(&lib, &mut source, length, 0.1);
            let count = circuit.count(&["t".into(), "tdg".into()], &lib);
            let target = Target::new(
                PartialMatrix::new(circuit.matrix().clone(), vec![true; 4], "random 1q").unwrap(),
                false,
                false,
            );
            let kernel = SingleQubitKernel::new(&target, &lib, count).unwrap();
            let options = SearchOptions::tuned();
            let mut solved = false;
            for _ in 0..1000 {
                let result = kernel.run(&target, &lib, &mut rng, &options);
                if result.found {
                    assert!(result.circuit.count(&["t".into(), "tdg".into()], &lib) <= count);
                    assert_eq!(
                        target.original.exact_cost(result.circuit.matrix(), 1e-6),
                        0.0
                    );
                    solved = true;
                    break;
                }
            }
            assert!(solved, "length {length} count {count}");
        }
    }
}
