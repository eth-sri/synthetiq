//! Simulated annealing over finite parity-phase gadgets, with local basis changes.
//! All proposals expand to the user's native Clifford+T gates. No target-specific
//! circuit templates or deterministic phase-polynomial solver are used.
use crate::{
    circuit::Circuit,
    complex::Complex,
    gates::GateLibrary,
    matrix::Matrix,
    partial::Target,
    rng::Rng,
    search::{SearchOptions, SearchResult},
};

fn roots() -> [Complex; 8] {
    let a = std::f64::consts::FRAC_1_SQRT_2;
    [
        Complex::ONE,
        Complex::new(a, a),
        Complex::I,
        Complex::new(-a, a),
        Complex::new(-1.0, 0.0),
        Complex::new(-a, -a),
        -Complex::I,
        Complex::new(a, -a),
    ]
}
fn same_phase(a: &Matrix, b: &Matrix) -> bool {
    let z = b.trace_conjugate_product(a);
    if z.abs() < 1e-12 {
        return false;
    }
    let phase = z / z.abs();
    a.data
        .iter()
        .zip(&b.data)
        .all(|(&x, &y)| (x - phase * y).norm_sqr() < 1e-20)
}
#[derive(Clone)]
struct Basis {
    matrix: Matrix,
    words: Vec<usize>,
    tcount: usize,
}
fn local_axes() -> Vec<Basis> {
    let a = std::f64::consts::FRAC_1_SQRT_2;
    let r = roots();
    let diagonal = |v| Matrix {
        n: 2,
        data: vec![Complex::ONE, Complex::ZERO, Complex::ZERO, v],
    };
    let generators = [
        Matrix {
            n: 2,
            data: vec![a.into(), a.into(), a.into(), (-a).into()],
        },
        diagonal(Complex::I),
        diagonal(-Complex::I),
        diagonal(r[1]),
        diagonal(r[7]),
    ];
    let mut states = vec![Basis {
        matrix: Matrix::identity(2),
        words: Vec::new(),
        tcount: 0,
    }];
    let mut next = 0;
    while next < states.len() {
        for (id, gate) in generators.iter().enumerate() {
            let tcount = states[next].tcount + usize::from(id >= 3);
            if tcount > 1 {
                continue;
            }
            let product = gate.mul(&states[next].matrix);
            if states.iter().any(|s| same_phase(&s.matrix, &product)) {
                continue;
            }
            assert!(states.len() < 96, "unexpected single-T Clifford closure");
            let mut words = states[next].words.clone();
            words.push(id);
            states.push(Basis {
                matrix: product,
                words,
                tcount,
            });
        }
        next += 1;
    }
    states.sort_by_key(|s| (s.tcount, s.words.len()));
    let z = diagonal((-1.0).into());
    let mut axes = Vec::<(Matrix, Basis)>::new();
    for s in states {
        let axis = s.matrix.mul(&z).mul(&s.matrix.adjoint());
        if !axes.iter().any(|(a, _)| same_phase(a, &axis)) {
            axes.push((axis, s));
        }
    }
    axes.into_iter().map(|(_, basis)| basis).collect()
}

/// Boolean finite differences characterize integer parity polynomials mod 8:
/// coefficients of degree d are divisible by 2^(d-1), and vanish for d >= 4.
fn parity_representable(phases: &[u8], qubits: usize) -> bool {
    let mut differences = phases.to_vec();
    for bit in 0..qubits {
        for mask in 0..phases.len() {
            if mask & (1 << bit) != 0 {
                differences[mask] = (differences[mask] + 8 - differences[mask ^ (1 << bit)]) % 8;
            }
        }
    }
    (1..phases.len()).all(|mask| {
        let divisor = 1usize << (mask.count_ones() as usize - 1).min(3);
        (differences[mask] as usize).is_multiple_of(divisor)
    })
}

/// Exact minimum number of odd coefficients in the parity-only representation,
/// for at most four qubits. This is an eligibility bound for this kernel, not
/// an optimality claim about unrestricted Clifford+T circuits.
fn parity_minimum_t(phases: &[u8], qubits: usize) -> Option<usize> {
    if qubits > 4 || !parity_representable(phases, qubits) {
        return None;
    }
    let n = phases.len();
    let mut differences = phases.to_vec();
    for bit in 0..qubits {
        for mask in 0..n {
            if mask & (1 << bit) != 0 {
                differences[mask] = (differences[mask] + 8 - differences[mask ^ (1 << bit)]) % 8;
            }
        }
    }
    // Degree-one, -two and -three differences determine sums of odd parity
    // coefficients over supersets. The degree-four coefficient is free.
    let mut minimum = usize::MAX;
    for top in 0..if qubits == 4 { 2 } else { 1 } {
        let mut odd = vec![0u8; n];
        if qubits == 4 {
            odd[n - 1] = top;
        }
        for mask in (1..n).rev() {
            let degree = mask.count_ones();
            if degree > 3 {
                continue;
            }
            let signature = (differences[mask] >> (degree - 1)) & 1;
            let supersets = ((mask + 1)..n)
                .filter(|&s| s & mask == mask)
                .fold(0, |sum, s| sum ^ odd[s]);
            odd[mask] = signature ^ supersets;
        }
        minimum = minimum.min(odd.iter().map(|&x| x as usize).sum());
    }
    Some(minimum)
}

#[derive(Clone)]
pub(crate) struct NativeClifford {
    pub native: Vec<[usize; 5]>,
    pub cx: Vec<Vec<Option<usize>>>,
}
impl NativeClifford {
    pub fn new(lib: &GateLibrary) -> Option<Self> {
        let q = lib.n_qubits;
        let n = 1 << q;
        let names = ["h", "s", "sdg", "t", "tdg"];
        let r = roots();
        let a = std::f64::consts::FRAC_1_SQRT_2;
        let diagonal = |z| Matrix {
            n: 2,
            data: vec![Complex::ONE, Complex::ZERO, Complex::ZERO, z],
        };
        let standard = [
            Matrix {
                n: 2,
                data: vec![a.into(), a.into(), a.into(), (-a).into()],
            },
            diagonal(Complex::I),
            diagonal(-Complex::I),
            diagonal(r[1]),
            diagonal(r[7]),
        ];
        let mut native = Vec::new();
        for bit in 0..q {
            let mut ids = [0; 5];
            let mut order: Vec<_> = (0..q).collect();
            order.swap(0, bit);
            for (k, name) in names.iter().enumerate() {
                ids[k] = lib.find(name, &[bit])?;
                if !same_phase(
                    &lib.gates[ids[k]].matrix,
                    &standard[k].kron_identity(q).permuted(&order),
                ) {
                    return None;
                }
            }
            native.push(ids);
        }
        let mut cx = vec![vec![None; q]; q];
        for (control, row) in cx.iter_mut().enumerate() {
            for (to, entry) in row.iter_mut().enumerate() {
                if control == to {
                    continue;
                }
                let id = lib.find("cx", &[control, to])?;
                let mut expected = Matrix::zero(n);
                for input in 0..n {
                    let output = input ^ (((input >> control) & 1) << to);
                    expected[(output, input)] = Complex::ONE;
                }
                if !same_phase(&lib.gates[id].matrix, &expected) {
                    return None;
                }
                *entry = Some(id);
            }
        }
        Some(Self { native, cx })
    }
}

struct NonlinearGadget {
    word: Vec<usize>,
    inverse: Vec<usize>,
    target: usize,
}

pub struct PhaseKernel {
    n_qubits: usize,
    phases: Vec<u8>,
    affected: Vec<Vec<usize>>,
    prefix: Vec<usize>,
    suffix: Vec<usize>,
    native: Vec<[usize; 5]>,
    cx: Vec<Vec<Option<usize>>>,
    inner_t_limit: usize,
    depth_limit: Option<usize>,
    nonlinear: Vec<NonlinearGadget>,
}
impl PhaseKernel {
    /// Only exact eighth-root diagonal functions and full specifications qualify.
    /// General synthesis remains available as a separate portfolio component.
    pub fn new(target: &Target, lib: &GateLibrary, t_limit: usize) -> Option<Self> {
        let p = &target.original;
        let q = lib.n_qubits;
        if q > 6 || p.cover.iter().any(|&v| !v) || target.inverses.iter().any(|&v| v) {
            return None;
        }
        let n = 1 << q;
        let NativeClifford { native, cx } = NativeClifford::new(lib)?;
        let diagonal = |z| Matrix {
            n: 2,
            data: vec![Complex::ONE, Complex::ZERO, Complex::ZERO, z],
        };
        let affine = crate::affine::reduce(&p.matrix, lib, &native, &cx);
        let working = affine.as_ref().map_or(&p.matrix, |a| &a.matrix);
        let mut suffix = Vec::new();
        let is_diagonal =
            (0..n).all(|i| (0..n).all(|j| i == j || working[(i, j)].norm_sqr() < 1e-20));
        if !is_diagonal {
            let z = diagonal((-1.0).into());
            let axes = local_axes();
            let mut local = true;
            for (bit, ids) in native.iter().enumerate() {
                let mut order: Vec<_> = (0..q).collect();
                order.swap(0, bit);
                let basis = axes.iter().find(|basis| {
                    let axis = basis
                        .matrix
                        .mul(&z)
                        .mul(&basis.matrix.adjoint())
                        .kron_identity(q)
                        .permuted(&order);
                    working.mul(&axis).max_abs_diff(&axis.mul(working)) < 1e-9
                });
                let Some(basis) = basis else {
                    local = false;
                    break;
                };
                suffix.extend(basis.words.iter().map(|&k| ids[k]));
            }
            if !local {
                let metadata = NativeClifford {
                    native: native.clone(),
                    cx: cx.clone(),
                };
                let word = crate::clifford_diagonal::basis_word(working, lib, &metadata)?;
                suffix = crate::affine::inverse_word(&word, lib)?;
            }
        }
        let mut prefix = crate::affine::inverse_word(&suffix, lib)?;
        let basis = Circuit::new(suffix.clone(), q, lib);
        let d = basis.matrix().adjoint().mul(working).mul(basis.matrix());
        if let Some(affine) = affine {
            let mut before = crate::affine::inverse_word(&affine.right, lib)?;
            before.extend(prefix);
            prefix = before;
            suffix.extend(crate::affine::inverse_word(&affine.left, lib)?);
        }
        let wrapper_t = suffix
            .iter()
            .chain(&prefix)
            .filter(|&&id| matches!(lib.gates[id].name.as_str(), "t" | "tdg"))
            .count();
        let inner_t_limit = t_limit.checked_sub(wrapper_t)?;
        if (0..n).any(|i| (0..n).any(|j| i != j && d[(i, j)].norm_sqr() > 1e-18)) {
            return None;
        }
        let diagonal: Vec<_> = (0..n).map(|x| d[(x, x)]).collect();
        let mut kernel =
            Self::from_diagonal(&diagonal, &NativeClifford { native, cx }, inner_t_limit)?;
        kernel.prefix = prefix;
        kernel.suffix = suffix;
        Some(kernel)
    }

    /// Reuse a validated native library when correcting many diagonal residuals.
    pub(crate) fn from_diagonal(
        diagonal: &[Complex],
        native: &NativeClifford,
        t_limit: usize,
    ) -> Option<Self> {
        Self::diagonal_basis(diagonal, native, t_limit, true)
    }
    fn diagonal_basis(
        diagonal: &[Complex],
        native: &NativeClifford,
        t_limit: usize,
        linear_only: bool,
    ) -> Option<Self> {
        let q = native.native.len();
        let n = 1 << q;
        if q > 6 || diagonal.len() != n {
            return None;
        }
        let r = roots();
        let mut phases = Vec::with_capacity(n);
        for &value in diagonal {
            if (value.norm_sqr() - 1.0).abs() > 1e-9 {
                return None;
            }
            let relative = value / diagonal[0];
            let index = r
                .iter()
                .position(|&root| (relative - root).norm_sqr() < 1e-18)?;
            phases.push(index as u8);
        }
        if linear_only
            && (!parity_representable(&phases, q)
                || parity_minimum_t(&phases, q).is_some_and(|minimum| minimum > t_limit))
        {
            return None;
        }
        let affected = (1usize..n)
            .map(|mask| {
                (0usize..n)
                    .filter(|&x| (mask & x).count_ones() % 2 == 1)
                    .collect()
            })
            .collect();
        Some(Self {
            n_qubits: q,
            phases,
            affected,
            prefix: Vec::new(),
            suffix: Vec::new(),
            native: native.native.clone(),
            cx: native.cx.clone(),
            inner_t_limit: t_limit,
            depth_limit: None,
            nonlinear: Vec::new(),
        })
    }

    /// Quadratic Boolean phases implemented by native relative-Toffoli
    /// conjugation expand the finite annealing basis beyond parity functions.
    pub fn new_nonlinear(target: &Target, lib: &GateLibrary, t_limit: usize) -> Option<Self> {
        let q = lib.n_qubits;
        let n = 1 << q;
        if !(3..=4).contains(&q)
            || t_limit < 9
            || target.original.cover.iter().any(|&b| !b)
            || target.inverses.iter().any(|&b| b)
        {
            return None;
        }
        let matrix = &target.original.matrix;
        if (0..n).any(|i| (0..n).any(|j| i != j && matrix[(i, j)].norm_sqr() > 1e-18)) {
            return None;
        }
        let native = NativeClifford::new(lib)?;
        let diagonal: Vec<_> = (0..n).map(|x| matrix[(x, x)]).collect();
        let mut kernel = Self::diagonal_basis(&diagonal, &native, t_limit, false)?;
        if parity_representable(&kernel.phases, q) {
            return None;
        }
        let layout = crate::phase_layers::PhaseLayout {
            native: &native.native,
            cx: &native.cx,
        };
        let mut seen = std::collections::HashSet::new();
        for u in 1usize..n {
            for v in 1usize..n {
                for w in v + 1..n {
                    if [0, v, w, v ^ w].contains(&u) {
                        continue;
                    }
                    let affected: Vec<_> = (0usize..n)
                        .filter(|&x| {
                            ((x & u).count_ones() % 2)
                                ^ (((x & v).count_ones() % 2) & ((x & w).count_ones() % 2))
                                != 0
                        })
                        .collect();
                    if !seen.insert(affected.clone()) {
                        continue;
                    }
                    // Below 18 T gates, at most one odd nonlinear gadget can occur.
                    // Reject bases whose residue cannot be completed by parity phases.
                    if t_limit < 18
                        && !(1u8..8).step_by(2).any(|coefficient| {
                            let mut residual = kernel.phases.clone();
                            for &x in &affected {
                                residual[x] = (residual[x] + 8 - coefficient) % 8;
                            }
                            parity_representable(&residual, q)
                        })
                    {
                        continue;
                    }
                    let mut rows = vec![v, w, u];
                    for bit in 0..q {
                        let unit = 1 << bit;
                        if !(0..1usize << rows.len()).any(|combination| {
                            rows.iter().enumerate().fold(0, |sum, (i, &row)| {
                                sum ^ if combination & (1 << i) != 0 { row } else { 0 }
                            }) == unit
                        }) {
                            rows.push(unit);
                        }
                    }
                    let mut word = layout.basis_word(rows);
                    let ids = native.native[2];
                    word.extend([
                        ids[0],
                        ids[3],
                        native.cx[1][2]?,
                        ids[4],
                        native.cx[0][2]?,
                        ids[3],
                        native.cx[1][2]?,
                        ids[4],
                        ids[0],
                    ]);
                    let inverse = crate::affine::inverse_word(&word, lib)?;
                    kernel.nonlinear.push(NonlinearGadget {
                        word,
                        inverse,
                        target: 2,
                    });
                    kernel.affected.push(affected);
                }
            }
        }
        if kernel.nonlinear.is_empty() {
            return None;
        }
        Some(kernel)
    }

    pub fn new_factored(target: &Target, lib: &GateLibrary, t_limit: usize) -> Option<Self> {
        if let Some(kernel) = Self::new(target, lib, t_limit) {
            return Some(kernel);
        }
        if target.original.cover.iter().any(|&b| !b) || target.inverses.iter().any(|&b| b) {
            return None;
        }
        let native = NativeClifford::new(lib)?;
        if let Some((input, output, matrix)) =
            crate::clifford_diagonal::two_sided_basis(&target.original.matrix, lib, &native)
        {
            let partial = crate::partial::PartialMatrix::new(
                matrix,
                target.original.cover.clone(),
                "Clifford factored",
            )
            .ok()?;
            if let Some(mut kernel) = Self::new(&Target::new(partial, false, false), lib, t_limit) {
                let mut prefix = input;
                prefix.extend(kernel.prefix);
                kernel.prefix = prefix;
                kernel
                    .suffix
                    .extend(crate::affine::inverse_word(&output, lib)?);
                return Some(kernel);
            }
        }
        let sparse = crate::affine::reduce_sparse_columns(
            &target.original.matrix,
            lib,
            &native.native,
            &native.cx,
        );
        let working = sparse
            .as_ref()
            .map_or(&target.original.matrix, |a| &a.matrix);
        let restore = if let Some(a) = &sparse {
            Some(crate::affine::inverse_word(&a.left, lib)?)
        } else {
            None
        };
        if sparse.is_some() {
            let partial = crate::partial::PartialMatrix::new(
                working.clone(),
                target.original.cover.clone(),
                "sparse affine factored",
            )
            .ok()?;
            if let Some(mut kernel) = Self::new(&Target::new(partial, false, false), lib, t_limit) {
                if let Some(restore) = restore {
                    kernel.suffix.extend(restore);
                }
                return Some(kernel);
            }
        }
        let inner_limit = t_limit.checked_sub(2)?;
        for (word, matrix) in crate::clifford_diagonal::near_pauli_bases(working, lib, &native) {
            let partial = crate::partial::PartialMatrix::new(
                matrix,
                target.original.cover.clone(),
                "one T basis",
            )
            .ok()?;
            if let Some(mut kernel) =
                Self::new(&Target::new(partial, false, false), lib, inner_limit)
            {
                let mut prefix = crate::affine::inverse_word(&word, lib)?;
                prefix.extend(kernel.prefix);
                kernel.prefix = prefix;
                kernel.suffix.extend(word);
                if let Some(restore) = restore {
                    kernel.suffix.extend(restore);
                }
                return Some(kernel);
            }
        }
        None
    }

    pub fn with_depth_limit(mut self, limit: usize) -> Self {
        self.depth_limit = Some(limit);
        self
    }

    fn scheduled_circuit(&self, coefficients: &[u8], lib: &GateLibrary, rng: &mut Rng) -> Circuit {
        if let Some(limit) = self.depth_limit.filter(|_| self.nonlinear.is_empty()) {
            let layout = crate::phase_layers::PhaseLayout {
                native: &self.native,
                cx: &self.cx,
            };
            let gates =
                layout.optimize(coefficients, (&self.prefix, &self.suffix), lib, rng, limit);
            Circuit::new(gates, self.n_qubits, lib)
        } else {
            self.circuit(coefficients, lib)
        }
    }

    fn circuit(&self, coefficients: &[u8], lib: &GateLibrary) -> Circuit {
        let mut gates = self.prefix.clone();
        // Gray order shares parity-computation CX gates and restores the basis
        // after each pivot group; zero phase terms need no visit.
        for pivot in 0..self.n_qubits {
            let mut previous = 0;
            for index in 0..1usize << pivot {
                let gray = index ^ (index >> 1);
                let coefficient = coefficients[((1 << pivot) | gray) - 1];
                if coefficient == 0 {
                    continue;
                }
                let difference = previous ^ gray;
                for bit in 0..pivot {
                    if difference & (1 << bit) != 0 {
                        gates.push(self.cx[bit][pivot].unwrap());
                    }
                }
                let ids = self.native[pivot];
                match coefficient {
                    1 => gates.push(ids[3]),
                    2 => gates.push(ids[1]),
                    3 => gates.extend([ids[1], ids[3]]),
                    4 => gates.extend([ids[1], ids[1]]),
                    5 => gates.extend([ids[2], ids[4]]),
                    6 => gates.push(ids[2]),
                    7 => gates.push(ids[4]),
                    _ => unreachable!(),
                }
                previous = gray;
            }
            for bit in 0..pivot {
                if previous & (1 << bit) != 0 {
                    gates.push(self.cx[bit][pivot].unwrap());
                }
            }
        }
        for (index, gadget) in self.nonlinear.iter().enumerate() {
            let coefficient = coefficients[self.phases.len() - 1 + index];
            if coefficient == 0 {
                continue;
            }
            gates.extend_from_slice(&gadget.word);
            let ids = self.native[gadget.target];
            match coefficient {
                1 => gates.push(ids[3]),
                2 => gates.push(ids[1]),
                3 => gates.extend([ids[1], ids[3]]),
                4 => gates.extend([ids[1], ids[1]]),
                5 => gates.extend([ids[2], ids[4]]),
                6 => gates.push(ids[2]),
                7 => gates.push(ids[4]),
                _ => unreachable!(),
            }
            gates.extend_from_slice(&gadget.inverse);
        }
        gates.extend(&self.suffix);
        Circuit::new(gates, self.n_qubits, lib)
    }

    pub fn run(
        &self,
        target: &Target,
        lib: &GateLibrary,
        rng: &mut Rng,
        options: &SearchOptions,
    ) -> SearchResult {
        let quality_weight =
            if matches!(options.quality_metric, crate::search::QualityMetric::TCount) {
                options.quality_weight
            } else {
                0.0
            };
        let n = self.phases.len();
        let slots = self.affected.len();
        let roots = roots();
        let mut coefficients: Vec<u8> = (0..n - 1).map(|_| rng.usize(8) as u8).collect();
        coefficients.resize(slots, 0);
        let mut nonlinear_cost = 0;
        if !self.nonlinear.is_empty() {
            let pos = n - 1 + rng.usize(self.nonlinear.len());
            let value = 1 + rng.usize(7) as u8;
            coefficients[pos] = value;
            nonlinear_cost = 8 + value as usize % 2;
        }
        let mut odd: Vec<_> = (0..n - 1).filter(|&p| coefficients[p] % 2 == 1).collect();
        while odd.len() + nonlinear_cost > self.inner_t_limit {
            let p = odd.swap_remove(rng.usize(odd.len()));
            coefficients[p] = (rng.usize(4) * 2) as u8;
        }
        let mut tcount = odd.len() + nonlinear_cost;
        let mut errors: Vec<u8> = self.phases.iter().map(|&p| (8 - p) % 8).collect();
        for (p, &coefficient) in coefficients.iter().enumerate() {
            for &x in &self.affected[p] {
                errors[x] = (errors[x] + coefficient) % 8;
            }
        }
        let mut histogram = [0usize; 8];
        for &error in &errors {
            histogram[error as usize] += 1;
        }
        let overlap =
            |h: &[usize; 8]| (0..8).fold(Complex::ZERO, |z, k| z + roots[k] * (h[k] as f64));
        let initial_eq2 = 2.0 * (1.0 - overlap(&histogram).abs() / n as f64).max(0.0);
        let mut best_energy = initial_eq2.powf(options.cost_power) + quality_weight * tcount as f64;
        let mut best = coefficients.clone();
        let mut last_improvement = 0;
        let mut steps = 0;
        let mut accepted = 0;
        let sweeps = (options.iterations_factor * self.n_qubits as f64).ceil() as usize;
        for sweep in 0..sweeps {
            if options.is_cancelled() {
                break;
            }
            let temperature = options.start_temp_base / (n as f64).sqrt()
                * options.sweep_cooling.powi(sweep as i32);
            for offset in 0..slots {
                steps += 1;
                let pos = if sweep % 2 == 0 {
                    offset
                } else {
                    slots - 1 - offset
                };
                let old = coefficients[pos];
                let overhead = if pos < n - 1 { 0 } else { 8 };
                let cost = |value: usize| value % 2 + overhead * usize::from(value != 0);
                let others = tcount - cost(old as usize);
                let mut affected = [0usize; 8];
                for &x in &self.affected[pos] {
                    affected[errors[x] as usize] += 1;
                }
                let total = overlap(&histogram);
                let changing = overlap(&affected);
                let mut energies = [f64::INFINITY; 8];
                for (value, energy) in energies.iter_mut().enumerate() {
                    let count = others + cost(value);
                    if count > self.inner_t_limit {
                        continue;
                    }
                    let delta = (value + 8 - old as usize) % 8;
                    if histogram[0] - affected[0] + affected[(8 - delta) % 8] == n {
                        coefficients[pos] = value as u8;
                        let circuit = self.scheduled_circuit(&coefficients, lib, rng);
                        let eq = target.original.cost(circuit.matrix(), options.simple_cost);
                        if target
                            .original
                            .exact_cost(circuit.matrix(), options.epsilon)
                            == 0.0
                        {
                            return SearchResult {
                                circuit,
                                best_eq: eq,
                                best_energy: eq.powf(2.0 * options.cost_power)
                                    + quality_weight * count as f64,
                                found: true,
                                steps,
                                accepted,
                            };
                        }
                        coefficients[pos] = old;
                    }
                    let trace = total + (roots[delta] - Complex::ONE) * changing;
                    let eq2 = 2.0 * (1.0 - trace.abs() / n as f64).max(0.0);
                    *energy = eq2.powf(options.cost_power) + quality_weight * count as f64;
                }
                let min = energies.iter().copied().fold(f64::INFINITY, f64::min);
                let mut weights = [0.0; 8];
                let mut sum = 0.0;
                for value in 0..8 {
                    let prior = if value == 0 {
                        options.pid
                    } else {
                        (1.0 - options.pid) / 7.0
                    };
                    weights[value] = (-(energies[value] - min) / temperature.max(1e-14)).exp()
                        * prior.max(1e-12).powf(options.gate_prior);
                    sum += weights[value];
                }
                let mut draw = rng.random01() * sum;
                let mut chosen = 7;
                for (value, &weight) in weights.iter().enumerate() {
                    draw -= weight;
                    if draw <= 0.0 {
                        chosen = value;
                        break;
                    }
                }
                let delta = (chosen + 8 - old as usize) % 8;
                for &x in &self.affected[pos] {
                    histogram[errors[x] as usize] -= 1;
                    errors[x] = (errors[x] + delta as u8) % 8;
                    histogram[errors[x] as usize] += 1;
                }
                coefficients[pos] = chosen as u8;
                tcount = others + cost(chosen);
                accepted += usize::from(delta != 0);
                if energies[chosen] + 1e-12 < best_energy {
                    best_energy = energies[chosen];
                    best.clone_from(&coefficients);
                    last_improvement = sweep;
                }
            }
            if options.stall_sweeps > 0 && sweep >= last_improvement + options.stall_sweeps {
                break;
            }
        }
        let mut circuit = self.circuit(&best, lib);
        if target
            .original
            .exact_cost(circuit.matrix(), options.epsilon)
            == 0.0
        {
            circuit = self.scheduled_circuit(&best, lib, rng);
        }
        let eq = target.original.cost(circuit.matrix(), options.simple_cost);
        let found = target
            .original
            .exact_cost(circuit.matrix(), options.epsilon)
            == 0.0;
        SearchResult {
            circuit,
            best_eq: eq,
            best_energy,
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
    fn library(n: usize) -> GateLibrary {
        GateLibrary::load(n, "data/gates/CliffordT", "").unwrap()
    }

    #[test]
    fn parity_budget_bound_matches_every_small_binary_coefficient_vector() {
        for q in 1..=4 {
            let n = 1usize << q;
            for vector in 0u32..(1 << (n - 1)) {
                let mut phases = vec![0u8; n];
                for (x, phase) in phases.iter_mut().enumerate() {
                    for mask in 1..n {
                        if (x & mask).count_ones() % 2 == 1 {
                            // Even coefficients exercise lifting from the binary
                            // signature back to all eight phase values.
                            let coefficient =
                                ((vector >> (mask - 1)) & 1) as u8 + 2 * (mask % 4) as u8;
                            *phase = (*phase + coefficient) % 8;
                        }
                    }
                }
                let count = vector.count_ones() as usize;
                let expected = if q == 4 {
                    count.min(n - 1 - count)
                } else {
                    count
                };
                assert_eq!(
                    parity_minimum_t(&phases, q),
                    Some(expected),
                    "q={q} vector={vector}"
                );
            }
        }
    }

    #[test]
    fn parity_kernel_rejects_impossible_t_budgets_but_accepts_clifford_phases() {
        let lib = library(3);
        let native = NativeClifford::new(&lib).unwrap();
        let r = roots();
        let mut diagonal = vec![Complex::ONE; 8];
        diagonal[7] = -Complex::ONE;
        assert!(PhaseKernel::from_diagonal(&diagonal, &native, 6).is_none());
        assert!(PhaseKernel::from_diagonal(&diagonal, &native, 7).is_some());
        for (x, value) in diagonal.iter_mut().enumerate() {
            *value = r[(2 * (x & 1) + 4 * usize::from(x & 6 == 6)) % 8];
        }
        assert!(PhaseKernel::from_diagonal(&diagonal, &native, 0).is_some());
    }

    #[test]
    fn cached_diagonal_initialization_preserves_seeded_search() {
        let mut source = Rng::new(68120);
        let r = roots();
        for q in 1..=4 {
            let lib = library(q);
            let native = NativeClifford::new(&lib).unwrap();
            let n = 1 << q;
            for sample in 0..12 {
                let coefficients: Vec<u8> = (1..n).map(|_| source.usize(8) as u8).collect();
                let diagonal: Vec<_> = (0usize..n)
                    .map(|x| {
                        let power: usize = coefficients
                            .iter()
                            .enumerate()
                            .filter(|(k, _)| ((k + 1) & x).count_ones() % 2 == 1)
                            .map(|(_, &c)| c as usize)
                            .sum();
                        r[(power + sample) % 8]
                    })
                    .collect();
                let mut matrix = Matrix::zero(n);
                for x in 0..n {
                    matrix[(x, x)] = diagonal[x];
                }
                let target = Target::new(
                    PartialMatrix::new(matrix, vec![true; n * n], "cached diagonal").unwrap(),
                    false,
                    false,
                );
                let original = PhaseKernel::new(&target, &lib, n).unwrap();
                let cached = PhaseKernel::from_diagonal(&diagonal, &native, n).unwrap();
                assert_eq!(original.phases, cached.phases);
                assert_eq!(original.affected, cached.affected);
                let options = SearchOptions::tuned();
                let mut a = Rng::new(sample as u64);
                let mut b = Rng::new(sample as u64);
                let first = original.run(&target, &lib, &mut a, &options);
                let second = cached.run(&target, &lib, &mut b, &options);
                assert_eq!(
                    (first.steps, first.accepted, first.found),
                    (second.steps, second.accepted, second.found)
                );
                assert_eq!(first.circuit.gates, second.circuit.gates);
                assert_eq!(a.random01(), b.random01());
            }
        }
    }

    #[test]
    fn two_sided_clifford_preprocessing_preserves_the_t_budget() {
        let mut rng = Rng::new(177625);
        for q in 1..=4 {
            let lib = library(q);
            let n = 1 << q;
            let identity = Target::new(
                PartialMatrix::new(Matrix::identity(n), vec![true; n * n], "id").unwrap(),
                false,
                false,
            );
            let factory = PhaseKernel::new(&identity, &lib, 4).unwrap();
            let ids: Vec<_> = lib
                .all
                .iter()
                .copied()
                .filter(|&id| matches!(lib.gates[id].name.as_str(), "h" | "s" | "sdg" | "cx"))
                .collect();
            for sample in 0..8 {
                let mut coeff = vec![0; n - 1];
                for _ in 0..4 {
                    let p = rng.usize(n - 1);
                    coeff[p] = (coeff[p] + 1) % 8;
                }
                let left: Circuit = Circuit::new(
                    (0..16).map(|_| ids[rng.usize(ids.len())]).collect(),
                    q,
                    &lib,
                );
                let right: Circuit = Circuit::new(
                    (0..16).map(|_| ids[rng.usize(ids.len())]).collect(),
                    q,
                    &lib,
                );
                let phase = factory.circuit(&coeff, &lib);
                let matrix = left.matrix().mul(phase.matrix()).mul(right.matrix());
                let target = Target::new(
                    PartialMatrix::new(matrix, vec![true; n * n], "two sided phase").unwrap(),
                    false,
                    false,
                );
                let kernel = PhaseKernel::new_factored(&target, &lib, 4).unwrap();
                let options = SearchOptions::tuned();
                let mut solved = false;
                for _ in 0..100 {
                    let result = kernel.run(&target, &lib, &mut rng, &options);
                    if result.found {
                        assert_eq!(
                            target.original.exact_cost(result.circuit.matrix(), 1e-6),
                            0.0
                        );
                        assert!(result.circuit.count(&["t".into(), "tdg".into()], &lib) <= 4);
                        solved = true;
                        break;
                    }
                }
                assert!(solved, "q{q} sample{sample}");
            }
        }
    }

    #[test]
    fn affine_wrapped_controlled_h_retains_the_nine_t_goal() {
        let mut rng = Rng::new(22671);
        let base = PartialMatrix::read("data/input/64/comparison/cch.txt").unwrap();
        for q in 3..=4 {
            let lib = library(q);
            let native = NativeClifford::new(&lib).unwrap();
            let base = base.matrix.kron_identity(q);
            for sample in 0..16 {
                let mut words = [Vec::new(), Vec::new()];
                for word in &mut words {
                    for _ in 0..12 {
                        let a = rng.usize(q);
                        let mut b = rng.usize(q - 1);
                        if b >= a {
                            b += 1;
                        }
                        word.push(native.cx[a][b].unwrap());
                        if rng.random01() < 0.2 {
                            let ids = native.native[b];
                            word.extend([ids[0], ids[1], ids[1], ids[0]]);
                        }
                    }
                }
                let left = Circuit::new(words[0].clone(), q, &lib);
                let right = Circuit::new(words[1].clone(), q, &lib);
                let matrix = left.matrix().mul(&base).mul(right.matrix());
                let target = Target::new(
                    PartialMatrix::new(matrix, vec![true; 1 << (2 * q)], "affine CCH").unwrap(),
                    false,
                    false,
                );
                let kernel = PhaseKernel::new_factored(&target, &lib, 9)
                    .unwrap_or_else(|| panic!("no basis q{q} sample{sample}"));
                let options = SearchOptions::tuned();
                let mut solved = false;
                for _ in 0..100 {
                    let result = kernel.run(&target, &lib, &mut rng, &options);
                    if result.found {
                        assert_eq!(
                            target.original.exact_cost(result.circuit.matrix(), 1e-6),
                            0.0
                        );
                        assert!(result.circuit.count(&["t".into(), "tdg".into()], &lib) <= 9);
                        solved = true;
                        break;
                    }
                }
                assert!(solved, "q{q} sample{sample}");
            }
        }
    }

    #[test]
    fn nonlinear_boolean_phase_gadgets_match_their_truth_tables() {
        for q in 3..=4 {
            let lib = library(q);
            let n = 1 << q;
            let source = PartialMatrix::read("data/input/65/ccrz_2.txt").unwrap();
            let diagonal: Vec<_> = (0..n).map(|x| source.matrix[(x % 8, x % 8)]).collect();
            let mut matrix = Matrix::zero(n);
            for x in 0..n {
                matrix[(x, x)] = diagonal[x];
            }
            let target = Target::new(
                PartialMatrix::new(matrix, vec![true; n * n], "nonlinear").unwrap(),
                false,
                false,
            );
            let kernel = PhaseKernel::new_nonlinear(&target, &lib, 12).unwrap();
            for index in 0..kernel.nonlinear.len() {
                for value in 1..8 {
                    let mut coefficients = vec![0; kernel.affected.len()];
                    coefficients[n - 1 + index] = value;
                    let circuit = kernel.circuit(&coefficients, &lib);
                    let mut expected = Matrix::identity(n);
                    for &x in &kernel.affected[n - 1 + index] {
                        expected[(x, x)] = roots()[value as usize];
                    }
                    assert!(same_phase(&expected, circuit.matrix()));
                    assert_eq!(
                        circuit.count(&["t".into(), "tdg".into()], &lib),
                        8 + value as usize % 2
                    );
                }
            }
        }
    }
    #[test]
    fn nonlinear_phase_annealing_reaches_the_controlled_rotation_bound() {
        let lib = library(3);
        let target = Target::new(
            PartialMatrix::read("data/input/65/ccrz_2.txt").unwrap(),
            false,
            false,
        );
        let kernel = PhaseKernel::new_nonlinear(&target, &lib, 12).unwrap();
        let options = SearchOptions::tuned();
        for seed in 0..16 {
            let mut rng = Rng::new(seed);
            let mut solved = false;
            for _ in 0..256 {
                let result = kernel.run(&target, &lib, &mut rng, &options);
                if result.found {
                    assert_eq!(
                        target.original.exact_cost(result.circuit.matrix(), 1e-6),
                        0.0
                    );
                    assert!(result.circuit.count(&["t".into(), "tdg".into()], &lib) <= 12);
                    solved = true;
                    break;
                }
            }
            assert!(solved, "seed {seed}");
        }
    }

    #[test]
    fn nonlinear_phase_bases_support_random_affine_control_parities() {
        let mut source = Rng::new(447928);
        let r = roots();
        for q in 3..=4 {
            let lib = library(q);
            let n = 1 << q;
            for sample in 0..12 {
                let mut rows: Vec<_> = (0..q).map(|i| 1usize << i).collect();
                for _ in 0..12 {
                    let a = source.usize(q);
                    let mut b = source.usize(q - 1);
                    if b >= a {
                        b += 1;
                    }
                    rows[b] ^= rows[a];
                }
                let mut matrix = Matrix::identity(n);
                for x in 0usize..n {
                    let bit = |v: usize| ((x & v).count_ones() % 2) as u8;
                    let f = bit(rows[2]) ^ (bit(rows[0]) & bit(rows[1]));
                    matrix[(x, x)] = r[f as usize];
                }
                let target = Target::new(
                    PartialMatrix::new(matrix, vec![true; n * n], "quadratic basis").unwrap(),
                    false,
                    false,
                );
                let kernel = PhaseKernel::new_nonlinear(&target, &lib, 9).unwrap();
                let mut rng = Rng::new(sample);
                let mut solved = false;
                for _ in 0..256 {
                    let result = kernel.run(&target, &lib, &mut rng, &SearchOptions::tuned());
                    if result.found {
                        assert_eq!(
                            target.original.exact_cost(result.circuit.matrix(), 1e-6),
                            0.0
                        );
                        assert!(result.circuit.count(&["t".into(), "tdg".into()], &lib) <= 9);
                        solved = true;
                        break;
                    }
                }
                assert!(solved, "q{q} sample{sample}");
            }
        }
    }

    #[test]
    fn local_basis_catalog_has_nine_unoriented_axes() {
        let axes = local_axes();
        assert_eq!(axes.len(), 9);
        assert_eq!(axes.iter().filter(|b| b.tcount == 0).count(), 3);
        for b in axes {
            assert!(b
                .matrix
                .adjoint()
                .mul(&b.matrix)
                .approximately_equal(&Matrix::identity(2), 1e-12));
        }
    }
    #[test]
    fn parity_criterion_rejects_phases_outside_the_integer_gadget_span() {
        assert!(parity_representable(&[0, 0, 0, 4], 2)); // CZ
        assert!(parity_representable(&[0, 0, 0, 0, 0, 0, 0, 4], 3)); // CCZ
        assert!(!parity_representable(&[0, 0, 0, 1], 2)); // controlled T
        let mut cccz = vec![0; 16];
        cccz[15] = 4;
        assert!(!parity_representable(&cccz, 4));
    }
    #[test]
    fn gray_order_native_circuits_match_integer_phase_polynomials() {
        for nq in 1..=4 {
            let lib = library(nq);
            let n = 1 << nq;
            let target = Target::new(
                PartialMatrix::new(Matrix::identity(n), vec![true; n * n], "id").unwrap(),
                false,
                false,
            );
            let kernel = PhaseKernel::new(&target, &lib, n).unwrap();
            let roots = roots();
            let mut rng = Rng::new(819);
            for _ in 0..40 {
                let coefficients: Vec<_> = (1..n).map(|_| rng.usize(8) as u8).collect();
                let c = kernel.circuit(&coefficients, &lib);
                let mut expected = Matrix::zero(n);
                for x in 0usize..n {
                    let phase = (1usize..n)
                        .filter(|&mask| (x & mask).count_ones() % 2 == 1)
                        .map(|mask| coefficients[mask - 1] as usize)
                        .sum::<usize>()
                        % 8;
                    expected[(x, x)] = roots[phase];
                }
                assert!(same_phase(c.matrix(), &expected));
                assert_eq!(
                    c.count(&["t".into(), "tdg".into()], &lib),
                    coefficients.iter().filter(|&&v| v % 2 == 1).count()
                );
            }
        }
    }
    #[test]
    fn basis_changed_phase_search_reaches_development_optimal_t_counts() {
        let lib = library(3);
        for (name, limit) in [("ccx", 7), ("cch", 9)] {
            let target = Target::new(
                PartialMatrix::read(format!("data/input/64/comparison/{name}.txt")).unwrap(),
                false,
                false,
            );
            let kernel = PhaseKernel::new(&target, &lib, limit).unwrap();
            let options = SearchOptions::tuned();
            for seed in 0..16 {
                let mut rng = Rng::new(seed);
                let mut found = false;
                for _ in 0..50 {
                    let result = kernel.run(&target, &lib, &mut rng, &options);
                    if result.found {
                        assert_eq!(
                            target
                                .original
                                .exact_cost(result.circuit.matrix(), options.epsilon),
                            0.0
                        );
                        assert!(result.circuit.count(&["t".into(), "tdg".into()], &lib) <= limit);
                        found = true;
                        break;
                    }
                }
                assert!(found, "{name}, seed{seed}");
            }
        }
    }
    #[test]
    fn affine_wrapped_nonlinear_permutations_retain_seven_t_solution() {
        let lib = library(3);
        let ccx = PartialMatrix::read("data/input/64/comparison/ccx.txt").unwrap();
        let mut rng = Rng::new(7891);
        let mut affine_gates = Vec::new();
        for a in 0..3 {
            for b in 0..3 {
                if a != b {
                    affine_gates.push(vec![lib.find("cx", &[a, b]).unwrap()]);
                }
            }
            affine_gates.push(
                ["h", "s", "s", "h"]
                    .map(|name| lib.find(name, &[a]).unwrap())
                    .to_vec(),
            );
        }
        for sample in 0..24 {
            let left = Circuit::new(
                (0..8)
                    .flat_map(|_| affine_gates[rng.usize(affine_gates.len())].clone())
                    .collect(),
                3,
                &lib,
            );
            let right = Circuit::new(
                (0..8)
                    .flat_map(|_| affine_gates[rng.usize(affine_gates.len())].clone())
                    .collect(),
                3,
                &lib,
            );
            let matrix = left.matrix().mul(&ccx.matrix).mul(right.matrix());
            let target = Target::new(
                PartialMatrix::new(matrix, vec![true; 64], "affine wrapped").unwrap(),
                false,
                false,
            );
            let kernel = PhaseKernel::new(&target, &lib, 7).unwrap();
            let options = SearchOptions::tuned();
            let found = (0..50)
                .find_map(|_| {
                    let result = kernel.run(&target, &lib, &mut rng, &options);
                    result.found.then_some(result.circuit)
                })
                .unwrap_or_else(|| panic!("sample {sample}"));
            assert_eq!(target.original.exact_cost(found.matrix(), 1e-6), 0.0);
            assert_eq!(found.count(&["t".into(), "tdg".into()], &lib), 7);
        }
    }

    #[test]
    fn phase_solutions_can_meet_optimal_t_depth() {
        let lib = library(3);
        for (name, depth) in [("ccx", 3), ("cch", 4)] {
            let target = Target::new(
                PartialMatrix::read(format!("data/input/64/comparison/{name}.txt")).unwrap(),
                false,
                false,
            );
            let kernel = PhaseKernel::new(&target, &lib, 3 * depth)
                .unwrap()
                .with_depth_limit(depth);
            let options = SearchOptions::tuned();
            for seed in 0..16 {
                let mut rng = Rng::new(seed);
                let mut best_depth = usize::MAX;
                for _ in 0..50 {
                    let result = kernel.run(&target, &lib, &mut rng, &options);
                    if result.found {
                        assert_eq!(
                            target.original.exact_cost(result.circuit.matrix(), 1e-6),
                            0.0
                        );
                        best_depth =
                            best_depth.min(result.circuit.depth(&["t".into(), "tdg".into()], &lib));
                        if best_depth <= depth {
                            break;
                        }
                    }
                }
                assert!(
                    best_depth <= depth,
                    "{name} seed {seed}: depth {best_depth}"
                );
            }
        }
    }

    #[test]
    fn four_qubit_adder_phase_search_reaches_count_and_depth_targets() {
        let lib = library(4);
        let target = Target::new(
            PartialMatrix::read("data/input/64/comparison/adder.txt").unwrap(),
            false,
            false,
        );
        for depth in [None, Some(2)] {
            let limit = if depth.is_some() { 8 } else { 7 };
            let mut kernel = PhaseKernel::new(&target, &lib, limit).unwrap();
            if let Some(depth) = depth {
                kernel = kernel.with_depth_limit(depth);
            }
            let options = SearchOptions::tuned();
            for seed in 0..8 {
                let mut rng = Rng::new(seed);
                let mut solved = false;
                for _ in 0..500 {
                    let result = kernel.run(&target, &lib, &mut rng, &options);
                    if !result.found {
                        continue;
                    }
                    assert_eq!(
                        target.original.exact_cost(result.circuit.matrix(), 1e-6),
                        0.0
                    );
                    assert!(result.circuit.count(&["t".into(), "tdg".into()], &lib) <= limit);
                    if depth.is_none_or(|depth| {
                        result.circuit.depth(&["t".into(), "tdg".into()], &lib) <= depth
                    }) {
                        solved = true;
                        break;
                    }
                }
                assert!(solved, "seed {seed}, depth {depth:?}");
            }
        }
    }

    #[test]
    fn entangled_commuting_phase_functions_use_native_clifford_wrappers() {
        let mut rng = Rng::new(99812);
        for q in 2..=4 {
            let lib = library(q);
            let cliffords: Vec<_> = lib
                .all
                .iter()
                .copied()
                .filter(|&id| matches!(lib.gates[id].name.as_str(), "h" | "s" | "sdg" | "cx"))
                .collect();
            for _ in 0..8 {
                let c = Circuit::new(
                    (0..18)
                        .map(|_| cliffords[rng.usize(cliffords.len())])
                        .collect(),
                    q,
                    &lib,
                );
                let d = Circuit::new(
                    vec![lib.find("t", &[0]).unwrap(), lib.find("tdg", &[1]).unwrap()],
                    q,
                    &lib,
                );
                let matrix = c.matrix().mul(d.matrix()).mul(&c.matrix().adjoint());
                let target = Target::new(
                    PartialMatrix::new(matrix, vec![true; 1 << (2 * q)], "entangled phase")
                        .unwrap(),
                    false,
                    false,
                );
                let kernel = PhaseKernel::new(&target, &lib, 2).unwrap();
                let options = SearchOptions::tuned();
                let solution = (0..500)
                    .find_map(|_| {
                        let r = kernel.run(&target, &lib, &mut rng, &options);
                        r.found.then_some(r.circuit)
                    })
                    .expect("phase synthesis");
                assert_eq!(target.original.exact_cost(solution.matrix(), 1e-6), 0.0);
                assert!(solution.count(&["t".into(), "tdg".into()], &lib) <= 2);
            }
        }
    }

    #[test]
    fn arbitrary_rotation_and_partial_specs_keep_general_search() {
        let lib = library(1);
        let t = Target::new(
            PartialMatrix::read("data/input/65/rz_4.txt").unwrap(),
            false,
            false,
        );
        assert!(PhaseKernel::new(&t, &lib, 7).is_none());
        let p = PartialMatrix::new(Matrix::identity(2), vec![true, false, true, false], "state")
            .unwrap();
        assert!(PhaseKernel::new(&Target::new(p, false, false), &lib, 0).is_none());
    }
}
