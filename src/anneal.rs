//! Coordinate heat-bath annealing with reusable contraction environments.
//! For a unitary circuit U=S G P, Tr(T† U)=Tr(P T† S G). Moving
//! to an adjacent site updates this environment with two sparse gate products.
#[path = "masked.rs"]
mod masked;

use crate::{
    circuit::Circuit,
    complex::Complex,
    gates::GateLibrary,
    matrix::Matrix,
    partial::{PartialMatrix, Target},
    rng::Rng,
    search::{random_circuit, QualityMetric, SearchOptions, SearchResult},
};

fn identity_trace(matrix: &Matrix) -> Complex {
    (0..matrix.n).fold(Complex::ZERO, |sum, k| sum + matrix[(k, k)])
}

#[derive(Clone)]
enum TraceTerms {
    Unit(Vec<usize>),
    Real(Vec<(usize, f64)>),
    Complex(Vec<(usize, Complex)>),
}
impl TraceTerms {
    fn new(entries: &[(usize, usize, Complex)], n: usize) -> Self {
        if entries.iter().all(|&(_, _, v)| v == Complex::ONE) {
            Self::Unit(entries.iter().map(|&(r, c, _)| c * n + r).collect())
        } else if entries.iter().all(|&(_, _, v)| v.im == 0.0) {
            Self::Real(entries.iter().map(|&(r, c, v)| (c * n + r, v.re)).collect())
        } else {
            Self::Complex(entries.iter().map(|&(r, c, v)| (c * n + r, v)).collect())
        }
    }
    fn constant_factor(&self) -> Option<(Vec<usize>, Complex)> {
        match self {
            Self::Unit(indices) if !indices.is_empty() => Some((indices.clone(), Complex::ONE)),
            Self::Real(entries) => {
                let factor = entries.first()?.1;
                entries
                    .iter()
                    .all(|&(_, v)| v == factor)
                    .then(|| (entries.iter().map(|&(i, _)| i).collect(), factor.into()))
            }
            Self::Complex(entries) => {
                let factor = entries.first()?.1;
                entries
                    .iter()
                    .all(|&(_, v)| v == factor)
                    .then(|| (entries.iter().map(|&(i, _)| i).collect(), factor))
            }
            _ => None,
        }
    }
    fn work(&self) -> usize {
        match self {
            Self::Unit(v) => 2 * v.len(),
            Self::Real(v) => 4 * v.len(),
            Self::Complex(v) => 8 * v.len(),
        }
    }
}
#[derive(Clone)]
enum ProductKind {
    Identity,
    Permutation {
        rows: Vec<usize>,
        columns: Vec<usize>,
    },
    Diagonal(Vec<Complex>),
    General,
}
#[derive(Clone)]
struct SparseGate {
    product: ProductKind,
    rows: Vec<Vec<(usize, Complex)>>,
    entries: Vec<(usize, usize, Complex)>,
    trace_terms: TraceTerms,
    trace_identity: bool,
}
impl SparseGate {
    fn new(m: &Matrix) -> Self {
        let mut rows = vec![Vec::new(); m.n];
        let mut entries = Vec::new();
        for r in 0..m.n {
            for c in 0..m.n {
                let z = m[(r, c)];
                if z != Complex::ZERO {
                    rows[r].push((c, z));
                    entries.push((r, c, z));
                }
            }
        }
        let product = if m.data.iter().enumerate().all(|(i, &z)| {
            z == if i / m.n == i % m.n {
                Complex::ONE
            } else {
                Complex::ZERO
            }
        }) {
            ProductKind::Identity
        } else if entries.iter().all(|&(r, c, _)| r == c) {
            ProductKind::Diagonal((0..m.n).map(|i| m[(i, i)]).collect())
        } else if rows
            .iter()
            .all(|row| row.len() == 1 && row[0].1 == Complex::ONE)
        {
            let row_map: Vec<_> = rows.iter().map(|row| row[0].0).collect();
            let mut columns = vec![usize::MAX; m.n];
            for (r, &c) in row_map.iter().enumerate() {
                columns[c] = r;
            }
            if columns.contains(&usize::MAX) {
                ProductKind::General
            } else {
                ProductKind::Permutation {
                    rows: row_map,
                    columns,
                }
            }
        } else {
            ProductKind::General
        };
        let full = TraceTerms::new(&entries, m.n);
        let mut delta = Vec::new();
        for r in 0..m.n {
            for c in 0..m.n {
                let value = m[(r, c)] - if r == c { Complex::ONE } else { Complex::ZERO };
                if value != Complex::ZERO {
                    delta.push((r, c, value));
                }
            }
        }
        let delta = TraceTerms::new(&delta, m.n);
        let trace_identity = delta.work() < full.work();
        let trace_terms = if trace_identity { delta } else { full };
        Self {
            product,
            rows,
            entries,
            trace_terms,
            trace_identity,
        }
    }
    fn trace(&self, e: &Matrix, identity: Complex) -> f64 {
        let mut z = if self.trace_identity {
            identity
        } else {
            Complex::ZERO
        };
        match &self.trace_terms {
            TraceTerms::Unit(indices) => {
                for &index in indices {
                    z += e.data[index];
                }
            }
            TraceTerms::Real(entries) => {
                for &(index, value) in entries {
                    z += e.data[index] * value;
                }
            }
            TraceTerms::Complex(entries) => {
                for &(index, value) in entries {
                    z += e.data[index] * value;
                }
            }
        }
        z.norm_sqr().sqrt()
    }

    fn left(&self, m: &Matrix, out: &mut Matrix) {
        let n = m.n;
        match &self.product {
            ProductKind::Identity => {
                out.data.copy_from_slice(&m.data);
                return;
            }
            ProductKind::Permutation { rows, .. } => {
                for (r, &k) in rows.iter().enumerate() {
                    out.data[r * n..(r + 1) * n].copy_from_slice(&m.data[k * n..(k + 1) * n]);
                }
                return;
            }
            ProductKind::Diagonal(phases) => {
                for (r, &a) in phases.iter().enumerate() {
                    let dst = &mut out.data[r * n..(r + 1) * n];
                    let src = &m.data[r * n..(r + 1) * n];
                    if a == Complex::ONE {
                        dst.copy_from_slice(src);
                    } else if a.im == 0.0 {
                        for (d, &b) in dst.iter_mut().zip(src) {
                            *d = b * a.re;
                        }
                    } else if a.re == 0.0 {
                        for (d, b) in dst.iter_mut().zip(src) {
                            d.re = -a.im * b.im;
                            d.im = a.im * b.re;
                        }
                    } else {
                        for (d, &b) in dst.iter_mut().zip(src) {
                            *d = a * b;
                        }
                    }
                }
                return;
            }
            ProductKind::General => (),
        }
        out.data.fill(Complex::ZERO);
        for (r, entries) in self.rows.iter().enumerate() {
            let dst = &mut out.data[r * n..(r + 1) * n];
            for &(k, a) in entries {
                let src = &m.data[k * n..(k + 1) * n];
                if a.im == 0.0 {
                    for (d, b) in dst.iter_mut().zip(src) {
                        d.re += a.re * b.re;
                        d.im += a.re * b.im;
                    }
                } else {
                    for (d, b) in dst.iter_mut().zip(src) {
                        d.re += a.re * b.re - a.im * b.im;
                        d.im += a.re * b.im + a.im * b.re;
                    }
                }
            }
        }
    }
    fn right(&self, m: &Matrix, out: &mut Matrix) {
        let n = m.n;
        match &self.product {
            ProductKind::Identity => {
                out.data.copy_from_slice(&m.data);
                return;
            }
            ProductKind::Permutation { columns, .. } => {
                for (src, dst) in m.data.chunks_exact(n).zip(out.data.chunks_exact_mut(n)) {
                    for (d, &k) in dst.iter_mut().zip(columns) {
                        *d = src[k];
                    }
                }
                return;
            }
            ProductKind::Diagonal(phases) => {
                for (src, dst) in m.data.chunks_exact(n).zip(out.data.chunks_exact_mut(n)) {
                    for ((&b, d), &a) in src.iter().zip(dst).zip(phases) {
                        *d = b * a;
                    }
                }
                return;
            }
            ProductKind::General => (),
        }
        out.data.fill(Complex::ZERO);
        for r in 0..n {
            let src = &m.data[r * n..(r + 1) * n];
            let dst = &mut out.data[r * n..(r + 1) * n];
            for &(k, c, a) in &self.entries {
                dst[c] += src[k] * a;
            }
        }
    }
}

/// Whole rows or columns of a unitary have constant total squared norm.
pub(crate) fn constant_norm(p: &PartialMatrix) -> Option<f64> {
    let n = p.matrix.n;
    let columns = (0..n).all(|c| (0..n).all(|r| p.cover[r * n + c] == p.cover[c]));
    let rows = p
        .cover
        .chunks_exact(n)
        .all(|row| row.iter().all(|&b| b == row[0]));
    if columns || rows {
        Some(p.n_constraints() as f64 / n as f64)
    } else {
        None
    }
}

#[derive(Clone, Copy)]
struct EqualityForm {
    squared_norm: f64,
    normalization: f64,
    simple: bool,
}
impl EqualityForm {
    fn new(p: &PartialMatrix, circuit_norm: f64, simple: bool) -> Self {
        Self {
            squared_norm: p.squared_norm + circuit_norm,
            normalization: (p.n_constraints() as f64).sqrt(),
            simple,
        }
    }
    fn evaluate(self, overlap: f64) -> f64 {
        if self.simple {
            2.0 * (1.0 - overlap / self.normalization).max(0.0)
        } else {
            (self.squared_norm - 2.0 * overlap).max(0.0) / self.normalization
        }
    }
}
#[cfg(test)]
fn equality_squared(p: &PartialMatrix, circuit_norm: f64, overlap: f64, simple: bool) -> f64 {
    EqualityForm::new(p, circuit_norm, simple).evaluate(overlap)
}

pub struct SweepKernel {
    trace_groups: Vec<Vec<usize>>,
    trace_reuse: Vec<Option<(usize, Complex)>>,
    sparse: Vec<SparseGate>,
    adjoints: Vec<SparseGate>,
    choices: Vec<usize>,
    priors: Vec<f64>,
    qualities: Vec<f64>,
    depth: Option<crate::quality::DepthModel>,
}
impl SweepKernel {
    pub fn new(lib: &GateLibrary, options: &SearchOptions) -> Self {
        let sparse: Vec<_> = lib
            .gates
            .iter()
            .map(|g| SparseGate::new(&g.matrix))
            .collect();
        let adjoints = lib
            .gates
            .iter()
            .map(|g| SparseGate::new(&g.matrix.adjoint()))
            .collect();
        let mut choices: Vec<usize> = Vec::new();
        for &id in &lib.all {
            if !choices.iter().any(|&other| {
                let g = &lib.gates[other];
                g.name == lib.gates[id].name && g.qubits == lib.gates[id].qubits
            }) {
                choices.push(id);
            }
        }
        // Phase gates on the same support share a sum of environment entries;
        // only their scalar coefficient differs. Reuse those linear forms.
        let mut groups: Vec<Vec<usize>> = Vec::new();
        let mut uses = Vec::new();
        let mut trace_reuse = vec![None; sparse.len()];
        for &id in &choices {
            if let Some((indices, factor)) = sparse[id].trace_terms.constant_factor() {
                let index = if let Some(index) = groups.iter().position(|g| *g == indices) {
                    index
                } else {
                    groups.push(indices);
                    uses.push(0);
                    groups.len() - 1
                };
                uses[index] += 1;
                trace_reuse[id] = Some((index, factor));
            }
        }
        let mut trace_groups = Vec::new();
        let mut remap = vec![None; groups.len()];
        for (i, indices) in groups.into_iter().enumerate() {
            if uses[i] >= 2 && indices.len() >= 2 {
                remap[i] = Some(trace_groups.len());
                trace_groups.push(indices);
            }
        }
        for reuse in &mut trace_reuse {
            *reuse = reuse.and_then(|(i, factor)| remap[i].map(|index| (index, factor)));
        }
        let priors = choices
            .iter()
            .map(|&id| {
                if options.gate_prior == 0.0 {
                    return 1.0;
                }
                let p = if id == lib.identity {
                    options.pid
                } else {
                    let basic = lib.basic.contains(&id);
                    let weight = if basic { 1.0 } else { options.pcomp };
                    let denom = lib.basic.len() as f64 + options.pcomp * lib.composite.len() as f64;
                    (1.0 - options.pid) * weight / denom
                };
                p.max(1e-12).powf(options.gate_prior)
            })
            .collect();
        let qualities = lib
            .gates
            .iter()
            .enumerate()
            .map(|(id, g)| {
                let ids = if options.quality_expand && !g.decomposition.is_empty() {
                    g.decomposition.clone()
                } else {
                    vec![id]
                };
                ids.iter()
                    .map(|&id| match options.quality_metric {
                        QualityMetric::TCount | QualityMetric::TDepth => {
                            f64::from(matches!(lib.gates[id].name.as_str(), "t" | "tdg"))
                        }
                        QualityMetric::GateCount => f64::from(id != lib.identity),
                        QualityMetric::WeightedCost => lib.gates[id].cost,
                    })
                    .sum()
            })
            .collect();
        Self {
            trace_groups,
            trace_reuse,
            sparse,
            adjoints,
            choices,
            priors,
            qualities,
            depth: (matches!(options.quality_metric, QualityMetric::TDepth)
                && (options.quality_weight > 0.0 || options.quality_slack.is_some()))
            .then(|| crate::quality::DepthModel::new(lib, options.quality_expand)),
        }
    }

    fn fill_trace_groups(&self, env: &Matrix, values: &mut [Complex]) {
        for (sum, indices) in values.iter_mut().zip(&self.trace_groups) {
            *sum = indices.iter().fold(Complex::ZERO, |z, &i| z + env.data[i]);
        }
    }

    fn reused_trace(&self, id: usize, env: &Matrix, identity: Complex, values: &[Complex]) -> f64 {
        if let Some((index, factor)) = self.trace_reuse[id] {
            let trace = values[index] * factor
                + if self.sparse[id].trace_identity {
                    identity
                } else {
                    Complex::ZERO
                };
            trace.norm_sqr().sqrt()
        } else {
            self.sparse[id].trace(env, identity)
        }
    }

    fn circuit_quality(&self, gates: &[usize]) -> f64 {
        self.depth.as_ref().map_or_else(
            || gates.iter().map(|&id| self.qualities[id]).sum(),
            |depth| depth.depth(gates) as f64,
        )
    }

    fn swap_environment(
        &self,
        e: &Matrix,
        out: &mut Matrix,
        scratch: &mut Matrix,
        neighbor: usize,
        forward: bool,
    ) {
        if forward {
            self.adjoints[neighbor].right(e, scratch);
            self.sparse[neighbor].left(scratch, out);
        } else {
            self.adjoints[neighbor].left(e, scratch);
            self.sparse[neighbor].right(scratch, out);
        }
    }

    pub fn supports(target: &Target, lib: &GateLibrary) -> bool {
        target.variants.iter().all(|p| constant_norm(p).is_some())
            && Self::supports_masked(target, lib)
    }
    pub fn supports_masked(target: &Target, lib: &GateLibrary) -> bool {
        !target.inverses.iter().any(|&b| b)
            && lib.gates.iter().all(|g| {
                g.matrix
                    .adjoint()
                    .mul(&g.matrix)
                    .approximately_equal(&Matrix::identity(g.matrix.n), 1e-10)
            })
    }
    pub fn run(
        &self,
        target: &Target,
        lib: &GateLibrary,
        rng: &mut Rng,
        n_gates: usize,
        options: &SearchOptions,
    ) -> Result<SearchResult, String> {
        if constant_norm(&target.original).is_none() {
            return self.run_masked(target, lib, rng, n_gates, options);
        }
        let mut initial = random_circuit(lib, rng, n_gates, options.pid);
        let quality_cap = options
            .quality_slack
            .zip(options.quality_goal)
            .map(|(slack, goal)| slack + goal);
        let mut quality: f64 = self.circuit_quality(&initial.gates);
        if let Some(cap) = quality_cap {
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
        if n_gates == 0 {
            return crate::search::search_from(target, lib, rng, initial, options);
        }
        let mut gates = initial.gates.clone();
        let mut best_gates = gates.clone();
        let specs: Vec<_> = if options.enable_permutations && options.sweep_permutations {
            target.variants.iter().collect()
        } else {
            vec![&target.original]
        };
        let forms: Vec<_> = specs
            .iter()
            .map(|p| EqualityForm::new(p, constant_norm(p).unwrap(), options.simple_cost))
            .collect();
        let mut env: Vec<_> = specs.iter().map(|p| p.matrix.adjoint()).collect();
        let mut swap_env = env.clone();
        let mut identity_traces = vec![Complex::ZERO; env.len()];
        let mut trace_values = vec![vec![Complex::ZERO; self.trace_groups.len()]; env.len()];
        let n = 1 << lib.n_qubits;
        let mut scratch = Matrix::zero(n);
        let mut best_eq = specs
            .iter()
            .map(|p| p.cost(initial.matrix(), options.simple_cost))
            .fold(f64::INFINITY, f64::min);
        let mut best_energy =
            best_eq.powf(2.0 * options.cost_power) + options.quality_weight * quality;
        let mut energies = vec![0.0; self.choices.len()];
        let mut equality = energies.clone();
        let mut weights = energies.clone();
        let mut proposal_quality = energies.clone();
        let mut depth_work = self.depth.as_ref().map(|d| d.workspace(n_gates));
        let sweeps = (options.iterations_factor * lib.n_qubits as f64).ceil() as usize;
        let mut steps = 0;
        let mut accepted = 0;
        let mut last_improvement = 0;
        for sweep in 0..sweeps {
            if options.is_cancelled() {
                break;
            }
            let forward = sweep % 2 == 0;
            if let (Some(model), Some(work)) = (&self.depth, &mut depth_work) {
                work.reset(model, &gates, forward);
            }
            // Periodic exact rebuild bounds drift from repeated inverse updates.
            if sweep % 8 == 0 {
                for (e, p) in env.iter_mut().zip(&specs) {
                    *e = p.matrix.adjoint();
                    for &g in gates[1..].iter().rev() {
                        self.sparse[g].right(e, &mut scratch);
                        std::mem::swap(e, &mut scratch);
                    }
                }
            }
            let temperature = options.start_temp_base / (n as f64).sqrt()
                * options.sweep_cooling.powi(sweep as i32);
            for offset in 0..n_gates {
                let position = if forward {
                    offset
                } else {
                    n_gates - 1 - offset
                };
                steps += 1;
                let others_quality = quality - self.qualities[gates[position]];
                for (trace, e) in identity_traces.iter_mut().zip(&env) {
                    *trace = (0..n).fold(Complex::ZERO, |sum, k| sum + e[(k, k)]);
                }
                for (e, values) in env.iter().zip(&mut trace_values) {
                    self.fill_trace_groups(e, values);
                }
                for (i, &id) in self.choices.iter().enumerate() {
                    proposal_quality[i] =
                        if let (Some(model), Some(work)) = (&self.depth, &depth_work) {
                            work.candidate(model, position, id) as f64
                        } else {
                            others_quality + self.qualities[id]
                        };
                    if quality_cap.is_some_and(|cap| proposal_quality[i] > cap + 1e-12) {
                        equality[i] = f64::INFINITY;
                        energies[i] = f64::INFINITY;
                        continue;
                    }
                    let mut eq2 = f64::INFINITY;
                    for (((e, form), &identity), values) in env
                        .iter()
                        .zip(&forms)
                        .zip(&identity_traces)
                        .zip(&trace_values)
                    {
                        let overlap = self.reused_trace(id, e, identity, values);
                        let value = form.evaluate(overlap);
                        eq2 = eq2.min(value);
                    }
                    equality[i] = eq2;
                    energies[i] = if options.cost_power == 0.5 {
                        eq2.sqrt()
                    } else if options.cost_power == 1.0 {
                        eq2
                    } else {
                        eq2.powf(options.cost_power)
                    } + options.quality_weight * proposal_quality[i];
                }
                let min = energies.iter().copied().fold(f64::INFINITY, f64::min);
                let best_equality = equality.iter().copied().fold(f64::INFINITY, f64::min);
                let best_choice = equality.iter().position(|&v| v == best_equality).unwrap();
                // Validate a possible solution using a freshly built circuit.
                if equality[best_choice]
                    <= (options.epsilon * std::f64::consts::SQRT_2 + 1e-7).powi(2)
                {
                    let old = gates[position];
                    gates[position] = self.choices[best_choice];
                    let candidate = Circuit::new(gates.clone(), lib.n_qubits, lib);
                    if let Some(circuit) = correct(candidate, target, lib, options.epsilon)? {
                        let eq = target.original.cost(circuit.matrix(), options.simple_cost);
                        let energy = eq.powf(2.0 * options.cost_power)
                            + options.quality_weight * self.circuit_quality(&circuit.gates);
                        return Ok(SearchResult {
                            circuit,
                            best_eq: eq,
                            best_energy: energy,
                            found: true,
                            steps,
                            accepted,
                        });
                    }
                    gates[position] = old;
                }
                let mut sum = 0.0;
                for ((w, &energy), &prior) in weights.iter_mut().zip(&energies).zip(&self.priors) {
                    *w = (-(energy - min) / temperature.max(1e-14)).exp() * prior;
                    sum += *w;
                }
                let mut draw = rng.random01() * sum;
                let mut chosen = weights.len() - 1;
                for (i, &w) in weights.iter().enumerate() {
                    draw -= w;
                    if draw <= 0.0 {
                        chosen = i;
                        break;
                    }
                }
                let new = self.choices[chosen];
                accepted += usize::from(gates[position] != new);
                gates[position] = new;
                quality = proposal_quality[chosen];
                let complete_energy = energies[chosen];
                if complete_energy + 1e-12 < best_energy {
                    best_energy = complete_energy;
                    best_eq = equality[chosen].sqrt();
                    best_gates.clone_from(&gates);
                    last_improvement = sweep;
                }
                if offset + 1 < n_gates {
                    let neighbor = if forward { position + 1 } else { position - 1 };
                    let b = gates[neighbor];
                    if options.swap_probability > 0.0
                        && b != new
                        && rng.random01() < options.swap_probability
                    {
                        let mut swap_eq2 = f64::INFINITY;
                        for ((e, out), form) in env.iter().zip(&mut swap_env).zip(&forms) {
                            self.swap_environment(e, out, &mut scratch, b, forward);
                            let overlap = self.sparse[new].trace(out, identity_trace(out));
                            swap_eq2 = swap_eq2.min(form.evaluate(overlap));
                        }
                        let swap_eq = swap_eq2.sqrt();
                        let swap_quality =
                            if let (Some(model), Some(work)) = (&self.depth, &mut depth_work) {
                                work.swapped(model, position, new, b) as f64
                            } else {
                                quality
                            };
                        let swap_energy = swap_eq2.powf(options.cost_power)
                            + options.quality_weight * swap_quality;
                        if !quality_cap.is_some_and(|cap| swap_quality > cap + 1e-12)
                            && crate::search::accept_mutation(
                                rng.random01(),
                                swap_energy,
                                complete_energy,
                                temperature,
                            )
                        {
                            gates.swap(position, neighbor);
                            quality = swap_quality;
                            if let (Some(model), Some(work)) = (&self.depth, &mut depth_work) {
                                work.advance(model, b);
                            }
                            accepted += 1;
                            std::mem::swap(&mut env, &mut swap_env);
                            if swap_energy + 1e-12 < best_energy {
                                best_energy = swap_energy;
                                best_eq = swap_eq;
                                best_gates.clone_from(&gates);
                                last_improvement = sweep;
                            }
                            if swap_eq <= options.epsilon * std::f64::consts::SQRT_2 + 1e-7 {
                                let candidate = Circuit::new(gates.clone(), lib.n_qubits, lib);
                                if let Some(circuit) =
                                    correct(candidate, target, lib, options.epsilon)?
                                {
                                    let eq =
                                        target.original.cost(circuit.matrix(), options.simple_cost);
                                    return Ok(SearchResult {
                                        circuit,
                                        best_eq: eq,
                                        best_energy: swap_energy,
                                        found: true,
                                        steps,
                                        accepted,
                                    });
                                }
                            }
                            // The swap environment already belongs to the next
                            // coordinate: the two gate costs are unchanged.
                            continue;
                        }
                    }
                    if let (Some(model), Some(work)) = (&self.depth, &mut depth_work) {
                        work.advance(model, new);
                    }
                    if forward {
                        let next = gates[position + 1];
                        for e in &mut env {
                            self.sparse[new].left(e, &mut scratch);
                            std::mem::swap(e, &mut scratch);
                            self.adjoints[next].right(e, &mut scratch);
                            std::mem::swap(e, &mut scratch);
                        }
                    } else {
                        let prev = gates[position - 1];
                        for e in &mut env {
                            self.adjoints[prev].left(e, &mut scratch);
                            std::mem::swap(e, &mut scratch);
                            self.sparse[new].right(e, &mut scratch);
                            std::mem::swap(e, &mut scratch);
                        }
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

fn correct(
    mut circuit: Circuit,
    target: &Target,
    lib: &GateLibrary,
    epsilon: f64,
) -> Result<Option<Circuit>, String> {
    for (i, p) in target.variants.iter().enumerate() {
        if !target.inverses[i] && p.exact_cost(circuit.matrix(), epsilon) == 0.0 {
            circuit.permute(&target.permutations[i], lib)?;
            return Ok(
                (target.original.exact_cost(circuit.matrix(), epsilon) == 0.0).then_some(circuit),
            );
        }
    }
    Ok(None)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn sparse_products_and_trace_match_dense() {
        for nq in 1..=4 {
            let lib = GateLibrary::load(nq, "data/gates/CliffordT", "").unwrap();
            let n = 1 << nq;
            let m = Matrix {
                n,
                data: (0..n * n)
                    .map(|i| Complex::new((i % 13) as f64 / 13.0, (i % 7) as f64 / 7.0))
                    .collect(),
            };
            for gate in &lib.gates {
                let sparse = SparseGate::new(&gate.matrix);
                let mut out = Matrix::zero(n);
                sparse.left(&m, &mut out);
                assert!(out.max_abs_diff(&gate.matrix.mul(&m)) < 1e-12);
                sparse.right(&m, &mut out);
                assert!(out.max_abs_diff(&m.mul(&gate.matrix)) < 1e-12);
                let prod = m.mul(&gate.matrix);
                let trace = (0..n).fold(Complex::ZERO, |a, i| a + prod[(i, i)]).abs();
                assert!((sparse.trace(&m, identity_trace(&m)) - trace).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn sparse_products_initialize_zero_rows_and_columns_and_match_dense_complex_products() {
        let mut rng = Rng::new(55129);
        for n in [2, 4, 8, 16] {
            for _ in 0..20 {
                let mut data = Matrix::zero(n);
                let mut gate = Matrix::zero(n);
                for r in 0..n {
                    for c in 0..n {
                        data[(r, c)] = Complex::new(rng.random01(), rng.random01());
                        if r > 0 && c > 0 && rng.random01() < 0.6 {
                            gate[(r, c)] = Complex::new(rng.random01(), rng.random01());
                        }
                    }
                }
                let sparse = SparseGate::new(&gate);
                let mut out = data.clone();
                sparse.left(&data, &mut out);
                assert!(out.max_abs_diff(&gate.mul(&data)) < 1e-12);
                out = data.clone();
                sparse.right(&data, &mut out);
                assert!(out.max_abs_diff(&data.mul(&gate)) < 1e-12);
            }
        }
    }

    #[test]
    fn shared_phase_traces_match_uncached_and_dense_products() {
        let mut rng = Rng::new(190447);
        for q in 1..=5 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let kernel = SweepKernel::new(&lib, &SearchOptions::tuned());
            if q > 1 {
                assert!(kernel.trace_reuse.iter().filter(|x| x.is_some()).count() >= 4 * q);
            }
            let n = 1 << q;
            for _ in 0..8 {
                let env = Matrix {
                    n,
                    data: (0..n * n)
                        .map(|_| Complex::new(rng.random01(), rng.random01()))
                        .collect(),
                };
                let mut values = vec![Complex::ZERO; kernel.trace_groups.len()];
                kernel.fill_trace_groups(&env, &mut values);
                let identity = identity_trace(&env);
                for &id in &kernel.choices {
                    let actual = kernel.reused_trace(id, &env, identity, &values);
                    let uncached = kernel.sparse[id].trace(&env, identity);
                    let dense = identity_trace(&env.mul(&lib.gates[id].matrix)).abs();
                    assert!((actual - uncached).abs() < 1e-12);
                    assert!((actual - dense).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn detects_only_constant_norm_covers() {
        let full = PartialMatrix::new(Matrix::identity(4), vec![true; 16], "full").unwrap();
        assert_eq!(constant_norm(&full), Some(4.0));
        let iso = PartialMatrix::new(
            Matrix::identity(4),
            (0..16).map(|i| i % 4 == 0).collect(),
            "iso",
        )
        .unwrap();
        assert_eq!(constant_norm(&iso), Some(1.0));
        let arbitrary = PartialMatrix::new(
            Matrix::identity(4),
            (0..16).map(|i| i == 0).collect(),
            "mask",
        )
        .unwrap();
        assert_eq!(constant_norm(&arbitrary), None);
    }
    #[test]
    fn contraction_environments_match_complete_circuits() {
        for nq in 1..=4 {
            let lib = GateLibrary::load(nq, "data/gates/CliffordT", "").unwrap();
            let options = SearchOptions::default();
            let kernel = SweepKernel::new(&lib, &options);
            let mut rng = Rng::new(123);
            let mut c = random_circuit(&lib, &mut rng, 17, 0.3);
            let target = random_circuit(&lib, &mut rng, 13, 0.0).matrix().clone();
            let mut env = target.adjoint();
            let mut scratch = Matrix::zero(env.n);
            for &g in c.gates[1..].iter().rev() {
                kernel.sparse[g].right(&env, &mut scratch);
                std::mem::swap(&mut env, &mut scratch);
            }
            for sweep in 0..10 {
                let forward = sweep % 2 == 0;
                for offset in 0..c.gates.len() {
                    let pos = if forward {
                        offset
                    } else {
                        c.gates.len() - 1 - offset
                    };
                    let id = lib.all[rng.usize(lib.all.len())];
                    c.replace(pos, id, &lib);
                    let estimate = kernel.sparse[id].trace(&env, identity_trace(&env));
                    let exact = target.trace_conjugate_product(c.matrix()).abs();
                    assert!(
                        (estimate - exact).abs() < 1e-10,
                        "nq={nq}, sweep={sweep}, pos={pos}"
                    );
                    if offset + 1 < c.gates.len() {
                        if forward {
                            kernel.sparse[id].left(&env, &mut scratch);
                            std::mem::swap(&mut env, &mut scratch);
                            kernel.adjoints[c.gates[pos + 1]].right(&env, &mut scratch);
                            std::mem::swap(&mut env, &mut scratch);
                        } else {
                            kernel.adjoints[c.gates[pos - 1]].left(&env, &mut scratch);
                            std::mem::swap(&mut env, &mut scratch);
                            kernel.sparse[id].right(&env, &mut scratch);
                            std::mem::swap(&mut env, &mut scratch);
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn sweep_solutions_obey_full_and_partial_specs() {
        let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
        let options = SearchOptions {
            cost_power: 1.0,
            gate_prior: 1.0,
            sweep_permutations: false,
            ..SearchOptions::default()
        };
        let kernel = SweepKernel::new(&lib, &options);
        let bell = Circuit::new(
            vec![
                lib.find("h", &[0]).unwrap(),
                lib.find("cx", &[0, 1]).unwrap(),
            ],
            2,
            &lib,
        );
        for partial in [false, true] {
            for seed in 0..8 {
                let target = Target::new(
                    PartialMatrix::new(
                        bell.matrix().clone(),
                        (0..16).map(|i| !partial || i % 4 == 0).collect(),
                        "bell",
                    )
                    .unwrap(),
                    true,
                    false,
                );
                let mut rng = Rng::new(seed);
                let mut found = false;
                for _ in 0..20 {
                    let result = kernel.run(&target, &lib, &mut rng, 8, &options).unwrap();
                    if result.found {
                        assert_eq!(
                            target.original.exact_cost(result.circuit.matrix(), 1e-6),
                            0.0
                        );
                        found = true;
                        break;
                    }
                }
                assert!(found, "seed={seed}, partial={partial}");
            }
        }
    }

    #[test]
    fn contraction_objective_matches_full_evaluation_for_rows_columns_and_full_masks() {
        for nq in 1..=4 {
            let lib = GateLibrary::load(nq, "data/gates/CliffordT", "").unwrap();
            let mut rng = Rng::new(918);
            let n = 1 << nq;
            for _ in 0..8 {
                let target_matrix = random_circuit(&lib, &mut rng, 11, 0.0).matrix().clone();
                let prefix = random_circuit(&lib, &mut rng, 7, 0.0).matrix().clone();
                let suffix = random_circuit(&lib, &mut rng, 9, 0.0).matrix().clone();
                for kind in 0..3 {
                    let mask = (0..n * n)
                        .map(|i| match kind {
                            0 => true,
                            1 => i % n < n / 2,
                            _ => i / n < n / 2,
                        })
                        .collect();
                    let p = PartialMatrix::new(target_matrix.clone(), mask, "fixture").unwrap();
                    let env = prefix.mul(&p.matrix.adjoint()).mul(&suffix);
                    for gate in &lib.gates {
                        let full = suffix.mul(&gate.matrix).mul(&prefix);
                        let overlap =
                            SparseGate::new(&gate.matrix).trace(&env, identity_trace(&env));
                        for simple in [false, true] {
                            let fast =
                                equality_squared(&p, constant_norm(&p).unwrap(), overlap, simple);
                            let reference = p.cost(&full, simple).powi(2);
                            assert!((fast - reference).abs() < 1e-11);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn raw_t_budget_is_respected_even_on_unsolved_chains() {
        let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
        let target_circuit = Circuit::new(
            vec![
                lib.find("h", &[0]).unwrap(),
                lib.find("cx", &[0, 1]).unwrap(),
            ],
            2,
            &lib,
        );
        let target = Target::new(
            PartialMatrix::new(target_circuit.matrix().clone(), vec![true; 16], "bell").unwrap(),
            false,
            false,
        );
        for cap in [0, 1, 2] {
            let options = SearchOptions {
                quality_metric: QualityMetric::TCount,
                quality_goal: Some(cap as f64),
                quality_slack: Some(0.0),
                quality_weight: 0.003,
                ..SearchOptions::tuned()
            };
            let kernel = SweepKernel::new(&lib, &options);
            for seed in 0..32 {
                let result = kernel
                    .run(&target, &lib, &mut Rng::new(seed), 13, &options)
                    .unwrap();
                assert!(result.circuit.count(&["t".into(), "tdg".into()], &lib) <= cap);
                assert!(result.best_energy.is_finite());
                if result.found {
                    assert_eq!(
                        target.original.exact_cost(result.circuit.matrix(), 1e-6),
                        0.0
                    );
                }
            }
        }
    }

    #[test]
    fn swapped_environments_match_independent_full_products_in_both_directions() {
        for nq in 1..=4 {
            let mut lib = GateLibrary::load(nq, "data/gates/CliffordT", "").unwrap();
            lib.add_clifford_macros();
            let kernel = SweepKernel::new(&lib, &SearchOptions::tuned());
            let mut rng = Rng::new(913);
            for _ in 0..16 {
                let target = random_circuit(&lib, &mut rng, 7, 0.3).matrix().clone();
                for forward in [false, true] {
                    let mut c = random_circuit(&lib, &mut rng, 9, 0.3);
                    let position = if forward {
                        rng.usize(8)
                    } else {
                        1 + rng.usize(8)
                    };
                    let neighbor = if forward { position + 1 } else { position - 1 };
                    let prefix = Circuit::new(c.gates[..position].to_vec(), nq, &lib);
                    let suffix = Circuit::new(c.gates[position + 1..].to_vec(), nq, &lib);
                    let env = prefix.matrix().mul(&target.adjoint()).mul(suffix.matrix());
                    let mut out = Matrix::zero(env.n);
                    let mut scratch = Matrix::zero(env.n);
                    kernel.swap_environment(
                        &env,
                        &mut out,
                        &mut scratch,
                        c.gates[neighbor],
                        forward,
                    );
                    let a = c.gates[position];
                    c.gates.swap(position, neighbor);
                    c.rebuild(&lib);
                    let estimate = kernel.sparse[a].trace(&out, identity_trace(&out));
                    let exact = target.trace_conjugate_product(c.matrix()).abs();
                    assert!((estimate - exact).abs() < 1e-10);
                    let prefix = Circuit::new(c.gates[..neighbor].to_vec(), nq, &lib);
                    let suffix = Circuit::new(c.gates[neighbor + 1..].to_vec(), nq, &lib);
                    let exact_env = prefix.matrix().mul(&target.adjoint()).mul(suffix.matrix());
                    assert!(out.max_abs_diff(&exact_env) < 1e-10);
                }
            }
        }
    }

    #[test]
    fn depth_capped_annealing_and_swaps_report_exact_depth_energy() {
        let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
        let c = Circuit::new(
            vec![
                lib.find("h", &[0]).unwrap(),
                lib.find("cx", &[0, 1]).unwrap(),
            ],
            2,
            &lib,
        );
        let target = Target::new(
            PartialMatrix::new(c.matrix().clone(), vec![true; 16], "bell").unwrap(),
            false,
            false,
        );
        for cap in [0, 1, 2] {
            let options = SearchOptions {
                quality_metric: QualityMetric::TDepth,
                quality_goal: Some(cap as f64),
                quality_slack: Some(0.0),
                quality_weight: 0.03,
                swap_probability: 0.5,
                ..SearchOptions::tuned()
            };
            let kernel = SweepKernel::new(&lib, &options);
            for seed in 0..16 {
                let result = kernel
                    .run(&target, &lib, &mut Rng::new(seed), 13, &options)
                    .unwrap();
                let depth = result.circuit.depth(&["t".into(), "tdg".into()], &lib);
                assert!(depth <= cap);
                let energy = target.original.cost(result.circuit.matrix(), false).powi(2)
                    + options.quality_weight * depth as f64;
                assert!((energy - result.best_energy).abs() < 1e-8);
                if result.found {
                    assert_eq!(
                        target.original.exact_cost(result.circuit.matrix(), 1e-6),
                        0.0
                    );
                }
            }
        }
    }

    #[test]
    fn rejects_nonunitary_gate_libraries() {
        let mut lib = GateLibrary::load(1, "data/gates/CliffordT", "").unwrap();
        let target = Target::new(
            PartialMatrix::new(Matrix::identity(2), vec![true; 4], "id").unwrap(),
            false,
            false,
        );
        assert!(SweepKernel::supports(&target, &lib));
        lib.gates[0].matrix[(0, 0)] = Complex::new(0.9, 0.0);
        assert!(!SweepKernel::supports(&target, &lib));
    }
}
