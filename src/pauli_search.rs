//! Annealed search over bounded Pauli-T words and an exact Clifford remainder.
//! The Clifford-distance objective removes free Clifford choices from the search.
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
use std::time::Instant;

/// (a + b sqrt(2)) / sqrt(2)^exponent, reduced in Z[sqrt(2)].
/// This guides proposals; the original matrix still decides acceptance.
#[derive(Clone, Copy, Debug)]
struct RootEntry {
    a: i64,
    b: i64,
    exponent: usize,
}
impl RootEntry {
    fn reduced(mut a: i64, mut b: i64, mut exponent: usize) -> Self {
        if a == 0 && b == 0 {
            exponent = 0;
        }
        while exponent > 0 && a % 2 == 0 {
            (a, b) = (b, a / 2);
            exponent -= 1;
        }
        Self { a, b, exponent }
    }
    fn raised(self, exponent: usize) -> (i64, i64) {
        let shift = exponent - self.exponent;
        let (a, b) = if shift % 2 == 1 {
            (2 * self.b, self.a)
        } else {
            (self.a, self.b)
        };
        (a << (shift / 2), b << (shift / 2))
    }
    fn rotated(self, other: Self, sign: i64) -> (Self, Self) {
        let exponent = self.exponent.max(other.exponent);
        let (a, b) = self.raised(exponent);
        let (c, d) = other.raised(exponent);
        (
            Self::reduced(a - sign * c, b - sign * d, exponent + 1),
            Self::reduced(sign * a + c, sign * b + d, exponent + 1),
        )
    }
    #[cfg(test)]
    fn value(self) -> f64 {
        (self.a as f64 + self.b as f64 * std::f64::consts::SQRT_2)
            / 2.0f64.powf(self.exponent as f64 / 2.0)
    }
}
#[derive(Clone)]
struct RingTransfer {
    data: Vec<RootEntry>,
    histogram: [usize; 64],
    support: Vec<Vec<usize>>,
    normalization: f64,
}
impl RingTransfer {
    fn decode(data: &[f64], limit: usize) -> Option<Self> {
        // Bound startup and intermediate integer coefficients. Other targets
        // retain the fourth-moment objective.
        if limit > 16 {
            return None;
        }
        let scale = 2.0f64.powf(limit as f64 / 2.0);
        let mut values = Vec::new();
        let b_bound = (scale / std::f64::consts::SQRT_2 + 1e-10).floor() as i64;
        for b in -b_bound..=b_bound {
            let bound =
                (scale - (b as f64 * std::f64::consts::SQRT_2).abs() + 1e-10).floor() as i64;
            for a in -bound..=bound {
                let value = (a as f64 + b as f64 * std::f64::consts::SQRT_2) / scale;
                values.push((value, RootEntry::reduced(a, b, limit)));
            }
        }
        values.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
        let decoded: Option<Vec<_>> = data
            .iter()
            .map(|&value| {
                let i = values.partition_point(|entry| entry.0 < value);
                let closest = [i.saturating_sub(1), i.min(values.len() - 1)]
                    .into_iter()
                    .min_by(|&a, &b| {
                        (values[a].0 - value)
                            .abs()
                            .total_cmp(&(values[b].0 - value).abs())
                    })?;
                ((values[closest].0 - value).abs() < 1e-9).then_some(values[closest].1)
            })
            .collect();
        let data = decoded?;
        let mut histogram = [0; 64];
        for entry in &data {
            histogram[entry.exponent] += 1;
        }
        let d = (data.len() as f64).sqrt() as usize;
        let support = data
            .chunks(d)
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter_map(|(c, x)| (x.a != 0 || x.b != 0).then_some(c))
                    .collect()
            })
            .collect();
        Some(Self {
            data,
            histogram,
            support,
            normalization: (limit + 1) as f64,
        })
    }
    fn energy_of(&self, histogram: &[usize; 64]) -> f64 {
        let max = histogram.iter().rposition(|&count| count != 0).unwrap_or(0);
        let mass: f64 = histogram[..=max]
            .iter()
            .enumerate()
            .map(|(k, &count)| count as f64 * 2.0f64.powi(k as i32 - max as i32))
            .sum();
        (max as f64 + 0.25 * mass / self.data.len() as f64) / self.normalization
    }
    fn rotated_energy(&self, d: usize, pairs: &[(usize, usize, f64)]) -> f64 {
        let mut histogram = self.histogram;
        for &(a, b, _) in pairs {
            let left = &self.support[a];
            let right = &self.support[b];
            let (mut i, mut j) = (0, 0);
            while i < left.len() || j < right.len() {
                let xcol = left.get(i).copied().unwrap_or(usize::MAX);
                let ycol = right.get(j).copied().unwrap_or(usize::MAX);
                let column = xcol.min(ycol);
                i += usize::from(xcol == column);
                j += usize::from(ycol == column);
                let x = self.data[a * d + column];
                let y = self.data[b * d + column];
                let (first, second) = x.rotated(y, 1);
                histogram[x.exponent] -= 1;
                histogram[y.exponent] -= 1;
                histogram[first.exponent] += 1;
                histogram[second.exponent] += 1;
            }
        }
        self.energy_of(&histogram)
    }
    fn rotate_into(&self, d: usize, pairs: &[(usize, usize, f64)], direction: f64, out: &mut Self) {
        out.data.copy_from_slice(&self.data);
        out.histogram = self.histogram;
        out.support.clone_from(&self.support);
        for &(a, b, sign) in pairs {
            out.support[a].clear();
            out.support[b].clear();
            for c in 0..d {
                let x = self.data[a * d + c];
                let y = self.data[b * d + c];
                let (first, second) = x.rotated(y, (sign * direction) as i64);
                out.data[a * d + c] = first;
                out.data[b * d + c] = second;
                if first.a != 0 || first.b != 0 {
                    out.support[a].push(c);
                }
                if second.a != 0 || second.b != 0 {
                    out.support[b].push(c);
                }
                out.histogram[x.exponent] -= 1;
                out.histogram[y.exponent] -= 1;
                out.histogram[first.exponent] += 1;
                out.histogram[second.exponent] += 1;
            }
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Pauli {
    x: usize,
    z: usize,
    negative: bool,
}
impl Pauli {
    fn axis(bits: usize, q: usize) -> Self {
        Self {
            x: bits & ((1 << q) - 1),
            z: bits >> q,
            negative: false,
        }
    }
    fn phase(self, column: usize) -> Complex {
        let value = match (self.x & self.z).count_ones() % 4 {
            0 => Complex::ONE,
            1 => Complex::I,
            2 => -Complex::ONE,
            _ => -Complex::I,
        };
        if ((self.z & column).count_ones() % 2 == 1) ^ self.negative {
            -value
        } else {
            value
        }
    }
    fn h(&mut self, bit: usize) {
        let x = (self.x >> bit) & 1;
        let z = (self.z >> bit) & 1;
        self.negative ^= x & z != 0;
        self.x ^= (x ^ z) << bit;
        self.z ^= (x ^ z) << bit;
    }
    fn s(&mut self, bit: usize, inverse: bool) {
        let x = (self.x >> bit) & 1;
        let z = (self.z >> bit) & 1;
        self.negative ^= x & (z ^ usize::from(inverse)) != 0;
        self.z ^= x << bit;
    }
    fn cx(&mut self, control: usize, target: usize) {
        let xc = (self.x >> control) & 1;
        let xt = (self.x >> target) & 1;
        let zc = (self.z >> control) & 1;
        let zt = (self.z >> target) & 1;
        self.negative ^= xc & zt & (xt ^ zc ^ 1) != 0;
        self.x ^= xc << target;
        self.z ^= zt << control;
    }
}

#[derive(Clone)]
struct Transfer {
    data: Vec<f64>,
    fourth: f64,
    ring: Option<RingTransfer>,
}
impl Transfer {
    fn from_matrix(matrix: &Matrix, q: usize) -> Self {
        let n = matrix.n;
        let d = n * n;
        let basis: Vec<_> = (0..d)
            .map(|axis| {
                let p = Pauli::axis(axis, q);
                (0..n)
                    .map(|c| ((c ^ p.x) * n + c, p.phase(c)))
                    .collect::<Vec<_>>()
            })
            .collect();
        let adjoint = matrix.adjoint();
        let mut shifted = Matrix::zero(n);
        let mut conjugated = Matrix::zero(n);
        let mut data = vec![0.0; d * d];
        for (column, input) in basis.iter().enumerate() {
            for &(index, phase) in input {
                let (r, c) = (index / n, index % n);
                for k in 0..n {
                    shifted[(r, k)] = phase * adjoint[(c, k)];
                }
            }
            matrix.mul_into(&shifted, &mut conjugated);
            for (row, output) in basis.iter().enumerate() {
                let value = output.iter().fold(Complex::ZERO, |sum, &(index, phase)| {
                    sum + phase.conj() * conjugated.data[index]
                });
                data[row * d + column] = value.re / n as f64;
            }
        }
        let fourth = data.iter().map(|x| (x * x) * (x * x)).sum();
        Self {
            data,
            fourth,
            ring: None,
        }
    }
    fn cost(&self, d: usize) -> f64 {
        (1.0 - (self.fourth - 1.0) / (d - 1) as f64).max(0.0)
    }
    fn rotated_fourth(&self, d: usize, pairs: &[(usize, usize, f64)]) -> f64 {
        let mut change = 0.0;
        for &(a, b, _) in pairs {
            for (&x, &y) in self.data[a * d..(a + 1) * d]
                .iter()
                .zip(&self.data[b * d..(b + 1) * d])
            {
                let x2 = x * x;
                let y2 = y * y;
                change += 3.0 * x2 * y2 - 0.5 * (x2 * x2 + y2 * y2);
            }
        }
        self.fourth + change
    }
    fn rotate_into(&self, d: usize, pairs: &[(usize, usize, f64)], direction: f64, out: &mut Self) {
        out.data.copy_from_slice(&self.data);
        let h = std::f64::consts::FRAC_1_SQRT_2;
        for &(a, b, sign) in pairs {
            let sign = sign * direction;
            for c in 0..d {
                let x = self.data[a * d + c];
                let y = self.data[b * d + c];
                out.data[a * d + c] = h * (x - sign * y);
                out.data[b * d + c] = h * (sign * x + y);
            }
        }
        out.fourth = self.rotated_fourth(d, pairs);
        if let (Some(ring), Some(next)) = (&self.ring, &mut out.ring) {
            ring.rotate_into(d, pairs, direction, next);
        }
    }
    fn clifford_images(&self, q: usize) -> Option<Vec<Pauli>> {
        let d = 1 << (2 * q);
        (0..2 * q)
            .map(|bit| {
                let column = 1 << bit;
                let row = (1..d).find(|&r| self.data[r * d + column].abs() > 1.0 - 1e-8)?;
                let mut pauli = Pauli::axis(row, q);
                pauli.negative = self.data[row * d + column] < 0.0;
                Some(pauli)
            })
            .collect()
    }
}

fn rotation_pairs(q: usize) -> Vec<Vec<(usize, usize, f64)>> {
    let d = 1usize << (2 * q);
    (0..d)
        .map(|axis| {
            let a = Pauli::axis(axis, q);
            (0..d)
                .filter_map(|bits| {
                    let other = bits ^ axis;
                    let b = Pauli::axis(bits, q);
                    if bits >= other
                        || ((a.x & b.z).count_ones() + (a.z & b.x).count_ones()).is_multiple_of(2)
                    {
                        return None;
                    }
                    let phase = ((a.x & a.z).count_ones() as i32 + (b.x & b.z).count_ones() as i32
                        - ((a.x ^ b.x) & (a.z ^ b.z)).count_ones() as i32
                        + 2 * (a.z & b.x).count_ones() as i32)
                        .rem_euclid(4);
                    Some((bits, other, if phase == 1 { 1.0 } else { -1.0 }))
                })
                .collect()
        })
        .collect()
}

fn clifford_word(
    mut images: Vec<Pauli>,
    lib: &GateLibrary,
    native: &NativeClifford,
) -> Option<Vec<usize>> {
    let q = lib.n_qubits;
    let mut elimination = Vec::new();
    let apply = |id: usize, images: &mut [Pauli], word: &mut Vec<usize>| {
        let gate = &lib.gates[id];
        for p in images {
            match gate.name.as_str() {
                "h" => p.h(gate.qubits[0]),
                "s" => p.s(gate.qubits[0], false),
                "sdg" => p.s(gate.qubits[0], true),
                "cx" => p.cx(gate.qubits[0], gate.qubits[1]),
                _ => unreachable!(),
            }
        }
        word.push(id);
    };
    for bit in 0..q {
        let pivot = (bit..q).find(|&j| ((images[bit].x | images[bit].z) >> j) & 1 != 0)?;
        if pivot != bit {
            for (a, b) in [(bit, pivot), (pivot, bit), (bit, pivot)] {
                apply(native.cx[a][b]?, &mut images, &mut elimination);
            }
        }
        for j in bit..q {
            let x = (images[bit].x >> j) & 1;
            let z = (images[bit].z >> j) & 1;
            if x == 0 && z == 1 {
                apply(native.native[j][0], &mut images, &mut elimination);
            } else if x == 1 && z == 1 {
                apply(native.native[j][2], &mut images, &mut elimination);
            }
            if j != bit && (images[bit].x >> j) & 1 != 0 {
                apply(native.cx[bit][j]?, &mut images, &mut elimination);
            }
        }
        if images[bit].negative {
            for _ in 0..2 {
                apply(native.native[bit][1], &mut images, &mut elimination);
            }
        }
        for j in bit + 1..q {
            let x = (images[q + bit].x >> j) & 1;
            let z = (images[q + bit].z >> j) & 1;
            if x == 1 {
                if z == 1 {
                    apply(native.native[j][2], &mut images, &mut elimination);
                }
                apply(native.native[j][0], &mut images, &mut elimination);
            }
            if (images[q + bit].z >> j) & 1 != 0 {
                apply(native.cx[j][bit]?, &mut images, &mut elimination);
            }
        }
        if (images[q + bit].x >> bit) & 1 != 0 {
            for k in [0, 1, 0] {
                apply(native.native[bit][k], &mut images, &mut elimination);
            }
        }
        if images[q + bit].negative {
            for k in [0, 1, 1, 0] {
                apply(native.native[bit][k], &mut images, &mut elimination);
            }
        }
    }
    if images.iter().enumerate().any(|(i, p)| {
        p.negative
            || p.x != if i < q { 1 << i } else { 0 }
            || p.z != if i >= q { 1 << (i - q) } else { 0 }
    }) {
        return None;
    }
    crate::affine::inverse_word(&elimination, lib)
}

fn rotation_word(axis: usize, lib: &GateLibrary, native: &NativeClifford) -> Vec<usize> {
    let p = Pauli::axis(axis, lib.n_qubits);
    let pivot = (p.x | p.z).trailing_zeros() as usize;
    let mut basis = Vec::new();
    for bit in 0..lib.n_qubits {
        if (p.x >> bit) & 1 != 0 {
            if (p.z >> bit) & 1 != 0 {
                basis.push(native.native[bit][2]);
            }
            basis.push(native.native[bit][0]);
        }
    }
    for bit in 0..lib.n_qubits {
        if bit != pivot && ((p.x | p.z) >> bit) & 1 != 0 {
            basis.push(native.cx[bit][pivot].unwrap());
        }
    }
    let inverse = crate::affine::inverse_word(&basis, lib).unwrap();
    basis.push(native.native[pivot][3]);
    basis.extend(inverse);
    basis
}

pub struct PauliKernel {
    q: usize,
    d: usize,
    limit: usize,
    initial: Transfer,
    pairs: Vec<Vec<(usize, usize, f64)>>,
    words: Vec<Vec<usize>>,
    native: NativeClifford,
}
impl PauliKernel {
    pub fn new(target: &Target, lib: &GateLibrary, limit: usize) -> Option<Self> {
        let q = lib.n_qubits;
        if q > 4
            || limit > 24
            || target.original.cover.iter().any(|&v| !v)
            || target.inverses.iter().any(|&v| v)
        {
            return None;
        }
        let native = NativeClifford::new(lib)?;
        let matrix = &target.original.matrix;
        if !matrix
            .mul(&matrix.adjoint())
            .approximately_equal(&Matrix::identity(matrix.n), 1e-9)
        {
            return None;
        }
        let d = 1 << (2 * q);
        let mut words = vec![Vec::new()];
        words.extend((1..d).map(|axis| rotation_word(axis, lib, &native)));
        let mut initial = Transfer::from_matrix(matrix, q);
        initial.ring = RingTransfer::decode(&initial.data, limit);
        Some(Self {
            q,
            d,
            limit,
            initial,
            pairs: rotation_pairs(q),
            words,
            native,
        })
    }
    fn candidate(
        &self,
        state: &Transfer,
        axes: &[usize],
        target: &Target,
        lib: &GateLibrary,
        options: &SearchOptions,
    ) -> Option<Circuit> {
        let mut word = clifford_word(state.clifford_images(self.q)?, lib, &self.native)?;
        for &axis in axes.iter().rev() {
            word.extend(&self.words[axis]);
        }
        let circuit = Circuit::new(word, self.q, lib);
        (target
            .original
            .exact_cost(circuit.matrix(), options.epsilon)
            == 0.0)
            .then_some(circuit)
    }
    pub fn run_until(
        &self,
        target: &Target,
        lib: &GateLibrary,
        rng: &mut Rng,
        options: &SearchOptions,
        deadline: Option<Instant>,
    ) -> SearchResult {
        self.run_impl::<true>(target, lib, rng, options, deadline)
    }
    fn run_impl<const CACHE: bool>(
        &self,
        target: &Target,
        lib: &GateLibrary,
        rng: &mut Rng,
        options: &SearchOptions,
        deadline: Option<Instant>,
    ) -> SearchResult {
        let mut states = vec![self.initial.clone(); self.limit + 1];
        let mut scratch = self.initial.clone();
        let mut axes = Vec::with_capacity(self.limit);
        // A prefix can be revisited many times, especially at the length bound.
        // Proposal energies depend only on that state and objective, not temperature.
        let mut fourth_energies: Vec<Option<Vec<f64>>> = vec![None; self.limit + 1];
        let mut ring_energies: Vec<Option<Vec<f64>>> = vec![None; self.limit + 1];
        let mut weights = vec![0.0; self.d];
        let iterations = (options.iterations_factor * (self.limit + 1) as f64).ceil() as usize;
        let mut steps = 0;
        let mut accepted = 0;
        let warmup = 2 * (self.limit + 1);
        for iteration in 0..iterations {
            // Try the inexpensive fourth moment first. If it stalls, restart
            // from the target and anneal the exact denominator structure.
            if self.initial.ring.is_some() && iteration == warmup {
                axes.clear();
            }
            let ring = (iteration >= warmup)
                .then_some(())
                .and(self.initial.ring.as_ref());
            if options.is_cancelled() || deadline.is_some_and(|end| Instant::now() >= end) {
                break;
            }
            if states[axes.len()].cost(self.d) < 1e-8 {
                if let Some(circuit) =
                    self.candidate(&states[axes.len()], &axes, target, lib, options)
                {
                    let eq = target.original.cost(circuit.matrix(), options.simple_cost);
                    return SearchResult {
                        circuit,
                        best_eq: eq,
                        best_energy: eq.powf(2.0 * options.cost_power),
                        found: true,
                        steps,
                        accepted,
                    };
                }
            }
            if self.limit == 0 {
                break;
            }
            // Append, replace the last rotation, or revisit an earlier prefix.
            // Heat-bath sampling anneals the same Clifford-distance objective.
            let prefix = if !axes.is_empty() && rng.random01() < 0.12 {
                rng.usize(axes.len())
            } else if axes.len() == self.limit {
                axes.len() - 1
            } else {
                axes.len()
            };
            let state = &states[prefix];
            let ring = ring.and(state.ring.as_ref());
            // A prefix already represents a variable-length word. An identity
            // proposal can trap it forever when every useful next T rotation
            // temporarily raises the denominator; prefix rewinds supply deletion.
            let cache = if ring.is_some() {
                &mut ring_energies
            } else {
                &mut fourth_energies
            };
            if !CACHE {
                cache[prefix] = None;
            }
            let energies = cache[prefix].get_or_insert_with(|| {
                let mut values = vec![f64::INFINITY; self.d];
                for (axis, energy) in values.iter_mut().enumerate().skip(1) {
                    *energy = ring.map_or_else(
                        || {
                            let fourth = state.rotated_fourth(self.d, &self.pairs[axis]);
                            (1.0 - (fourth - 1.0) / (self.d - 1) as f64).max(0.0)
                        },
                        |ring| ring.rotated_energy(self.d, &self.pairs[axis]),
                    );
                }
                values
            });
            let min = energies.iter().copied().fold(f64::INFINITY, f64::min);
            let temperature = (options.start_temp_base
                * 0.01
                * options.sweep_cooling.powf(iteration as f64 / 8.0))
            .max(1e-6);
            let mut sum = 0.0;
            for axis in 0..self.d {
                let prior = if axis == 0 {
                    options.pid
                } else {
                    (1.0 - options.pid) / (self.d - 1) as f64
                };
                weights[axis] = if axis == 0 {
                    0.0
                } else {
                    (-(energies[axis] - min) / temperature).exp()
                        * prior.max(1e-12).powf(options.gate_prior)
                };
                sum += weights[axis];
            }
            let mut draw = rng.random01() * sum;
            let mut chosen = self.d - 1;
            for (axis, &weight) in weights.iter().enumerate() {
                draw -= weight;
                if draw <= 0.0 {
                    chosen = axis;
                    break;
                }
            }
            axes.truncate(prefix);
            if chosen != 0 {
                states[prefix].rotate_into(self.d, &self.pairs[chosen], -1.0, &mut scratch);
                std::mem::swap(&mut states[prefix + 1], &mut scratch);
                fourth_energies[prefix + 1] = None;
                ring_energies[prefix + 1] = None;
                axes.push(chosen);
            }
            steps += 1;
            accepted += 1;
        }
        if states[axes.len()].cost(self.d) < 1e-8 {
            if let Some(circuit) = self.candidate(&states[axes.len()], &axes, target, lib, options)
            {
                let eq = target.original.cost(circuit.matrix(), options.simple_cost);
                return SearchResult {
                    circuit,
                    best_eq: eq,
                    best_energy: eq.powf(2.0 * options.cost_power),
                    found: true,
                    steps,
                    accepted,
                };
            }
        }
        let circuit = Circuit::new(Vec::new(), self.q, lib);
        let eq = target.original.cost(circuit.matrix(), options.simple_cost);
        SearchResult {
            circuit,
            best_eq: eq,
            best_energy: eq.powf(2.0 * options.cost_power),
            found: false,
            steps,
            accepted,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partial::PartialMatrix;
    fn same(a: &Matrix, b: &Matrix) -> bool {
        let z = a.trace_conjugate_product(b);
        let phase = z / z.abs();
        a.data
            .iter()
            .zip(&b.data)
            .all(|(&x, &y)| (y - phase * x).norm_sqr() < 1e-18)
    }
    #[test]
    fn exact_ring_arithmetic_matches_independent_real_values() {
        for exponent in 0..=8 {
            for a in -12..=12 {
                for b in -12..=12 {
                    let entry = RootEntry::reduced(a, b, exponent);
                    let expected = (a as f64 + b as f64 * std::f64::consts::SQRT_2)
                        / 2.0f64.powf(exponent as f64 / 2.0);
                    assert!((entry.value() - expected).abs() < 1e-12);
                    let other = RootEntry::reduced(b - a, a + b, 8 - exponent);
                    for sign in [-1, 1] {
                        let (first, second) = entry.rotated(other, sign);
                        let h = std::f64::consts::FRAC_1_SQRT_2;
                        assert!(
                            (first.value() - h * (entry.value() - sign as f64 * other.value()))
                                .abs()
                                < 1e-11
                        );
                        assert!(
                            (second.value() - h * (sign as f64 * entry.value() + other.value()))
                                .abs()
                                < 1e-11
                        );
                    }
                }
            }
        }
        assert!(RingTransfer::decode(&[0.1234567], 8).is_none());
        assert!(RingTransfer::decode(&[1.0], 17).is_none());
        for exponent in 0..=16 {
            let ring = RingTransfer::decode(&[0.0, 1.0, -1.0], exponent).unwrap();
            assert!(ring.data.iter().all(|entry| entry.exponent == 0));
        }
    }

    #[test]
    fn clifford_tableaux_reconstruct_random_native_circuits() {
        let mut rng = Rng::new(41537);
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let native = NativeClifford::new(&lib).unwrap();
            for _ in 0..24 {
                let mut word = Vec::new();
                for _ in 0..32 {
                    let bit = rng.usize(q);
                    word.push(native.native[bit][rng.usize(3)]);
                    if q > 1 && rng.random01() < 0.4 {
                        let a = rng.usize(q);
                        let mut b = rng.usize(q - 1);
                        if b >= a {
                            b += 1;
                        }
                        word.push(native.cx[a][b].unwrap());
                    }
                }
                let source = Circuit::new(word, q, &lib);
                let ptm = Transfer::from_matrix(source.matrix(), q);
                assert!(ptm.cost(1 << (2 * q)) < 1e-10);
                let word = clifford_word(ptm.clifford_images(q).unwrap(), &lib, &native).unwrap();
                let result = Circuit::new(word, q, &lib);
                assert!(same(source.matrix(), result.matrix()));
                assert_eq!(result.count(&["t".into(), "tdg".into()], &lib), 0);
            }
        }
    }
    #[test]
    fn transfer_rotations_and_collision_scores_match_dense_native_operators() {
        let mut rng = Rng::new(8162);
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let native = NativeClifford::new(&lib).unwrap();
            let d = 1 << (2 * q);
            let word = (0..14)
                .map(|_| lib.basic[rng.usize(lib.basic.len())])
                .collect();
            let source = Circuit::new(word, q, &lib);
            let mut ptm = Transfer::from_matrix(source.matrix(), q);
            let count = source.count(&["t".into(), "tdg".into()], &lib);
            ptm.ring = RingTransfer::decode(&ptm.data, count);
            assert!(ptm.ring.is_some());
            let pairs = rotation_pairs(q);
            for _ in 0..8 {
                let axis = 1 + rng.usize(d - 1);
                let word = rotation_word(axis, &lib, &native);
                let rotation = Circuit::new(word, q, &lib);
                assert_eq!(rotation.count(&["t".into(), "tdg".into()], &lib), 1);
                for direction in [-1.0, 1.0] {
                    let matrix = if direction < 0.0 {
                        rotation.matrix().adjoint()
                    } else {
                        rotation.matrix().clone()
                    };
                    let expected = Transfer::from_matrix(&matrix.mul(source.matrix()), q);
                    let mut actual = ptm.clone();
                    ptm.rotate_into(d, &pairs[axis], direction, &mut actual);
                    assert!(actual
                        .data
                        .iter()
                        .zip(&expected.data)
                        .all(|(a, b)| (a - b).abs() < 1e-11));
                    assert!((actual.fourth - expected.fourth).abs() < 1e-9);
                    let ring = actual.ring.as_ref().unwrap();
                    assert!(ring
                        .data
                        .iter()
                        .zip(&expected.data)
                        .all(|(&a, &b)| (a.value() - b).abs() < 1e-10));
                    let mut histogram = [0; 64];
                    for entry in &ring.data {
                        histogram[entry.exponent] += 1;
                    }
                    assert_eq!(ring.histogram, histogram);
                    let predicted = ptm.ring.as_ref().unwrap().rotated_energy(d, &pairs[axis]);
                    assert!((predicted - ring.energy_of(&histogram)).abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn every_two_qubit_clifford_remainder_reconstructs_exactly() {
        let q = 2;
        let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
        let native = NativeClifford::new(&lib).unwrap();
        let initial: Vec<_> = (0..2 * q)
            .map(|bit| Pauli {
                x: if bit < q { 1 << bit } else { 0 },
                z: if bit >= q { 1 << (bit - q) } else { 0 },
                negative: false,
            })
            .collect();
        let key = |images: &[Pauli]| {
            images.iter().enumerate().fold(0usize, |key, (i, p)| {
                key | ((p.x | (p.z << q) | (usize::from(p.negative) << (2 * q)))
                    << (i * (2 * q + 1)))
            })
        };
        let mut visited = vec![false; 1 << (2 * q * (2 * q + 1))];
        visited[key(&initial)] = true;
        let mut queue = vec![(initial, Matrix::identity(1 << q))];
        let generators: Vec<_> = native
            .native
            .iter()
            .flat_map(|ids| ids[..3].iter().copied())
            .chain([native.cx[0][1].unwrap(), native.cx[1][0].unwrap()])
            .collect();
        let mut head = 0;
        while head < queue.len() {
            let (images, matrix) = queue[head].clone();
            head += 1;
            for &id in &generators {
                let gate = &lib.gates[id];
                let mut next = images.clone();
                for p in &mut next {
                    match gate.name.as_str() {
                        "h" => p.h(gate.qubits[0]),
                        "s" => p.s(gate.qubits[0], false),
                        "sdg" => p.s(gate.qubits[0], true),
                        "cx" => p.cx(gate.qubits[0], gate.qubits[1]),
                        _ => unreachable!(),
                    }
                }
                let index = key(&next);
                if visited[index] {
                    continue;
                }
                visited[index] = true;
                let expected = gate.matrix.mul(&matrix);
                let transfer = Transfer::from_matrix(&expected, q);
                let observed = transfer.clifford_images(q).unwrap();
                assert_eq!(key(&observed), index);
                let word = clifford_word(next.clone(), &lib, &native).unwrap();
                let actual = Circuit::new(word, q, &lib);
                assert!(same(&expected, actual.matrix()));
                assert_eq!(actual.count(&["t".into(), "tdg".into()], &lib), 0);
                queue.push((next, expected));
            }
        }
        assert_eq!(queue.len(), 11520);
    }
    #[test]
    fn denominator_channels_remain_consistent_through_long_rotation_chains() {
        let mut rng = Rng::new(903527);
        for q in 1..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let native = NativeClifford::new(&lib).unwrap();
            let mut word = Vec::new();
            for _ in 0..6 {
                let bit = rng.usize(q);
                word.extend([native.native[bit][0], native.native[bit][3]]);
            }
            let source = Circuit::new(word, q, &lib);
            let mut matrix = source.matrix().clone();
            let mut state = Transfer::from_matrix(&matrix, q);
            state.ring = RingTransfer::decode(&state.data, 6);
            assert!(state.ring.is_some());
            let pairs = rotation_pairs(q);
            let d = 1 << (2 * q);
            for step in 0..24 {
                let axis = 1 + rng.usize(d - 1);
                let direction = if rng.usize(2) == 0 { -1.0 } else { 1.0 };
                let rotation = Circuit::new(rotation_word(axis, &lib, &native), q, &lib);
                let operator = if direction < 0.0 {
                    rotation.matrix().adjoint()
                } else {
                    rotation.matrix().clone()
                };
                matrix = operator.mul(&matrix);
                let mut next = state.clone();
                state.rotate_into(d, &pairs[axis], direction, &mut next);
                let expected = Transfer::from_matrix(&matrix, q);
                let ring = next.ring.as_ref().unwrap();
                assert!(
                    ring.data
                        .iter()
                        .zip(&expected.data)
                        .all(|(&a, &b)| (a.value() - b).abs() < 1e-9),
                    "q={q}, step={step}"
                );
                assert!(ring.data.iter().all(|entry| entry.exponent <= 6 + step + 1));
                for (row, support) in ring.support.iter().enumerate() {
                    let expected: Vec<_> = ring.data[row * d..(row + 1) * d]
                        .iter()
                        .enumerate()
                        .filter_map(|(c, x)| (x.a != 0 || x.b != 0).then_some(c))
                        .collect();
                    assert_eq!(*support, expected);
                }
                let predicted = state.ring.as_ref().unwrap().rotated_energy(d, &pairs[axis]);
                assert!((predicted - ring.energy_of(&ring.histogram)).abs() < 1e-12);
                state = next;
            }
        }
    }
    #[test]
    fn cached_prefix_costs_preserve_complete_annealing_trajectories() {
        let mut source_rng = Rng::new(7352061);
        for q in 2..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let native = NativeClifford::new(&lib).unwrap();
            for count in [4, 8, 12] {
                let mut word = Vec::new();
                for _ in 0..count {
                    for _ in 0..2 {
                        let a = source_rng.usize(q);
                        let mut b = source_rng.usize(q - 1);
                        if b >= a {
                            b += 1;
                        }
                        word.extend([
                            native.native[a][0],
                            native.cx[a][b].unwrap(),
                            native.native[b][1],
                        ]);
                    }
                    word.push(native.native[source_rng.usize(q)][3]);
                }
                let circuit = Circuit::new(word, q, &lib);
                let target = Target::new(
                    PartialMatrix::new(
                        circuit.matrix().clone(),
                        vec![true; 1 << (2 * q)],
                        "cache-replay",
                    )
                    .unwrap(),
                    false,
                    false,
                );
                let mut kernel = PauliKernel::new(&target, &lib, count).unwrap();
                assert!(kernel.initial.ring.is_some());
                for exact in [true, false] {
                    if !exact {
                        kernel.initial.ring = None;
                    }
                    for seed in [5167, 871023] {
                        let options = SearchOptions {
                            iterations_factor: 8.0,
                            ..SearchOptions::tuned()
                        };
                        let mut cached_rng = Rng::new(seed);
                        let mut direct_rng = cached_rng.clone();
                        let cached =
                            kernel.run_impl::<true>(&target, &lib, &mut cached_rng, &options, None);
                        let direct = kernel.run_impl::<false>(
                            &target,
                            &lib,
                            &mut direct_rng,
                            &options,
                            None,
                        );
                        assert_eq!(
                            cached.circuit.gates, direct.circuit.gates,
                            "q={q}, T={count}, exact={exact}, seed={seed}"
                        );
                        assert_eq!(cached.found, direct.found);
                        assert_eq!(cached.steps, direct.steps);
                        assert_eq!(cached.accepted, direct.accepted);
                        assert_eq!(cached.best_eq.to_bits(), direct.best_eq.to_bits());
                        assert_eq!(cached.best_energy.to_bits(), direct.best_energy.to_bits());
                        for _ in 0..16 {
                            assert_eq!(cached_rng.next_u32(), direct_rng.next_u32());
                        }
                    }
                }
            }
        }
    }
    #[test]
    fn annealed_pauli_words_reproduce_unstructured_small_t_targets() {
        let mut source_rng = Rng::new(62738);
        for q in 2..=4 {
            let lib = GateLibrary::load(q, "data/gates/CliffordT", "").unwrap();
            let native = NativeClifford::new(&lib).unwrap();
            for count in 1..=4 {
                let mut word = Vec::new();
                for _ in 0..count {
                    for _ in 0..3 {
                        let bit = source_rng.usize(q);
                        word.push(native.native[bit][source_rng.usize(3)]);
                        let a = source_rng.usize(q);
                        let mut b = source_rng.usize(q - 1);
                        if b >= a {
                            b += 1;
                        }
                        word.push(native.cx[a][b].unwrap());
                    }
                    word.push(native.native[source_rng.usize(q)][3]);
                }
                let circuit = Circuit::new(word, q, &lib);
                let target = Target::new(
                    PartialMatrix::new(
                        circuit.matrix().clone(),
                        vec![true; 1 << (2 * q)],
                        "random",
                    )
                    .unwrap(),
                    false,
                    false,
                );
                let kernel = PauliKernel::new(&target, &lib, count).unwrap();
                let options = SearchOptions {
                    iterations_factor: 4.0,
                    ..SearchOptions::tuned()
                };
                let mut solved = false;
                for seed in 0..8 {
                    let result =
                        kernel.run_until(&target, &lib, &mut Rng::new(9153 + seed), &options, None);
                    if result.found {
                        assert!(same(circuit.matrix(), result.circuit.matrix()));
                        assert!(result.circuit.count(&["t".into(), "tdg".into()], &lib) <= count);
                        solved = true;
                        break;
                    }
                }
                assert!(solved, "q={q}, T={count}");
            }
        }
    }
}
