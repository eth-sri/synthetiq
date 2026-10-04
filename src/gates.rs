//! Gate-file loading, qubit embedding, and composite decomposition.
use crate::complex::Complex;
use crate::matrix::{next_permutation, Matrix};
use std::collections::HashMap;
use std::fs;
use std::path::{Path, PathBuf};

#[derive(Clone, Debug)]
pub struct Gate {
    pub name: String,
    pub matrix: Matrix,
    pub qubits: Vec<usize>,
    pub cost: f64,
    /// Empty for basic gates; composite entries refer to basic gate IDs.
    pub decomposition: Vec<usize>,
}
impl Gate {
    pub fn is_basic(&self) -> bool {
        self.decomposition.is_empty()
    }
}

#[derive(Clone, Debug)]
pub struct GateLibrary {
    pub n_qubits: usize,
    pub gates: Vec<Gate>,
    pub basic: Vec<usize>,
    pub composite: Vec<usize>,
    pub basic_by_name: Vec<Vec<usize>>,
    pub composite_by_name: Vec<Vec<usize>>,
    pub all: Vec<usize>,
    pub readable: Vec<usize>,
    pub identity: usize,
    pub max_cost_basic: f64,
    pub max_cost_all: f64,
    lookup: HashMap<(String, Vec<usize>), usize>,
}
impl GateLibrary {
    pub fn load(n_qubits: usize, basic: &str, composite: &str) -> Result<Self, String> {
        if n_qubits == 0 || n_qubits > 16 {
            return Err("The gate library requires between 1 and 16 qubits".into());
        }
        let identity = Gate {
            name: "id".into(),
            matrix: Matrix::identity(1 << n_qubits),
            qubits: vec![0],
            cost: 0.0,
            decomposition: Vec::new(),
        };
        let mut lib = Self {
            n_qubits,
            gates: vec![identity],
            basic: Vec::new(),
            composite: Vec::new(),
            basic_by_name: Vec::new(),
            composite_by_name: Vec::new(),
            all: vec![0],
            readable: vec![0],
            identity: 0,
            max_cost_basic: 0.0,
            max_cost_all: 0.0,
            lookup: HashMap::new(),
        };
        lib.lookup.insert(("id".into(), vec![0]), 0);
        for path in folder_files(basic)? {
            lib.load_basic(&path)?;
        }
        lib.max_cost_all = lib.max_cost_basic;
        for path in folder_files(composite)? {
            lib.load_composite(&path, false)?;
        }
        Ok(lib)
    }

    /// Add the finite one-qubit Clifford closure of the supplied H/S/SDG gates
    /// as search macros. Every macro stores an exact native-gate decomposition.
    /// A non-Clifford/infinite custom closure is rejected before mutating the library.
    pub fn add_clifford_macros(&mut self) -> usize {
        let generators: Vec<_> = ["h", "s", "sdg"]
            .iter()
            .filter_map(|name| self.find(name, &[0]))
            .filter(|&id| self.gates[id].is_basic())
            .collect();
        if generators.len() < 2 {
            return 0;
        }
        let mut states = vec![(Matrix::identity(1 << self.n_qubits), Vec::<usize>::new())];
        let mut next = 0;
        while next < states.len() {
            for &id in &generators {
                let product = self.gates[id].matrix.mul(&states[next].0);
                if states.iter().any(|(m, _)| projectively_equal(m, &product)) {
                    continue;
                }
                if states.len() == 24 {
                    return 0;
                }
                let mut path = states[next].1.clone();
                path.push(id);
                states.push((product, path));
            }
            next += 1;
        }
        if states.len() != 24 {
            return 0;
        }
        let native: Vec<_> = self
            .all
            .iter()
            .map(|&id| self.gates[id].matrix.clone())
            .collect();
        let mut added = 0;
        for (index, (matrix, path)) in states.into_iter().enumerate() {
            if path.is_empty() || native.iter().any(|m| projectively_equal(m, &matrix)) {
                continue;
            }
            let mut group = Vec::new();
            for q in 0..self.n_qubits {
                let mut order: Vec<_> = (0..self.n_qubits).collect();
                order.swap(0, q);
                let decomposition: Vec<_> = path
                    .iter()
                    .map(|&id| self.find(&self.gates[id].name, &[q]).unwrap())
                    .collect();
                let cost = decomposition.iter().map(|&id| self.gates[id].cost).sum();
                let id = self.insert(Gate {
                    name: format!("__clifford_{index}"),
                    matrix: matrix.permuted(&order),
                    qubits: vec![q],
                    cost,
                    decomposition,
                });
                self.composite.push(id);
                self.all.push(id);
                group.push(id);
                added += 1;
            }
            self.composite_by_name.push(group);
        }
        self.max_cost_all = self.gates.iter().map(|g| g.cost).fold(0.0, f64::max);
        added
    }

    pub fn add_composite_folder(&mut self, folder: &str) -> Result<(), String> {
        for path in folder_files(folder)? {
            self.load_composite(&path, false)?;
        }
        self.max_cost_all = self.gates.iter().map(|g| g.cost).fold(0.0, f64::max);
        Ok(())
    }
    pub fn add_read_folder(&mut self, folder: &str) -> Result<(), String> {
        for path in folder_files(folder)? {
            self.load_composite(&path, true)?;
        }
        Ok(())
    }
    pub fn find(&self, name: &str, qubits: &[usize]) -> Option<usize> {
        self.lookup
            .get(&(name.to_owned(), qubits.to_vec()))
            .copied()
    }
    pub fn find_all(&self, name: &str, qubits: &[usize]) -> Option<usize> {
        self.all
            .iter()
            .copied()
            .find(|&id| self.gates[id].name == name && self.gates[id].qubits == qubits)
    }
    /// Preserve CircuitHelper::invertGate's exact elementwise conjugate match.
    /// This is intentionally distinct from finding a mathematical adjoint.
    pub fn inverse_gate(&self, id: usize) -> Option<usize> {
        let matrix = &self.gates[id].matrix;
        self.all.iter().copied().find(|&other| {
            matrix
                .data
                .iter()
                .zip(&self.gates[other].matrix.data)
                .all(|(&a, &b)| a == b.conj())
        })
    }
    fn insert(&mut self, gate: Gate) -> usize {
        let id = self.gates.len();
        self.lookup
            .entry((gate.name.clone(), gate.qubits.clone()))
            .or_insert(id);
        self.gates.push(gate);
        self.readable.push(id);
        id
    }
    fn load_basic(&mut self, path: &Path) -> Result<(), String> {
        let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let mut lines = text.lines();
        let name = next_nonempty(&mut lines, path, "gate name")?
            .trim()
            .to_owned();
        let n_qubits = parse_usize(next_nonempty(&mut lines, path, "qubit count")?, path)?;
        if n_qubits > self.n_qubits {
            return Ok(());
        }
        let cost = next_nonempty(&mut lines, path, "gate cost")?
            .trim()
            .parse::<f64>()
            .map_err(|e| format!("{}: invalid cost: {e}", path.display()))?;
        if !cost.is_finite() || cost < 0.0 {
            return Err(format!(
                "{}: gate cost must be finite and nonnegative",
                path.display()
            ));
        }
        let acting = lines
            .next()
            .ok_or_else(|| format!("{}: missing acting qubits", path.display()))?;
        let qubits = parse_plain_qubits(acting)?;
        validate_qubits(&qubits, n_qubits, path)?;
        let data = parse_complex_values(&lines.collect::<Vec<_>>().join("\n"))?;
        let matrix = Matrix::from_data(1 << n_qubits, data)?.kron_identity(self.n_qubits);
        self.max_cost_basic = self.max_cost_basic.max(cost);
        let mut group = Vec::new();
        // Enumerate the first full permutation for each distinct placement.
        // The reference enumerates n! permutations and discards the repeats.
        for order in distinct_placements(self.n_qubits, &qubits) {
            let mapped: Vec<_> = qubits.iter().map(|&q| order[q]).collect();
            if self.find(&name, &mapped).is_some() {
                continue;
            }
            let id = self.insert(Gate {
                name: name.clone(),
                matrix: matrix.permuted(&order),
                qubits: mapped,
                cost,
                decomposition: Vec::new(),
            });
            self.basic.push(id);
            self.all.push(id);
            group.push(id);
        }
        self.basic_by_name.push(group);
        Ok(())
    }
    fn load_composite(&mut self, path: &Path, read_only: bool) -> Result<(), String> {
        let text = fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        let (name, required, qubits, body) = if path.extension().is_some_and(|x| x == "qasm") {
            let n = qasm_qubit_count(&text)?;
            // Most composite libraries use .txt. Keep the legacy path-derived
            // name for .qasm composites to preserve accepted input names.
            let path_text = path.to_string_lossy();
            let start = path_text.rfind('/').unwrap_or(0);
            let length = path_text.rfind('.').unwrap_or(path_text.len());
            let end = (start + length).min(path_text.len());
            (
                path_text[start..end].to_owned(),
                n,
                (0..n).collect(),
                text.lines()
                    .filter(|line| !is_qasm_header(line))
                    .map(str::to_owned)
                    .collect::<Vec<_>>(),
            )
        } else {
            let mut lines = text.lines();
            let name = next_nonempty(&mut lines, path, "gate name")?
                .trim()
                .to_owned();
            let n = parse_usize(next_nonempty(&mut lines, path, "qubit count")?, path)?;
            if n > self.n_qubits {
                return Ok(());
            }
            let qubits = parse_plain_qubits(
                lines
                    .next()
                    .ok_or_else(|| format!("{}: missing acting qubits", path.display()))?,
            )?;
            (name, n, qubits, lines.map(str::to_owned).collect())
        };
        if required > self.n_qubits {
            return Ok(());
        }
        if self.basic.is_empty() {
            return Err(format!(
                "{}: composite gate needs a basic gate library",
                path.display()
            ));
        }
        validate_qubits(&qubits, required, path)?;
        let mut decomposition = Vec::new();
        let mut cost = 0.0;
        for line in &body {
            if let Some((name, acting)) = parse_gate_line(line)? {
                let id = self
                    .basic
                    .iter()
                    .copied()
                    .find(|&id| self.gates[id].name == name && self.gates[id].qubits == acting)
                    .ok_or_else(|| {
                        format!("{}: gate cannot be constructed: {line}", path.display())
                    })?;
                cost += self.gates[id].cost;
                decomposition.push(id);
            }
        }
        if !read_only {
            self.max_cost_all = self.max_cost_all.max(cost);
        }
        let mut group = Vec::new();
        let mut order: Vec<_> = (0..self.n_qubits).collect();
        loop {
            // Match the reference loader's candidate multiplicity: its presence
            // check repeats the first acting qubit. Multi-qubit composites
            // consequently retain all full permutations, including duplicates.
            let check = vec![order[qubits[0]]; qubits.len()];
            if self.find(&name, &check).is_none() {
                let mapped: Vec<_> = qubits.iter().map(|&q| order[q]).collect();
                let mut mapped_decomp = Vec::with_capacity(decomposition.len());
                let mut matrix = Matrix::identity(1 << self.n_qubits);
                let mut scratch = Matrix::zero(matrix.n);
                for &id in &decomposition {
                    let gate = &self.gates[id];
                    let acted: Vec<_> = gate.qubits.iter().map(|&q| order[q]).collect();
                    let mapped_id = self
                        .basic
                        .iter()
                        .copied()
                        .find(|&id| {
                            self.gates[id].name == gate.name && self.gates[id].qubits == acted
                        })
                        .ok_or_else(|| {
                            format!("{}: composite gate cannot be permuted", path.display())
                        })?;
                    self.gates[mapped_id].matrix.mul_into(&matrix, &mut scratch);
                    std::mem::swap(&mut matrix, &mut scratch);
                    mapped_decomp.push(mapped_id);
                }
                let id = self.insert(Gate {
                    name: name.clone(),
                    matrix,
                    qubits: mapped,
                    cost,
                    decomposition: mapped_decomp,
                });
                if !read_only {
                    self.composite.push(id);
                    self.all.push(id);
                    group.push(id);
                }
            }
            if !next_permutation(&mut order) {
                break;
            }
        }
        if !read_only {
            self.composite_by_name.push(group);
        }
        Ok(())
    }
}
fn projectively_equal(a: &Matrix, b: &Matrix) -> bool {
    let overlap = b.trace_conjugate_product(a);
    let length = overlap.abs();
    if length < 1e-12 {
        return false;
    }
    let phase = overlap / length;
    a.data
        .iter()
        .zip(&b.data)
        .all(|(&x, &y)| (x - phase * y).norm_sqr() < 1e-20)
}

fn folder_files(folder: &str) -> Result<Vec<PathBuf>, String> {
    if folder.is_empty() || !Path::new(folder).is_dir() {
        return Ok(Vec::new());
    }
    fs::read_dir(folder)
        .map_err(|e| format!("{folder}: {e}"))?
        .filter_map(|e| match e {
            Ok(entry) if entry.path().is_file() => Some(Ok(entry.path())),
            Ok(_) => None,
            Err(e) => Some(Err(format!("{folder}: {e}"))),
        })
        .collect()
}
fn next_nonempty<'a>(
    lines: &mut impl Iterator<Item = &'a str>,
    path: &Path,
    field: &str,
) -> Result<&'a str, String> {
    lines
        .find(|x| !x.trim().is_empty())
        .ok_or_else(|| format!("{}: missing {field}", path.display()))
}
fn parse_usize(s: &str, path: &Path) -> Result<usize, String> {
    s.trim()
        .parse()
        .map_err(|e| format!("{}: invalid integer: {e}", path.display()))
}
fn parse_plain_qubits(s: &str) -> Result<Vec<usize>, String> {
    s.split_whitespace()
        .map(|x| x.parse().map_err(|_| format!("Invalid qubit index: {x}")))
        .collect()
}
fn validate_qubits(qubits: &[usize], n: usize, path: &Path) -> Result<(), String> {
    if qubits.is_empty() || qubits.iter().any(|&q| q >= n) {
        return Err(format!(
            "{}: acting qubits must be nonempty and within the gate width",
            path.display()
        ));
    }
    for (i, q) in qubits.iter().enumerate() {
        if qubits[..i].contains(q) {
            return Err(format!("{}: repeated acting qubit", path.display()));
        }
    }
    Ok(())
}
/// All distinct acting-qubit placements in the same first-encounter order as
/// lexicographic full permutations, without factorial work on idle qubits.
fn distinct_placements(n: usize, qubits: &[usize]) -> Vec<Vec<usize>> {
    fn visit(
        pos: usize,
        n: usize,
        relevant: &[usize],
        order: &mut [usize],
        used: &mut [bool],
        out: &mut Vec<Vec<usize>>,
    ) {
        if pos == n {
            out.push(order.to_vec());
            return;
        }
        if relevant.contains(&pos) {
            for value in 0..n {
                if !used[value] {
                    used[value] = true;
                    order[pos] = value;
                    visit(pos + 1, n, relevant, order, used, out);
                    used[value] = false;
                }
            }
        } else {
            // Idle positions before an active one still determine permutation
            // order, so enumerate then deduplicate outside for those positions.
            if relevant.iter().any(|&q| q > pos) {
                for value in 0..n {
                    if !used[value] {
                        used[value] = true;
                        order[pos] = value;
                        visit(pos + 1, n, relevant, order, used, out);
                        used[value] = false;
                    }
                }
            } else {
                let remaining: Vec<_> = (0..n).filter(|&v| !used[v]).collect();
                order[pos..].copy_from_slice(&remaining);
                out.push(order.to_vec());
            }
        }
    }
    let mut permutations = Vec::new();
    visit(
        0,
        n,
        qubits,
        &mut vec![0; n],
        &mut vec![false; n],
        &mut permutations,
    );
    let mut seen = std::collections::HashSet::new();
    permutations.retain(|order| seen.insert(qubits.iter().map(|&q| order[q]).collect::<Vec<_>>()));
    permutations
}

/// Parse the stream format accepted by std::complex<double>, including spaces
/// inside parentheses and real-only entries.
pub fn parse_complex_values(text: &str) -> Result<Vec<Complex>, String> {
    let mut values = Vec::new();
    let mut rest = text.trim();
    while !rest.is_empty() {
        let end = if rest.starts_with('(') {
            rest.find(')')
                .map(|x| x + 1)
                .ok_or("Unclosed complex number")?
        } else {
            rest.find(char::is_whitespace).unwrap_or(rest.len())
        };
        let value: Complex = rest[..end].parse()?;
        if !value.re.is_finite() || !value.im.is_finite() {
            return Err("Gate matrices must contain finite numbers".into());
        }
        values.push(value);
        rest = rest[end..].trim_start();
    }
    Ok(values)
}
pub fn is_qasm_header(line: &str) -> bool {
    let line = line.trim();
    line.starts_with("OPENQASM") || line.starts_with("include") || line.starts_with("qreg")
}
pub fn qasm_qubit_count(text: &str) -> Result<usize, String> {
    let line = text
        .lines()
        .map(str::trim)
        .find(|line| line.starts_with("qreg "))
        .ok_or("Missing QASM qreg declaration")?;
    let start = line.find('[').ok_or("Invalid QASM qreg declaration")? + 1;
    let end = line[start..]
        .find(']')
        .ok_or("Invalid QASM qreg declaration")?
        + start;
    line[start..end]
        .trim()
        .parse()
        .map_err(|_| "Invalid QASM qubit count".into())
}
pub fn parse_gate_line(line: &str) -> Result<Option<(String, Vec<usize>)>, String> {
    let line = line.split("//").next().unwrap_or("").trim();
    if line.is_empty() {
        return Ok(None);
    }
    let end = line.find(char::is_whitespace).unwrap_or(line.len());
    let name = line[..end].to_owned();
    let args = line[end..].trim().trim_end_matches(';').trim();
    if args.is_empty() {
        return Err(format!("Missing acting qubits: {line}"));
    }
    let mut qubits = Vec::new();
    if args.contains('[') {
        let mut rest = args;
        while let Some(start) = rest.find('[') {
            rest = &rest[start + 1..];
            let end = rest
                .find(']')
                .ok_or_else(|| format!("Invalid qubit reference: {line}"))?;
            qubits.push(
                rest[..end]
                    .trim()
                    .parse()
                    .map_err(|_| format!("Invalid qubit index: {line}"))?,
            );
            rest = &rest[end + 1..];
        }
    } else {
        for token in args
            .split(|c: char| c.is_whitespace() || c == ',')
            .filter(|s| !s.is_empty())
        {
            qubits.push(
                token
                    .parse()
                    .map_err(|_| format!("Invalid qubit index: {line}"))?,
            );
        }
    }
    if qubits.is_empty() {
        return Err(format!("Missing acting qubits: {line}"));
    }
    Ok(Some((name, qubits)))
}

#[cfg(test)]
mod tests {
    use super::*;
    pub fn basic_lib(n: usize) -> GateLibrary {
        GateLibrary::load(n, "data/gates/basic_gates", "").unwrap()
    }
    #[test]
    fn basic_gate_counts_and_placement() {
        for n in 1..=5 {
            let lib = basic_lib(n);
            assert_eq!(lib.basic.len(), 5 * n + n * (n - 1));
            assert_eq!(lib.gates.len(), lib.basic.len() + 1);
            for target in 0..n {
                let h = &lib.gates[lib.find("h", &[target]).unwrap()];
                assert!(
                    h.matrix
                        .mul(&h.matrix)
                        .max_abs_diff(&Matrix::identity(1 << n))
                        < 1e-14
                );
            }
        }
    }
    #[test]
    fn cx_little_endian_truth_tables() {
        let lib = basic_lib(4);
        for control in 0..4 {
            for target in 0..4 {
                if target == control {
                    continue;
                }
                let gate = &lib.gates[lib.find("cx", &[control, target]).unwrap()];
                for col in 0..16 {
                    let row = if (col >> control) & 1 == 1 {
                        col ^ (1 << target)
                    } else {
                        col
                    };
                    assert_eq!(gate.matrix[(row, col)], Complex::ONE);
                }
            }
        }
    }
    #[test]
    fn exact_reference_inverses() {
        let lib = basic_lib(3);
        for &id in &lib.all {
            let other = lib.inverse_gate(id).unwrap();
            assert_eq!(lib.gates[id].matrix, lib.gates[other].matrix.conjugate());
            assert!(
                lib.gates[id]
                    .matrix
                    .mul(&lib.gates[other].matrix)
                    .max_abs_diff(&Matrix::identity(8))
                    < 1e-14
            );
        }
    }
    #[test]
    fn composite_decomposition_and_permutations() {
        let lib =
            GateLibrary::load(3, "data/gates/basic_gates", "data/gates/composite_ccx").unwrap();
        assert_eq!(lib.composite.len(), 6);
        for &id in &lib.composite {
            let gate = &lib.gates[id];
            let mut matrix = Matrix::identity(8);
            let mut cost = 0.0;
            for &basic in &gate.decomposition {
                matrix = lib.gates[basic].matrix.mul(&matrix);
                cost += lib.gates[basic].cost;
            }
            assert_eq!(matrix, gate.matrix);
            assert_eq!(cost, gate.cost);
        }
    }
    #[test]
    fn composite_retains_legacy_duplicate_candidates() {
        let lib =
            GateLibrary::load(4, "data/gates/basic_gates", "data/gates/composite_rccx").unwrap();
        assert_eq!(lib.composite.len(), 24);
    }
    #[test]
    fn parser_handles_both_gate_formats() {
        assert_eq!(
            parse_gate_line(" cx q[0], q[2]; // comment").unwrap(),
            Some(("cx".into(), vec![0, 2]))
        );
        assert_eq!(
            parse_gate_line("cx 0 2").unwrap(),
            Some(("cx".into(), vec![0, 2]))
        );
        assert_eq!(parse_gate_line("  // empty").unwrap(), None);
        assert!(parse_gate_line("cx q[x]").is_err());
    }
    #[test]
    fn complex_stream_with_spaces() {
        assert_eq!(
            parse_complex_values("(1, 2)\n -3.5 (0,-1e-2)").unwrap(),
            vec![
                Complex::new(1.0, 2.0),
                Complex::new(-3.5, 0.0),
                Complex::new(0.0, -0.01)
            ]
        );
    }
    #[test]
    fn placement_order_matches_exhaustive_reference() {
        for n in 1..=5 {
            for qubits in [vec![0], vec![n - 1], (0..n.min(3)).rev().collect()] {
                let mut expected = Vec::new();
                let mut seen = std::collections::HashSet::new();
                let mut p: Vec<_> = (0..n).collect();
                loop {
                    if seen.insert(qubits.iter().map(|&q| p[q]).collect::<Vec<_>>()) {
                        expected.push(p.clone());
                    }
                    if !next_permutation(&mut p) {
                        break;
                    }
                }
                assert_eq!(distinct_placements(n, &qubits), expected);
            }
        }
    }
}

#[cfg(test)]
mod clifford_macro_tests {
    use super::*;
    use crate::circuit::Circuit;

    #[test]
    fn every_clifford_macro_expands_exactly_on_each_qubit() {
        for nq in 1..=4 {
            let mut lib = GateLibrary::load(nq, "data/gates/CliffordT", "").unwrap();
            let first = lib.gates.len();
            assert_eq!(lib.add_clifford_macros(), 20 * nq);
            for id in first..lib.gates.len() {
                let gate = &lib.gates[id];
                assert!(!gate.decomposition.is_empty());
                let mut circuit = Circuit::new(vec![id], nq, &lib);
                let before = circuit.matrix().clone();
                circuit.expand(&lib);
                assert!(before.max_abs_diff(circuit.matrix()) < 1e-12);
                assert_eq!(circuit.count(&["t".into(), "tdg".into()], &lib), 0);
                assert!((circuit.cost() - gate.cost).abs() < 1e-12);
            }
        }
    }
}
