//! Command-line compatibility with the original `bin/main`.
use crate::{
    circuit::Circuit,
    gates::GateLibrary,
    partial::{PartialMatrix, Target},
    resynthesis::Resynth,
    rng::Rng,
    search::{self, GateScheme, SearchOptions},
};
use std::{
    fs,
    io::Write,
    sync::{Mutex, OnceLock},
    time::Instant,
};

#[derive(Clone, Debug)]
pub struct Config {
    pub input: String,
    pub auto_search: bool,
    pub general_only: bool,
    pub general_every: usize,
    pub structured_time_slice: f64,
    pub bounded_structure: bool,
    pub legacy_fallback_after: f64,
    pub resource_fallback_after: f64,
    specified: Vec<String>,
    pub output: String,
    pub format: String,
    pub basic: String,
    pub composite: String,
    pub times_file: String,
    pub threads: usize,
    pub circuits: usize,
    pub ancillas: usize,
    pub dirty: usize,
    pub tcount: i64,
    pub tdepth: i64,
    pub gatecount: i64,
    pub cost_required: f64,
    pub time: f64,
    pub save_times: bool,
    pub expand: bool,
    pub clifford_macros: bool,
    pub phase_restarts: usize,
    pub phase_fold: bool,
    pub pauli_fold: bool,
    pub pauli_restarts: usize,
    pub clifford_fold: bool,
    pub phase_factor: bool,
    pub clifford_factor: bool,
    pub permutation_restarts: usize,
    pub single_qubit_restarts: usize,
    pub nonlinear_phase_restarts: usize,
    pub affine_search: bool,
    pub rccx_macros: bool,
    pub independent: bool,
    pub inverse: bool,
    pub save_all: bool,
    pub save_circuits: bool,
    pub resynth: bool,
    pub update_scheme: bool,
    pub optimization: usize,
    pub optimize_depth: bool,
    pub depth_gates: Vec<String>,
    pub seed: u64,
    pub search: SearchOptions,
    pub scheme: GateScheme,
}

impl Config {
    pub fn parse(args: &[String]) -> Result<Self, String> {
        let input = args
            .first()
            .ok_or("missing input specification (use --help)")?
            .clone();
        if input.starts_with('-') {
            return Err(format!("expected input file, got {input}"));
        }
        let (stem, ext) = input.rsplit_once('.').unwrap_or((&input, "txt"));
        let mut c = Self {
            input: input.clone(),
            auto_search: true,
            general_only: false,
            general_every: 0,
            structured_time_slice: 0.0,
            bounded_structure: false,
            legacy_fallback_after: 0.0,
            resource_fallback_after: 0.0,
            specified: args.to_vec(),
            output: format!("{stem}/"),
            format: ext.into(),
            basic: "CliffordT".into(),
            composite: "composite_gates".into(),
            times_file: "data/times.csv".into(),
            threads: 1,
            circuits: 10,
            ancillas: 0,
            dirty: 0,
            tcount: -1,
            tdepth: -1,
            gatecount: -1,
            cost_required: -1.0,
            time: 100.0,
            save_times: false,
            expand: false,
            clifford_macros: false,
            phase_restarts: 0,
            phase_fold: false,
            pauli_fold: false,
            pauli_restarts: 0,
            clifford_fold: false,
            phase_factor: false,
            clifford_factor: false,
            permutation_restarts: 0,
            single_qubit_restarts: 0,
            nonlinear_phase_restarts: 0,
            affine_search: false,
            rccx_macros: false,
            independent: true,
            inverse: false,
            save_all: false,
            save_circuits: true,
            resynth: true,
            update_scheme: true,
            optimization: 12,
            optimize_depth: true,
            depth_gates: vec!["t".into(), "tdg".into()],
            seed: 0,
            search: SearchOptions::tuned(),
            scheme: GateScheme::default(),
        };
        let mut quality_objective = None;
        let mut absolute_input = false;
        let mut absolute_output = false;
        let mut absolute_gates = false;
        let mut i = 1;
        while i < args.len() {
            let arg = args[i].as_str();
            macro_rules! value {
                () => {{
                    i += 1;
                    args.get(i)
                        .ok_or_else(|| format!("{arg} requires a value"))?
                }};
            }
            macro_rules! number {
                () => {
                    value!()
                        .parse()
                        .map_err(|_| format!("invalid value for {arg}"))?
                };
            }
            match arg {
                "--threads" | "-h" => c.threads = number!(),
                "--save" => c.save_times = true,
                "--time" | "-t" => c.time = number!(),
                "--circuits" | "-c" => c.circuits = number!(),
                "--gates" | "-g" => {
                    let n: usize = number!();
                    c.scheme
                        .set_min_start_gates((c.scheme.min_factor * n as f64) as usize);
                    c.scheme
                        .set_max_start_gates((c.scheme.max_factor * n as f64) as usize);
                    c.scheme.set_start_best_gates(n);
                    c.scheme.reset();
                }
                "--output" | "-o" => c.output = format!("{}/", value!()),
                "--format" | "-f" => c.format = value!().clone(),
                "--ancilla" | "-a" => c.ancillas = number!(),
                "--dirty" | "-d" => c.dirty = number!(),
                "--expand" | "-e" => c.expand = true,
                "--clifford-macros" => c.clifford_macros = true,
                "--general-every" => c.general_every = number!(),
                "--structured-time-slice" => c.structured_time_slice = number!(),
                "--bounded-structure" => c.bounded_structure = true,
                "--legacy-fallback-after" => c.legacy_fallback_after = number!(),
                "--resource-fallback-after" => c.resource_fallback_after = number!(),
                "--phase-gadgets" => c.phase_restarts = 16,
                "--phase-fold" => c.phase_fold = true,
                "--pauli-fold" => c.pauli_fold = true,
                "--pauli-restarts" => c.pauli_restarts = number!(),
                "--clifford-fold" => c.clifford_fold = true,
                "--phase-factor" => c.phase_factor = true,
                "--clifford-factor" => c.clifford_factor = true,
                "--no-clifford-factor" => c.clifford_factor = false,
                "--no-phase-fold" => c.phase_fold = false,
                "--no-pauli-fold" => c.pauli_fold = false,
                "--affine-search" => c.affine_search = true,
                "--rccx-macros" => c.rccx_macros = true,
                "--phase-restarts" => c.phase_restarts = number!(),
                "--single-qubit-restarts" => c.single_qubit_restarts = number!(),
                "--nonlinear-phase-restarts" => c.nonlinear_phase_restarts = number!(),
                "--permutation-gadgets" => c.permutation_restarts = 4,
                "--permutation-restarts" => c.permutation_restarts = number!(),
                "--qubit_independence" | "-q" => c.independent = false,
                "--tcount" | "-tc" => c.tcount = number!(),
                "--tdepth" | "-td" => c.tdepth = number!(),
                "--gate-count" | "-gc" => c.gatecount = number!(),
                "--cost-required" | "-cr" => c.cost_required = number!(),
                "--gate-set" | "-gs" => c.basic = value!().clone(),
                "--composite-gates" | "-cg" => c.composite = value!().clone(),
                "--save-all" | "-sa" => c.save_all = true,
                "--start-temp" | "-st" => c.search.start_temp_base = number!(),
                "--epsilon" | "-eps" => c.search.epsilon = number!(),
                "--beta" => c.scheme.beta = number!(),
                "--fmin" => c.scheme.min_factor = number!(),
                "--fmax" => c.scheme.max_factor = number!(),
                "--pcomp" => c.search.pcomp = number!(),
                "--pid" => c.search.pid = number!(),
                "--no-perms" => c.search.enable_permutations = false,
                "--no-resynth" => c.resynth = false,
                "--simple" => c.search.simple_cost = true,
                "--n-norm" => c.search.n_norm = number!(),
                "--iterations-factor" => c.search.iterations_factor = number!(),
                "--gs-no-update" => c.update_scheme = false,
                "--n-start-gates" => {
                    let n: usize = number!();
                    c.scheme.set_min_start_gates(n);
                    c.scheme
                        .set_max_start_gates(n.checked_add(1).ok_or("gate count overflow")?);
                    c.update_scheme = false;
                }
                "--times-file" => c.times_file = value!().clone(),
                "--optimization-number" => c.optimization = number!(),
                "--no-save-circuits" => c.save_circuits = false,
                "--no-optimize-depth" => c.optimize_depth = false,
                "--inverse-independent" => c.inverse = true,
                "--depth-gates" => c.depth_gates = value!().split(',').map(String::from).collect(),
                "--absolute-input" => absolute_input = true,
                "--absolute-output" => absolute_output = true,
                "--absolute-gates" => absolute_gates = true,
                "--seed" => c.seed = number!(),
                "--search-mode" => {
                    let mode = value!().as_str();
                    c.auto_search = matches!(mode, "auto" | "general");
                    c.general_only = mode == "general";
                    c.search.optimized = match mode {
                        "legacy" => false,
                        "sweep" | "auto" | "general" => true,
                        _ => {
                            return Err("search-mode must be legacy, sweep, general, or auto".into())
                        }
                    }
                }
                "--cost-power" => c.search.cost_power = number!(),
                "--sweep-cooling" => c.search.sweep_cooling = number!(),
                "--stall-sweeps" => c.search.stall_sweeps = number!(),
                "--gate-prior" => c.search.gate_prior = number!(),
                "--sweep-original" => c.search.sweep_permutations = false,
                "--sweep-permutations" => c.search.sweep_permutations = true,
                "--no-cache-rollback" => c.search.cache_rollback = false,
                "--improved-post" => c.search.improved_post = true,
                "--quality-objective" => {
                    quality_objective = Some(match value!().as_str() {
                        "tcount" => search::QualityMetric::TCount,
                        "tdepth" => search::QualityMetric::TDepth,
                        "gatecount" => search::QualityMetric::GateCount,
                        "cost" => search::QualityMetric::WeightedCost,
                        _ => {
                            return Err(
                                "quality-objective must be tcount, tdepth, gatecount, or cost"
                                    .into(),
                            )
                        }
                    });
                }
                "--quality-weight" => c.search.quality_weight = number!(),
                "--swap-probability" => c.search.swap_probability = number!(),
                "--masked-sweep" => c.search.masked_sweep = true,
                "--masked-proposals" => c.search.masked_proposals = number!(),
                "--quality-slack" => c.search.quality_slack = Some(number!()),
                _ => return Err(format!("unknown option: {arg}")),
            }
            i += 1;
        }
        if !absolute_input {
            c.input = format!("data/input/{}", c.input);
        }
        if !absolute_output {
            c.output = format!("data/output/{}", c.output);
        }
        if !absolute_gates {
            c.basic = format!("data/gates/{}", c.basic);
            c.composite = format!("data/gates/{}", c.composite);
        }
        if c.tcount >= 0 {
            c.search.quality_goal = Some(c.tcount as f64);
            c.search.quality_metric = search::QualityMetric::TCount;
        } else if c.gatecount >= 0 {
            c.search.quality_goal = Some(c.gatecount as f64);
            c.search.quality_metric = search::QualityMetric::GateCount;
        } else if c.cost_required >= 0.0 {
            c.search.quality_goal = Some(c.cost_required);
        } else if c.tdepth >= 0 {
            c.search.quality_metric = search::QualityMetric::TDepth;
            c.search.quality_goal = Some(c.tdepth as f64);
        }
        if let Some(metric) = quality_objective {
            let goal = match metric {
                search::QualityMetric::TCount => c.tcount as f64,
                search::QualityMetric::TDepth => c.tdepth as f64,
                search::QualityMetric::GateCount => c.gatecount as f64,
                search::QualityMetric::WeightedCost => c.cost_required,
            };
            if goal < 0.0 {
                return Err("quality-objective requires the corresponding resource limit".into());
            }
            c.search.quality_metric = metric;
            c.search.quality_goal = Some(goal);
        }
        if c.clifford_macros
            || c.rccx_macros
            || c.permutation_restarts > 0
            || c.phase_fold
            || c.pauli_fold
        {
            c.expand = true;
            c.search.improved_post = true;
        }
        c.search.quality_expand = c.expand;
        c.validate()?;
        Ok(c)
    }

    fn validate(&self) -> Result<(), String> {
        if self.general_only
            && (self.clifford_macros
                || self.rccx_macros
                || self.affine_search
                || self.phase_factor
                || self.clifford_factor
                || self.phase_restarts > 0
                || self.nonlinear_phase_restarts > 0
                || self.single_qubit_restarts > 0
                || self.permutation_restarts > 0
                || self.pauli_restarts > 0
                || self.phase_fold
                || self.pauli_fold
                || self.clifford_fold)
        {
            return Err("general search excludes structural kernels, automatic macros, and Clifford/phase folding; select the gate library explicitly".into());
        }
        if self.threads == 0 {
            return Err("--threads must be positive".into());
        }
        for (name, value) in [
            ("time", self.time),
            ("structured-time-slice", self.structured_time_slice),
            ("legacy-fallback-after", self.legacy_fallback_after),
            ("resource-fallback-after", self.resource_fallback_after),
            ("quality-weight", self.search.quality_weight),
            ("quality-slack", self.search.quality_slack.unwrap_or(0.0)),
            ("pcomp", self.search.pcomp),
            ("epsilon", self.search.epsilon),
            ("start-temp", self.search.start_temp_base),
            ("n-norm", self.search.n_norm),
            ("iterations-factor", self.search.iterations_factor),
            ("beta", self.scheme.beta),
            ("fmin", self.scheme.min_factor),
            ("fmax", self.scheme.max_factor),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(format!("--{name} must be finite and nonnegative"));
            }
        }
        if self.resource_fallback_after > 0.0 && self.cost_required < 0.0 {
            return Err("--resource-fallback-after requires --cost-required".into());
        }
        if !(1..=64).contains(&self.search.masked_proposals) {
            return Err("--masked-proposals must be in 1..=64".into());
        }
        if self.search.n_norm == 0.0 || self.search.start_temp_base == 0.0 {
            return Err("temperature and n-norm must be positive".into());
        }
        if !(0.0..=1.0).contains(&self.search.swap_probability) {
            return Err("--swap-probability must be in [0,1]".into());
        }
        if !(0.0..=1.0).contains(&self.search.pid) {
            return Err("--pid must be in [0,1]".into());
        }
        if !self.search.cost_power.is_finite() || self.search.cost_power <= 0.0 {
            return Err("cost-power must be finite and positive".into());
        }
        if !(0.0..=1.0).contains(&self.search.sweep_cooling) || self.search.sweep_cooling == 0.0 {
            return Err("sweep-cooling must be in (0,1]".into());
        }
        if !self.search.gate_prior.is_finite() || self.search.gate_prior < 0.0 {
            return Err("gate-prior must be finite and nonnegative".into());
        }
        if self.scheme.min_factor > self.scheme.max_factor {
            return Err("fmin must not exceed fmax".into());
        }
        if !self.cost_required.is_finite() {
            return Err("cost-required must be finite".into());
        }
        if self.ancillas.saturating_add(self.dirty) > 10 {
            return Err("too many auxiliary qubits for dense synthesis".into());
        }
        if !["qasm", "txt"].contains(&self.format.as_str()) {
            return Err("format must be txt or qasm".into());
        }
        Ok(())
    }

    pub fn print(&self) {
        macro_rules! p {
            ($name:expr,$value:expr) => {
                println!("{}: {}", $name, $value);
            };
        }
        macro_rules! b {
            ($name:expr,$value:expr) => {
                p!($name, u8::from($value));
            };
        }
        p!(
            "Search mode",
            if self.general_only {
                "general (simulated annealing, supplied gate library only)"
            } else if self.auto_search {
                "auto (simulated annealing)"
            } else if self.search.optimized {
                "sweep (legacy fallback for unsupported covers)"
            } else {
                "legacy"
            }
        );
        p!("Input file", self.input);
        p!("Output folder", self.output);
        p!("Format", self.format);
        p!("Gate set", self.basic);
        p!("Composite gate folder", self.composite);
        p!("Number of threads", self.threads);
        p!("Number of circuits to find", self.circuits);
        p!("Number of ancillas", self.ancillas);
        p!("Number of dirty qubits", self.dirty);
        p!("Optimal t-count", self.tcount);
        p!("Optimal t-depth", self.tdepth);
        p!("Optimal gate count", self.gatecount);
        p!("Cost required", self.cost_required);
        p!("Time allowed", self.time);
        p!("Start temp base", self.search.start_temp_base);
        p!("Epsilon", self.search.epsilon);
        p!("Beta", self.scheme.beta);
        p!("Min factor", self.scheme.min_factor);
        p!("Max factor", self.scheme.max_factor);
        p!("Pcomp", self.search.pcomp);
        p!("Pid", self.search.pid);
        b!("Enable permutations", self.search.enable_permutations);
        b!("Do resynth", self.resynth);
        b!("Clifford basis preprocessing", self.clifford_factor);
        b!("Phase folding", self.phase_fold);
        b!("Pauli folding", self.pauli_fold);
        b!("Simple cost", self.search.simple_cost);
        p!("N norm", self.search.n_norm);
        p!("Iterations factor", self.search.iterations_factor);
        b!("Update gate scheme", self.update_scheme);
        b!("Save times", self.save_times);
        b!("Expand composite", self.expand);
        b!("Qubit independent", self.independent);
        b!("Save all circuits", self.save_all);
        b!("Save any circuit", self.save_circuits);
        p!("Times file", self.times_file);
        p!("Optimization number", self.optimization);
        b!("Optimize depth", self.optimize_depth);
        b!("Inverse independent", self.inverse);
        p!("Depth gates", format!("{} ", self.depth_gates.join(" ")));
    }

    fn unspecified(&self, flags: &[&str]) -> bool {
        !self
            .specified
            .iter()
            .any(|arg| flags.contains(&arg.as_str()))
    }

    /// Settings use operator/gate properties only, never an input path or name.
    fn resolve_search(&self, target: &Target, lib: &GateLibrary, unitary: bool) -> Self {
        let mut result = self.clone();
        if !self.auto_search || !self.search.optimized || !unitary {
            return result;
        }
        let full = target.original.cover.iter().all(|&b| b);
        let arbitrary_mask = crate::anneal::constant_norm(&target.original).is_none();
        let exact = self.search.epsilon <= 1e-6;
        if self.unspecified(&["--general-every"]) {
            result.general_every = 4;
        }
        if self.unspecified(&["--structured-time-slice"]) {
            result.structured_time_slice = 0.02;
        }
        // Opaque user composites keep their original resource accounting.
        let native = lib.n_qubits <= 6
            && (self.expand || lib.composite.is_empty())
            && crate::phase::NativeClifford::new(lib).is_some();
        if exact {
            if self.unspecified(&["--quality-slack"]) {
                result.search.quality_slack = Some(0.0);
            }
            if self.unspecified(&["--swap-probability"]) {
                result.search.swap_probability = 0.25;
            }
        }
        if arbitrary_mask {
            if self.unspecified(&["--masked-sweep"]) {
                result.search.masked_sweep = true;
            }
            if self.unspecified(&["--masked-proposals"]) {
                result.search.masked_proposals = 4;
            }
            if self.unspecified(&["--start-temp", "-st"]) {
                result.search.start_temp_base = 0.3;
            }
        }
        if native
            && exact
            && self.resynth
            && self.tcount > 0
            && self.gatecount >= 0
            && self.cost_required >= 0.0
            && self.unspecified(&[
                "--resource-fallback-after",
                "--quality-objective",
                "--quality-weight",
                "--quality-slack",
            ])
        {
            result.resource_fallback_after = 0.3;
        }
        // General mode retains only representation-independent annealing settings.
        // Return before any structural recognition, macro injection, or new rewrite pass.
        if self.general_only {
            result.search.quality_expand = result.expand;
            return result;
        }
        // Budget-free defaults do not add gates or invoke specialized synthesis.
        // Budgeted kernels and permutation/RCCX methods require explicit flags.
        result.clifford_fold |= native && self.resynth;
        if self.unspecified(&["--clifford-factor", "--no-clifford-factor"]) {
            result.clifford_factor = native && full && lib.n_qubits <= 4 && !self.inverse;
        }
        if native && self.resynth && lib.n_qubits <= 4 {
            if self.unspecified(&["--phase-fold", "--no-phase-fold"]) {
                result.phase_fold = true;
            }
            if self.unspecified(&["--pauli-fold", "--no-pauli-fold"]) {
                result.pauli_fold = true;
            }
        }
        if result.clifford_macros || result.rccx_macros || result.permutation_restarts > 0 {
            result.expand = true;
            result.search.improved_post = true;
        }
        result.search.quality_expand = result.expand;
        result
    }

    fn accepts_quality(&self, circuit: &Circuit, lib: &GateLibrary) -> bool {
        let names = ["t".to_owned(), "tdg".to_owned()];
        (self.tcount < 0 || circuit.count(&names, lib) as i64 <= self.tcount)
            && (self.tdepth < 0 || circuit.depth(&names, lib) as i64 <= self.tdepth)
            && (self.gatecount < 0 || circuit.non_identity() as i64 <= self.gatecount)
            && (self.cost_required < 0.0 || circuit.cost() <= self.cost_required)
    }

    /// Preserve a feasible native candidate across rewrites with other objectives.
    fn feasible_backup(&self, circuit: &Circuit, lib: &GateLibrary) -> Option<Circuit> {
        if !self.search.optimized {
            return None;
        }
        let mut native = circuit.clone();
        if self.expand {
            native.expand(lib);
        }
        self.accepts_quality(&native, lib).then_some(native)
    }

    pub fn load_circuit(&self) -> Result<(GateLibrary, Circuit), String> {
        let text = fs::read_to_string(&self.input).map_err(|e| format!("{}: {e}", self.input))?;
        let n = qasm_qubits(&text)?;
        let mut lib = GateLibrary::load(n, &self.basic, &self.composite)?;
        lib.add_read_folder("data/gates/read_gates")?;
        let circuit = Circuit::from_qasm(&text, &lib)?;
        Ok((lib, circuit))
    }

    pub fn load_target(&self) -> Result<Target, String> {
        let original = if self.format == "qasm" {
            let (_, circ) = self.load_circuit()?;
            PartialMatrix::new(
                circ.matrix().clone(),
                vec![true; circ.matrix().data.len()],
                "matrix",
            )?
        } else {
            PartialMatrix::read(&self.input)?
        };
        let total = original
            .n_qubits
            .saturating_add(self.ancillas)
            .saturating_add(self.dirty);
        if total > 10 {
            return Err(
                "dense synthesis is limited to 10 qubits to prevent unbounded allocation".into(),
            );
        }
        let original = original
            .add_ancillas(self.ancillas)?
            .add_dirty(self.dirty)?;
        if original.n_constraints() == 0 {
            return Err("specification has no covered entries".into());
        }
        let independent = self.independent && !self.original_sweep(&original);
        Ok(Target::new(original, independent, self.inverse))
    }

    fn original_sweep(&self, original: &PartialMatrix) -> bool {
        self.search.optimized
            && !self.search.sweep_permutations
            && !self.inverse
            && (self.auto_search
                || self.search.masked_sweep
                || crate::anneal::constant_norm(original).is_some())
    }
}

pub fn qasm_qubits(text: &str) -> Result<usize, String> {
    for line in text.lines() {
        let line = line.split("//").next().unwrap().trim();
        for statement in line.split(';') {
            let s = statement.trim();
            if s.starts_with("qreg ") {
                let n = s
                    .split_once('[')
                    .and_then(|(_, v)| v.split_once(']'))
                    .ok_or("invalid qreg declaration")?
                    .0;
                return n.parse().map_err(|_| "invalid qubit count".into());
            }
        }
    }
    Err("missing qreg declaration".into())
}

#[derive(Default, Debug)]
pub struct RunStats {
    pub runs: usize,
    pub successful: usize,
    pub optimal: usize,
    pub found: usize,
    pub times: Vec<f64>,
    pub best_tcount: Option<f64>,
    pub best_tdepth: Option<f64>,
    pub best_gatecount: Option<f64>,
    pub best_cost: Option<f64>,
}
struct Shared {
    stats: RunStats,
    scheme: GateScheme,
    stop: bool,
    last: Instant,
}

pub fn run(config: &Config) -> Result<RunStats, String> {
    run_target(config, config.load_target()?)
}
pub fn run_cli(config: &Config) -> Result<RunStats, String> {
    run_target_inner(config, config.load_target()?, true)
}
pub fn run_target(config: &Config, target: Target) -> Result<RunStats, String> {
    run_target_inner(config, target, false)
}
fn run_target_inner(
    config: &Config,
    target: Target,
    print_config: bool,
) -> Result<RunStats, String> {
    let start = Instant::now();
    let deadline = std::time::Duration::try_from_secs_f64(config.time)
        .ok()
        .and_then(|duration| start.checked_add(duration));
    fs::create_dir_all(&config.output).map_err(|e| format!("{}: {e}", config.output))?;
    let mut lib = GateLibrary::load(target.original.n_qubits, &config.basic, &config.composite)?;
    let needs_unitary = config.auto_search
        || config.search.masked_sweep
        || target
            .variants
            .iter()
            .all(|p| crate::anneal::constant_norm(p).is_some());
    let basic_unitary = config.search.optimized
        && needs_unitary
        && crate::anneal::SweepKernel::supports_masked(&target, &lib);
    let mut resolved = config.resolve_search(&target, &lib, basic_unitary);
    resolved.search.cancellation = (resolved.threads > 1)
        .then(|| std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false)));
    let config = &resolved;
    if print_config {
        config.print();
    }

    let fallback_library = (config.auto_search
        && (config.gatecount >= 0 || config.cost_required >= 0.0)
        && (config.clifford_macros || config.rccx_macros))
        .then(|| lib.clone());
    if config.clifford_macros {
        lib.add_clifford_macros();
    }
    if config.rccx_macros && crate::phase::NativeClifford::new(&lib).is_some() {
        lib.add_composite_folder("data/gates/composite_rccx")?;
    }
    let original = target.original.clone();
    let affine = if config.affine_search
        && config.search.optimized
        && (config.tcount >= 0 || config.tdepth >= 0)
        && (config.bounded_structure || (config.gatecount < 0 && config.cost_required < 0.0))
        && original.cover.iter().all(|&v| v)
        && !target.inverses.iter().any(|&v| v)
    {
        crate::phase::NativeClifford::new(&lib).and_then(|native| {
            crate::affine::reduce_full(&original.matrix, &lib, &native.native, &native.cx)
        })
    } else {
        None
    };
    let (target, restore) = if let Some(affine) = affine {
        let prefix =
            crate::affine::inverse_word(&affine.right, &lib).ok_or("invalid affine prefix")?;
        let suffix =
            crate::affine::inverse_word(&affine.left, &lib).ok_or("invalid affine suffix")?;
        let spec = PartialMatrix::new(affine.matrix, original.cover.clone(), &original.name)?;
        (Target::new(spec, false, false), Some((prefix, suffix)))
    } else {
        (target, None)
    };
    let basis = if restore.is_none()
        && config.clifford_factor
        && target.original.cover.iter().all(|&b| b)
        && !target.inverses.iter().any(|&b| b)
    {
        crate::clifford_diagonal::reduce(&target.original.matrix, &lib)
    } else {
        None
    };
    let basis_factored = basis.is_some();
    let (target, restore) = if let Some(basis) = basis {
        let spec = PartialMatrix::new(basis.matrix, original.cover.clone(), &original.name)?;
        (
            Target::new(spec, false, false),
            Some((basis.prefix, basis.suffix)),
        )
    } else {
        (target, restore)
    };
    let phase_kernel = if config.search.optimized
        && config.phase_restarts > 0
        && (config.tcount >= 0 || config.tdepth >= 0)
        && (config.bounded_structure || (config.gatecount < 0 && config.cost_required < 0.0))
    {
        let t_limit = if config.tcount >= 0 {
            config.tcount as usize
        } else {
            lib.n_qubits.saturating_mul(config.tdepth as usize)
        };
        let phase = if config.phase_factor {
            crate::phase::PhaseKernel::new_factored(&target, &lib, t_limit)
        } else {
            crate::phase::PhaseKernel::new(&target, &lib, t_limit)
        };
        phase.map(|phase| {
            if config.tdepth >= 0 {
                phase.with_depth_limit(config.tdepth as usize)
            } else {
                phase
            }
        })
    } else {
        None
    };
    let nonlinear_kernel = if config.search.optimized
        && config.nonlinear_phase_restarts > 0
        && config.tcount >= 9
        && (config.bounded_structure || (config.gatecount < 0 && config.cost_required < 0.0))
    {
        crate::phase::PhaseKernel::new_nonlinear(&target, &lib, config.tcount as usize)
    } else {
        None
    };
    let single_kernel = if config.search.optimized
        && config.single_qubit_restarts > 0
        && (config.tcount >= 0 || config.tdepth >= 0)
        && (config.bounded_structure || (config.gatecount < 0 && config.cost_required < 0.0))
    {
        let limit = if config.tcount >= 0 {
            config.tcount as usize
        } else {
            config.tdepth as usize
        };
        crate::single_qubit::SingleQubitKernel::new(&target, &lib, limit)
    } else {
        None
    };
    let permutation_kernel = if config.search.optimized
        && config.permutation_restarts > 0
        && (config.tcount >= 4 || config.tdepth >= 0)
        && (config.bounded_structure || (config.gatecount < 0 && config.cost_required < 0.0))
    {
        let limit = if config.tcount >= 0 {
            config.tcount as usize
        } else {
            lib.n_qubits.saturating_mul(config.tdepth as usize)
        };
        crate::permutation::PermutationKernel::new(&target, &lib, limit).map(|kernel| {
            if config.tcount < 0 {
                kernel.vary_t_budget()
            } else {
                kernel
            }
        })
    } else {
        None
    };
    let pauli_kernel = if config.search.optimized
        && config.pauli_restarts > 0
        && config.tcount >= 0
        && (config.bounded_structure || (config.gatecount < 0 && config.cost_required < 0.0))
        && phase_kernel.is_none()
        && nonlinear_kernel.is_none()
        && single_kernel.is_none()
        && permutation_kernel.is_none()
    {
        crate::pauli_search::PauliKernel::new(&target, &lib, config.tcount as usize)
    } else {
        None
    };
    let unitary = basic_unitary
        && (!(config.clifford_macros || config.rccx_macros)
            || crate::anneal::SweepKernel::supports_masked(&target, &lib));
    let constant = target
        .variants
        .iter()
        .all(|p| crate::anneal::constant_norm(p).is_some());
    let kernel = (config.search.optimized && unitary && (constant || config.search.masked_sweep))
        .then(|| crate::anneal::SweepKernel::new(&lib, &config.search));
    // A nonunitary gate library needs the general objective and legacy variants.
    let target =
        if kernel.is_none() && config.independent && config.original_sweep(&target.original) {
            Target::new(target.original, true, config.inverse)
        } else {
            target
        };
    let fallback_kernel = fallback_library.as_ref().and_then(|library| {
        (config.search.optimized && basic_unitary && (constant || config.search.masked_sweep))
            .then(|| crate::anneal::SweepKernel::new(library, &config.search))
    });
    let original_fallback = (basis_factored
        || ((restore.is_some() || fallback_library.is_some())
            && (config.gatecount >= 0 || config.cost_required >= 0.0)))
        .then(|| {
            Target::new(
                original.clone(),
                kernel.is_none() && config.independent,
                config.inverse,
            )
        });
    let legacy_target = OnceLock::new();
    let resource_options = (config.resource_fallback_after > 0.0 && config.cost_required >= 0.0)
        .then(|| {
            let mut options = config.search.clone();
            options.quality_metric = search::QualityMetric::WeightedCost;
            options.quality_goal = Some(config.cost_required);
            options
        });
    let resource_kernel = OnceLock::new();
    let shared = Mutex::new(Shared {
        stats: RunStats::default(),
        scheme: config.scheme.clone(),
        stop: false,
        last: Instant::now(),
    });
    let mut round = 0u64;
    while start.elapsed().as_secs_f64() < config.time
        && shared.lock().unwrap().stats.found < config.circuits
    {
        {
            let mut state = shared.lock().unwrap();
            state.stop = false;
            if let Some(flag) = &config.search.cancellation {
                flag.store(false, std::sync::atomic::Ordering::Relaxed);
            }
            state.last = Instant::now();
            state.scheme.reset();
        }
        std::thread::scope(|scope| -> Result<(), String> {
            let mut workers = Vec::new();
            for id in 0..config.threads {
                let shared = &shared;
                let lib = &lib;
                let target = &target;
                let kernel = &kernel;
                let phase_kernel = &phase_kernel;
                let permutation_kernel = &permutation_kernel;
                let single_kernel = &single_kernel;
                let nonlinear_kernel = &nonlinear_kernel;
                let pauli_kernel = &pauli_kernel;
                let original = &original;
                let restore = &restore;
                let original_fallback = &original_fallback;
                let fallback_library = &fallback_library;
                let fallback_kernel = &fallback_kernel;
                let legacy_target = &legacy_target;
                let resource_options = &resource_options;
                let resource_kernel = &resource_kernel;
                workers.push(scope.spawn(move || -> Result<(), String> {
                    let mut rng = Rng::new(
                        config
                            .seed
                            .wrapping_add(id as u64)
                            .wrapping_add(round.wrapping_mul(config.threads as u64)),
                    );
                    let resynth = Resynth::new(
                        config.optimization,
                        config.optimize_depth,
                        vec!["t".into(), "tdg".into()],
                    );
                    let tnames = vec!["t".into(), "tdg".into()];
                    let mut legacy_options = config.search.clone();
                    legacy_options.optimized = false;
                    let mut chain = 0usize;
                    let mut general_chain = 0usize;
                    loop {
                        {
                            let mut state = shared.lock().unwrap();
                            if state.stop
                                || state.stats.found >= config.circuits
                                || start.elapsed().as_secs_f64() >= config.time
                            {
                                break;
                            }
                            state.stats.runs += 1;
                            chain += 1;
                        }
                        let mut phase_solution = None;
                        // Reserve general annealing attempts even when a specialized
                        // kernel repeatedly finds equivalents that miss another goal.
                        if config.general_every == 0 || !chain.is_multiple_of(config.general_every)
                        {
                            let specialized_deadline = if config.structured_time_slice > 0.0 {
                                Instant::now()
                                    .checked_add(std::time::Duration::from_secs_f64(
                                        config.structured_time_slice,
                                    ))
                                    .map(|slice| deadline.map_or(slice, |end| end.min(slice)))
                            } else {
                                deadline
                            };
                            if let Some(phase) = phase_kernel {
                                for _ in 0..config.phase_restarts {
                                    let result = phase.run(target, lib, &mut rng, &config.search);
                                    if result.found {
                                        phase_solution = Some(result);
                                        break;
                                    }
                                }
                            }
                            if phase_solution.is_none() {
                                if let Some(nonlinear) = nonlinear_kernel {
                                    for _ in 0..config.nonlinear_phase_restarts {
                                        if start.elapsed().as_secs_f64() >= config.time {
                                            break;
                                        }
                                        let result =
                                            nonlinear.run(target, lib, &mut rng, &config.search);
                                        if result.found {
                                            phase_solution = Some(result);
                                            break;
                                        }
                                    }
                                }
                            }
                            if phase_solution.is_none() {
                                if let Some(single) = single_kernel {
                                    for _ in 0..config.single_qubit_restarts {
                                        if start.elapsed().as_secs_f64() >= config.time {
                                            break;
                                        }
                                        let result =
                                            single.run(target, lib, &mut rng, &config.search);
                                        if result.found {
                                            phase_solution = Some(result);
                                            break;
                                        }
                                    }
                                }
                            }
                            if phase_solution.is_none() {
                                if let Some(permutation) = permutation_kernel {
                                    for _ in 0..config.permutation_restarts {
                                        if start.elapsed().as_secs_f64() >= config.time {
                                            break;
                                        }
                                        let result = permutation.run_until(
                                            target,
                                            lib,
                                            &mut rng,
                                            &config.search,
                                            specialized_deadline,
                                        );
                                        if result.found {
                                            phase_solution = Some(result);
                                            break;
                                        }
                                    }
                                }
                            }
                            if phase_solution.is_none() {
                                if let Some(pauli) = pauli_kernel {
                                    let pauli_deadline =
                                        if config.unspecified(&["--structured-time-slice"]) {
                                            Instant::now()
                                                .checked_add(std::time::Duration::from_secs(1))
                                                .map(|slice| {
                                                    deadline.map_or(slice, |end| end.min(slice))
                                                })
                                        } else {
                                            specialized_deadline
                                        };
                                    for _ in 0..config.pauli_restarts {
                                        let result = pauli.run_until(
                                            target,
                                            lib,
                                            &mut rng,
                                            &config.search,
                                            pauli_deadline,
                                        );
                                        if result.found {
                                            phase_solution = Some(result);
                                            break;
                                        }
                                        if pauli_deadline.is_some_and(|end| Instant::now() >= end) {
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                        if start.elapsed().as_secs_f64() >= config.time {
                            break;
                        }
                        // Count only general attempts when scheduling fallback searches.
                        let use_legacy = phase_solution.is_none()
                            && config.search.optimized
                            && config.legacy_fallback_after > 0.0
                            && (general_chain + 1).is_multiple_of(3)
                            && start.elapsed().as_secs_f64() >= config.legacy_fallback_after;
                        let use_resource = phase_solution.is_none()
                            && !use_legacy
                            && resource_options.is_some()
                            && (general_chain + 1).is_multiple_of(3)
                            && start.elapsed().as_secs_f64() >= config.resource_fallback_after;
                        let n_gates = if phase_solution.is_none() {
                            general_chain += 1;
                            shared.lock().unwrap().scheme.get_start_gates(&mut rng)
                        } else {
                            0
                        };
                        let use_original = use_legacy
                            || (phase_solution.is_none()
                                && original_fallback.is_some()
                                && (!basis_factored || general_chain.is_multiple_of(4)));
                        let search_target = if use_legacy {
                            legacy_target.get_or_init(|| {
                                Target::new(
                                    original.clone(),
                                    config.independent && original.n_qubits <= 4,
                                    config.inverse,
                                )
                            })
                        } else if use_original {
                            original_fallback.as_ref().unwrap()
                        } else {
                            target
                        };
                        let active_restore = if use_original { None } else { restore.as_ref() };
                        let general_lib = fallback_library.as_ref().unwrap_or(lib);
                        let mut general_kernel = if fallback_library.is_some() {
                            fallback_kernel.as_ref()
                        } else {
                            kernel.as_ref()
                        };
                        let search_options = if use_resource {
                            resource_options.as_ref().unwrap()
                        } else {
                            &config.search
                        };
                        if use_resource && general_kernel.is_some() {
                            general_kernel = Some(resource_kernel.get_or_init(|| {
                                crate::anneal::SweepKernel::new(general_lib, search_options)
                            }));
                        }
                        let result = if let Some(solution) = phase_solution {
                            solution
                        } else if use_legacy {
                            search::search(
                                search_target,
                                general_lib,
                                &mut rng,
                                n_gates,
                                &legacy_options,
                            )?
                        } else if let Some(kernel) = general_kernel {
                            kernel.run(
                                search_target,
                                general_lib,
                                &mut rng,
                                n_gates,
                                search_options,
                            )?
                        } else {
                            search::search(
                                search_target,
                                general_lib,
                                &mut rng,
                                n_gates,
                                search_options,
                            )?
                        };
                        if config.search.is_cancelled() {
                            break;
                        }
                        let mut best = result.circuit;
                        if search_target.exact_cost(best.matrix(), config.search.epsilon) != 0.0 {
                            continue;
                        }
                        let core_count = best.non_identity();
                        if let Some((prefix, suffix)) = active_restore {
                            let mut gates = prefix.clone();
                            gates.extend(&best.gates);
                            gates.extend(suffix);
                            best = crate::circuit::Circuit::new(gates, lib.n_qubits, lib);
                        }
                        let before = (
                            best.non_identity(),
                            best.count(&tnames, lib),
                            best.depth(&tnames, lib),
                        );
                        let mut feasible = config.feasible_backup(&best, lib);
                        if config.resynth {
                            if config.search.improved_post {
                                resynth.run_optimized(&mut best, lib);
                            } else {
                                resynth.run(&mut best, lib);
                            }
                        }
                        let scheme_count = if active_restore.is_some() {
                            core_count
                        } else {
                            best.non_identity()
                        };
                        if config.expand {
                            best.expand(lib);
                            if config.resynth {
                                if config.search.improved_post {
                                    resynth.run_optimized(&mut best, lib);
                                } else {
                                    resynth.run(&mut best, lib);
                                }
                            }
                        }

                        // Rewriting may create an internal Clifford macro after
                        // the first expansion. Internal names must never escape
                        // into saved QASM or native-gate quality accounting.
                        if config.clifford_macros
                            || config.rccx_macros
                            || config.permutation_restarts > 0
                        {
                            best.expand(lib);
                        }
                        if let Some(candidate) = config.feasible_backup(&best, lib) {
                            feasible = Some(candidate);
                        }
                        if config.phase_fold {
                            let metric = if config.tdepth >= 0 {
                                search::QualityMetric::TDepth
                            } else {
                                config.search.quality_metric
                            };
                            crate::phase_fold::run(&mut best, lib, &mut rng, metric);
                            if let Some(candidate) = config.feasible_backup(&best, lib) {
                                feasible = Some(candidate);
                            }
                        }
                        if config.pauli_fold {
                            crate::pauli_fold::run_with_metric(
                                &mut best,
                                lib,
                                &mut rng,
                                config.search.quality_metric,
                                (config.tdepth >= 0).then_some(config.tdepth as usize),
                            );
                        }
                        if config.clifford_fold {
                            crate::linear_clifford::fold(&mut best, lib);
                        }
                        if !config.accepts_quality(&best, lib) {
                            if let Some(candidate) = feasible {
                                best = candidate;
                            }
                        }
                        if config.search.optimized
                            && original.exact_cost(best.matrix(), config.search.epsilon) != 0.0
                        {
                            continue;
                        }
                        let after = (
                            best.non_identity(),
                            best.count(&tnames, lib),
                            best.depth(&tnames, lib),
                        );
                        let mut state = shared.lock().unwrap();
                        state.stats.successful += 1;
                        state.stats.best_tcount = Some(
                            state
                                .stats
                                .best_tcount
                                .unwrap_or(f64::INFINITY)
                                .min(after.1 as f64),
                        );
                        state.stats.best_tdepth = Some(
                            state
                                .stats
                                .best_tdepth
                                .unwrap_or(f64::INFINITY)
                                .min(after.2 as f64),
                        );
                        state.stats.best_gatecount = Some(
                            state
                                .stats
                                .best_gatecount
                                .unwrap_or(f64::INFINITY)
                                .min(after.0 as f64),
                        );
                        state.stats.best_cost = Some(
                            state
                                .stats
                                .best_cost
                                .unwrap_or(f64::INFINITY)
                                .min(best.cost()),
                        );
                        if config.update_scheme {
                            state.scheme.update(scheme_count);
                        }
                        let has_goal = config.tcount >= 0
                            || config.tdepth >= 0
                            || config.gatecount >= 0
                            || config.cost_required >= 0.0;
                        let meets = config.accepts_quality(&best, lib);
                        let save = meets && !state.stop && state.stats.found < config.circuits;
                        if save {
                            if has_goal {
                                state.stop = true;
                                state.stats.optimal += 1;
                            }
                            let elapsed = state.last.elapsed().as_secs_f64();
                            state.stats.times.push(elapsed);
                            state.stats.found += 1;
                            if state.stop || state.stats.found >= config.circuits {
                                if let Some(flag) = &config.search.cancellation {
                                    flag.store(true, std::sync::atomic::Ordering::Relaxed);
                                }
                            }
                            println!(
                                "{} {} {}\nGates: {} -> {}\nT-count: {} -> {}\nT-depth: {} -> {}",
                                config.output,
                                elapsed,
                                state.stats.found,
                                before.0,
                                after.0,
                                before.1,
                                after.1,
                                before.2,
                                after.2
                            );
                            state.last = Instant::now();
                        }
                        if (save || config.save_all) && config.save_circuits {
                            let file = format!(
                                "{}{:.6}-{}-{}-{}-{}.qasm",
                                config.output,
                                best.cost(),
                                best.count(&config.depth_gates, lib),
                                best.depth(&config.depth_gates, lib),
                                id,
                                state.stats.found
                            );
                            fs::write(&file, best.qasm(lib)).map_err(|e| format!("{file}: {e}"))?;
                        }
                    }
                    Ok(())
                }));
            }
            let mut error = None;
            for worker in workers {
                match worker.join() {
                    Ok(Ok(())) => {}
                    Ok(Err(e)) => {
                        error.get_or_insert(e);
                        shared.lock().unwrap().stop = true;
                    }
                    Err(_) => {
                        error.get_or_insert("synthesis worker panicked".into());
                        shared.lock().unwrap().stop = true;
                    }
                }
            }
            error.map_or(Ok(()), Err)
        })?;
        round += 1;
    }
    let stats = shared.into_inner().unwrap().stats;
    if config.save_times {
        let mean = stats.times.iter().sum::<f64>() / stats.times.len() as f64;
        let variance =
            stats.times.iter().map(|v| v * v).sum::<f64>() / stats.times.len() as f64 - mean * mean;
        let mut file = fs::OpenOptions::new()
            .create(true)
            .append(true)
            .open(&config.times_file)
            .map_err(|e| format!("{}: {e}", config.times_file))?;
        let stdev = if stats.times.is_empty() {
            f64::NAN
        } else {
            variance.max(0.0).sqrt()
        };
        writeln!(
            file,
            "{}\t{}\t{}\t{}\t{}\t{}",
            config.output, mean, stdev, stats.runs, stats.successful, stats.optimal
        )
        .map_err(|e| e.to_string())?;
    }
    Ok(stats)
}

pub const HELP:&str="Synthetiq — Rust quantum circuit synthesis\nUsage: synthetiq INPUT [OPTIONS]\nInputs default to data/input; gates default to data/gates.\n  --threads, -h N             parallel synthesis workers (default 1)\n  --time, -t SECONDS          time budget (default 100)\n  --circuits, -c N            stop after N accepted circuits (default 10)\n  --output, -o DIR            output folder under data/output\n  --gate-set, -gs DIR         finite basic gate set (default CliffordT)\n  --composite-gates, -cg DIR  composite gates (default composite_gates)\n  --gates, -g N               initial gate-count estimate\n  --n-start-gates N           fixed initial length (legacy N+1)\n  --tcount/-tc, --tdepth/-td, --gate-count/-gc, --cost-required/-cr N\n  --ancilla/-a N, --dirty/-d N, --epsilon/-eps E, --expand/-e\n  --start-temp/-st T, --beta B, --fmin F, --fmax F, --pid P, --pcomp P\n  --n-norm N, --iterations-factor N, --optimization-number N\n  --no-perms, --no-resynth, --simple, --gs-no-update, --no-optimize-depth\n  --qubit_independence/-q     disable qubit-independent target variants\n  --inverse-independent, --depth-gates t,tdg\n  --save, --times-file FILE, --save-all/-sa, --no-save-circuits\n  --absolute-input, --absolute-output, --absolute-gates\n  --search-mode legacy|sweep|general|auto (default auto), --cost-power P, --gate-prior W\n  --sweep-cooling RATE, --stall-sweeps N, --sweep-original, --sweep-permutations\n  --no-cache-rollback, --improved-post\n  --quality-weight W, --quality-slack N, --clifford-macros\n  --quality-objective tcount|tdepth|gatecount|cost\n  --swap-probability P, --phase-gadgets, --phase-restarts N, --affine-search, --rccx-macros\n  --masked-sweep, --masked-proposals N\n  --permutation-gadgets, --permutation-restarts N\n  --phase-fold, --phase-factor, --clifford-factor\n  --no-phase-fold, --no-pauli-fold, --no-clifford-factor\n  --single-qubit-restarts N, --nonlinear-phase-restarts N, --pauli-fold, --clifford-fold, --pauli-restarts N\n  --general-every N, --structured-time-slice SECONDS, --bounded-structure, --legacy-fallback-after SECONDS\n  --resource-fallback-after SECONDS (requires --cost-required)\n  --format/-f txt|qasm, --seed N, --help\n";

#[cfg(test)]
mod tests {
    use super::*;
    fn parse(s: &str) -> Result<Config, String> {
        Config::parse(&s.split_whitespace().map(String::from).collect::<Vec<_>>())
    }
    #[test]
    fn automatic_budget_free_methods_never_enable_extra_gates_or_budgeted_search() {
        for input in [
            "64/comparison/U1.txt",
            "64/comparison/ccx.txt",
            "65/rz_4.txt",
        ] {
            let original = PartialMatrix::read(format!("data/input/{input}")).unwrap();
            let lib = GateLibrary::load(original.n_qubits, "data/gates/CliffordT", "").unwrap();
            let target = Target::new(original, false, false);
            for goal in ["", "--tcount 12", "--tdepth 5"] {
                let c = parse(&format!("{input} {goal}"))
                    .unwrap()
                    .resolve_search(&target, &lib, true);
                assert!(c.phase_fold && c.pauli_fold && c.clifford_fold && c.clifford_factor);
                assert!(
                    !c.rccx_macros && !c.clifford_macros && !c.affine_search && !c.phase_factor
                );
                assert_eq!(
                    (
                        c.phase_restarts,
                        c.pauli_restarts,
                        c.permutation_restarts,
                        c.nonlinear_phase_restarts,
                        c.single_qubit_restarts
                    ),
                    (0, 0, 0, 0, 0)
                );
                assert!(!c.expand && !c.search.improved_post);
            }
            let c = parse(&format!(
                "{input} --no-phase-fold --no-pauli-fold --no-clifford-factor"
            ))
            .unwrap()
            .resolve_search(&target, &lib, true);
            assert!(!c.phase_fold && !c.pauli_fold && !c.clifford_factor);
            let c = parse(&format!("{input} --no-resynth"))
                .unwrap()
                .resolve_search(&target, &lib, true);
            assert!(!c.phase_fold && !c.pauli_fold && !c.clifford_fold);
        }
    }

    #[test]
    fn general_search_keeps_supplied_library_and_disables_every_structural_path() {
        for (input, goal) in [
            ("64/comparison/U1.txt", "--tcount 11"),
            ("64/comparison/ccx.txt", "--tdepth 3"),
            ("65/rz_4.txt", "--tcount 16 --epsilon .01"),
        ] {
            let original = PartialMatrix::read(format!("data/input/{input}")).unwrap();
            let lib = GateLibrary::load(original.n_qubits, "data/gates/CliffordT", "").unwrap();
            let target = Target::new(original, false, false);
            let c = parse(&format!("{input} --search-mode general {goal}")).unwrap();
            let resolved = c.resolve_search(&target, &lib, true);
            assert!(resolved.general_only && resolved.search.optimized);
            assert!(!resolved.clifford_macros && !resolved.rccx_macros);
            assert!(!resolved.affine_search && !resolved.phase_factor && !resolved.clifford_factor);
            assert!(!resolved.phase_fold && !resolved.pauli_fold && !resolved.clifford_fold);
            assert!(!resolved.expand && !resolved.search.improved_post);
            assert_eq!(
                (
                    resolved.phase_restarts,
                    resolved.nonlinear_phase_restarts,
                    resolved.single_qubit_restarts,
                    resolved.permutation_restarts,
                    resolved.pauli_restarts
                ),
                (0, 0, 0, 0, 0)
            );
            assert_eq!(resolved.basic, c.basic);
            assert_eq!(resolved.composite, c.composite);
        }
    }

    #[test]
    fn general_search_rejects_hidden_macros_but_allows_explicit_shared_libraries() {
        for flag in [
            "--clifford-macros",
            "--rccx-macros",
            "--affine-search",
            "--phase-factor",
            "--clifford-factor",
            "--phase-gadgets",
            "--permutation-gadgets",
            "--single-qubit-restarts 1",
            "--nonlinear-phase-restarts 1",
            "--pauli-restarts 1",
            "--phase-fold",
            "--pauli-fold",
            "--clifford-fold",
        ] {
            assert!(
                parse(&format!("cx.txt --search-mode general {flag}")).is_err(),
                "{flag}"
            );
        }
        let c = parse("cx.txt --search-mode general --composite-gates composite_rccx --expand")
            .unwrap();
        assert!(c.expand && !c.rccx_macros);
        assert!(c.composite.ends_with("composite_rccx"));
    }

    #[test]
    fn general_masked_search_retains_generic_updates_and_explicit_overrides() {
        let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
        let mut cover = vec![true; 16];
        cover[1] = false;
        let target = Target::new(
            PartialMatrix::new(crate::matrix::Matrix::identity(4), cover, "mask").unwrap(),
            false,
            false,
        );
        let c = parse("cx.txt --search-mode general --tcount 0 --start-temp .2 --gates 7").unwrap();
        let r = c.resolve_search(&target, &lib, true);
        assert!(r.search.masked_sweep && !r.clifford_fold);
        assert_eq!(r.search.masked_proposals, 4);
        assert_eq!(r.search.start_temp_base, 0.2);
        assert_eq!((r.scheme.min_gates, r.scheme.max_gates), (17, 24));
        assert!(!r.rccx_macros && !r.affine_search && r.phase_restarts == 0);
    }

    #[test]
    fn automatic_native_rewrites_preserve_opaque_user_composite_accounting() {
        let lib =
            GateLibrary::load(3, "data/gates/CliffordT", "data/gates/composite_rccx").unwrap();
        let target = Target::new(
            PartialMatrix::read("data/input/64/comparison/ccx.txt").unwrap(),
            false,
            false,
        );
        let c = parse("target.txt --search-mode auto --tcount 7").unwrap();
        let resolved = c.resolve_search(&target, &lib, true);
        assert!(!resolved.expand);
        assert!(!resolved.rccx_macros);
        assert!(!resolved.affine_search);
        assert_eq!(resolved.phase_restarts, 0);
        let expanded = parse("target.txt --search-mode auto --tcount 7 --expand")
            .unwrap()
            .resolve_search(&target, &lib, true);
        assert!(expanded.phase_fold && expanded.pauli_fold && expanded.clifford_factor);
        assert_eq!(expanded.phase_restarts, 0);
    }

    #[test]
    fn automatic_macro_additions_preserve_original_gate_indices() {
        let original = GateLibrary::load(3, "data/gates/CliffordT", "").unwrap();
        let mut extended = original.clone();
        extended.add_clifford_macros();
        extended
            .add_composite_folder("data/gates/composite_rccx")
            .unwrap();
        for (a, b) in original.gates.iter().zip(&extended.gates) {
            assert_eq!(
                (&a.name, &a.qubits, a.cost, &a.decomposition),
                (&b.name, &b.qubits, b.cost, &b.decomposition)
            );
            assert_eq!(a.matrix, b.matrix);
        }
        assert_eq!(original.identity, extended.identity);
    }

    #[test]
    fn all_modes_preserve_original_adaptive_lengths_for_full_and_masked_targets() {
        let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
        for mode in ["legacy", "sweep", "general", "auto"] {
            for masked in [false, true] {
                let mut cover = vec![true; 16];
                cover[1] = !masked;
                let target = Target::new(
                    PartialMatrix::new(crate::matrix::Matrix::identity(4), cover, "target")
                        .unwrap(),
                    false,
                    false,
                );
                for goal in ["", "--tcount 0", "--gate-count 2", "--epsilon .01"] {
                    let c = parse(&format!("target.txt --search-mode {mode} {goal}")).unwrap();
                    let mut resolved = c.resolve_search(&target, &lib, true);
                    let mut original = GateScheme::default();
                    let mut actual_rng = Rng::new(94);
                    let mut expected_rng = Rng::new(94);
                    for found in [None, Some(20), Some(4)] {
                        if let Some(n) = found {
                            resolved.scheme.update(n);
                            original.update(n);
                        }
                        for _ in 0..128 {
                            assert_eq!(
                                resolved.scheme.get_start_gates(&mut actual_rng),
                                original.get_start_gates(&mut expected_rng),
                                "{mode}, masked={masked}, goal={goal}, found={found:?}"
                            );
                        }
                    }
                    resolved.scheme.reset();
                    assert_eq!(
                        (resolved.scheme.min_gates, resolved.scheme.max_gates),
                        (30, 120)
                    );
                }
            }
        }
    }

    #[test]
    fn explicit_lengths_survive_resolution_for_full_and_masked_targets() {
        let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
        for mode in ["legacy", "sweep", "general", "auto"] {
            for masked in [false, true] {
                let mut cover = vec![true; 16];
                cover[1] = !masked;
                let target = Target::new(
                    PartialMatrix::new(crate::matrix::Matrix::identity(4), cover, "target")
                        .unwrap(),
                    false,
                    false,
                );
                for (flag, bounds, adaptive) in [
                    ("--gates 9", (22, 31), true),
                    ("--n-start-gates 50", (50, 51), false),
                ] {
                    let c = parse(&format!(
                        "target.txt --search-mode {mode} --gate-count 2 {flag}"
                    ))
                    .unwrap();
                    let resolved = c.resolve_search(&target, &lib, true);
                    assert_eq!(
                        (resolved.scheme.min_gates, resolved.scheme.max_gates),
                        bounds
                    );
                    assert_eq!(resolved.update_scheme, adaptive);
                }
            }
        }
    }

    #[test]
    fn removed_length_shortcuts_are_rejected() {
        for flag in ["--short-start 5", "--long-start-every 4", "--goal-length"] {
            assert!(parse(&format!("cx.txt {flag}"))
                .unwrap_err()
                .contains("unknown option"));
        }
    }

    #[test]
    fn quality_backup_accounts_for_native_expansion_and_all_goals() {
        let lib =
            GateLibrary::load(3, "data/gates/CliffordT", "data/gates/composite_rccx").unwrap();
        let circuit = Circuit::new(vec![lib.find("rccx", &[0, 1, 2]).unwrap()], 3, &lib);
        let mut native = circuit.clone();
        native.expand(&lib);
        let names = ["t".into(), "tdg".into()];
        let mut c = parse("target.txt --search-mode auto --expand --tcount 4").unwrap();
        c.tdepth = native.depth(&names, &lib) as i64;
        c.gatecount = native.non_identity() as i64;
        c.cost_required = native.cost();
        let backup = c.feasible_backup(&circuit, &lib).unwrap();
        assert_eq!(backup.gates, native.gates);
        for constraint in 0..4 {
            let mut tighter = c.clone();
            match constraint {
                0 => tighter.tcount -= 1,
                1 => tighter.tdepth -= 1,
                2 => tighter.gatecount -= 1,
                _ => tighter.cost_required -= 0.001,
            }
            assert!(tighter.feasible_backup(&circuit, &lib).is_none());
        }
        c.search.optimized = false;
        assert!(c.feasible_backup(&circuit, &lib).is_none());
    }

    #[test]
    fn joint_native_limits_enable_delayed_cost_search_without_overriding_explicit_options() {
        let lib = GateLibrary::load(3, "data/gates/CliffordT", "").unwrap();
        let target = Target::new(
            PartialMatrix::read("data/input/64/comparison/ccx.txt").unwrap(),
            false,
            false,
        );
        let c =
            parse("arbitrary.txt --tcount 7 --tdepth 3 --gate-count 25 --cost-required 9").unwrap();
        let resolved = c.resolve_search(&target, &lib, true);
        assert_eq!(resolved.resource_fallback_after, 0.3);
        assert!(matches!(
            resolved.search.quality_metric,
            search::QualityMetric::TCount
        ));
        for flag in [
            "--resource-fallback-after 0",
            "--quality-objective tcount",
            "--quality-weight .01",
            "--quality-slack 2",
        ] {
            let explicit = parse(&format!(
                "arbitrary.txt --tcount 7 --gate-count 25 --cost-required 9 {flag}"
            ))
            .unwrap();
            assert_eq!(
                explicit
                    .resolve_search(&target, &lib, true)
                    .resource_fallback_after,
                0.0
            );
        }
        let zero = parse("arbitrary.txt --tcount 0 --gate-count 25 --cost-required 9").unwrap();
        assert_eq!(
            zero.resolve_search(&target, &lib, true)
                .resource_fallback_after,
            0.0
        );
    }

    #[test]
    fn automatic_joint_bounds_use_original_native_search_by_default() {
        let lib = GateLibrary::load(3, "data/gates/CliffordT", "").unwrap();
        let c = parse("unnamed.txt --search-mode auto --tcount 7 --gate-count 21").unwrap();
        let target = Target::new(
            PartialMatrix::read("data/input/64/comparison/ccx.txt").unwrap(),
            false,
            false,
        );
        let resolved = c.resolve_search(&target, &lib, true);
        assert!(!resolved.bounded_structure);
        assert_eq!(resolved.pauli_restarts, 0);
        assert!(!resolved.affine_search);
        assert!(!resolved.rccx_macros);
        assert_eq!(resolved.phase_restarts, 0);
        assert!(resolved.clifford_fold);
    }

    #[test]
    fn explicit_joint_structure_retains_general_search() {
        let lib = GateLibrary::load(3, "data/gates/CliffordT", "").unwrap();
        let c = parse("unnamed.txt --search-mode auto --bounded-structure --phase-restarts 4 --phase-factor --tcount 7 --tdepth 3 --gate-count 25 --cost-required 9").unwrap();
        let target = Target::new(
            PartialMatrix::read("data/input/64/comparison/ccx.txt").unwrap(),
            false,
            false,
        );
        let resolved = c.resolve_search(&target, &lib, true);
        assert!(resolved.bounded_structure);
        assert!(resolved.phase_restarts > 0);
        assert_eq!(resolved.general_every, 4);
        assert!(resolved.structured_time_slice > 0.0);
    }

    #[test]
    fn automatic_settings_follow_features_and_preserve_explicit_options() {
        let lib = GateLibrary::load(3, "data/gates/CliffordT", "").unwrap();
        let c=parse("arbitrary-name.txt --search-mode auto --tcount 7 --phase-restarts 0 --pauli-restarts 1 --start-temp .17 --swap-probability 0 --quality-slack 2 --gates 9").unwrap();
        let target = Target::new(
            PartialMatrix::read("data/input/64/comparison/ccx.txt").unwrap(),
            false,
            false,
        );
        let resolved = c.resolve_search(&target, &lib, true);
        assert_eq!(resolved.phase_restarts, 0);
        assert_eq!(resolved.pauli_restarts, 1);
        let mut disabled = c.clone();
        disabled.pauli_restarts = 0;
        assert_eq!(
            disabled.resolve_search(&target, &lib, true).pauli_restarts,
            0
        );
        assert_eq!(resolved.search.start_temp_base, 0.17);
        assert_eq!(resolved.search.swap_probability, 0.0);
        assert_eq!(resolved.search.quality_slack, Some(2.0));
        assert_eq!(resolved.scheme.min_start_gates, c.scheme.min_start_gates);
        assert!(!resolved.affine_search && !resolved.rccx_macros);
        assert_eq!(resolved.permutation_restarts, 0);
        assert!(resolved.clifford_factor && resolved.phase_fold && resolved.pauli_fold);
        let mut other = c.clone();
        other.input = "unrelated-path.dat".into();
        let second = other.resolve_search(&target, &lib, true);
        assert_eq!(
            (
                resolved.phase_restarts,
                resolved.permutation_restarts,
                resolved.affine_search,
                resolved.rccx_macros
            ),
            (
                second.phase_restarts,
                second.permutation_restarts,
                second.affine_search,
                second.rccx_macros
            )
        );
        let fallback = c.resolve_search(&target, &lib, false);
        assert_eq!(fallback.permutation_restarts, 0);
        assert!(!fallback.rccx_macros);
    }
    #[test]
    fn automatic_masked_search_skips_unused_target_variants() {
        let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
        let mut cover = vec![true; 16];
        cover[1] = false;
        let target = Target::new(
            PartialMatrix::new(crate::matrix::Matrix::identity(4), cover, "mask").unwrap(),
            false,
            false,
        );
        let c = parse("anything.txt --search-mode auto --tcount 0").unwrap();
        assert!(c.original_sweep(&target.original));
        let resolved = c.resolve_search(&target, &lib, true);
        assert!(resolved.search.masked_sweep);
        assert!(!resolved.clifford_factor);
        assert!(resolved.phase_fold && resolved.pauli_fold);
        assert_eq!(resolved.search.masked_proposals, 4);
        assert_eq!(resolved.search.start_temp_base, 0.3);
        assert_eq!(
            (resolved.scheme.min_gates, resolved.scheme.max_gates),
            (30, 120)
        );
        assert_eq!(resolved.phase_restarts, 0);
        assert_eq!(resolved.search.quality_slack, Some(0.0));
    }
    #[test]
    fn automatic_approximation_keeps_uncapped_search_and_budgeted_kernels_opt_in() {
        let lib = GateLibrary::load(1, "data/gates/CliffordT", "").unwrap();
        let c = parse("anything.txt --search-mode auto --tcount 16 --epsilon .01").unwrap();
        let target = Target::new(
            PartialMatrix::new(crate::matrix::Matrix::identity(2), vec![true; 4], "id").unwrap(),
            false,
            false,
        );
        let resolved = c.resolve_search(&target, &lib, true);
        assert!(resolved.search.quality_slack.is_none());
        assert_eq!(resolved.search.swap_probability, 0.0);
        assert_eq!(resolved.single_qubit_restarts, 0);
        assert!(resolved.phase_fold && resolved.pauli_fold && resolved.clifford_factor);
    }

    #[test]
    fn defaults_use_tuned_search_and_original_limits() {
        let c = parse("cx.txt").unwrap();
        assert_eq!(c.input, "data/input/cx.txt");
        assert_eq!(c.output, "data/output/cx/");
        assert_eq!(c.threads, 1);
        assert_eq!(c.circuits, 10);
        assert_eq!(c.time, 100.0);
        assert!(c.independent);
        assert!(c.search.optimized);
        assert!(c.auto_search);
        assert!(!parse("cx.txt --search-mode sweep").unwrap().auto_search);
        assert!(!parse("cx.txt --search-mode legacy").unwrap().auto_search);
        assert_eq!(c.search.cost_power, 1.0);
        assert_eq!(c.search.stall_sweeps, 12);
        assert!(!c.search.sweep_permutations);
        assert!(
            parse("cx.txt --search-mode legacy")
                .unwrap()
                .load_target()
                .unwrap()
                .variants
                .len()
                > 1
        );
        assert_eq!(c.load_target().unwrap().variants.len(), 1);
    }
    #[test]
    fn joint_search_objective_can_change_without_relaxing_acceptance_limits() {
        let c = parse("cx.txt --tcount 8 --tdepth 4 --gate-count 35 --quality-objective tdepth")
            .unwrap();
        assert!(matches!(
            c.search.quality_metric,
            search::QualityMetric::TDepth
        ));
        assert_eq!(c.search.quality_goal, Some(4.0));
        assert_eq!((c.tcount, c.tdepth, c.gatecount), (8, 4, 35));
        assert!(parse("cx.txt --quality-objective tdepth --tcount 8").is_err());
        assert!(parse("cx.txt --quality-objective invalid").is_err());
        assert!(parse("cx.txt --legacy-fallback-after NaN").is_err());
        assert!(parse("cx.txt --resource-fallback-after 1").is_err());
        assert!(parse("cx.txt --resource-fallback-after 1 --cost-required 10").is_ok());
    }
    #[test]
    fn aliases_and_absolute_paths() {
        let c=parse("/tmp/test.qasm --absolute-input --absolute-output -o /tmp/out -h 4 -tc 7 -td 3 -eps 0.01 -q --seed 42").unwrap();
        assert_eq!(c.input, "/tmp/test.qasm");
        assert_eq!(c.output, "/tmp/out/");
        assert_eq!(c.tcount, 7);
        assert_eq!(c.tdepth, 3);
        assert_eq!(c.seed, 42);
        assert!(!c.independent);
    }
    #[test]
    fn rejects_invalid_arguments() {
        for s in [
            "cx.txt --threads 0",
            "cx.txt --time NaN",
            "cx.txt --time -1",
            "cx.txt --pid 2",
            "cx.txt --epsilon",
            "cx.txt --unknown",
        ] {
            assert!(parse(s).is_err(), "{s}");
        }
    }
    #[test]
    fn composite_weight_can_exceed_one() {
        assert_eq!(parse("cx.txt --pcomp 2.5").unwrap().search.pcomp, 2.5);
        assert!(parse("cx.txt --pcomp -1").is_err());
        assert!(parse("cx.txt --pcomp NaN").is_err());
    }
    #[test]
    fn parse_qasm_qubits() {
        assert_eq!(qasm_qubits("OPENQASM 2.0;\nqreg q[3];\n"), Ok(3));
    }
}
