use std::{
    fs,
    path::PathBuf,
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};
use synthetiq::{circuit::Circuit, gates::GateLibrary, partial::PartialMatrix};

struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "synthetiq-cli-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::create_dir(&p).unwrap();
        Self(p)
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

fn synthesize(threads: usize, extras: &[&str]) -> (Temp, String) {
    let temp = Temp::new();
    let output = temp.0.join("circuits");
    let times = temp.0.join("times.tsv");
    let result = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
        .args(["cx.txt", "--absolute-output", "--output"])
        .arg(&output)
        .args([
            "--threads",
            &threads.to_string(),
            "--circuits",
            "5",
            "--time",
            "10",
            "--n-start-gates",
            "0",
            "--qubit_independence",
            "--seed",
            "123",
            "--save",
            "--times-file",
        ])
        .arg(&times)
        .args(extras)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let stdout = String::from_utf8(result.stdout).unwrap();
    assert_eq!(stdout.matches("Gates: ").count(), 5, "{stdout}");
    let target = PartialMatrix::read("data/input/cx.txt").unwrap();
    let lib = GateLibrary::load(2, "data/gates/CliffordT", "data/gates/composite_gates").unwrap();
    if !extras.contains(&"--no-save-circuits") {
        let paths: Vec<_> = fs::read_dir(&output)
            .unwrap()
            .map(|e| e.unwrap().path())
            .collect();
        assert_eq!(
            paths.len(),
            5,
            "parallel search must neither lose nor overshoot results"
        );
        for path in paths {
            let circuit = Circuit::from_qasm(&fs::read_to_string(&path).unwrap(), &lib).unwrap();
            assert_eq!(
                target.exact_cost(circuit.matrix(), 1e-6),
                0.0,
                "{}",
                path.display()
            );
            assert_eq!(circuit.non_identity(), 1);
            let fields: Vec<_> = path
                .file_stem()
                .unwrap()
                .to_str()
                .unwrap()
                .split('-')
                .map(String::from)
                .collect();
            assert_eq!(fields.len(), 5);
            assert_eq!(fields[0], format!("{:.6}", circuit.cost()));
        }
    } else {
        assert_eq!(fs::read_dir(output).unwrap().count(), 0);
    }
    let timing = fs::read_to_string(times).unwrap();
    let fields: Vec<_> = timing.trim().split('\t').collect();
    assert_eq!(fields.len(), 6);
    assert!(fields[1].parse::<f64>().unwrap().is_finite());
    assert!(fields[3].parse::<usize>().unwrap() >= 5);
    (temp, stdout)
}

#[test]
fn single_worker_outputs_satisfy_target() {
    synthesize(1, &[]);
}
#[test]
fn parallel_workers_outputs_satisfy_target() {
    synthesize(4, &[]);
}
#[test]
fn constrained_runs_reset_and_count_optima() {
    let (temp, _) = synthesize(
        2,
        &[
            "--tcount",
            "0",
            "--tdepth",
            "0",
            "--gate-count",
            "1",
            "--cost-required",
            "10",
        ],
    );
    let text = fs::read_to_string(temp.0.join("times.tsv")).unwrap();
    assert_eq!(text.trim().split('\t').next_back(), Some("5"));
}
#[test]
fn legacy_cost_and_no_save_flags() {
    synthesize(
        1,
        &[
            "--simple",
            "--no-perms",
            "--no-resynth",
            "--no-save-circuits",
        ],
    );
}
#[test]
fn qasm_input_and_resynthesis_executable() {
    let temp = Temp::new();
    let input = temp.0.join("h.qasm");
    fs::write(
        &input,
        "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg qubits[1];\nh qubits[0];\n",
    )
    .unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
        .arg(&input)
        .args(["--absolute-input", "--absolute-output", "--output"])
        .arg(temp.0.join("out"))
        .args([
            "--circuits",
            "2",
            "--time",
            "10",
            "--n-start-gates",
            "0",
            "--seed",
            "5",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&result.stdout)
            .matches("Gates: ")
            .count(),
        2
    );
    let result = Command::new(env!("CARGO_BIN_EXE_main_resynth"))
        .args(["example_circuit_to_simplify.qasm"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stdout).contains("Gates:"));
}
#[test]
fn invalid_input_fails_cleanly() {
    for args in [
        vec![],
        vec!["missing.txt"],
        vec!["cx.txt", "--threads", "0"],
        vec!["cx.txt", "--pid", "NaN"],
        vec!["cx.txt", "--unknown"],
    ] {
        let output = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
            .args(args)
            .output()
            .unwrap();
        assert!(!output.status.success());
        assert!(!String::from_utf8_lossy(&output.stderr).contains("panicked"));
    }
}

#[test]
fn legacy_recomputed_rollback_remains_available() {
    synthesize(1, &["--search-mode", "legacy", "--no-cache-rollback"]);
}
#[test]
fn legacy_cached_rollback_remains_available() {
    synthesize(1, &["--search-mode", "legacy"]);
}
#[test]
fn sweep_permutation_scoring_remains_available() {
    synthesize(1, &["--sweep-permutations"]);
}

#[test]
fn sweep_quality_cap_and_penalty_preserve_output_constraints() {
    synthesize(
        1,
        &[
            "--quality-slack",
            "0",
            "--quality-weight",
            "0.003",
            "--tcount",
            "0",
        ],
    );
}

#[test]
fn clifford_macros_are_expanded_before_saving() {
    synthesize(1, &["--clifford-macros", "--tcount", "0"]);
}

#[test]
fn a_single_macro_slot_produces_a_valid_native_clifford_circuit() {
    let temp = Temp::new();
    let input = temp.0.join("x.qasm");
    fs::write(&input, "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg qubits[1];\nh qubits[0];\ns qubits[0];\ns qubits[0];\nh qubits[0];\n").unwrap();
    let output = temp.0.join("out");
    let result = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
        .arg(&input)
        .args(["--absolute-input", "--absolute-output", "--output"])
        .arg(&output)
        .args([
            "--clifford-macros",
            "--tcount",
            "0",
            "--n-start-gates",
            "0",
            "--circuits",
            "1",
            "--time",
            "5",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let paths: Vec<_> = fs::read_dir(output)
        .unwrap()
        .map(|p| p.unwrap().path())
        .collect();
    assert_eq!(paths.len(), 1);
    let qasm = fs::read_to_string(&paths[0]).unwrap();
    assert!(!qasm.contains("__clifford"));
    let lib = GateLibrary::load(1, "data/gates/CliffordT", "").unwrap();
    let actual = Circuit::from_qasm(&qasm, &lib).unwrap();
    let expected = Circuit::from_qasm(&fs::read_to_string(input).unwrap(), &lib).unwrap();
    let target = PartialMatrix::new(expected.matrix().clone(), vec![true; 4], "x").unwrap();
    assert_eq!(target.exact_cost(actual.matrix(), 1e-6), 0.0);
    assert_eq!(actual.count(&["t".into(), "tdg".into()], &lib), 0);
}

#[test]
fn macro_postprocessing_cannot_reintroduce_internal_output_gates() {
    let temp = Temp::new();
    let output = temp.0.join("out");
    let result = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
        .args(["61/2qbs/10_0_2_0_5.txt", "--absolute-output", "--output"])
        .arg(&output)
        .args([
            "--clifford-macros",
            "--tcount",
            "2",
            "--seed",
            "202",
            "--circuits",
            "3",
            "--time",
            "10",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let paths: Vec<_> = fs::read_dir(output)
        .unwrap()
        .map(|p| p.unwrap().path())
        .collect();
    assert_eq!(paths.len(), 3);
    let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
    let target = PartialMatrix::read("data/input/61/2qbs/10_0_2_0_5.txt").unwrap();
    for path in paths {
        let qasm = fs::read_to_string(path).unwrap();
        assert!(!qasm.contains("__clifford"));
        let circuit = Circuit::from_qasm(&qasm, &lib).unwrap();
        assert_eq!(target.exact_cost(circuit.matrix(), 1e-6), 0.0);
        assert!(circuit.count(&["t".into(), "tdg".into()], &lib) <= 2);
    }
}

#[test]
fn phase_gadget_synthesis_restores_original_basis_and_optimal_t_count() {
    let temp = Temp::new();
    let output = temp.0.join("out");
    let result = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
        .args(["64/comparison/cch.txt", "--absolute-output", "--output"])
        .arg(&output)
        .args([
            "--phase-gadgets",
            "--tcount",
            "9",
            "--seed",
            "201",
            "--circuits",
            "4",
            "--time",
            "10",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let paths: Vec<_> = fs::read_dir(output)
        .unwrap()
        .map(|p| p.unwrap().path())
        .collect();
    assert_eq!(paths.len(), 4);
    let lib = GateLibrary::load(3, "data/gates/CliffordT", "").unwrap();
    let target = PartialMatrix::read("data/input/64/comparison/cch.txt").unwrap();
    for path in paths {
        let qasm = fs::read_to_string(path).unwrap();
        let circuit = Circuit::from_qasm(&qasm, &lib).unwrap();
        assert_eq!(target.exact_cost(circuit.matrix(), 1e-6), 0.0);
        assert!(circuit.count(&["t".into(), "tdg".into()], &lib) <= 9);
    }
}

#[test]
fn phase_layout_outputs_meet_original_t_depth_constraints() {
    for (name, qubits, depth) in [("ccx", 3, 3), ("cch", 3, 4), ("adder", 4, 2)] {
        let temp = Temp::new();
        let output = temp.0.join("circuits");
        let result = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
            .arg(format!("64/comparison/{name}.txt"))
            .args(["--absolute-output", "--output"])
            .arg(&output)
            .args([
                "--phase-gadgets",
                "--tdepth",
                &depth.to_string(),
                "--seed",
                "201",
                "--circuits",
                "4",
                "--time",
                "10",
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let paths: Vec<_> = fs::read_dir(output)
            .unwrap()
            .map(|p| p.unwrap().path())
            .collect();
        assert_eq!(paths.len(), 4, "{name}");
        let lib = GateLibrary::load(qubits, "data/gates/CliffordT", "").unwrap();
        let target = PartialMatrix::read(format!("data/input/64/comparison/{name}.txt")).unwrap();
        for path in paths {
            let circuit = Circuit::from_qasm(&fs::read_to_string(path).unwrap(), &lib).unwrap();
            assert_eq!(target.exact_cost(circuit.matrix(), 1e-6), 0.0);
            assert!(circuit.depth(&["t".into(), "tdg".into()], &lib) <= depth);
        }
    }
}

#[test]
fn arbitrary_mask_coordinate_sweeps_preserve_original_constraints() {
    let temp = Temp::new();
    let output = temp.0.join("circuits");
    let result = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
        .args(["61/2qbs/10_2_2_0_2.txt", "--absolute-output", "--output"])
        .arg(&output)
        .args([
            "--masked-sweep",
            "--masked-proposals",
            "4",
            "--quality-slack",
            "0",
            "--tcount",
            "2",
            "--seed",
            "201",
            "--circuits",
            "4",
            "--time",
            "10",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let paths: Vec<_> = fs::read_dir(output)
        .unwrap()
        .map(|p| p.unwrap().path())
        .collect();
    assert_eq!(paths.len(), 4);
    let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
    let target = PartialMatrix::read("data/input/61/2qbs/10_2_2_0_2.txt").unwrap();
    for path in paths {
        let circuit = Circuit::from_qasm(&fs::read_to_string(path).unwrap(), &lib).unwrap();
        assert_eq!(target.exact_cost(circuit.matrix(), 1e-6), 0.0);
        assert!(circuit.count(&["t".into(), "tdg".into()], &lib) <= 2);
    }
}

#[test]
fn phase_fold_follows_final_macro_expansion() {
    let temp = Temp::new();
    let input = temp.0.join("relative.qasm");
    let output = temp.0.join("circuits");
    let lib = GateLibrary::load(3, "data/gates/CliffordT", "data/gates/composite_rccx").unwrap();
    let mut source = Circuit::new(vec![lib.find("rccx", &[0, 1, 2]).unwrap()], 3, &lib);
    source.expand(&lib);
    fs::write(&input, source.qasm(&lib)).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
        .arg(&input)
        .args(["--absolute-input", "--absolute-output", "--output"])
        .arg(&output)
        .args([
            "--rccx-macros",
            "--phase-fold",
            "--tdepth",
            "2",
            "--n-start-gates",
            "0",
            "--circuits",
            "3",
            "--time",
            "5",
            "--seed",
            "0",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let paths: Vec<_> = fs::read_dir(output)
        .unwrap()
        .map(|p| p.unwrap().path())
        .collect();
    assert_eq!(paths.len(), 3);
    let spec = PartialMatrix::new(source.matrix().clone(), vec![true; 64], "relative").unwrap();
    for path in paths {
        let candidate = Circuit::from_qasm(&fs::read_to_string(path).unwrap(), &lib).unwrap();
        assert_eq!(spec.exact_cost(candidate.matrix(), 1e-6), 0.0);
        assert_eq!(candidate.count(&["t".into(), "tdg".into()], &lib), 4);
        assert_eq!(candidate.depth(&["t".into(), "tdg".into()], &lib), 2);
    }
}

#[test]
fn explicit_start_length_synthesizes_a_two_gate_circuit() {
    let temp = Temp::new();
    let input = temp.0.join("two.qasm");
    let output = temp.0.join("circuits");
    fs::write(&input,"OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg qubits[2];\nh qubits[0];\ncx qubits[0], qubits[1];\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
        .arg(&input)
        .args(["--absolute-input", "--absolute-output", "--output"])
        .arg(&output)
        .args([
            "--n-start-gates",
            "1",
            "--circuits",
            "1",
            "--time",
            "5",
            "--seed",
            "0",
            "--gate-count",
            "2",
        ])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let paths: Vec<_> = fs::read_dir(output)
        .unwrap()
        .map(|p| p.unwrap().path())
        .collect();
    assert_eq!(paths.len(), 1);
    let lib = GateLibrary::load(2, "data/gates/CliffordT", "").unwrap();
    let original = Circuit::from_qasm(&fs::read_to_string(input).unwrap(), &lib).unwrap();
    let found = Circuit::from_qasm(&fs::read_to_string(&paths[0]).unwrap(), &lib).unwrap();
    let target = PartialMatrix::new(original.matrix().clone(), vec![true; 16], "two").unwrap();
    assert_eq!(target.exact_cost(found.matrix(), 1e-6), 0.0);
    assert_eq!(found.non_identity(), 2);
}

#[test]
fn automatic_mode_preserves_native_parallel_constraints() {
    synthesize(
        4,
        &[
            "--search-mode",
            "auto",
            "--tcount",
            "0",
            "--gate-count",
            "1",
        ],
    );
}

#[test]
fn explicit_pauli_search_reaches_t_bounds_on_generic_and_barrier_targets() {
    for (fixture, qubits, bound) in [
        ("tests/fixtures/pauli-generic.qasm", 3, 8),
        ("tests/fixtures/pauli-barrier.qasm", 2, 12),
    ] {
        let temp = Temp::new();
        let output = temp.0.join("out");
        let result = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
            .args([fixture, "--absolute-input", "--absolute-output", "--output"])
            .arg(&output)
            .args([
                "--pauli-restarts",
                "1",
                "--tcount",
                &bound.to_string(),
                "--circuits",
                "1",
                "--time",
                "3",
                "--seed",
                "423",
            ])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let paths: Vec<_> = fs::read_dir(&output)
            .unwrap()
            .map(|p| p.unwrap().path())
            .collect();
        assert_eq!(
            paths.len(),
            1,
            "{}",
            String::from_utf8_lossy(&result.stdout)
        );
        let lib = GateLibrary::load(qubits, "data/gates/CliffordT", "").unwrap();
        let source = Circuit::from_qasm(&fs::read_to_string(fixture).unwrap(), &lib).unwrap();
        let found = Circuit::from_qasm(&fs::read_to_string(&paths[0]).unwrap(), &lib).unwrap();
        let target = PartialMatrix::new(
            source.matrix().clone(),
            vec![true; 1 << (2 * qubits)],
            "fixture",
        )
        .unwrap();
        assert_eq!(target.exact_cost(found.matrix(), 1e-6), 0.0);
        assert!(found.count(&["t".into(), "tdg".into()], &lib) <= bound);
    }
}

#[test]
fn default_budget_free_basis_search_restores_the_original_operator() {
    let temp = Temp::new();
    let input = temp.0.join("basis.qasm");
    let output = temp.0.join("circuits");
    fs::write(&input, "OPENQASM 2.0;\ninclude \"qelib1.inc\";\nqreg qubits[1];\nh qubits[0];\nt qubits[0];\nh qubits[0];\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_synthetiq"))
        .arg(&input)
        .args(["--absolute-input", "--absolute-output", "--output"])
        .arg(&output)
        .args(["--circuits", "1", "--time", "5", "--seed", "172"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let paths: Vec<_> = fs::read_dir(&output)
        .unwrap()
        .map(|p| p.unwrap().path())
        .collect();
    assert_eq!(
        paths.len(),
        1,
        "{}",
        String::from_utf8_lossy(&result.stdout)
    );
    let lib = GateLibrary::load(1, "data/gates/CliffordT", "").unwrap();
    let original = Circuit::from_qasm(&fs::read_to_string(input).unwrap(), &lib).unwrap();
    let found = Circuit::from_qasm(&fs::read_to_string(&paths[0]).unwrap(), &lib).unwrap();
    let spec = PartialMatrix::new(original.matrix().clone(), vec![true; 4], "basis").unwrap();
    assert_eq!(spec.exact_cost(found.matrix(), 1e-6), 0.0);
    assert!(found
        .gates
        .iter()
        .all(|&g| lib.gates[g].decomposition.is_empty()));
}
