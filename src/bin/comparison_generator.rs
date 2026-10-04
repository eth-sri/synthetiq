//! Compatibility port of comparison_generator.cpp, with optional size limits.
use std::{env, fs, path::Path};
use synthetiq::{
    cli::{self, Config},
    gates::GateLibrary,
    generator::{generate_circuit, generate_matrix, SpecificationKind},
    partial::Target,
    rng::Rng,
};

fn main() {
    if let Err(error) = run() {
        eprintln!("comparison_generator: {error}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let mut qubits = vec![2, 3, 4];
    let mut kinds = vec![
        SpecificationKind::Full,
        SpecificationKind::Isometry,
        SpecificationKind::Random,
    ];
    let mut gates: usize = 10;
    let mut extra: usize = 0;
    let mut number: usize = 10;
    let mut seed: u64 = 42;
    let mut threads: usize = 128;
    let mut seconds: f64 = 100.0;
    let mut gate_set = "CliffordT".to_string();
    let mut output: Option<String> = None;
    let mut max_attempts: Option<usize> = None;
    let args: Vec<String> = env::args().skip(1).collect();
    let mut i = 0;
    while i < args.len() {
        let flag = args[i].as_str();
        if flag == "--help" || flag == "-h" {
            println!(
                "Usage: comparison_generator [--qubits N] [--gates N] [--type 0|1|2]\n\
                [--count N] [--extra-gates N] [--seed N] [--threads N] [--time SECONDS]\n\
                [--gate-set NAME] [--output DIRECTORY] [--max-attempts N]\n\n\
                Defaults reproduce the reference experiment: 2, 3, 4 qubits; 10 gates;\n\
                all specification types; 10 files each; 128 threads; 100 s per candidate.\n\
                This can take a long time. The depth field in filenames retains the\n\
                reference's legacy value of zero."
            );
            return Ok(());
        }
        i += 1;
        let value = args
            .get(i)
            .ok_or_else(|| format!("{flag} requires a value"))?;
        macro_rules! number {
            () => {
                value
                    .parse()
                    .map_err(|_| format!("invalid value for {flag}"))?
            };
        }
        match flag {
            "--qubits" => qubits = vec![number!()],
            "--gates" => gates = number!(),
            "--extra-gates" => extra = number!(),
            "--type" => {
                kinds = vec![SpecificationKind::try_from(
                    value.parse::<usize>().map_err(|_| "invalid type")?,
                )?]
            }
            "--count" => number = number!(),
            "--seed" => seed = number!(),
            "--threads" => threads = number!(),
            "--time" => seconds = number!(),
            "--gate-set" => gate_set = value.clone(),
            "--output" => output = Some(value.clone()),
            "--max-attempts" => max_attempts = Some(number!()),
            _ => return Err(format!("unknown option: {flag}")),
        }
        i += 1;
    }
    if threads == 0
        || !seconds.is_finite()
        || seconds < 0.0
        || qubits.iter().any(|&n| n == 0 || n > 10)
    {
        return Err(
            "threads must be positive, time finite and nonnegative, and qubits in 1..=10".into(),
        );
    }
    if gates.checked_add(extra).is_none() || gates > i64::MAX as usize {
        return Err("gate count is too large".into());
    }
    for n_qubits in qubits {
        for &kind in &kinds {
            let directory = output
                .clone()
                .unwrap_or_else(|| format!("data/input/61{n_qubits}qbs"));
            fs::create_dir_all(&directory).map_err(|e| format!("{directory}: {e}"))?;
            let needle = format!("{gates}_{}_", kind as u8);
            let mut count = 0;
            for entry in fs::read_dir(&directory).map_err(|e| e.to_string())? {
                if entry
                    .map_err(|e| e.to_string())?
                    .path()
                    .to_string_lossy()
                    .contains(&needle)
                {
                    count += 1;
                }
            }
            println!("Working on qbs {n_qubits}\nWorking on set {gate_set}\nWorking on gates {gates}\nWorking on type {}\n{count}", kind as u8);
            let lib = GateLibrary::load(
                n_qubits,
                &format!("data/gates/{gate_set}"),
                "data/gates/composite_gates",
            )?;
            let mut rng = Rng::new(seed);
            let mut attempts = 0;
            while count < number && max_attempts.is_none_or(|maximum| attempts < maximum) {
                attempts += 1;
                let original = generate_circuit(&lib, &mut rng, gates, extra, 0.0, true);
                let partial = generate_matrix(&original, kind, &mut rng)?;
                let mut config = Config::parse(&["comparison.txt".into()])?;
                config.basic = format!("data/gates/{gate_set}");
                config.threads = threads;
                config.circuits = 1;
                config.gatecount = gates as i64 - 1;
                config.time = seconds;
                config.save_circuits = false;
                let target = Target::new(partial, true, false);
                println!("starting run");
                let result = cli::run_target(&config, target.clone())?;
                // Preserve the original acceptance criterion, even though the
                // synthesis run's requested maximum is one fewer gate.
                if let (Some(tcount), Some(gatecount)) = (result.best_tcount, result.best_gatecount)
                {
                    if gatecount >= gates as f64 {
                        let best_tcount = original
                            .count(&["t".into(), "tdg".into()], &lib)
                            .min(tcount as usize);
                        // C++ reads result["depth"], an absent std::map key, so
                        // the filename depth is always min(original_depth, 0).
                        let path = Path::new(&directory).join(format!(
                            "{gates}_{}_{best_tcount}_0_{count}.txt",
                            kind as u8
                        ));
                        println!("{}", path.display());
                        target.original.write(path)?;
                        count += 1;
                    }
                }
            }
        }
    }
    Ok(())
}
