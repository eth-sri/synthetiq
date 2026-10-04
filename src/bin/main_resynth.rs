fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a == "--help") {
        print!("{}", synthetiq::cli::HELP);
        return;
    }
    let result = synthetiq::cli::Config::parse(&args).and_then(|config| {
        config.print();
        let (lib, mut circuit) = config.load_circuit()?;
        let names = vec!["t".into(), "tdg".into()];
        let before = (
            circuit.non_identity(),
            circuit.count(&names, &lib),
            circuit.depth(&names, &lib),
        );
        // main_resynth in the original fixes these parameters, independently of CLI flags.
        synthetiq::resynthesis::Resynth::new(12, true, names.clone()).run(&mut circuit, &lib);
        println!(
            "Gates: {} -> {}\nT-count: {} -> {}\nT-depth: {} -> {}",
            before.0,
            circuit.non_identity(),
            before.1,
            circuit.count(&names, &lib),
            before.2,
            circuit.depth(&names, &lib)
        );
        Ok::<(), String>(())
    });
    if let Err(error) = result {
        eprintln!("main_resynth: {error}");
        std::process::exit(1);
    }
}
