fn main() {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a == "--help") {
        print!("{}", synthetiq::cli::HELP);
        return;
    }
    let result = synthetiq::cli::Config::parse(&args)
        .and_then(|config| synthetiq::cli::run_cli(&config).map(|_| ()));
    if let Err(error) = result {
        eprintln!("synthetiq: {error}");
        std::process::exit(1);
    }
}
