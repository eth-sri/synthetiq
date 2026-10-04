FROM rust:1.96-bookworm AS engine
WORKDIR /usr/synthetiq
COPY Cargo.toml Cargo.lock build.rs ./
COPY src ./src
RUN cargo build --release --bin synthetiq --bin main_resynth --bin comparison_generator

RUN apt-get update && apt-get install -y --no-install-recommends g++ make \
    && rm -rf /var/lib/apt/lists/*
COPY makefile ./
COPY synthetiq ./synthetiq
COPY include/eigen-3.3.9 ./include/eigen-3.3.9
RUN make cpp

FROM python:3.10-bookworm
WORKDIR /usr/synthetiq
RUN apt-get update && apt-get install -y --no-install-recommends libgomp1 \
    && rm -rf /var/lib/apt/lists/*
COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
COPY --from=engine /usr/synthetiq/bin/cpp/main ./bin/main
COPY --from=engine /usr/synthetiq/bin/cpp/main_resynth ./bin/main_resynth
COPY --from=engine /usr/synthetiq/bin/cpp/comparison_generator ./bin/comparison_generator
COPY --from=engine /usr/synthetiq/target/release/synthetiq ./bin/rust
COPY --from=engine /usr/synthetiq/target/release/main_resynth ./bin/rust_resynth
COPY --from=engine /usr/synthetiq/target/release/comparison_generator ./bin/rust_comparison_generator
