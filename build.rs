//! Precompute a target-independent Clifford rewrite table at build time.
use std::{env, fs, path::PathBuf};

fn main() {
    println!("cargo:rerun-if-changed=build.rs");
    let q = 4;
    let identity = (0..q).fold(0, |state, bit| state | (1 << (bit * q + bit)));
    let mut parent = vec![u16::MAX; 1 << (q * q)];
    let mut step = vec![0u8; parent.len()];
    let mut queue = vec![identity];
    parent[identity] = identity as u16;
    let mut head = 0;
    while head < queue.len() {
        let state = queue[head];
        head += 1;
        for control in 0..q {
            let source = (state >> (control * q)) & ((1 << q) - 1);
            for target in 0..q {
                if control == target {
                    continue;
                }
                let next = state ^ (source << (target * q));
                if parent[next] == u16::MAX {
                    parent[next] = state as u16;
                    step[next] = (control * q + target) as u8;
                    queue.push(next);
                }
            }
        }
    }
    assert_eq!(queue.len(), 20160);
    let directory = PathBuf::from(env::var_os("OUT_DIR").expect("Cargo output directory"));
    let bytes: Vec<_> = parent.into_iter().flat_map(u16::to_le_bytes).collect();
    fs::write(directory.join("linear4-parent.bin"), bytes).expect("write CNOT parents");
    fs::write(directory.join("linear4-step.bin"), step).expect("write CNOT steps");
}
