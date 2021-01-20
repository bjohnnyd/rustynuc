use criterion::{criterion_group, criterion_main, Criterion};
use log::error;

fn is_s(seq: &[u8], pos: usize) -> bool {
    if pos >= seq.len() {
        error!("Reference sequence is shorter than BAM alignment positions");
        std::process::exit(1)
    } else {
        match seq[pos] {
            b'G' | b'C' | b'g' | b'c' => true,
            _ => false,
        }
    }
}

fn criterion_benchmark(c: &mut Criterion) {
    let seq = b"acGTCttACGaTA";
    c.bench_function("is_s G", |b| b.iter(|| is_s(seq, 3)));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
