
use criterion::{criterion_group, criterion_main, Criterion};
use rand_python::MersenneTwister;
use std::slice;

fn make_rand() -> MersenneTwister {
    let mut rand = MersenneTwister::new();
    let s: u32 = 63245986;
    rand.init_by_array(slice::from_ref(&s));

    rand
}

pub fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("genrand_res53 implementation");

    let mut rand1 = make_rand();
    group.bench_function("genrand_res53", |b| b.iter(|| rand1.genrand_res53() ));

    // let mut rand2 = make_rand();
    // group.bench_function("genrand_res53_64", |b| b.iter(|| rand2.genrand_res53_64() ));

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);