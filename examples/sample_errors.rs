//! Error sampling for UnifOMR RLWEenc, Parameter Set 1 — the Rust mirror
//! of `script/sample_errors.sage`. Run with:
//!
//! ```sh
//! cargo run --release --example sample_errors
//! ```

use detect2::error::{sample_error, sample_secret, sample_uniform};
use detect2::field_alg::Rq;
use detect2::param;
use rand::rngs::ThreadRng;

fn main() {
    let mut rng = ThreadRng::default();

    println!(
        "ring                  : Z_{}[X]/(X^{} + 1)",
        param::Q,
        param::N_RING
    );
    println!(
        "defining relation     : Phi_{}(X) = X^{} + 1",
        param::M,
        param::N_RING
    );

    let num_samples = 2000;
    let mut flat: Vec<i64> = Vec::with_capacity(num_samples * param::N_RING);
    for _ in 0..num_samples {
        let e = sample_error(&mut rng);
        for i in 0..param::N_RING {
            flat.push(e.coeff_centered(i));
        }
    }

    let m = flat.len() as f64;
    let mean = flat.iter().sum::<i64>() as f64 / m;
    let var = flat.iter().map(|k| (*k as f64 - mean).powi(2)).sum::<f64>() / m;
    let std = var.sqrt();
    let count = |v: i64| flat.iter().filter(|k| **k == v).count() as f64;
    let maxabs = flat.iter().map(|k| k.abs()).max().unwrap();

    println!("samples               : {} coefficients", flat.len());
    println!("empirical mean        : {mean:+.5}");
    println!(
        "empirical std         : {std:.5}   (target {})",
        param::SIGMA
    );
    println!("empirical P(0)        : {:.5}", count(0) / m);
    println!("empirical P(+1) = P(-1): {:.5}", count(1) / m);
    println!(
        "max |coeff|           : {maxabs}          (bound r = {})",
        param::R
    );
    assert!(maxabs < param::R);

    let eff = param::effective_noise_std();
    let tail = param::noise_tail_probability(param::R);
    println!("effective noise std   : {eff:.4} = sqrt(2h+1)*sigma");
    println!(
        "P(|noise| > r)        : {tail:.3e}   (need <= eps_n/ell = {:.3e})",
        param::EPS_N / param::ELL as f64
    );

    let s: Rq = sample_secret(&mut rng);
    let alpha = sample_uniform(&mut rng);
    let e = sample_error(&mut rng);
    let pk1 = &(&alpha * &s) + &e;
    let noise = &pk1 - &(&alpha * &s);
    println!("demo pk noise max|c|  : {}", noise.max_centered_abs());
}
