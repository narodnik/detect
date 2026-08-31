//! The random distributions of RLWEenc (Definition 4.6):
//!
//! - $\chi_\sigma$: the error distribution, i.i.d. discrete Gaussian
//!   coefficients $D_{\mathbb{Z}, \sigma}$ with
//!   $\Pr[k] \propto e^{-k^2/2\sigma^2}$;
//! - $\mathcal{D}$: the secret distribution, ternary
//!   $\{-1, 0, 1\}^{n}$ of Hamming weight $h$;
//! - uniform over $R_q$ for the public $\alpha$.

use std::sync::{Mutex, OnceLock};

use rand::Rng;
use rand::rngs::ThreadRng;
use rand::seq::SliceRandom;

use crate::field_alg::{PolyFp, Rq};
use crate::param;

/// Exact cumulative distribution table for one $D_{\mathbb{Z}, \sigma}$
/// draw: the support with its (unnormalized) Gaussian weights and their
/// running sum for inversion sampling.
struct DiscreteGaussian {
    support: Vec<i64>,
    cumsum: Vec<f64>,
    total: f64,
}

fn build_discrete_gaussian(sigma: f64) -> DiscreteGaussian {
    let cutoff = (12.0 * sigma).ceil() as i64;
    let support: Vec<i64> = (-cutoff..=cutoff).collect();
    let mut cumsum = Vec::with_capacity(support.len());
    let mut acc = 0.0f64;
    for k in &support {
        acc += f64::exp(-((*k as f64) * (*k as f64)) / (2.0 * sigma * sigma));
        cumsum.push(acc);
    }
    DiscreteGaussian {
        support,
        cumsum,
        total: acc,
    }
}

/// The table for `sigma`, built on first use and cached for the
/// process lifetime (tables are keyed by the bit pattern of `sigma`).
fn discrete_gaussian(sigma: f64) -> &'static DiscreteGaussian {
    static CACHE: OnceLock<Mutex<Vec<(u64, &'static DiscreteGaussian)>>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(Vec::new()));
    let key = sigma.to_bits();
    let mut tables = cache.lock().unwrap();
    if let Some((_, table)) = tables.iter().find(|(k, _)| *k == key) {
        return table;
    }
    let table: &'static DiscreteGaussian = Box::leak(Box::new(build_discrete_gaussian(sigma)));
    tables.push((key, table));
    table
}

/// One coefficient draw from $\chi_\sigma$, i.e. $k \leftarrow
/// D_{\mathbb{Z}, \sigma}$: a discrete Gaussian over the integers with
/// standard deviation $\sigma = 0.6$.
///
/// The support is cut at $\pm\lceil 12\sigma \rceil = \pm 8$; the mass
/// beyond is $e^{-72} \sim 10^{-31}$.
pub fn chi_coeff<R: Rng>(rng: &mut R) -> i64 {
    chi_coeff_sigma(rng, param::SIGMA)
}

/// Generic-$\sigma$ version of [`chi_coeff`], e.g. the BFV error
/// distribution $\chi_{\sigma_{\mathrm{BFV}}}$ with $\sigma_{\mathrm{BFV}} = 3.2$
/// (support cut at $\pm 38$).
pub fn chi_coeff_sigma<R: Rng>(rng: &mut R, sigma: f64) -> i64 {
    let table = discrete_gaussian(sigma);
    let u = rng.gen_range(0.0..table.total);
    let idx = table.cumsum.partition_point(|c| *c <= u);
    table.support[idx.min(table.support.len() - 1)]
}

/// The probability mass function of $\chi_\sigma$ at `k`, for tests and
/// diagnostics.
pub fn chi_probability(k: i64) -> f64 {
    let cutoff = (12.0 * param::SIGMA).ceil() as i64;
    if k.abs() > cutoff {
        return 0.0;
    }
    let point = f64::exp(-((k * k) as f64) / (2.0 * param::SIGMA * param::SIGMA));
    let total: f64 = (-cutoff..=cutoff)
        .map(|j| f64::exp(-((j * j) as f64) / (2.0 * param::SIGMA * param::SIGMA)))
        .sum();
    point / total
}

/// An error polynomial $x \leftarrow \chi_\sigma$ in $R_q$: $n'$
/// i.i.d. discrete Gaussian coefficients.
pub fn sample_error<R: Rng>(rng: &mut R) -> Rq {
    sample_error_ring(rng)
}

/// Generic-ring version of [`sample_error`].
pub fn sample_error_ring<const P: u64, const N: usize, R: Rng>(rng: &mut R) -> PolyFp<P, N> {
    sample_error_sigma::<P, N, R>(rng, param::SIGMA)
}

/// Generic-ring error sampling at an arbitrary $\sigma$, e.g. the BFV
/// layer's $\chi_{\sigma_{\mathrm{BFV}}}$.
pub fn sample_error_sigma<const P: u64, const N: usize, R: Rng>(
    rng: &mut R,
    sigma: f64,
) -> PolyFp<P, N> {
    let coeffs: Vec<i64> = (0..N).map(|_| chi_coeff_sigma(rng, sigma)).collect();
    PolyFp::from_i64_coeffs(&coeffs)
}

/// A secret key $s \leftarrow \mathcal{D}$: exactly $h$ of the first $n$
/// coefficients are $\pm1$, the rest zero, zero-padded to $n'$.
pub fn sample_secret<R: Rng>(rng: &mut R) -> Rq {
    sample_secret_ring::<{ param::Q }, { param::N_RING }, { param::H }, R>(rng, param::N)
}

/// Generic-ring version of [`sample_secret`]: Hamming weight `H` over
/// the first `active` positions (zero-padded up to `N`).
pub fn sample_secret_ring<const P: u64, const N: usize, const H: usize, R: Rng>(
    rng: &mut R,
    active: usize,
) -> PolyFp<P, N> {
    debug_assert!(H <= active && active <= N);
    let mut positions: Vec<usize> = (0..active).collect();
    positions.partial_shuffle(rng, H);
    let mut coeffs = vec![0i64; N];
    for pos in positions.iter().take(H) {
        coeffs[*pos] = if rng.gen_bool(0.5) { 1 } else { -1 };
    }
    PolyFp::from_i64_coeffs(&coeffs)
}

/// A uniform ring element, e.g. the public $\alpha \leftarrow\$ R_q$.
pub fn sample_uniform<R: Rng>(rng: &mut R) -> Rq {
    sample_uniform_ring(rng)
}

/// Generic-ring version of [`sample_uniform`].
pub fn sample_uniform_ring<const P: u64, const N: usize, R: Rng>(rng: &mut R) -> PolyFp<P, N> {
    let mut coeffs = [0u64; N];
    for c in coeffs.iter_mut() {
        *c = rng.gen_range(0..P);
    }
    PolyFp::from_raw(coeffs)
}

/// Convenience wrapper drawing from the thread-local RNG.
pub fn sample_error_thread_rng() -> Rq {
    let mut rng = ThreadRng::default();
    sample_error(&mut rng)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chi_probability_matches_theory() {
        assert!((chi_probability(0) - 0.6638150448693395).abs() < 1e-12);
        assert!((chi_probability(1) - 0.16552374765776978).abs() < 1e-12);
        assert!((chi_probability(2) - 0.002566255950845424).abs() < 1e-12);
        assert_eq!(chi_probability(9), 0.0);
        let total: f64 = discrete_gaussian(param::SIGMA)
            .support
            .iter()
            .map(|k| chi_probability(*k))
            .sum();
        assert!((total - 1.0).abs() < 1e-12);
    }

    #[test]
    fn chi_statistics_converge() {
        let mut rng = ThreadRng::default();
        let n = 200_000;
        let draws: Vec<i64> = (0..n).map(|_| chi_coeff(&mut rng)).collect();
        let mean = draws.iter().sum::<i64>() as f64 / n as f64;
        let var = draws
            .iter()
            .map(|k| (*k as f64 - mean) * (*k as f64 - mean))
            .sum::<f64>()
            / n as f64;
        let std = var.sqrt();
        assert!(mean.abs() < 0.01);
        assert!((std - 0.5931).abs() < 0.005);
        let p0 = draws.iter().filter(|k| **k == 0).count() as f64 / n as f64;
        let p1 = draws.iter().filter(|k| **k == 1).count() as f64 / n as f64;
        assert!((p0 - 0.6638).abs() < 0.003);
        assert!((p1 - 0.1655).abs() < 0.003);
        assert!(draws.iter().all(|k| k.abs() <= 5));
    }

    #[test]
    fn secret_is_ternary_with_fixed_weight() {
        let mut rng = ThreadRng::default();
        let s = sample_secret(&mut rng);
        let mut weight = 0;
        let mut padded_zero = true;
        for i in 0..param::N_RING {
            let c = s.coeff_centered(i);
            assert!(c == 0 || c == -1 || c == 1);
            if c != 0 {
                weight += 1;
            }
            if i >= param::N && c != 0 {
                padded_zero = false;
            }
        }
        assert_eq!(weight, param::H);
        assert!(padded_zero);
    }

    #[test]
    fn keygen_noise_round_trip() {
        let mut rng = ThreadRng::default();
        let s = sample_secret(&mut rng);
        let alpha = sample_uniform(&mut rng);
        let e = sample_error(&mut rng);
        let pk1 = &(&alpha * &s) + &e;
        let noise = &pk1 - &(&alpha * &s);
        let maxabs = noise.max_centered_abs();
        assert!(maxabs <= 3);
    }
}
