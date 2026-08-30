//! Parameter Set 1 of UnifOMR (Table 1, Section 7.1 of the paper).
//!
//! Every constant carries its math symbol from the paper. Constants that
//! are *derived* rather than freely chosen are re-verified by unit tests
//! at the bottom of this file.

/// RLWEenc secret key dimension ($n$).
///
/// Zero-padded up to [`N_RING`] per footnote 10 of the paper.
pub const N: usize = 900;

/// Ring dimension, $n$ padded to the next power of two ($n'$).
///
/// $R_q = \mathbb{Z}_q[X]/(X^{n'} + 1)$.
pub const N_RING: usize = 1024;

/// Cyclotomic index ($m = 2n'$), so $\Phi_m(X) = X^{n'} + 1$.
pub const M: usize = 2 * N_RING;

/// RLWEenc ciphertext modulus ($q$), a prime.
///
/// $\log_2 q \approx 22$; an element of $R_q$ costs
/// $n' \cdot 22$ bits = 2816 bytes packed.
pub const Q: u64 = 4169729;

/// RLWEenc error standard deviation ($\sigma$), per coefficient of
/// $\chi_\sigma = D_{\mathbb{Z}, \sigma}$.
pub const SIGMA: f64 = 0.6;

/// RLWEenc secret key Hamming weight ($h$).
pub const H: usize = 80;

/// RLWEenc error range bound ($r$).
///
/// Derived: the minimum integer with
/// $\operatorname{erfc}\!\big(r / (\sqrt2\,\sqrt{2h+1}\,\sigma)\big)
/// \le \epsilon_n / \ell$ is 47; the table's 48 adds safety margin.
pub const R: i64 = 48;

/// Bits encrypted per RLWEenc ciphertext ($\ell$).
pub const ELL: usize = 1;

/// False negative rate ($\epsilon_n = 2^{-30}$).
pub const EPS_N: f64 = 9.313_225_746_154_785e-10;

/// False positive rate ($\epsilon_p = 2^{-15}$).
pub const EPS_P: f64 = 3.051_757_812_5e-05;

/// BFV plaintext modulus ($t$).
///
/// Algorithm 1 sets $t = q$, doubling the RLWEenc modulus as the BFV
/// plaintext modulus so partial decryption lands directly in
/// $\mathbb{Z}_t$.
pub const T: u64 = Q;

/// BFV ring dimension and SIMD slot count ($D$).
///
/// Power of two for the NTT; $D \ge n$ per Algorithm 1; $t \equiv 1
/// \pmod{2D}$ for the slot decomposition.
pub const D_BFV: usize = 2048;

/// BFV initial ciphertext modulus ($Q \approx 2^{60}$).
pub const Q_BFV: u64 = 1 << 60;

/// BFV secret key Hamming weight ($h_{\mathrm{BFV}}$).
pub const H_BFV: usize = 400;

/// BFV error standard deviation ($\sigma_{\mathrm{BFV}}$).
pub const SIGMA_BFV: f64 = 3.2;

/// BFV ciphertext modulus after modulus switching ($Q'$), a prime with
/// $Q' = 2^{26} - 2^{12} + 1$.
pub const Q_PRIME: u64 = 67104769;

/// Final error bound after merging via modulus switching ($r'$).
///
/// Derived (Algorithm 1):
/// $r' = (Q'/t)\,r + \mathrm{eModSW} + \lceil Q'/t \rceil$ with
/// $\mathrm{eModSW} \le h_{\mathrm{BFV}}/2$.
pub const R_PRIME: i64 = 992;

/// Standard deviation of the accumulated RLWEenc decryption noise,
/// $\sqrt{2h+1}\,\sigma$ (Definition 4.6).
pub fn effective_noise_std() -> f64 {
    ((2 * H + 1) as f64).sqrt() * SIGMA
}

/// Probability that the accumulated decryption noise exceeds `bound`,
/// i.e. $\Pr[|\text{noise}| > r] =
/// \operatorname{erfc}\!\big(r / (\sqrt2\,\sqrt{2h+1}\,\sigma)\big)$
/// for Gaussian noise of standard deviation
/// $\sqrt{2h+1}\,\sigma$ (Definition 4.6's $1 - \operatorname{erf}$
/// condition).
///
/// Uses [`erfc`], a Chebyshev fit with near machine-precision relative
/// accuracy even deep in the tail.
pub fn noise_tail_probability(bound: i64) -> f64 {
    erfc(bound as f64 / (std::f64::consts::SQRT_2 * effective_noise_std()))
}

/// Complementary error function
/// $\operatorname{erfc}(x) = \frac{2}{\sqrt\pi}\int_x^\infty e^{-t^2}dt$,
/// valid for $x \ge 0$ with relative accuracy $\sim 10^{-15}$ even deep
/// in the tail.
///
/// Two branches: the Maclaurin series of $\operatorname{erf}$ for small
/// $x$, and the classical continued fraction (Abramowitz & Stegun
/// 7.1.14, via the modified Lentz algorithm) for large $x$, where the
/// series would cancel catastrophically.
pub fn erfc(x: f64) -> f64 {
    if x < 2.0 {
        1.0 - erf_series(x)
    } else {
        erfc_continued_fraction(x)
    }
}

fn erf_series(x: f64) -> f64 {
    let mut sum = 0.0f64;
    let mut term = x;
    let mut n = 0u32;
    let xx = -x * x;
    while term.abs() > 1e-18 * sum.abs().max(1e-300) && n < 200 {
        sum += term / (2 * n + 1) as f64;
        n += 1;
        term *= xx / n as f64;
    }
    2.0 / std::f64::consts::PI.sqrt() * sum
}

fn erfc_continued_fraction(z: f64) -> f64 {
    let tiny = 1e-300;
    let mut f = if z == 0.0 { tiny } else { z };
    let mut c = f;
    let mut d = 0.0f64;
    for n in 1..300u32 {
        let a_n = n as f64 / 2.0;
        d = z + a_n * d;
        if d == 0.0 {
            d = tiny;
        }
        c = z + a_n / c;
        if c == 0.0 {
            c = tiny;
        }
        d = 1.0 / d;
        let delta = c * d;
        f *= delta;
        if (delta - 1.0).abs() < 1e-17 {
            break;
        }
    }
    (-z * z).exp() / std::f64::consts::PI.sqrt() / f
}

/// Deterministic Miller-Rabin primality test for `u64` (valid below
/// $3.3 \times 10^{24}$ with these bases).
#[cfg(test)]
fn is_prime_u64(n: u64) -> bool {
    const BASES: [u64; 12] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
    if n < 2 {
        return false;
    }
    for p in BASES {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let mut d = n - 1;
    let mut s = 0u32;
    while d.is_multiple_of(2) {
        d /= 2;
        s += 1;
    }
    'witness: for a in BASES {
        let mut x = modpow_u64(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 0..s - 1 {
            x = ((x as u128 * x as u128) % n as u128) as u64;
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

#[cfg(test)]
fn modpow_u64(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    let mut result = 1u64;
    while exp > 0 {
        if !exp.is_multiple_of(2) {
            result = ((result as u128 * base as u128) % modulus as u128) as u64;
        }
        base = ((base as u128 * base as u128) % modulus as u128) as u64;
        exp >>= 1;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erfc_matches_reference_values() {
        assert!((erfc(0.5) - 0.4795001221869535).abs() < 1e-14);
        assert!((erfc(1.0) - 0.15729920705028513).abs() < 1e-14);
        assert!((erfc(2.0) - 0.004677734981047265).abs() < 1e-16);
        let deep = erfc(4.4593);
        assert!(deep > 0.0 && (deep - 2.8559567523976677e-10).abs() / 2.8559567523976677e-10 < 1e-12);
        assert!((erfc(5.0) - 1.5374597944280351e-12).abs() / 1.5374597944280351e-12 < 1e-12);
    }

    #[test]
    fn q_is_prime_and_ntt_friendly() {
        assert!(is_prime_u64(Q));
        assert_eq!((Q - 1) % (2 * M as u64), 0);
    }

    #[test]
    fn q_prime_shape_and_primality() {
        assert_eq!(Q_PRIME, (1 << 26) - (1 << 12) + 1);
        assert!(is_prime_u64(Q_PRIME));
        assert_eq!((Q_PRIME - 1) % (2 * D_BFV as u64), 0);
    }

    #[test]
    fn effective_noise_std_value() {
        let s = effective_noise_std();
        assert!((s - 7.613146524269712).abs() < 1e-9);
    }

    #[test]
    fn r_meets_the_error_budget() {
        assert!(noise_tail_probability(R) <= EPS_N / ELL as f64);
        let mut r_min = 1i64;
        while noise_tail_probability(r_min) > EPS_N {
            r_min += 1;
        }
        assert_eq!(r_min, 47);
        assert!(R >= r_min);
    }

    #[test]
    fn r_prime_matches_its_derivation() {
        let q_prime_over_t = Q_PRIME as f64 / T as f64;
        let calc = q_prime_over_t * R as f64 + (H_BFV / 2) as f64 + q_prime_over_t.ceil();
        assert!((R_PRIME as f64 - calc).abs() <= 3.0);
    }

    #[test]
    fn eps_p_budget_reproduces() {
        let k = (T * R_PRIME as u64) / Q_PRIME;
        assert_eq!(k, 61);
        let eps_p = (2 * k + 1) as f64 / T as f64;
        assert!(eps_p <= EPS_P);
    }
}
