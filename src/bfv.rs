//! A minimal BFV homomorphic-encryption layer — exactly the operations
//! UnifOMR's detection circuit (Algorithm 1) needs and nothing more:
//! fresh encryption of *constant* plaintexts, plaintext-by-ciphertext
//! multiplication, ciphertext addition/negation, plaintext addition,
//! and modulus switching. Ciphertext-ciphertext multiplication,
//! relinearization and the SIMD slot encoding are deliberately absent:
//! the circuit is a linear form in encrypted constants,
//! $\sum_w p_w(X) \cdot \mathrm{ct}_{s_w}$, so no two ciphertexts are
//! ever multiplied, and the paper packs its scalars into *coefficients*
//! of $p_w$, not slots.
//!
//! $$\begin{aligned}
//! \textsf{KeyGen}: \quad & \mathrm{pk} = (a,\; a s + e) \\
//! \textsf{Enc}(m): \quad & c_0 = \mathrm{pk}[1]\, u + e_1 + \lfloor P/t \rfloor m, \quad
//!    c_1 = a u + e_2 \\
//! \textsf{Dec}(\mathrm{ct}): \quad & c_0 - c_1 s
//!    \;=\; \lfloor P/t \rfloor m + \underbrace{(e u + e_1 - e_2 s)}_{\text{noise}}
//! \end{aligned}$$
//!
//! with $s, u$ ternary of Hamming weight $h_{\mathrm{BFV}}$, errors
//! $\leftarrow \chi_{\sigma_{\mathrm{BFV}}}$, and plaintext modulus
//! $t = q$ (Algorithm 1 sets the RLWEenc modulus equal to the BFV
//! plaintext modulus), so the RLWEenc partial decryption lands directly
//! in $\mathbb{Z}_t$.
//!
//! # Noise accounting at Parameter Set 1
//!
//! A fresh ciphertext's noise $e u + e_1 - e_2 s$ has per-coefficient
//! std $\sigma_{\mathrm{BFV}} \sqrt{2 h_{\mathrm{BFV}} + 1} \approx 90.6$.
//! One plain-multiplication multiplies that noise by a packed
//! polynomial (coefficients $< t$), growing it by
//! $\sqrt{D} \cdot 90.6 \cdot t/\sqrt{12} \approx 2^{32}$; summing
//! $n = 900$ independent terms leaves $\approx 2^{37}$ before the
//! switch, i.e. $\approx 9$ per coefficient after scaling by
//! $Q'/Q = 2^{34}$. Together with the modulus-switch rounding
//! $\mathrm{eModSW} \le h_{\mathrm{BFV}}/2 = 200$ and the RLWEenc term
//! $(Q'/t) r \le 772$, pertinent clues decrypt well inside the client
//! threshold $r' Q'/t = 15963$ — verified empirically by the tests in
//! [`crate::detection`].
//!
//! # Moduli and arithmetic
//!
//! The ciphertext modulus $Q = 2^{60}$ is *not* prime (and not
//! NTT-friendly); all products go through the structured
//! multiplications of [`crate::field_alg::PolyFp`] —
//! [`PolyFp::mul_ternary`] for the sparse-secret products and
//! [`PolyFp::mul_plain`] for plaintext-by-ciphertext products — both
//! safe at $P = 2^{60}$, $N = 2048$ where generic schoolbook
//! [`PolyFp::mul`] could overflow its `i128` accumulator.

use std::ops::{Add, Neg};

use rand::Rng;

use crate::error::{sample_error_sigma, sample_secret_ring, sample_uniform_ring};
use crate::field_alg::{Fp, PolyFp};
use crate::param;

/// A BFV secret key $s$, ternary of Hamming weight $h_{\mathrm{BFV}}$.
#[derive(Clone, Debug)]
pub struct SecretKey<const P: u64, const N: usize>(PolyFp<P, N>);

/// A BFV public key $\mathrm{pk} = (a,\; a s + e)$.
#[derive(Clone, Debug)]
pub struct PublicKey<const P: u64, const N: usize> {
    /// The uniform base $a = \mathrm{pk}[0]$.
    pub a: PolyFp<P, N>,
    /// The noisy product $\mathrm{pk}[1] = a s + e$.
    pub b: PolyFp<P, N>,
}

/// A BFV ciphertext $(c_0, c_1) \in R_P^2$.
#[derive(Clone, Debug)]
pub struct Ciphertext<const P: u64, const N: usize> {
    pub c0: PolyFp<P, N>,
    pub c1: PolyFp<P, N>,
}

impl<const P: u64, const N: usize> Ciphertext<P, N> {
    /// The encryption of the zero plaintext — the additive identity,
    /// and the accumulator for the detection circuit's homomorphic sum.
    pub fn zero() -> Self {
        Ciphertext {
            c0: PolyFp::zero(),
            c1: PolyFp::zero(),
        }
    }
}

impl<const P: u64, const N: usize> SecretKey<P, N> {
    pub fn new(s: PolyFp<P, N>) -> Self {
        SecretKey(s)
    }

    pub fn as_poly(&self) -> &PolyFp<P, N> {
        &self.0
    }
}

impl<const P: u64, const N: usize> Add for &Ciphertext<P, N> {
    type Output = Ciphertext<P, N>;
    fn add(self, rhs: Self) -> Ciphertext<P, N> {
        Ciphertext {
            c0: &self.c0 + &rhs.c0,
            c1: &self.c1 + &rhs.c1,
        }
    }
}

impl<const P: u64, const N: usize> Neg for &Ciphertext<P, N> {
    type Output = Ciphertext<P, N>;
    fn neg(self) -> Ciphertext<P, N> {
        Ciphertext {
            c0: -&self.c0,
            c1: -&self.c1,
        }
    }
}

/// $\textsf{KeyGen}$: draw a ternary secret $s$ (weight `H` over the
/// first `active` positions), a uniform $a$, error $e \leftarrow
/// \chi_\sigma$, and output $\mathrm{pk} = (a, a s + e)$.
pub fn keygen<const P: u64, const N: usize, const H: usize, R: Rng>(
    rng: &mut R,
    active: usize,
    sigma: f64,
) -> (SecretKey<P, N>, PublicKey<P, N>) {
    let s = sample_secret_ring::<P, N, H, R>(rng, active);
    let a = sample_uniform_ring::<P, N, R>(rng);
    let e = sample_error_sigma::<P, N, R>(rng, sigma);
    let sk = SecretKey(s.clone());
    let pk = PublicKey {
        b: &a.mul_ternary(&s) + &e,
        a,
    };
    (sk, pk)
}

/// Deterministic encryption core: given the ephemeral ternary $u$ and
/// fresh errors $e_1, e_2$, encrypt the *constant* $m$ (degree 0) —
/// the shape the detection key uses, one ciphertext per secret
/// coefficient $s_w \in \{-1, 0, 1\}$.
///
/// $$c_0 = \mathrm{pk}[1] u + e_1 + \lfloor P/t \rfloor m, \quad
///   c_1 = a u + e_2 .$$
pub fn encrypt_with<const P: u64, const N: usize>(
    pk: &PublicKey<P, N>,
    u: &PolyFp<P, N>,
    e1: &PolyFp<P, N>,
    e2: &PolyFp<P, N>,
    t: u64,
    m: i64,
) -> Ciphertext<P, N> {
    let delta = (P / t) as i128;
    let encoded = (delta * m as i128).rem_euclid(P as i128) as i64;
    let c0 = &(&pk.b.mul_ternary(u) + e1) + &PolyFp::<P, N>::from_i64_coeffs(&[encoded]);
    let c1 = &pk.a.mul_ternary(u) + e2;
    Ciphertext { c0, c1 }
}

/// $\textsf{Enc}$ of the constant $m$: draw $u$ ternary (weight `H`)
/// and $e_1, e_2 \leftarrow \chi_\sigma$, then run [`encrypt_with`].
pub fn encrypt<const P: u64, const N: usize, const H: usize, R: Rng>(
    rng: &mut R,
    pk: &PublicKey<P, N>,
    active: usize,
    sigma: f64,
    t: u64,
    m: i64,
) -> Ciphertext<P, N> {
    let u = sample_secret_ring::<P, N, H, R>(rng, active);
    let e1 = sample_error_sigma::<P, N, R>(rng, sigma);
    let e2 = sample_error_sigma::<P, N, R>(rng, sigma);
    encrypt_with(pk, &u, &e1, &e2, t, m)
}

/// Plaintext-by-ciphertext multiplication $(c_0 p,\; c_1 p)$: because
/// the plaintext $p$ multiplies both components, the decryption phase
/// becomes $c_0 p - c_1 p\, s$ — the noise and the message are both
/// scaled by $p$, and the ciphertext stays a fresh 2-component one (no
/// relinearization needed).
pub fn plain_mul<const P: u64, const N: usize>(
    ct: &Ciphertext<P, N>,
    p: &PolyFp<P, N>,
) -> Ciphertext<P, N> {
    Ciphertext {
        c0: ct.c0.mul_plain(p),
        c1: ct.c1.mul_plain(p),
    }
}

/// Add a plaintext at message scale: $c_0 \mathrel{+}= \lfloor P/t
/// \rfloor \cdot p$, leaving $c_1$ alone. Used to fold the packed
/// public $b$-term into the detection circuit's result.
pub fn add_plaintext<const P: u64, const N: usize>(
    ct: &Ciphertext<P, N>,
    p: &PolyFp<P, N>,
    t: u64,
) -> Ciphertext<P, N> {
    let delta = Fp::<P>::new(P / t);
    Ciphertext {
        c0: &ct.c0 + &p.scalar_mul(delta),
        c1: ct.c1.clone(),
    }
}

/// Modulus switching $P \to P_2$ (Algorithm 1's ModSW): scale every
/// coefficient of both components by $P_2/P$ and round,
/// $c' = \lfloor P_2/P \cdot c \rceil \bmod P_2$. Each rounding is off
/// by at most $1/2$, which the switched secret product accumulates to
/// at most $h_{\mathrm{BFV}}/2$ — the $\mathrm{eModSW}$ term of $r'$.
pub fn mod_switch<const P: u64, const P2: u64, const N: usize>(
    ct: &Ciphertext<P, N>,
) -> Ciphertext<P2, N> {
    let switch = |poly: &PolyFp<P, N>| -> PolyFp<P2, N> {
        let mut out = [0u64; N];
        for (o, i) in out.iter_mut().zip(0..N) {
            let c = poly.coeff(i).raw() as u128;
            let scaled = (c * (P2 as u128) + (P as u128) / 2) / (P as u128);
            *o = (scaled % (P2 as u128)) as u64;
        }
        PolyFp::from_raw(out)
    };
    Ciphertext {
        c0: switch(&ct.c0),
        c1: switch(&ct.c1),
    }
}

/// The decryption phase $c_0 - c_1 s \in R_P$: message at scale
/// $\lfloor P/t \rfloor$ plus noise, for noise diagnostics and tests.
pub fn decrypt_phase<const P: u64, const N: usize>(
    sk: &SecretKey<P, N>,
    ct: &Ciphertext<P, N>,
) -> PolyFp<P, N> {
    &ct.c0 - &ct.c1.mul_ternary(sk.as_poly())
}

/// The BFV ring $R_Q = \mathbb{Z}_Q[X]/(X^D + 1)$ at Parameter Set 1
/// ($Q = 2^{60}$, $D = 2048$).
pub type BfvRing = PolyFp<{ param::Q_BFV }, { param::D_BFV }>;

/// The modulus-switched ring $R_{Q'}$ at Parameter Set 1
/// ($Q' = 2^{26} - 2^{12} + 1$).
pub type BfvRingQp = PolyFp<{ param::Q_PRIME }, { param::D_BFV }>;

/// BFV secret key at Parameter Set 1.
pub type SkParam1 = SecretKey<{ param::Q_BFV }, { param::D_BFV }>;

/// BFV public key at Parameter Set 1.
pub type PkParam1 = PublicKey<{ param::Q_BFV }, { param::D_BFV }>;

/// BFV ciphertext at Parameter Set 1.
pub type CtParam1 = Ciphertext<{ param::Q_BFV }, { param::D_BFV }>;

/// Modulus-switched BFV ciphertext at Parameter Set 1 — the server's
/// reply format.
pub type CtQpParam1 = Ciphertext<{ param::Q_PRIME }, { param::D_BFV }>;

/// The message-scale factor $\Delta = \lfloor Q/t \rfloor \approx 2.77
/// \times 10^{14}$ at Parameter Set 1.
pub const DELTA: u64 = param::Q_BFV / param::T;

/// $\textsf{KeyGen}$ at Parameter Set 1 ($h_{\mathrm{BFV}} = 400$,
/// $\sigma_{\mathrm{BFV}} = 3.2$, secret dense over all $D$ positions).
pub fn keygen_param1<R: Rng>(rng: &mut R) -> (SkParam1, PkParam1) {
    keygen::<{ param::Q_BFV }, { param::D_BFV }, { param::H_BFV }, R>(
        rng,
        param::D_BFV,
        param::SIGMA_BFV,
    )
}

/// $\textsf{Enc}$ of the constant $m$ at Parameter Set 1.
pub fn encrypt_param1<R: Rng>(rng: &mut R, pk: &PkParam1, m: i64) -> CtParam1 {
    encrypt::<{ param::Q_BFV }, { param::D_BFV }, { param::H_BFV }, R>(
        rng,
        pk,
        param::D_BFV,
        param::SIGMA_BFV,
        param::T,
        m,
    )
}

/// Modulus switching $Q \to Q'$ at Parameter Set 1.
pub fn mod_switch_param1(ct: &CtParam1) -> CtQpParam1 {
    mod_switch::<{ param::Q_BFV }, { param::Q_PRIME }, { param::D_BFV }>(ct)
}

/// The decryption phase of a modulus-switched reply: $c_0' - c_1' s$
/// in $R_{Q'}$, with the ternary secret re-embedded from $R_Q$ (its
/// $\pm 1$ coefficients are canonical in both rings).
///
/// Coefficient $j$ of the result is $(Q'/t)\, d_j + \text{merged
/// error}$ for clue $j$ of the window — the client range-checks it
/// against $r' Q'/t$ ([`crate::detection::client_decode`]).
pub fn decrypt_phase_param1(sk: &SkParam1, ct: &CtQpParam1) -> BfvRingQp {
    let sk_qp = SecretKey::new(sk.as_poly().lift_to::<{ param::Q_PRIME }>());
    decrypt_phase::<{ param::Q_PRIME }, { param::D_BFV }>(&sk_qp, ct)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::ThreadRng;

    /// Toy BFV: plaintext modulus $t = 17$ (the toy RLWEenc $q$),
    /// ciphertext modulus $P = 4352 = 256 \cdot 17$ so $\Delta = 256$,
    /// switched modulus $P_2 = 257$.
    type ToySk = SecretKey<4352, 4>;
    type ToyPk = PublicKey<4352, 4>;

    #[test]
    fn toy_encrypt_decrypt_phase() {
        let mut rng = ThreadRng::default();
        let (sk, pk): (ToySk, ToyPk) = keygen::<4352, 4, 2, _>(&mut rng, 4, 0.6);
        let ct = encrypt::<4352, 4, 2, _>(&mut rng, &pk, 4, 0.6, 17, 3);
        let v = decrypt_phase(&sk, &ct);
        // noise = e*u + e1 - e2*s, toy std sqrt(2*2+1)*0.6 ~ 1.3
        assert!((v.coeff_centered(0) - 3 * 256).abs() <= 8);
        for i in 1..4 {
            assert!(v.coeff_centered(i).abs() <= 8);
        }
    }

    #[test]
    fn toy_plain_mul_and_plaintext_add() {
        let mut rng = ThreadRng::default();
        let (sk, pk): (ToySk, ToyPk) = keygen::<4352, 4, 2, _>(&mut rng, 4, 0.6);
        let ct = encrypt::<4352, 4, 2, _>(&mut rng, &pk, 4, 0.6, 17, 2);
        let p = PolyFp::<4352, 4>::from_i64_coeffs(&[1, -3, 0, 4]);
        let b = PolyFp::<4352, 4>::from_i64_coeffs(&[5, 1, 0, -2]);

        let scaled = plain_mul(&ct, &p);
        let with_b = add_plaintext(&scaled, &b, 17);
        let v = decrypt_phase(&sk, &with_b);
        for i in 0..4 {
            let expected = 256 * (2 * p.coeff_centered(i) + b.coeff_centered(i));
            assert!((v.coeff_centered(i) - expected).abs() <= 30, "coeff {i}");
        }
    }

    #[test]
    fn toy_modulus_switch() {
        let mut rng = ThreadRng::default();
        let (sk, pk): (ToySk, ToyPk) = keygen::<4352, 4, 2, _>(&mut rng, 4, 0.6);
        let ct = encrypt::<4352, 4, 2, _>(&mut rng, &pk, 4, 0.6, 17, 5);
        let switched: Ciphertext<257, 4> = mod_switch::<4352, 257, 4>(&ct);
        // re-embed the ternary secret at 257 and decrypt there
        let sk_qp = SecretKey::new(sk.as_poly().lift_to::<257>());
        let v = decrypt_phase::<257, 4>(&sk_qp, &switched);
        // (257/4352) * 256 * 5 = 75.53, plus switch rounding <= h/2 = 1
        // and the scaled BFV noise (< 0.1)
        assert!((v.coeff_centered(0) as f64 - 257.0 / 17.0 * 5.0).abs() <= 2.0);
        assert!(v.max_centered_abs() <= 80);
    }

    #[test]
    fn param1_encrypt_round_trip() {
        let mut rng = ThreadRng::default();
        let (sk, pk) = keygen_param1(&mut rng);
        for m in [-1i64, 0, 1] {
            let ct = encrypt_param1(&mut rng, &pk, m);
            let v = decrypt_phase(&sk, &ct);
            // fresh noise std ~ 90.6; 6 sigma ~ 550, far below Delta/2
            let err = (v.coeff_centered(0) as i128 - m as i128 * DELTA as i128).abs();
            assert!(err < 600, "m = {m}, err = {err}");
            for i in 1..param::D_BFV {
                assert!(v.coeff_centered(i).abs() < 600, "coeff {i}");
            }
        }
    }

    #[test]
    fn param1_plain_mul_noise_stays_small() {
        let mut rng = ThreadRng::default();
        let (sk, pk) = keygen_param1(&mut rng);
        let ct = encrypt_param1(&mut rng, &pk, 1);
        // a packed-looking plaintext: coefficients < q = t
        let p = PolyFp::<{ param::Q_BFV }, { param::D_BFV }>::from_raw(std::array::from_fn(|_| {
            rng.gen_range(0..param::Q)
        }));
        let scaled = plain_mul(&ct, &p);
        let v = decrypt_phase(&sk, &scaled);
        for i in 0..param::D_BFV {
            let expected =
                (DELTA as i128 * p.coeff(i).raw() as i128).rem_euclid(param::Q_BFV as i128);
            let expected_centered = if expected > param::Q_BFV as i128 / 2 {
                expected - param::Q_BFV as i128
            } else {
                expected
            };
            let err = (v.coeff_centered(i) as i128 - expected_centered).abs();
            // one plain-mult noise: sqrt(D)*90.6*t/sqrt(12) ~ 4.9e9, 5 sigma ~ 2.5e10
            assert!(err < 5 * (1 << 33), "coeff {i}, err {err}");
        }
    }
}
