//! The RLWEenc public-key encryption scheme (Definition 4.6 of the
//! UnifOMR paper) over the negacyclic ring
//! $R_q = \mathbb{F}_q[X]/(X^{N} + 1)$.
//!
//! $$\begin{aligned}
//! \textsf{KeyGen}: \quad & \mathrm{pk} = (\alpha,\; \alpha s + x) \\
//! \textsf{Enc}(\vec m): \quad & t = \sum_i \tfrac{q}{2} m_i X^i, \quad
//!    a = \alpha y + x', \quad b = \mathrm{pk}[1]\, y + t + x'' \\
//! \textsf{Dec}: \quad & d_i = b_i - (a s)_i, \quad
//!    m_i = \mathbb 1[\,|d_i| > r\,]
//! \end{aligned}$$
//!
//! with $s, y \leftarrow \mathcal{D}$ (ternary, weight $h$), $x, x',
//! x'' \leftarrow \chi_\sigma$, $\alpha$ uniform. The ciphertext
//! transmits the ring element $a$ in full but truncates $b$ to the
//! $\ell$ payload coefficients, matching the type
//! $(a, \vec b) \in R_q \times \mathbb{Z}_q^{\ell}$.

use rand::Rng;

use crate::error::{
    sample_error_ring,
    sample_secret_ring,
    sample_uniform_ring,
};
use crate::field_alg::{Fp, PolyFp};
use crate::param;

/// A secret key $s \leftarrow \mathcal{D}$.
#[derive(Clone, Debug)]
pub struct SecretKey<const P: u64, const N: usize>(PolyFp<P, N>);

/// A public key $\mathrm{pk} = (\alpha, \alpha s + x)$.
#[derive(Clone, Debug)]
pub struct PublicKey<const P: u64, const N: usize> {
    /// The public base $\alpha = \mathrm{pk}[0]$.
    pub alpha: PolyFp<P, N>,
    /// The noisy product $\mathrm{pk}[1] = \alpha s + x$.
    pub beta: PolyFp<P, N>,
}

/// A ciphertext $(a, \vec b) \in R_q \times \mathbb{Z}_q^{\ell}$: the
/// full ring element $a$, plus only the $\ell$ payload coefficients of
/// $b$.
#[derive(Clone, Debug)]
pub struct Ciphertext<const P: u64, const N: usize> {
    pub a: PolyFp<P, N>,
    pub b: Vec<Fp<P>>,
}

impl<const P: u64, const N: usize> SecretKey<P, N> {
    pub fn new(s: PolyFp<P, N>) -> Self {
        SecretKey(s)
    }

    pub fn as_poly(&self) -> &PolyFp<P, N> {
        &self.0
    }
}

/// $\textsf{KeyGen}$: draw $s \leftarrow \mathcal{D}$ (weight `H` over
/// the first `active` positions), $\alpha$ uniform, $x \leftarrow
/// \chi_\sigma$, and output $\mathrm{pk} = (\alpha, \alpha s + x)$.
pub fn keygen<const P: u64, const N: usize, const H: usize, R: Rng>(
    rng: &mut R,
    active: usize,
) -> (SecretKey<P, N>, PublicKey<P, N>) {
    let s = sample_secret_ring::<P, N, H, R>(rng, active);
    let alpha = sample_uniform_ring::<P, N, R>(rng);
    let x = sample_error_ring::<P, N, R>(rng);
    let sk = SecretKey(s.clone());
    let pk = PublicKey { beta: &(&alpha * &s) + &x, alpha };
    (sk, pk)
}

/// Deterministic encryption core: given the ephemeral key $y$ and fresh
/// noises $x', x''$, compute the ciphertext of the payload bits
/// $\vec m \in \{0,1\}^{\ell}$ (with the encoding
/// $t = \sum_i \frac{q}{2} m_i X^i$).
pub fn encrypt_with<const P: u64, const N: usize>(
    pk: &PublicKey<P, N>,
    y: &PolyFp<P, N>,
    x1: &PolyFp<P, N>,
    x2: &PolyFp<P, N>,
    message: &[bool],
) -> Ciphertext<P, N> {
    debug_assert!(message.len() <= N);
    let t = PolyFp::<P, N>::from_i64_coeffs(
        &message.iter().map(|m| if *m { (P / 2) as i64 } else { 0 }).collect::<Vec<_>>(),
    );
    let a = &(&pk.alpha * y) + x1;
    let b_full = &(&(&pk.beta * y) + &t) + x2;
    let b = (0..message.len()).map(|i| b_full.coeff(i)).collect();
    Ciphertext { a, b }
}

/// $\textsf{Enc}$: draw $y \leftarrow \mathcal{D}$ (weight `H` over the
/// first `active` positions) and $x', x'' \leftarrow \chi_\sigma$, then
/// run [`encrypt_with`].
pub fn encrypt<const P: u64, const N: usize, const H: usize, R: Rng>(
    rng: &mut R,
    pk: &PublicKey<P, N>,
    active: usize,
    message: &[bool],
) -> Ciphertext<P, N> {
    let y = sample_secret_ring::<P, N, H, R>(rng, active);
    let x1 = sample_error_ring::<P, N, R>(rng);
    let x2 = sample_error_ring::<P, N, R>(rng);
    encrypt_with(pk, &y, &x1, &x2, message)
}

/// $\textsf{Dec}$: compute $d_i = b_i - (a s)_i$ for each payload
/// position and decode $m_i = \mathbb 1[|d_i| > r]$ (with $d_i$ read in
/// $(-q/2, q/2]$).
pub fn decrypt<const P: u64, const N: usize>(
    sk: &SecretKey<P, N>,
    ct: &Ciphertext<P, N>,
    r: i64,
) -> Vec<bool> {
    let as_poly = &ct.a * sk.as_poly();
    ct.b
        .iter()
        .enumerate()
        .map(|(i, b)| (*b - as_poly.coeff(i)).lift_centered().abs() > r)
        .collect()
}

/// The centered decryption values $d_i = b_i - (a s)_i$ (message
/// encoding still included), for noise diagnostics.
pub fn decryption_values<const P: u64, const N: usize>(
    sk: &SecretKey<P, N>,
    ct: &Ciphertext<P, N>,
) -> Vec<i64> {
    let as_poly = &ct.a * sk.as_poly();
    ct.b.iter().enumerate().map(|(i, b)| (*b - as_poly.coeff(i)).lift_centered()).collect()
}

/// RLWEenc secret key at Parameter Set 1.
pub type SkParam1 = SecretKey<{ param::Q }, { param::N_RING }>;

/// RLWEenc public key at Parameter Set 1.
pub type PkParam1 = PublicKey<{ param::Q }, { param::N_RING }>;

/// RLWEenc ciphertext at Parameter Set 1.
pub type CtParam1 = Ciphertext<{ param::Q }, { param::N_RING }>;

/// $\textsf{KeyGen}$ at Parameter Set 1 ($n = 900$ active positions,
/// $h = 80$).
pub fn keygen_param1<R: Rng>(rng: &mut R) -> (SkParam1, PkParam1) {
    keygen::<{ param::Q }, { param::N_RING }, { param::H }, R>(rng, param::N)
}

/// $\textsf{Enc}$ at Parameter Set 1.
pub fn encrypt_param1<R: Rng>(rng: &mut R, pk: &PkParam1, message: &[bool]) -> CtParam1 {
    encrypt::<{ param::Q }, { param::N_RING }, { param::H }, R>(rng, pk, param::N, message)
}

/// $\textsf{Dec}$ at Parameter Set 1 ($r = 48$).
pub fn decrypt_param1(sk: &SkParam1, ct: &CtParam1) -> Vec<bool> {
    decrypt(sk, ct, param::R)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::ThreadRng;

    type ToyRing = PolyFp<17, 4>;
    type ToySk = SecretKey<17, 4>;
    type ToyPk = PublicKey<17, 4>;

    #[test]
    fn worked_toy_keygen_enc_dec() {
        let alpha = ToyRing::from_i64_coeffs(&[3, 5, -2, 1]);
        let s = ToyRing::from_i64_coeffs(&[1, 0, 0, -1]);
        let x = ToyRing::from_i64_coeffs(&[1, -1]);
        let pk = PublicKey { beta: &(&alpha * &s) + &x, alpha: alpha.clone() };
        assert_eq!(pk.beta.coeff_centered(0), -8);

        let y = ToyRing::from_i64_coeffs(&[-1, 0, 1, 0]);
        let x1 = ToyRing::from_i64_coeffs(&[1]);
        let x2 = ToyRing::from_i64_coeffs(&[-1]);
        let ct = encrypt_with(&pk, &y, &x1, &x2, &[true]);
        assert_eq!(ct.a.max_centered_abs(), 6);
        assert_eq!(ct.b[0].lift_centered(), -1);

        let sk: ToySk = SecretKey(s);
        let d = decryption_values(&sk, &ct);
        assert_eq!(d, vec![5]);
        assert_eq!(decrypt(&sk, &ct, 3), vec![true]);
    }

    #[test]
    fn toy_round_trip() {
        let mut rng = ThreadRng::default();
        let mut correct = 0;
        for _ in 0..200 {
            let (sk, pk): (ToySk, ToyPk) = keygen::<17, 4, 2, _>(&mut rng, 4);
            let m = vec![rng.gen_bool(0.5)];
            let ct = encrypt::<17, 4, 2, _>(&mut rng, &pk, 4, &m);
            if decrypt(&sk, &ct, 3) == m {
                correct += 1;
            }
        }
        // the toy budget: noise std sqrt(2h+1)sigma ~ 1.34 against r = 3
        // leaves a ~2% false-negative tail, like script/packing.sage's trials
        assert!(correct >= 190);
    }

    #[test]
    fn wrong_key_outputs_one() {
        // with a wrong key, a*s' is uniform over Z_q, so d lands at a
        // uniform point and |d| > r with probability 1 - (2r+1)/q; the
        // toy ring makes that only 10/17, so test at Parameter Set 1
        // where it holds with probability 1 - 97/4169729
        let mut rng = ThreadRng::default();
        let (sk, pk) = keygen_param1(&mut rng);
        let (wrong_sk, _) = keygen_param1(&mut rng);
        let ct = encrypt_param1(&mut rng, &pk, &[false]);
        let d = decryption_values(&wrong_sk, &ct);
        assert!(d[0].abs() > param::R);
        assert_eq!(decrypt_param1(&wrong_sk, &ct), vec![true]);
        assert_eq!(decrypt_param1(&sk, &ct), vec![false]);
    }

    #[test]
    fn param1_round_trip() {
        let mut rng = ThreadRng::default();
        let (sk, pk) = keygen_param1(&mut rng);
        for m in [vec![false], vec![true]] {
            let ct = encrypt_param1(&mut rng, &pk, &m);
            let d = decryption_values(&sk, &ct);
            // d = (q/2)m + noise with noise std sqrt(2h+1)*sigma ~ 7.6;
            // the centered lift maps q/2 + n to -q/2 + n when n > 0, so
            // measure the distance to the nearer of +-q/2 (or to 0)
            let q2 = (param::Q / 2) as i64;
            let dist = if m[0] {
                (d[0] - q2).abs().min((d[0] + q2).abs())
            } else {
                d[0].abs()
            };
            assert!(dist < 40);
            assert_eq!(decrypt_param1(&sk, &ct), m);
        }
    }

    #[test]
    fn ciphertext_b_is_truncated() {
        let mut rng = ThreadRng::default();
        let (_, pk) = keygen_param1(&mut rng);
        let ct = encrypt_param1(&mut rng, &pk, &[true]);
        assert_eq!(ct.b.len(), param::ELL);
    }
}
