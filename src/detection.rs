//! The detector side of UnifOMR (UnifOMD): SIMD partial decryption of
//! many clues via coefficient packing, plus the pertinence decision.
//!
//! Mirrors §3.4 of `crypto.md` / `script/packing.sage`. The detector's
//! plaintext ring is $\mathbb{Z}_t[X]/(X^{D}+1)$ with $t = q$ and
//! $D = 2048$ slots; the homomorphic layer is *simulated on plaintexts*
//! (constant $\times$ polynomial), which preserves the algorithm's cost
//! shape: $n$ scalar multiplications serve $D$ clues at once.
//!
//! Pipeline for decoded bit $k$ and clues $j = 0, \dots$:
//!
//! 1. **Unfold.** The $k$-th coefficient of each negacyclic product is a
//!    linear form in the secret coefficients,
//!    $(a^{(j)} s)_k = \sum_i A_j[k][i]\, s_i$ with
//!    $A[k][i] = \pm a_{(k-i) \bmod N}$ (sign flip exactly when the
//!    product wraps past $X^N = -1$).
//! 2. **Pack.** $p_i(X) = \sum_j A_j[k][i]\, X^j$: clue $j$'s scalar in
//!    coefficient $j$ of one plaintext polynomial.
//! 3. **Multiply.** $\sum_i p_i \cdot s_i$, where $s_i$ enters as an
//!    encrypted constant; coefficient $j$ of the result is
//!    $(a^{(j)} s)_k$, all clues at once.
//! 4. **Decide.** Per clue, $d_j = b^{(j)}[k] - (a^{(j)} s)_k$; the
//!    message is *pertinent* iff every decoded bit is 0, i.e.
//!    $|d_j| \le r$ for all payload positions.

use crate::field_alg::{Fp, PolyFp, Rq};
use crate::param;
use crate::rlwenc::{Ciphertext, SecretKey, SkParam1};

/// The detector's plaintext ring $\mathbb{Z}_t[X]/(X^{D}+1)$ with
/// $t = q$ and $D = 2048$ (Algorithm 1 sets $t = q$).
pub type DetectorRing = PolyFp<{ param::Q }, { param::D_BFV }>;

/// The unfolding scalar $A[k][i]$: position $i$ of the secret pairs with
/// coefficient $(k-i) \bmod N$ of $a$, negated when the product wraps
/// negacyclically ($i > k$).
pub fn lin_form_scalar<const P: u64, const N: usize>(
    a: &PolyFp<P, N>,
    k: usize,
    i: usize,
) -> Fp<P> {
    debug_assert!(k < N && i < N);
    let j = (k + N - i) % N;
    let val = a.coeff(j);
    if i <= k { val } else { -val }
}

/// Pack the $D$-clue scalar column for secret position $i$ and decoded
/// bit $k$ into one plaintext polynomial
/// $p_i(X) = \sum_j A_j[k][i]\, X^j$.
///
/// Requires `a_list.len() <= D` (the packing capacity is one clue per
/// slot).
pub fn pack_scalar_column<const P: u64, const N: usize, const D: usize>(
    a_list: &[PolyFp<P, N>],
    k: usize,
    i: usize,
) -> PolyFp<P, D> {
    debug_assert!(a_list.len() <= D);
    let mut coeffs = [0u64; D];
    for (j, a) in a_list.iter().enumerate() {
        coeffs[j] = lin_form_scalar(a, k, i).raw();
    }
    PolyFp::from_raw(coeffs)
}

/// SIMD partial decryption of decoded bit $k$ for every clue:
/// $\sum_{i < \text{active}} p_i \cdot s_i$, returning the centered
/// value of $(a^{(j)} s)_k$ per clue $j$ (coefficient $j$ of the packed
/// result).
///
/// The homomorphic layer is simulated: each $p_i$ is multiplied by the
/// plaintext constant $s_i$, i.e. one scalar-broadcast per secret
/// position — $n$ scalar multiplications for up to $D$ clues.
pub fn packed_partial_decrypt<const P: u64, const N: usize, const D: usize>(
    a_list: &[PolyFp<P, N>],
    sk: &SecretKey<P, N>,
    k: usize,
    active: usize,
) -> Vec<i64> {
    debug_assert!(a_list.len() <= D);
    let mut total = PolyFp::<P, D>::zero();
    for i in 0..active {
        let p_i = pack_scalar_column::<P, N, D>(a_list, k, i);
        let s_i = sk.as_poly().coeff(i);
        total = &total + &p_i.scalar_mul(s_i);
    }
    (0..a_list.len()).map(|j| total.coeff_centered(j)).collect()
}

/// Honest per-clue partial decryption
/// $(a^{(j)} s)_k$ by direct ring multiplication, for cross-checking
/// the packed path (they agree exactly: packing is exact algebra).
pub fn direct_partial_decrypt<const P: u64, const N: usize>(
    a_list: &[PolyFp<P, N>],
    sk: &SecretKey<P, N>,
    k: usize,
) -> Vec<i64> {
    a_list.iter().map(|a| (a * sk.as_poly()).coeff_centered(k)).collect()
}

/// The pertinence decision from the centered decryption values of all
/// payload bits: a message is pertinent iff every bit decodes to 0,
/// i.e. $|d| \le r$ throughout (Definition 4.6 decoding, $\ell$ bits).
pub fn is_pertinent(d_values: &[i64], r: i64) -> bool {
    d_values.iter().all(|d| d.abs() <= r)
}

/// The client's detection pass over a server reply: for each clue,
/// finish the partial decryption,
/// $d_j = b^{(j)}[k] - \text{reply}[j]$, read the value centered in
/// $(-q/2, q/2]$, and mark the clue pertinent iff $|d_j| \le r$.
/// Returns the indices of the pertinent clues, in order.
///
/// This is the recipient-side step of UnifOMD: the server's packed
/// SIMD pass ([`packed_partial_decrypt`]) already computed the
/// $(a^{(j)} s)_k$ half of every clue at once; the client only adds
/// the public $b^{(j)}[k]$ term and range-checks — one subtraction, one
/// comparison per clue. Recall the semantics: a clue decrypting to
/// $0^\ell$ is the pertinence signal (senders encrypt $m = 0$ to flag
/// "this is for you"; under any other key $d_j$ is pseudorandom and
/// fails the check with probability $1 - \epsilon_p$).
pub fn client_detect<const P: u64, const N: usize>(
    reply: &[i64],
    cts: &[Ciphertext<P, N>],
    k: usize,
    r: i64,
) -> Vec<usize> {
    debug_assert!(reply.len() == cts.len());
    let mut pertinent = Vec::new();
    for (j, ct) in cts.iter().enumerate() {
        let b_k = Fp::<P>::new(ct.b[k].raw());
        let d = (b_k - Fp::<P>::from_i64(reply[j])).lift_centered();
        if d.abs() <= r {
            pertinent.push(j);
        }
    }
    pertinent
}

/// End-to-end detection of decoded bit $k$ for a batch of clues: packed
/// SIMD partial decryption, then the client's pass
/// ([`client_detect`]) over the reply — $d_j = b^{(j)}[k] -
/// (a^{(j)} s)_k$ and the range check $|d_j| \le r$, per clue.
pub fn detect<const P: u64, const N: usize, const D: usize>(
    cts: &[Ciphertext<P, N>],
    sk: &SecretKey<P, N>,
    k: usize,
    active: usize,
    r: i64,
) -> Vec<bool> {
    let a_list: Vec<PolyFp<P, N>> = cts.iter().map(|ct| ct.a.clone()).collect();
    let packed = packed_partial_decrypt::<P, N, D>(&a_list, sk, k, active);
    let pertinent = client_detect(&packed, cts, k, r);
    let mut out = vec![false; cts.len()];
    for j in pertinent {
        out[j] = true;
    }
    out
}

/// The server's packed reply for one recipient at Parameter Set 1:
/// SIMD partial decryption of decoded bit $k = 0$ for every clue in
/// the window — $n = 900$ scalar broadcasts of dimension $D = 2048$,
/// one clue per slot — yielding the centered $(a^{(j)} s)_0$ value per
/// clue.
///
/// In the real scheme these values leave the server as $\ell$ elements
/// of $\mathbb{Z}_{Q'}$ per message (4 B each after modulus
/// switching); the recipient finishes with [`client_detect`].
pub fn packed_reply(a_list: &[Rq], sk: &SkParam1) -> Vec<i64> {
    packed_partial_decrypt::<{ param::Q }, { param::N_RING }, { param::D_BFV }>(
        a_list, sk, 0, param::N,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rlwenc::{decrypt, keygen, keygen_param1, encrypt_param1};
    use rand::rngs::ThreadRng;
    use rand::Rng;

    type ToyRing = PolyFp<17, 4>;
    type ToySk = SecretKey<17, 4>;

    #[test]
    fn worked_toy_packing() {
        let clues = [
            ToyRing::from_i64_coeffs(&[1, 2, 0, 5]),
            ToyRing::from_i64_coeffs(&[4, -3, 2, 0]),
            ToyRing::from_i64_coeffs(&[-2, 1, 0, -4]),
            ToyRing::from_i64_coeffs(&[6, -1, 3, 2]),
        ];
        let p0 = pack_scalar_column::<17, 4, 4>(&clues, 0, 0);
        let p1 = pack_scalar_column::<17, 4, 4>(&clues, 0, 1);
        let p2 = pack_scalar_column::<17, 4, 4>(&clues, 0, 2);
        let p3 = pack_scalar_column::<17, 4, 4>(&clues, 0, 3);
        assert_eq!(
            (0..4).map(|i| p0.coeff_centered(i)).collect::<Vec<_>>(),
            vec![1, 4, -2, 6]
        );
        assert_eq!(
            (0..4).map(|i| p1.coeff_centered(i)).collect::<Vec<_>>(),
            vec![-5, 0, 4, -2]
        );
        assert_eq!(
            (0..4).map(|i| p2.coeff_centered(i)).collect::<Vec<_>>(),
            vec![0, -2, 0, -3]
        );
        assert_eq!(
            (0..4).map(|i| p3.coeff_centered(i)).collect::<Vec<_>>(),
            vec![-2, 3, -1, 1]
        );

        let s = ToyRing::from_i64_coeffs(&[1, 0, 0, -1]);
        let sk: ToySk = SecretKey::new(s);
        let packed = packed_partial_decrypt::<17, 4, 4>(&clues, &sk, 0, 4);
        assert_eq!(packed, vec![3, 1, -1, 5]);

        let direct = direct_partial_decrypt(&clues, &sk, 0);
        assert_eq!(packed, direct);
    }

    #[test]
    fn packed_matches_direct_on_random_clues() {
        let mut rng = ThreadRng::default();
        for _ in 0..50 {
            let (sk, _): (ToySk, _) = keygen::<17, 4, 2, _>(&mut rng, 4);
            let clues: Vec<ToyRing> =
                (0..4).map(|_| ToyRing::from_raw(std::array::from_fn(|_| rng.gen_range(0..17)))).collect();
            let packed = packed_partial_decrypt::<17, 4, 4>(&clues, &sk, 0, 4);
            let direct = direct_partial_decrypt(&clues, &sk, 0);
            assert_eq!(packed, direct);
        }
    }

    #[test]
    fn pertinence_boundary() {
        assert!(is_pertinent(&[0, 3, -3], 3));
        assert!(!is_pertinent(&[4], 3));
        assert!(!is_pertinent(&[0, -4], 3));
    }

    #[test]
    fn client_detect_returns_indices() {
        let mut rng = ThreadRng::default();
        let (sk, pk) = keygen_param1(&mut rng);
        let (_, other_pk) = keygen_param1(&mut rng);

        let own_ct = encrypt_param1(&mut rng, &pk, &[false]);
        let other_ct = encrypt_param1(&mut rng, &other_pk, &[false]);
        let cts = [own_ct, other_ct];

        let reply = direct_partial_decrypt(
            &cts.iter().map(|ct| ct.a.clone()).collect::<Vec<_>>(),
            &sk,
            0,
        );
        assert_eq!(client_detect(&reply, &cts, 0, param::R), vec![0]);

        let reply = packed_reply(&cts.iter().map(|ct| ct.a.clone()).collect::<Vec<_>>(), &sk);
        assert_eq!(reply.len(), 2);
        assert_eq!(client_detect(&reply, &cts, 0, param::R), vec![0]);
    }

    #[test]
    fn param1_detection_mix() {
        let mut rng = ThreadRng::default();
        let (sk, pk) = keygen_param1(&mut rng);
        let (_, other_pk) = keygen_param1(&mut rng);

        let mut cts = Vec::new();
        for _ in 0..3 {
            cts.push(encrypt_param1(&mut rng, &pk, &[false]));
        }
        for _ in 0..3 {
            cts.push(encrypt_param1(&mut rng, &other_pk, &[false]));
        }

        let decisions = detect::<{ param::Q }, { param::N_RING }, { param::D_BFV }>(
            &cts,
            &sk,
            0,
            param::N,
            param::R,
        );
        assert_eq!(decisions, vec![true, true, true, false, false, false]);

        let honest: Vec<bool> = cts
            .iter()
            .map(|ct| decrypt(&sk, ct, param::R) == vec![false])
            .collect();
        assert_eq!(decisions, honest);

        let pertinent = &cts[0];
        assert_eq!(decrypt(&sk, pertinent, param::R), vec![false]);
        let impertinent = &cts[3];
        assert_eq!(decrypt(&sk, impertinent, param::R), vec![true]);
    }
}
