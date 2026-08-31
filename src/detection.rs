//! The detector side of UnifOMR (UnifOMD): packed partial decryption of
//! many clues against an *encrypted* detection key, plus the pertinence
//! decision.
//!
//! Mirrors §3.4 of `crypto.md` / `script/packing.sage` and Algorithm 1
//! of the paper. The BFV plaintext ring is $\mathbb{Z}_t[X]/(X^{D}+1)$
//! with $t = q$ and $D = 2048$ slots; the homomorphic layer is a real
//! (minimal) BFV scheme ([`crate::bfv`]), so the server only ever
//! touches ciphertexts.
//!
//! Pipeline for decoded bit $k$ and clues $j = 0, \dots$:
//!
//! 1. **Unfold.** The $k$-th coefficient of each negacyclic product is a
//!    linear form in the secret coefficients,
//!    $(a^{(j)} s)_k = \sum_i A_j[k][i]\, s_i$ with
//!    $A[k][i] = \pm a_{(k-i) \bmod N}$ (sign flip exactly when the
//!    product wraps past $X^N = -1$).
//! 2. **Detect.** The recipient uploads a detection key
//!    ([`detection_key`]): every secret coefficient $s_w$ encrypted as
//!    a BFV constant under the recipient's own BFV key.
//! 3. **Pack.** $p_w(X) = \sum_j A_j[k][w]\, X^j$: clue $j$'s scalar in
//!    coefficient $j$ of one plaintext polynomial.
//! 4. **Multiply.** $\sum_w p_w \cdot \mathrm{ct}_{s_w}$, homomorphically:
//!    coefficient $j$ of the still-encrypted result is
//!    $(a^{(j)} s)_k$, all clues at once; the packed public
//!    $b^{(j)}[k]$ term is folded in at message scale and the
//!    ciphertext is modulus-switched to $Q'$ ([`packed_reply`]).
//! 5. **Decide.** The client decrypts with the BFV secret key and
//!    range-checks $|v_j| \le r' Q'/t$ ([`client_decode`]) — the
//!    message is *pertinent* iff every decoded bit is 0.
//!
//! [`packed_partial_decrypt`] / [`detect`] evaluate the same algebra on
//! the *raw* secret key — a plaintext simulation kept as a cross-check
//! oracle for the encrypted path (and as the bench's arithmetic-only
//! microbenchmark).

use rand::Rng;

use crate::bfv;
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
    a_list
        .iter()
        .map(|a| (a * sk.as_poly()).coeff_centered(k))
        .collect()
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

/// The recipient's public detection key (Algorithm 1, KeyGen):
/// $\mathrm{pk}_{det} = (\mathrm{pk}_{BFV}, \mathrm{ct}_{s_0}, \dots,
/// \mathrm{ct}_{s_{n-1}})$ — one BFV ciphertext per RLWEenc secret-key
/// coefficient, each encrypting $s_w \in \{-1, 0, 1\}$ as a *constant*
/// polynomial under the recipient's own BFV key.
///
/// The server evaluates the packed partial decryption purely against
/// these ciphertexts ([`packed_reply`]); the matching BFV secret key
/// never leaves the client, so the server cannot finish any decryption
/// — including the pertinence range check — for itself.
pub struct DetectionKey {
    /// The BFV public key under which the secret coefficients are
    /// encrypted (part of the paper's $\mathrm{pk}_{det}$; the
    /// detection circuit itself only reads `cts`).
    pub pk: bfv::PkParam1,
    /// The encrypted secret coefficients $\mathrm{ct}_{s_w}$, $w \in [n]$.
    pub cts: Vec<bfv::CtParam1>,
}

/// Recipient-side detection-key generation (Algorithm 1, KeyGen): draw
/// a BFV key pair and encrypt every RLWEenc secret coefficient as a
/// constant. Returns the uploadable public detection key together with
/// the client-held BFV secret key ([`bfv::SkParam1`]).
///
/// At Parameter Set 1 this is $n = 900$ BFV ciphertexts of
/// $2 \times D$ coefficients at 60 bits — 27,648,000 B packed.
pub fn detection_key<R: Rng>(rng: &mut R, sk: &SkParam1) -> (DetectionKey, bfv::SkParam1) {
    let (bfv_sk, bfv_pk) = bfv::keygen_param1(rng);
    let cts = (0..param::N)
        .map(|w| bfv::encrypt_param1(rng, &bfv_pk, sk.as_poly().coeff(w).lift_centered()))
        .collect();
    (DetectionKey { pk: bfv_pk, cts }, bfv_sk)
}

/// The server's packed reply for one recipient over a window of clues
/// (Algorithm 1, Retrieve0): the homomorphic SIMD partial decryption of
/// decoded bit $k = 0$ for every clue, evaluated against the encrypted
/// detection key only — the server never sees a secret key or a
/// decrypted value.
///
/// 1. For each secret position $w \in [n]$, pack the linear-form
///    scalars $p_w(X) = \sum_j A_j[0][w] X^j$ (clue $j$'s scalar in
///    coefficient $j$) and plaintext-multiply the detection-key
///    ciphertext $\mathrm{ct}_{s_w}$ by it — $n = 900$ plain-mults of
///    dimension $D = 2048$, one clue per slot, all clues at once.
/// 2. Sum over $w$, negate, and homomorphically add the packed public
///    $b$-term: the ciphertext now encrypts the decryption value
///    $d_j = b^{(j)}[0] - (a^{(j)} s)_0$ itself, at message scale.
/// 3. Modulus-switch $Q \to Q'$.
///
/// The reply is one BFV ciphertext in $R_{Q'}^2$: $2D$ coefficients of
/// 26 bits, 13,312 B per window of up to $D = 2048$ clues (6.5 B per
/// message) — and it is opaque to everyone but the holder of the BFV
/// secret key, who finishes with [`client_decode`].
pub fn packed_reply(
    cts: &[Ciphertext<{ param::Q }, { param::N_RING }>],
    dk: &DetectionKey,
) -> bfv::CtQpParam1 {
    debug_assert!(cts.len() <= param::D_BFV);
    let a_list: Vec<Rq> = cts.iter().map(|ct| ct.a.clone()).collect();
    let mut acc = bfv::CtParam1::zero();
    for w in 0..param::N {
        let p_w =
            pack_scalar_column::<{ param::Q }, { param::N_RING }, { param::D_BFV }>(&a_list, 0, w)
                .lift_to::<{ param::Q_BFV }>();
        acc = &acc + &bfv::plain_mul(&dk.cts[w], &p_w);
    }
    let mut b_coeffs = [0u64; param::D_BFV];
    for (j, ct) in cts.iter().enumerate() {
        b_coeffs[j] = ct.b[0].raw();
    }
    let b_packed = PolyFp::<{ param::Q_BFV }, { param::D_BFV }>::from_raw(b_coeffs);
    let with_b = bfv::add_plaintext(&(-&acc), &b_packed, param::T);
    bfv::mod_switch_param1(&with_b)
}

/// The client's detection pass over an encrypted reply (Algorithm 1,
/// Decode0): decrypt $c_0' - c_1' s_{BFV}$ in $R_{Q'}$, read
/// coefficient $j$ centered for each clue $j$ of the window, and mark
/// the clue pertinent iff $|v_j| \le r' \cdot Q'/t$.
///
/// $v_j$ is the decryption value $d_j$ scaled by $Q'/t \approx 16.1$
/// plus the merged BFV/modulus-switch error, so the check is the
/// encrypted-flow counterpart of the plaintext range check $|d_j| \le
/// r$: a clue addressed to this recipient decrypts to $0^\ell$ (the
/// pertinence signal, senders encrypt $m = 0$), while under anyone
/// else's key $d_j$ is pseudorandom over $\mathbb{Z}_t$ and fails with
/// probability $\approx 1 - 2r'/Q'$.
///
/// Returns the indices of the pertinent clues, in order. Only the
/// holder of `bfv_sk` can perform this step — the server, holding just
/// the detection key, cannot.
pub fn client_decode(
    reply: &bfv::CtQpParam1,
    bfv_sk: &bfv::SkParam1,
    num_clues: usize,
) -> Vec<usize> {
    debug_assert!(num_clues <= param::D_BFV);
    let v = bfv::decrypt_phase_param1(bfv_sk, reply);
    (0..num_clues)
        .filter(|j| v.coeff_centered(*j).abs() <= param::REPLY_RANGE)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rlwenc::{decrypt, encrypt_param1, keygen, keygen_param1};
    use rand::Rng;
    use rand::rngs::ThreadRng;

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
            let clues: Vec<ToyRing> = (0..4)
                .map(|_| ToyRing::from_raw(std::array::from_fn(|_| rng.gen_range(0..17))))
                .collect();
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

        let reply = packed_partial_decrypt::<{ param::Q }, { param::N_RING }, { param::D_BFV }>(
            &cts.iter().map(|ct| ct.a.clone()).collect::<Vec<_>>(),
            &sk,
            0,
            param::N,
        );
        assert_eq!(reply.len(), 2);
        assert_eq!(client_detect(&reply, &cts, 0, param::R), vec![0]);
    }

    #[test]
    fn param1_encrypted_detection_flow() {
        let mut rng = ThreadRng::default();
        let (sk, pk) = keygen_param1(&mut rng);
        let (_, other_pk) = keygen_param1(&mut rng);

        let cts = [
            encrypt_param1(&mut rng, &pk, &[false]),
            encrypt_param1(&mut rng, &other_pk, &[false]),
            encrypt_param1(&mut rng, &pk, &[false]),
        ];

        let (dk, bfv_sk) = detection_key(&mut rng, &sk);
        let reply = packed_reply(&cts, &dk);
        assert_eq!(client_decode(&reply, &bfv_sk, cts.len()), vec![0, 2]);

        // the encrypted path agrees with the plaintext simulation
        let sim = detect::<{ param::Q }, { param::N_RING }, { param::D_BFV }>(
            &cts,
            &sk,
            0,
            param::N,
            param::R,
        );
        assert_eq!(sim, vec![true, false, true]);

        // a third recipient's detection key finds nothing: their BFV
        // decryption is exact, but every d_j is pseudorandom under a
        // clue key that is not theirs
        let (third_sk, _) = keygen_param1(&mut rng);
        let (third_dk, third_bfv_sk) = detection_key(&mut rng, &third_sk);
        let third_reply = packed_reply(&cts, &third_dk);
        assert_eq!(
            client_decode(&third_reply, &third_bfv_sk, cts.len()),
            Vec::<usize>::new()
        );
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
