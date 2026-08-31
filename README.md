# detect — UnifOMR Detection

A Rust library (plus SageMath scripts and a writeup) studying the
RLWEenc + BFV detection flow of **UnifOMR** — oblivious message
detection/retrieval.

## Overview

1. Each wallet has a 5.6k public key for receiving payments from senders.
2. The wallet uploads a 27Mb detection key to the server: the secret key
   encrypted coefficient-by-coefficient (n = 900 BFV ciphertexts) under
   a fresh BFV key whose secret half never leaves the wallet. ~1.3 s,
   once.
3. Senders attach 2.8k clue key which takes 1ms to generate to their
   payments.
4. Servers scan a window of 2048 payments for each wallet with one
   packed homomorphic pass (900 plaintext-by-ciphertext multiplications,
   ~5.5 s per window in this schoolbook implementation, ~270 s for
   100k payments) and reply with a single 13kB modulus-switched BFV
   ciphertext — which the server itself cannot decrypt (~650 kB of
   replies for 100k payments).
5. Wallets decrypt the packed reply (~1 ms per window, ~30 ms for 100k
   payments) and get a list of indices indicating which payments in the
   window belong to them.

See `example/full_usage.rs` for more info.

The crypto is quite simple although there are some nuances. The security
parameters are taken as-is from the UnifOMR paper.

Check **Algorithm 1** from the UnifOMR paper on page 30. We do up to step
`Decode0` then stop.

## Parameters

Using **Parameter Set 1** (Table 1, Section 7.1 of the paper):

| symbol | value | meaning |
|---|---|---|
| $n$ | 900 | secret key dimension (padded to $n' = 1024$) |
| $q$ | 4169729 | ciphertext modulus (prime, $q - 1 = 2^{13}\cdot 509$) |
| $\sigma$ | 0.6 | error std (discrete Gaussian $D_{\mathbb{Z},\sigma}$) |
| $h$ | 80 | secret Hamming weight (ternary) |
| $r$ | 48 | error range bound |
| $\ell$ | 1 | payload bits per ciphertext |
| $D$ | 2048 | BFV ring dimension / SIMD slots |
| $Q'$ | 67104769 | BFV modulus after switching ($2^{26}-2^{12}+1$) |

The negacyclic ring is $R_q = \mathbb{F}_q[X]/(X^{1024}+1)$; an element
costs 1024 x 22 bits = **2816 bytes** packed.

## Layout

```
src/
├── param.rs       Parameter Set 1 constants + unit tests verifying every
│                  derived parameter (r via the erfc budget, q/Q' primality
│                  and NTT conditions, r', eps_p, reply range; machine-
│                  precision erfc)
├── field_alg.rs   F_p field + PolyFp<P, N>: the negacyclic polynomial
│                  ring F_p[X]/(X^N + 1) with std::ops traits and the
│                  structured products (ternary x dense, dense x small)
│                  that stay exact at the 60-bit BFV modulus
├── error.rs       chi_sigma (exact discrete-Gaussian CDT, any sigma),
│                  ternary secret distribution D, uniform sampling —
│                  generic over the ring
├── rlwenc.rs      RLWEenc PKE (Definition 4.6): KeyGen / Enc / Dec,
│                  ciphertext (a, b) in R_q x Z_q^ell (b truncated)
├── bfv.rs         Minimal BFV: fresh encryption of constants,
│                  plaintext-by-ciphertext multiplication, modulus
│                  switching — the exact homomorphic vocabulary of
│                  Algorithm 1's detection circuit (no ct x ct mult,
│                  no relinearization, no NTT)
└── detection.rs   Detector side: the encrypted detection key
│                   (sk coefficient-wise under the client's BFV key),
│                   linear-form unfolding, coefficient packing
│                   p_w(X) = sum_j A_j[k][w] X^j, the server's packed
│                   homomorphic reply, and the client's decrypt +
│                   |v| <= r'*Q'/t pertinence decision

script/           SageMath companions: error sampling, the rounded-normal
                   variant, and the toy-ring packing walkthrough
examples/         sample_errors.rs — Rust mirror of the sampler script
                  full_usage.rs — 4 recipients, 4 clues, 1 server, fully
                  encrypted
benches/          bench_param1.rs — full-window Param-1 benchmark
```

The full flow is real: RLWEenc clues, BFV detection keys, the
homomorphic packed evaluation, the Q -> Q' modulus switch and the
client-side decryption. Ring multiplication is schoolbook (no NTT —
Q = 2^60 is not NTT-friendly), so the server pass costs ~6 s per
2048-clue window here; an NTT-backed implementation would be orders of
magnitude faster. The plaintext-simulation pass
(`packed_partial_decrypt`) is kept as a cross-check oracle.

## Commands

```sh
# tests: 35 unit tests incl. parameter-derivation checks and the
# end-to-end encrypted detection flow
cargo test

# error-sampling verification (mirror of script/sample_errors.sage)
cargo run --release --example sample_errors

# full end-to-end encrypted flow: 4 recipients, 4 clues, 1 server
cargo run --release --example full_usage

# benchmark: full detector window of D = 2048 clues
cargo bench --bench bench_param1

# sage companions
sage script/sample_errors.sage
sage script/packing.sage
```

## Benchmark (Parameter Set 1, window of D = 2048 clues)

Data sizes (packed-minimal unless noted):

| flow | size |
|---|---|
| public key (recipient) | 5632 B |
| detection key (client -> server) | 27,648,000 B (900 real BFV cts @ 2x2048x60 bits) |
| one clue (sender -> server) | 2819 B (a: 2816 B + b: 3 B) |
| clues for 100k messages | 281,900,000 B (268.8 MiB, 49 windows) |
| packed reply (server -> client) | 13,312 B/window (1 switched BFV ct @ 2x2048x26 bits) -> 6.5 B/message, 663,305 B for 100k |

Timings (release build, schoolbook ring arithmetic, no NTT):

| phase | per window | for 100k clues |
|---|---|---|
| clue generation (senders) | ~1.7 s (813 us/clue) | ~81 s |
| detection key gen (client, once) | ~1.3 s | — |
| packed reply, real BFV (server) | ~5.5 s | ~270 s |
| packed partial decrypt, sim (server) | ~16 ms (7.4 us/clue) | ~0.7 s |
| decode + range check (client) | ~1 ms (324 ns/clue) | ~30 ms |

The reply compression (2819 B/clue in, 6.5 B/message out) is the
SIMD-packing win of UnifOMR made concrete; the ~300x gap between the
real BFV pass and the plaintext simulation is the price of this
implementation's NTT-free schoolbook arithmetic, not of the scheme.
Sanity checks: all pertinent clues detected over the encrypted path,
cross-checked against the plaintext simulation.

## Sage scripts and writeup

- `script/sample_errors.sage` — error sampling + full verification (stats,
  erfc budget, KeyGen demo). `script/sample_errors2.sage` is the rounded-normal
  variant (which, note, misses Param 1's false-negative budget by ~14x).
- `script/packing.sage` — toy-ring (Z_17[X]/(X^4+1)) walkthrough of
  KeyGen/Enc/Dec, four-clue packing, and randomized checks.
- `crypto.md` — the full writeup: scheme, ElGamal structure, the $q/2$
  encoding, packing step by step, worked examples, the erf/erfc
  correctness budget, and derivations of every Table 1 parameter
  (build the PDF with `pandoc crypto.md -o crypto.pdf --toc
  --pdf-engine=pdflatex`).
