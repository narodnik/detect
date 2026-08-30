# detect2 — UnifOMR RLWEenc study & simulation

A Rust library (plus SageMath scripts and a writeup) studying and
simulating the RLWEenc layer of **UnifOMR** — oblivious message
detection/retrieval — instantiated at **Parameter Set 1** (Table 1,
Section 7.1 of the paper):

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
│                  and NTT conditions, r', eps_p; machine-precision erfc)
├── field_alg.rs   F_p field + PolyFp<P, N>: the negacyclic polynomial
│                  ring F_p[X]/(X^N + 1) with std::ops traits
├── error.rs       chi_sigma (exact discrete-Gaussian CDT), ternary secret
│                  distribution D, uniform sampling — generic over the ring
├── rlwenc.rs      RLWEenc PKE (Definition 4.6): KeyGen / Enc / Dec,
│                  ciphertext (a, b) in R_q x Z_q^ell (b truncated)
└── detection.rs   Detector side: linear-form unfolding, coefficient
                   packing p_i(X) = sum_j A_j[k][i] X^j, SIMD partial
                   decryption (n scalar broadcasts for D clues), and the
                   pertinence decision |d| <= r

script/           SageMath companions: error sampling, the rounded-normal
                   variant, and the toy-ring packing walkthrough
examples/         sample_errors.rs — Rust mirror of the sampler script
benches/          bench_param1.rs — full-window Param-1 benchmark
```

The homomorphic BFV layer is simulated on plaintexts (constant x
polynomial), preserving the real scheme's cost shape and communication
profile.

## Commands

```sh
# tests: 24 unit tests incl. parameter-derivation checks
cargo test

# error-sampling verification (mirror of script/sample_errors.sage)
cargo run --release --example sample_errors

# full end-to-end flow: 4 recipients, 4 clues, 1 server
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
| detection key (client -> server) | 27,648,000 B (900 BFV cts @ 2x2048x60 bits; analytic, BFV simulated) |
| one clue (sender -> server) | 2819 B (a: 2816 B + b: 3 B) |
| clues for 100k messages | 281,900,000 B (268.8 MiB, 49 windows) |
| packed reply (server -> client) | 4 B/message -> 400,000 B for 100k messages |

Timings (min of 20 windows, release build):

| phase | per window | for 100k clues |
|---|---|---|
| clue generation (senders) | ~2.1 s (1050 us/clue) | ~105 s |
| packed partial decrypt (server) | ~15 ms (7.1 us/clue) | ~0.7 s |
| detection pass (client) | ~6 us (3 ns/clue) | ~0.3 ms |

The ~150x gap between per-clue encryption (two ring multiplications
each) and the server's amortized packed pass, and the 2819 B -> 4 B
reply compression, are the SIMD-packing win of UnifOMR made concrete.
Sanity checks: all pertinent clues detected, packed path cross-checked
against direct decryption.

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
