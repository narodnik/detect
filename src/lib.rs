//! Study crate for the RLWEenc layer of UnifOMR (oblivious message
//! retrieval), Parameter Set 1.
//!
//! The modules mirror the accompanying writeup `crypto.md`:
//!
//! - [`param`]: the hardcoded Parameter Set 1 values of Table 1, with unit
//!   tests verifying every derived parameter against its defining
//!   constraint.
//! - [`field_alg`]: the finite field $\mathbb{F}_q = \mathbb{Z}/q\mathbb{Z}$
//!   and the negacyclic polynomial ring
//!   $R_q = \mathbb{F}_q[X]/(X^{n'} + 1)$.
//! - [`rlwenc`]: the RLWEenc PKE itself (KeyGen, Enc, Dec).
//! - [`bfv`]: a minimal BFV layer — fresh encryption of constants,
//!   plaintext-by-ciphertext multiplication, and modulus switching —
//!   the exact homomorphic vocabulary of the detection circuit.
//! - [`detection`]: the detector side — the encrypted detection key,
//!   the server's packed homomorphic partial decryption, and the
//!   client's decryption + pertinence decision.
//! - [`error`]: the error distribution $\chi_\sigma$ (discrete Gaussian),
//!   the ternary secret distribution $\mathcal{D}$, and uniform sampling.
//!
//! The `sample_errors` example reproduces the verification output of
//! `script/sample_errors.sage`; the `full_usage` example runs the full
//! encrypted detection flow (4 recipients, 4 clues, 1 server).

pub mod bfv;
pub mod detection;
pub mod error;
pub mod field_alg;
pub mod param;
pub mod rlwenc;

pub use detection::{client_decode, client_detect, detection_key, packed_reply};
