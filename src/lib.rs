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
//! - [`detection`]: the detector side — SIMD partial decryption via
//!   coefficient packing and the pertinence decision.
//! - [`error`]: the error distribution $\chi_\sigma$ (discrete Gaussian),
//!   the ternary secret distribution $\mathcal{D}$, and uniform sampling.
//!
//! The `sample_errors` example reproduces the verification output of
//! `script/sample_errors.sage`.

pub mod detection;
pub mod error;
pub mod field_alg;
pub mod param;
pub mod rlwenc;

pub use detection::{client_detect, packed_reply};
