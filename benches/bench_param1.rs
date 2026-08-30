//! Parameter Set 1 benchmark over a full detector window of
//! $D = 2048$ clues, mirroring the data flow of UnifOMR:
//!
//! ```text
//! client --detection key--> server   (n BFV ciphertexts, simulated)
//! senders --2048 clues----> server   (RLWEenc ciphertexts)
//! server --packed reply---> client   (ell elements of Z_Q' per message)
//! ```
//!
//! Reports the byte sizes of each flow (packed-minimal for the
//! ciphertexts we implement; analytic for the BFV layer we simulate)
//! and wall-clock times, each phase repeated over several windows:
//! clue generation, the server's packed partial decryption, and the
//! client's detection pass.
//!
//! Run with:
//!
//! ```sh
//! cargo bench --bench bench_param1
//! ```

use std::time::Instant;

use detect2::detection::{detect, is_pertinent, packed_partial_decrypt};
use detect2::param;
use detect2::rlwenc::{encrypt_param1, keygen_param1, Ciphertext};
use rand::rngs::ThreadRng;

const CLUE_WINDOWS: usize = 3;
const PACK_ITERATIONS: usize = 20;
const CLIENT_ITERATIONS: usize = 20;

fn main() {
    let mut rng = ThreadRng::default();

    let num_clues = param::D_BFV;
    let num_pertinent = num_clues / 2;

    println!("=== UnifOMR Param-1 benchmark: window of D = {num_clues} clues ===");
    println!();
    data_sizes(num_clues);

    let (sk, pk) = keygen_param1(&mut rng);
    let (_, other_pk) = keygen_param1(&mut rng);

    let mut clue_times = Vec::with_capacity(CLUE_WINDOWS);
    let mut cts = Vec::with_capacity(num_clues);
    for _ in 0..CLUE_WINDOWS {
        cts.clear();
        let t0 = Instant::now();
        for i in 0..num_clues {
            let target = if i < num_pertinent { &pk } else { &other_pk };
            cts.push(encrypt_param1(&mut rng, target, &[false]));
        }
        clue_times.push(t0.elapsed());
    }

    let a_list: Vec<_> = cts.iter().map(|ct| ct.a.clone()).collect();

    let mut pack_times = Vec::with_capacity(PACK_ITERATIONS);
    let mut packed = Vec::new();
    for _ in 0..PACK_ITERATIONS {
        let t1 = Instant::now();
        packed = packed_partial_decrypt::<{ param::Q }, { param::N_RING }, { param::D_BFV }>(
            &a_list,
            &sk,
            0,
            param::N,
        );
        pack_times.push(t1.elapsed());
    }

    let mut client_times = Vec::with_capacity(CLIENT_ITERATIONS);
    let mut decisions = Vec::new();
    for _ in 0..CLIENT_ITERATIONS {
        let t2 = Instant::now();
        decisions = client_pass(&cts, &packed);
        client_times.push(t2.elapsed());
    }

    let detected = decisions.iter().filter(|d| **d).count();
    assert_eq!(detected, num_pertinent);
    let cross = detect::<{ param::Q }, { param::N_RING }, { param::D_BFV }>(
        &cts[..16],
        &sk,
        0,
        param::N,
        param::R,
    );
    assert_eq!(&cross[..8], &decisions[..8]);

    println!(
        "--- timings ({} windows each; server window = n = {} scalar broadcasts of dim {}, serving D = {} clues) ---",
        CLUE_WINDOWS.max(PACK_ITERATIONS),
        param::N,
        param::D_BFV,
        param::D_BFV
    );
    report("clue generation (senders)", &clue_times, num_clues, "us", 1e6);
    report("packed partial decrypt (server)", &pack_times, num_clues, "us", 1e6);
    report("detection pass (client)", &client_times, num_clues, "ns", 1e9);
    println!();
    println!("sanity: {detected}/{} pertinent detected, cross-check vs detect() ok", num_pertinent);
}

fn client_pass(cts: &[Ciphertext<{ param::Q }, { param::N_RING }>], packed: &[i64]) -> Vec<bool> {
    let q = param::Q as i64;
    cts.iter()
        .enumerate()
        .map(|(j, ct)| {
            let d = (ct.b[0].raw() as i64 + q - packed[j].rem_euclid(q)) % q;
            let centered = if d > (param::Q / 2) as i64 { d - q } else { d };
            is_pertinent(&[centered], param::R)
        })
        .collect()
}

const CLUES_100K: f64 = 100_000.0;

fn report(name: &str, times: &[std::time::Duration], per: usize, unit: &str, scale: f64) {
    let mut sorted: Vec<f64> = times.iter().map(|t| t.as_secs_f64()).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let min = sorted[0];
    let mean = sorted.iter().sum::<f64>() / sorted.len() as f64;
    println!(
        "{name:28}: mean {mean:8.3} s/window   min {min:8.3} s/window   {per_clue:8.1} {unit}/clue   for 100k clues: {big:8.1} s",
        per_clue = min / per as f64 * scale,
        big = min / per as f64 * CLUES_100K
    );
}

fn data_sizes(num_clues: usize) {
    let q_bits = 22u64;
    let q_prime_bits = 26u64;
    let ring_elem = (param::N_RING as u64 * q_bits).div_ceil(8);
    let b_bytes = (param::ELL as u64 * q_bits).div_ceil(8);
    let clue = ring_elem + b_bytes;
    let reply_per = (param::ELL as u64 * q_prime_bits).div_ceil(8);
    let bfv_ct = 2 * param::D_BFV as u64 * 60 / 8;

    println!("--- data sizes (packed-minimal unless noted) ---");
    println!(
        "public key (recipient)    : {pk:9} B   (2 ring elements: alpha, alpha*s + x)",
        pk = 2 * ring_elem
    );
    println!(
        "detection key (client->server): {dk:9} B   ({n} BFV ciphertexts of 2 x {d} coeffs @ 60 bits; analytic, BFV simulated)",
        dk = param::N as u64 * bfv_ct,
        n = param::N,
        d = param::D_BFV
    );
    println!(
        "one clue (sender->server) : {clue:9} B   (a: {a} B + b: {b} B for ell = {ell})",
        a = ring_elem,
        b = b_bytes,
        ell = param::ELL
    );
    println!(
        "window of {num_clues} clues          : {win:9} B   ({mib:6.2} MiB)",
        win = clue * num_clues as u64,
        mib = (clue * num_clues as u64) as f64 / 2f64.powi(20)
    );
    println!(
        "clues for 100k messages   : {win:9} B   ({mib:6.2} MiB, {n} windows of D)",
        win = clue * 100_000,
        mib = clue as f64 * 100_000.0 / 2f64.powi(20),
        n = (100_000_f64 / num_clues as f64).ceil()
    );
    println!(
        "reply for 100k messages   : {win:9} B   ({kib:6.2} KiB, {per} B/message)",
        win = reply_per * 100_000,
        kib = reply_per as f64 * 100_000.0 / 2f64.powi(10),
        per = reply_per
    );
    println!(
        "packed reply (server->client): {re:9} B   ({per} B/message: ell elements of Z_Q', Q' = {qp})",
        re = reply_per * num_clues as u64,
        per = reply_per,
        qp = param::Q_PRIME
    );    println!(
        "in-memory a storage       : {mem:9} B   (u64 per coefficient, this simulation)",
        mem = param::N_RING as u64 * 8 * num_clues as u64
    );
    println!();
}
