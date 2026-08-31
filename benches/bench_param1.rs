//! Parameter Set 1 benchmark over a full detector window of
//! $D = 2048$ clues, mirroring the data flow of UnifOMR:
//!
//! ```text
//! client --detection key--> server   (n = 900 real BFV ciphertexts)
//! senders --2048 clues----> server   (RLWEenc ciphertexts)
//! server --packed reply---> client   (1 modulus-switched BFV ct)
//! client --decode---------> indices  (decrypt + range check)
//! ```
//!
//! Reports the byte sizes of each flow (packed-minimal for the
//! ciphertexts we implement) and wall-clock times for each phase:
//! clue generation, detection-key generation, the server's packed
//! homomorphic reply (real BFV; one iteration, since the packed
//! plaintexts are dense over a full window and this implementation is
//! NTT-free schoolbook), the plaintext-simulation arithmetic pass for
//! comparison, and the client's decode.
//!
//! Run with:
//!
//! ```sh
//! cargo bench --bench bench_param1
//! ```

use std::time::Instant;

use detect2::client_decode;
use detect2::detection::{
    detect, detection_key, is_pertinent, packed_partial_decrypt, packed_reply,
};
use detect2::param;
use detect2::rlwenc::{Ciphertext, encrypt_param1, keygen_param1};
use rand::rngs::ThreadRng;

const CLUE_WINDOWS: usize = 3;
const SIM_ITERATIONS: usize = 20;
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

    let t_dk = Instant::now();
    let (dk, bfv_sk) = detection_key(&mut rng, &sk);
    let dk_time = t_dk.elapsed();

    let a_list: Vec<_> = cts.iter().map(|ct| ct.a.clone()).collect();

    let mut sim_times = Vec::with_capacity(SIM_ITERATIONS);
    let mut packed = Vec::new();
    for _ in 0..SIM_ITERATIONS {
        let t1 = Instant::now();
        packed = packed_partial_decrypt::<{ param::Q }, { param::N_RING }, { param::D_BFV }>(
            &a_list,
            &sk,
            0,
            param::N,
        );
        sim_times.push(t1.elapsed());
    }

    println!("--- timings (server window = n = 900 plain-mults of dim D, serving D clues) ---");
    report(
        "clue generation (senders)",
        &clue_times,
        num_clues,
        "us",
        1e6,
    );
    println!(
        "{name:28}: {secs:8.3} s/one-off   ({n} BFV encryptions of a constant)",
        name = "detection key gen (client)",
        secs = dk_time.as_secs_f64(),
        n = param::N
    );

    let t2 = Instant::now();
    let reply = packed_reply(&cts, &dk);
    let reply_time = t2.elapsed();
    println!(
        "{name:28}: {secs:8.3} s/window   (real BFV: n plain-mults + mod switch, schoolbook, 1 iteration)",
        name = "packed reply (server)",
        secs = reply_time.as_secs_f64()
    );
    report(
        "packed partial decrypt (sim)",
        &sim_times,
        num_clues,
        "us",
        1e6,
    );

    let mut client_times = Vec::with_capacity(CLIENT_ITERATIONS);
    let mut detections = Vec::new();
    for _ in 0..CLIENT_ITERATIONS {
        let t3 = Instant::now();
        detections = client_decode(&reply, &bfv_sk, num_clues);
        client_times.push(t3.elapsed());
    }
    report(
        "decode + range check (client)",
        &client_times,
        num_clues,
        "ns",
        1e9,
    );

    let detected = detections.len();
    assert_eq!(detected, num_pertinent);

    let decisions = sim_pass(&cts, &packed);
    let sim_detected = decisions.iter().filter(|d| **d).count();
    assert_eq!(sim_detected, num_pertinent);
    for d in &detections {
        assert!(
            decisions[*d],
            "clue {d} flagged by decode but not by the sim"
        );
    }
    let cross = detect::<{ param::Q }, { param::N_RING }, { param::D_BFV }>(
        &cts[..16],
        &sk,
        0,
        param::N,
        param::R,
    );
    assert_eq!(&cross[..8], &decisions[..8]);

    println!();
    println!(
        "sanity: {detected}/{} pertinent detected over the encrypted path, sim agrees, cross-check vs detect() ok",
        num_pertinent
    );
}

fn sim_pass(cts: &[Ciphertext<{ param::Q }, { param::N_RING }>], packed: &[i64]) -> Vec<bool> {
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
    let bfv_ct = 2 * param::D_BFV as u64 * 60 / 8;
    let reply_ct = 2 * param::D_BFV as u64 * q_prime_bits / 8;

    println!("--- data sizes (packed-minimal unless noted) ---");
    println!(
        "public key (recipient)    : {pk:9} B   (2 ring elements: alpha, alpha*s + x)",
        pk = 2 * ring_elem
    );
    println!(
        "detection key (client->server): {dk:9} B   ({n} real BFV ciphertexts of 2 x {d} coeffs @ 60 bits)",
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
        "packed reply (server->client): {re:9} B   (1 switched BFV ct: 2 x {d} coeffs @ {qpb} bits = {per:8.2} B/message, Q' = {qp})",
        re = reply_ct,
        d = param::D_BFV,
        qpb = q_prime_bits,
        per = reply_ct as f64 / num_clues as f64,
        qp = param::Q_PRIME
    );
    println!(
        "replies for 100k messages : {win:9} B   ({kib:6.2} KiB)",
        win = reply_ct * (100_000 + num_clues as u64 - 1) / num_clues as u64,
        kib = reply_ct as f64 * (100_000.0 / num_clues as f64).ceil() / 2f64.powi(10)
    );
    println!(
        "in-memory a storage       : {mem:9} B   (u64 per coefficient, this simulation)",
        mem = param::N_RING as u64 * 8 * num_clues as u64
    );
    println!();
}
