//! Full end-to-end usage of the UnifOMR detection flow at Parameter
//! Set 1: four recipients, four clues, one server — fully unrolled, no
//! loops, every value named.
//!
//! ```text
//! setup    : 4 recipients keygen, publish pk, upload detection keys
//! senders  : 4 clues, each encrypted to one recipient's pk -> server
//! server   : for each recipient, packed SIMD partial decryption of all
//!            clues, sends a tiny per-clue reply
//! clients  : detection pass over the reply -> pertinent clue indices
//! ```
//!
//! The server's homomorphic layer is simulated on plaintexts (it holds
//! the secret coefficients directly); everything else — key sizes,
//! clue sizes, reply sizes, the packed evaluation, the client's range
//! checks — is the real Parameter Set 1 arithmetic.
//!
//! Run with:
//!
//! ```sh
//! cargo run --release --example full_usage
//! ```

use detect2::client_detect;
use detect2::packed_reply;
use detect2::param;
use detect2::rlwenc::{encrypt_param1, keygen_param1};
use rand::rngs::ThreadRng;

fn main() {
    let mut rng = ThreadRng::default();

    println!("=== UnifOMR full usage: 4 recipients, 4 clues, 1 server (Param 1) ===");

    // ------------------------------------------------------------------
    // setup: each recipient generates an RLWEenc key pair. The public
    // key is published for senders; the secret key's coefficients are
    // what the server receives (as n = 900 BFV ciphertexts in the real
    // scheme, 27,648,000 B; stored directly in this simulation)
    // ------------------------------------------------------------------
    let (receiver_0_sk, receiver_0_pk) = keygen_param1(&mut rng);
    let (receiver_1_sk, receiver_1_pk) = keygen_param1(&mut rng);

    // no sender targets receiver 2, so their pk is never used
    let (receiver_2_sk, _receiver_2_pk) = keygen_param1(&mut rng);

    let (receiver_3_sk, receiver_3_pk) = keygen_param1(&mut rng);

    println!("recipient 0: pk published (5632 B), detection key uploaded (27,648,000 B analytic)");
    println!("recipient 1: pk published (5632 B), detection key uploaded (27,648,000 B analytic)");
    println!("recipient 2: pk published (5632 B), detection key uploaded (27,648,000 B analytic)");
    println!("recipient 3: pk published (5632 B), detection key uploaded (27,648,000 B analytic)");

    // ------------------------------------------------------------------
    // senders: encrypt one clue each, addressed to a chosen recipient
    //
    // why [false]? In UnifOMR the plaintext bit is inverted from
    // intuition: a clue that decrypts to 0^ell IS the pertinence
    // signal. Under the recipient's own detection key, d = (q/2)m +
    // noise passes the range check |d| <= r only when m = 0 (a 1-bit
    // sits at q/2 ~ 2.1e6, far outside r = 48); under anyone else's
    // key, d is pseudorandom and fails. So senders flag "this is for
    // you" by encrypting 0, and the bit carries no content -- the real
    // message is later fetched through the separate (batch)PIR channel
    // ------------------------------------------------------------------
    let clue_0 = encrypt_param1(&mut rng, &receiver_3_pk, &[false]);
    let clue_1 = encrypt_param1(&mut rng, &receiver_0_pk, &[false]);
    let clue_2 = encrypt_param1(&mut rng, &receiver_3_pk, &[false]);
    let clue_3 = encrypt_param1(&mut rng, &receiver_1_pk, &[false]);

    // receiver 2 receives no clue at all: nobody sent them anything
    println!("sender 0: clue_0 -> recipient 3 (2819 B uploaded, m = [false])");
    println!("sender 1: clue_1 -> recipient 0 (2819 B uploaded, m = [false])");
    println!("sender 2: clue_2 -> recipient 3 (2819 B uploaded, m = [false])");
    println!("sender 3: clue_3 -> recipient 1 (2819 B uploaded, m = [false])");

    let clues = [clue_0, clue_1, clue_2, clue_3];
    let a_list: Vec<_> = clues.iter().map(|ct| ct.a.clone()).collect();

    // ------------------------------------------------------------------
    // server: one packed SIMD pass per registered recipient -- n = 900
    // scalar broadcasts of dim D = 2048 serve all 4 clues at once. The
    // reply is ell elements of Z_Q' per message (4 B each, 16 B total)
    // ------------------------------------------------------------------
    let receiver_0_reply = packed_reply(&a_list, &receiver_0_sk);
    let receiver_1_reply = packed_reply(&a_list, &receiver_1_sk);
    let receiver_2_reply = packed_reply(&a_list, &receiver_2_sk);
    let receiver_3_reply = packed_reply(&a_list, &receiver_3_sk);

    println!("server: 4 clues scanned per recipient in 1 SIMD pass, reply 16 B each");

    // ------------------------------------------------------------------
    // clients: each recipient range-checks their reply -- the library
    // client_detect computes d_j = b_j[0] - reply_j, reads it centered,
    // and collects every clue with |d_j| <= r
    // ------------------------------------------------------------------
    let receiver_0_detections = client_detect(&receiver_0_reply, &clues, 0, param::R);
    let receiver_1_detections = client_detect(&receiver_1_reply, &clues, 0, param::R);
    let receiver_2_detections = client_detect(&receiver_2_reply, &clues, 0, param::R);
    let receiver_3_detections = client_detect(&receiver_3_reply, &clues, 0, param::R);

    println!("recipient 0: detects pertinent clues {receiver_0_detections:?} (expected [1])");
    println!("recipient 1: detects pertinent clues {receiver_1_detections:?} (expected [3])");
    println!("recipient 2: detects pertinent clues {receiver_2_detections:?} (expected [])");
    println!("recipient 3: detects pertinent clues {receiver_3_detections:?} (expected [0, 2])");

    // ------------------------------------------------------------------
    // sanity: exact expectations, verified
    // ------------------------------------------------------------------
    assert_eq!(receiver_0_detections, vec![1]);
    assert_eq!(receiver_1_detections, vec![3]);
    assert_eq!(receiver_2_detections, Vec::<usize>::new());
    assert_eq!(receiver_3_detections, vec![0, 2]);

    // receiver 2's empty result is the key-privacy property: clues
    // encrypted to other recipients' keys decrypt to pseudorandom
    // values under receiver 2's key and never pass the range check
    println!("sanity: all detections verified; recipient 2 (no clues sent to them) detects nothing");
}
