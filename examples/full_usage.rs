//! Full end-to-end usage of the UnifOMR detection flow at Parameter
//! Set 1: four recipients, four clues, one server — fully unrolled, no
//! loops, every value named.
//!
//! ```text
//! setup    : 4 recipients keygen, publish pk, upload detection keys
//!            (the sk encrypted coefficient-wise under their own BFV key)
//! senders  : 4 clues, each encrypted to one recipient's pk -> server
//! server   : for each recipient, one packed homomorphic pass over all
//!            clues against the encrypted detection key -> reply ct
//! clients  : decrypt the reply with the BFV sk, range-check -> indices
//! ```
//!
//! Everything is the real Parameter Set 1 arithmetic: the RLWEenc clue
//! layer, the BFV detection keys, the homomorphic packed evaluation,
//! the modulus switch to $Q'$ and the client-side decryption with the
//! $r' Q'/t$ range check. The server never touches a plaintext secret
//! key — it cannot learn any recipient's pertinent set.
//!
//! Run with:
//!
//! ```sh
//! cargo run --release --example full_usage
//! ```

use detect2::client_decode;
use detect2::detection_key;
use detect2::packed_reply;
use detect2::rlwenc::{encrypt_param1, keygen_param1};
use rand::rngs::ThreadRng;

fn main() {
    let mut rng = ThreadRng::default();

    println!("=== UnifOMR full usage: 4 recipients, 4 clues, 1 server (Param 1) ===");

    // ------------------------------------------------------------------
    // setup: each recipient generates an RLWEenc clue key pair (pk
    // published for senders) and a BFV detection key pair. The
    // detection key uploads n = 900 BFV ciphertexts, one per secret
    // coefficient s_w encrypted as a constant: 27,648,000 B at 60
    // bits. The BFV secret key stays client-side -- without it the
    // server cannot open anything it computes.
    // ------------------------------------------------------------------
    let (receiver_0_sk, receiver_0_pk) = keygen_param1(&mut rng);
    let (receiver_0_dk, receiver_0_bfv_sk) = detection_key(&mut rng, &receiver_0_sk);
    println!(
        "recipient 0: pk published (5632 B), detection key uploaded (27,648,000 B: 900 BFV cts)"
    );

    let (receiver_1_sk, receiver_1_pk) = keygen_param1(&mut rng);
    let (receiver_1_dk, receiver_1_bfv_sk) = detection_key(&mut rng, &receiver_1_sk);
    println!("recipient 1: pk published (5632 B), detection key uploaded (27,648,000 B)");

    // no sender targets receiver 2, so their pk is never used
    let (receiver_2_sk, _receiver_2_pk) = keygen_param1(&mut rng);
    let (receiver_2_dk, receiver_2_bfv_sk) = detection_key(&mut rng, &receiver_2_sk);
    println!("recipient 2: pk published (5632 B), detection key uploaded (27,648,000 B)");

    let (receiver_3_sk, receiver_3_pk) = keygen_param1(&mut rng);
    let (receiver_3_dk, receiver_3_bfv_sk) = detection_key(&mut rng, &receiver_3_sk);
    println!("recipient 3: pk published (5632 B), detection key uploaded (27,648,000 B)");

    // ------------------------------------------------------------------
    // senders: encrypt one clue each, addressed to a chosen recipient
    //
    // why [false]? In UnifOMR the plaintext bit is inverted from
    // intuition: a clue that decrypts to 0^ell IS the pertinence
    // signal. Under the recipient's own key, d = (q/2)m + noise is
    // small only when m = 0 (a 1-bit sits at q/2 ~ 2.1e6, far outside
    // the error range); under anyone else's key, d is pseudorandom.
    // So senders flag "this is for you" by encrypting 0, and the bit
    // carries no content -- the real message is later fetched through
    // the separate (batch)PIR channel
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

    // ------------------------------------------------------------------
    // server: one packed homomorphic pass per registered recipient --
    // n = 900 plaintext-by-ciphertext multiplications of dim D = 2048
    // serve all 4 clues at once (one clue per coefficient). The reply
    // folds in the packed public b-term and modulus-switches Q -> Q':
    // one ciphertext of 2 x 2048 coeffs @ 26 bits = 13,312 B per
    // window, ~6.5 B/message amortized -- and opaque to the server
    // itself, since it only ever held the encrypted detection key
    // ------------------------------------------------------------------
    let receiver_0_reply = packed_reply(&clues, &receiver_0_dk);
    let receiver_1_reply = packed_reply(&clues, &receiver_1_dk);
    let receiver_2_reply = packed_reply(&clues, &receiver_2_dk);
    let receiver_3_reply = packed_reply(&clues, &receiver_3_dk);

    println!("server: 4 clues scanned per recipient in 1 packed BFV pass (900 plain-mults)");
    println!("        reply: 1 switched ciphertext per recipient, 13,312 B (6.5 B/clue amortized)");

    // ------------------------------------------------------------------
    // clients: each recipient decrypts the reply with their BFV secret
    // key -- v_j = (Q'/t) d_j + merged error, centered -- and collects
    // every clue with |v_j| <= r' * Q'/t = 15,963
    // ------------------------------------------------------------------
    let num_clues = clues.len();
    let receiver_0_detections = client_decode(&receiver_0_reply, &receiver_0_bfv_sk, num_clues);
    let receiver_1_detections = client_decode(&receiver_1_reply, &receiver_1_bfv_sk, num_clues);
    let receiver_2_detections = client_decode(&receiver_2_reply, &receiver_2_bfv_sk, num_clues);
    let receiver_3_detections = client_decode(&receiver_3_reply, &receiver_3_bfv_sk, num_clues);

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

    // privacy, twice over: clues encrypted to other recipients' keys
    // decrypt to pseudorandom values and never pass the range check
    // (receiver 2's empty result), and the server could not have run
    // the check itself -- the reply is a BFV ciphertext under the
    // client's key, and all the server ever saw was the detection key
    println!(
        "sanity: all detections verified client-side; recipient 2 (no clues sent to them) detects nothing"
    );
    println!(
        "privacy: the server held only encrypted detection keys -- every pertinence decision above was made by the client"
    );
}
