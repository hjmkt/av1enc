extern crate av1enc;
extern crate rand;

use av1enc::bool_coder::*;
use std::collections::VecDeque;

#[test]
fn enc_dec_equal() {
    const TEST_LENGTH: usize = 1024;
    const PROBABILITY: i32 = 32;

    fn vec_to_queue(v: &Vec<u8>) -> VecDeque<u8> {
        let mut q: VecDeque<u8> = VecDeque::new();
        for bit in v { q.push_back(*bit); }
        q
    }

    let mut original: Vec<u8> = vec![];
    for _ in 0..TEST_LENGTH { original.push(if rand::random::<u8>()>=PROBABILITY as u8 {1} else {0}); }
    let mut coded: Vec<u8> = vec![];

    let mut coder = BoolCoder::new();
    coder.init_encoder();
    coder.encode_symbols(&original, &mut coded, PROBABILITY);
    coder.exit_encoder(&mut coded);

    let mut decoded: Vec<u8> = vec![];
    let mut q = vec_to_queue(&mut coded);
    coder.init_decoder(&mut q, (coded.len()/8) as u32);
    coder.decode_symbols(&mut q, &mut decoded, TEST_LENGTH as i32, PROBABILITY);
    coder.exit_decoder(&mut q);

    assert!(original.iter().eq(decoded.iter()));
}