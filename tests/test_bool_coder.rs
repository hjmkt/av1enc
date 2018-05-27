extern crate av1enc;
extern crate rand;

use av1enc::bool_coder::*;
use std::collections::VecDeque;

fn vec_to_queue(v: &Vec<u8>) -> VecDeque<u8> {
    let mut q: VecDeque<u8> = VecDeque::new();
    for bit in v { q.push_back(*bit); }
    q
}

#[test]
fn enc_dec_equal() {
    const TEST_LENGTH: usize = 1024;
    const PROBABILITY: i32 = 32;

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

#[test]
fn uvlc_enc_dec_equal() {
    let mut coder = BoolCoder::new();
    let mut t = |n: u32| {
        let mut out_bits: Vec<u8> = vec![];
        coder.encode_uvlc(&mut out_bits, n);
        println!("{:?}", out_bits);
        let mut q = vec_to_queue(&out_bits);
        let decoded = coder.decode_uvlc(&mut q);
        println!("{}, {}", n, decoded);
        assert_eq!(n, decoded);
    };

    let test_cases : Vec<u32> = vec![
        0b0, 0b1,
        0b10, 0b11,
        0b100, 0b101, 0b111,
        0b1000, 0b1001, 0b1011,
        0b10000, 0b10001, 0b10101,
        0b100000, 0b100001, 0b100101,
        0b1000000, 0b1000001, 0b1001001,
        0b10000000, 0b10000001, 0b10001001,
        0b100000000, 0b100000001, 0b100010001,
        0b1000000000, 0b1000000001, 0b1000010001,
        0b10000000000, 0b10000000001, 0b10000010001,
        0b100000000000, 0b100000000001, 0b100000100001,
        0b1000000000000, 0b1000000000001, 0b1000000100001,
        0b10000000000000, 0b10000000000001, 0b10000001000001,
        0b100000000000000, 0b100000000000001, 0b100000001000001,
        0b1000000000000000, 0b1000000000000001, 0b1000000010000001,
        0b10000000000000000, 0b10000000000000001, 0b10000000010000001,
        0b100000000000000000, 0b100000000000000001, 0b100000000100000001,
        0b1000000000000000000, 0b1000000000000000001, 0b1000000000100000001,
        0b10000000000000000000, 0b10000000000000000001, 0b10000000001000000001,
        0b100000000000000000000, 0b100000000000000000001, 0b100000000001000000001,
        0b1000000000000000000000, 0b1000000000000000000001, 0b1000000000010000000001,
        0b10000000000000000000000, 0b10000000000000000000001, 0b10000000000100000000001,
        0b100000000000000000000000, 0b100000000000000000000001, 0b100000000000100000000001,
        0b1000000000000000000000000, 0b1000000000000000000000001, 0b1000000000001000000000001,
        0b10000000000000000000000000, 0b10000000000000000000000001, 0b10000000000001000000000001,
        0b100000000000000000000000000, 0b100000000000000000000000001, 0b100000000000010000000000001,
        0b1000000000000000000000000000, 0b1000000000000000000000000001, 0b1000000000000010000000000001,
        0b10000000000000000000000000000, 0b10000000000000000000000000001, 0b10000000000000100000000000001,
        0b100000000000000000000000000000, 0b100000000000000000000000000001, 0b100000000000000100000000000001,
        0b1000000000000000000000000000000, 0b1000000000000000000000000000001, 0b1000000000000001000000000000001,
        0b10000000000000000000000000000000, 0b10000000000000000000000000000001, 0b10000000000000001000000000000001,
    ];

    for n in test_cases { t(n); }
}

#[test]
fn le_enc_dec_equal() {
    const TEST_LENGTH: usize = 8;
    let mut coder = BoolCoder::new();

    for n in 1..9 {
        for _ in 0..TEST_LENGTH {
            let mut original: u64 = 0;
            for _ in 0..n {
                original = (original<<1) + rand::random::<u8>() as u64;
            }
            let mut coded: Vec<u8> = vec![];
            coder.encode_le(&mut coded, original, n);
            let mut q = vec_to_queue(&coded);
            let decoded: u64 = coder.decode_le(&mut q, n);
            assert_eq!(original, decoded);
        }
    }
}

#[test]
fn leb128_enc_dec_equal() {
    const TEST_LENGTH: usize = 8;
    let mut coder = BoolCoder::new();

    for n in 1..8 {
        for _ in 0..TEST_LENGTH {
            let mut original: u64 = 0;
            for _ in 0..n {
                original = (original<<1) + rand::random::<u8>() as u64;
            }
            let mut coded: Vec<u8> = vec![];
            coder.encode_leb128(&mut coded, original);
            let mut q = vec_to_queue(&coded);
            let decoded: u64 = coder.decode_leb128(&mut q);
            assert_eq!(original, decoded);
        }
    }
}

#[test]
fn su_enc_dec_equal() {
    const TEST_LENGTH: usize = 8;
    let mut coder = BoolCoder::new();

    for n in 1..32 {
        for _ in 0..TEST_LENGTH {
            let mut original: i64 = (rand::random::<u32>() % (1u64<<(n+1)-1) as u32) as i64 + (-1i64 << (n-1));
            let mut coded: Vec<u8> = vec![];
            coder.encode_su(&mut coded, original, n);
            let mut q = vec_to_queue(&coded);
            let decoded: i64 = coder.decode_su(&mut q, n);
            assert_eq!(original, decoded);
        }
    }
}

#[test]
fn ns_enc_dec_equal() {
    const TEST_LENGTH: usize = 1024;
    let mut coder = BoolCoder::new();

    for _ in 0..TEST_LENGTH {
        let mut original: u64 = rand::random::<u32>() as u64;
        let mut n: u64 = original + rand::random::<u32>() as u64;
        let mut coded: Vec<u8> = vec![];
        coder.encode_ns(&mut coded, original, n);
        let mut q = vec_to_queue(&coded);
        let decoded: u64 = coder.decode_ns(&mut q, n);
        assert_eq!(original, decoded);
    }
}