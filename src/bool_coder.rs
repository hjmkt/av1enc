#![allow(dead_code)]
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use std::collections::VecDeque;
use util::*;
use constants::*;
use cdf::*;
use common::*;
use constants::BlockSize::*;
use constants::TxSize::*;
use constants::YMode::*;
use constants::Partition::*;
use constants::RefFrame::*;
use frame::*;
use constants::InterpFilter::*;

#[derive(Copy, Clone)]
pub struct BoolCoder{
    bool_value: i32,
    bool_range: i32,
    bool_max_bits: i32,
    bool_value_last_zero: bool,
    bool_value_zero_count: i32,
    bool_first: bool,
    bit_count: u64,
    leb_128_bytes: u64,
//    CDF tile_cdf;
//    std::array<CDF, 8> cdfs;
}

impl BoolCoder {

    pub fn new() -> BoolCoder {
        BoolCoder {
            bool_value: 0,
            bool_range: 0,
            bool_max_bits: 0,
            bool_value_last_zero: false,
            bool_value_zero_count: 0,
            bool_first: true,
            bit_count: 0,
            leb_128_bytes: 0,
        }
    }

    pub fn push_bit(&self, out_bits: &mut Vec<u8>, bit: u8) {
        out_bits.push(bit);
    }

    pub fn push_bits_with_size(&self, out_bits: &mut Vec<u8>, bits: u32, size: u8) {
        for i in (0..size).rev() { out_bits.push(((bits>>i)&1) as u8); }
    }

    pub fn push_bits(&self, out_bits: &mut Vec<u8>, bits: &Vec<u8>) {
        for bit in bits { out_bits.push(*bit); }
    }

    pub fn push_bits_align(&self, out_bits: &mut Vec<u8>) {
        let size = (8-out_bits.len()%8)%8;
        for _ in 0..size { out_bits.push(0); }
    }

    pub fn init_decoder(&mut self, in_bits: &mut VecDeque<u8>, sz: u32) {
        let mut buf: i32 = 0;
        for _ in 0..15 {
            if let Some(bit) = in_bits.pop_front() { buf = (buf << 1) + bit as i32; }
            else { assert!(false); }
        }
        self.bool_value = (((1i32 << 15) - 1) ^ buf) as i32;
        self.bool_range = 1<<15;
        self.bool_max_bits = 8 * sz as i32 - 15;
    }

    pub fn decode_symbol(&mut self, in_bits: &mut VecDeque<u8>, cdf: &mut [i32]) -> u64 {
        let N: i32 = (cdf.len() - 1) as i32;
        let mut cur: i32 = self.bool_range;
        let mut symbol: i32 = -1;
        let mut prev: i32;
        while {
            symbol += 1;
            prev = cur;
            let f = (1 << 15) - cdf[symbol as usize];
            cur = ((self.bool_range >> 8) * f) >> 7;
            self.bool_value < cur
        } {};

        self.bool_range = prev - cur;
        self.bool_value -= cur;
        let bits: u8 = 15 - msb16(self.bool_range as u16);
        self.bool_range <<= bits;
        let mut new_data: i32 = 0;
        for _ in 0..(bits as i32) {
            if let Some(bit) = in_bits.pop_front() {
                new_data = (new_data << 1) + bit as i32;
                self.bool_max_bits -= 1;
            }
            else { assert!(false); }
        }
        self.bool_value = new_data ^ (((self.bool_value + 1) << bits) - 1);
        self.bool_max_bits -= bits as i32;

        // update cdf
        let rate: i32 = 4 + (cdf[N as usize] > 31) as i32 + msb16(N as u16) as i32;
        let rate2: i32 = 5;
        let tmp0: i32 = 1 << rate2;
        let mut tmp: i32 = tmp0;
        let diff: i32 = (((1 << 15) - (N << rate2)) >> rate) << rate;
        for i in 0..(N-1) {
            tmp += if (i as i32)==symbol { diff } else { 0 };
            cdf[i as usize] -= (cdf[i as usize] - tmp) >> rate;
            tmp += tmp0;
        }
        cdf[N as usize] += (cdf[N as usize]<32) as i32;

        symbol as u64
    }

    pub fn decode_symbolp(&mut self, in_bits: &mut VecDeque<u8>, p: i32) -> u64 {
        let mut cdf: Vec<i32> = vec![ ((p << 15) + 256 - p) >> 8, 1<<15, 0 ];
        self.decode_symbol(in_bits, &mut cdf)
    }

    pub fn decode_symbols(&mut self, in_bits: &mut VecDeque<u8>, out_bits: &mut Vec<u8>, n: i32, p: i32) {
        let mut cdf: Vec<i32> = vec![ ((p << 15) + 256 - p) >> 8, 1<<15, 0 ];
        for _ in 0..n { out_bits.push(self.decode_symbol(in_bits, &mut cdf) as u8); }
    }

    pub fn exit_decoder(&mut self, in_bits: &mut VecDeque<u8>) {
        for _ in 0..self.bool_max_bits { in_bits.pop_front(); }
        self.bool_max_bits = 0;
    }

    pub fn init_encoder(&mut self) {
        self.bool_value = 0;
        self.bool_range = 1<<15;
        self.bool_value_last_zero = false;
        self.bool_value_zero_count = 0;
        self.bool_first = true;
    }

    pub fn encode_symbol(&mut self, symbol: i32, out_bits: &mut Vec<u8>, cdf: &mut [i32]) {
        let N: i32 = (cdf.len() - 1) as i32;
        let fl: i32 = if symbol>0 { (32768-cdf[(symbol-1) as usize]) } else { 0 };
        let fr: i32 = 32768-cdf[symbol as usize];
        let mut left: i32 = ((self.bool_range >> 8) * fl) >> 7;
        let right: i32 = ((self.bool_range >> 8) * fr) >> 7;
        if symbol==0 { left = self.bool_range };

        self.bool_value += self.bool_range - left;
        self.bool_range = left - right;

        if self.bool_value>=65536 {
            self.bool_value -= 65536;
            if self.bool_value_zero_count>0 {
                out_bits.push(1);
                self.bit_count += 1;
                for _ in 0..(self.bool_value_zero_count-1) {
                    out_bits.push(0);
                    self.bit_count += 1;
                }
                self.bool_value_zero_count = 0;
                self.bool_value_last_zero = true;
            }
            else {
                self.bool_value_zero_count = 1;
                self.bool_value_last_zero = false;
            }
        }

        while self.bool_range<32768 {
            let bit = (self.bool_value>>15) & 1;
            if !self.bool_first {
                if bit==0 {
                    if self.bool_value_last_zero {
                        out_bits.push(0);
                        self.bit_count += 1;
                    }
                    for _ in 0..self.bool_value_zero_count {
                        out_bits.push(1);
                        self.bit_count += 1;
                    }
                    self.bool_value_zero_count = 0;
                    self.bool_value_last_zero = true;
                }
                else {
                    self.bool_value_zero_count += 1;
                }
            }
            self.bool_first = false;
            self.bool_value = (self.bool_value<<1) & 0xffff;
            self.bool_range <<= 1;
        }

        // update cdf
        let rate: i32 = 4 + ((cdf[N as usize]>31) as i32) + msb16(N as u16) as i32;
        let rate2: i32 = 5;
        let tmp0: i32 = 1 << rate2;
        let mut tmp: i32 = tmp0;
        let diff: i32 = (((1 << 15) - (N << rate2)) >> rate) << rate;
        for i in 0..(N-1) {
            tmp += if i==symbol { diff } else { 0 };
            cdf[i as usize] -= (cdf[i as usize] - tmp) >> rate;
            tmp += tmp0;
        }
        cdf[N as usize] += (cdf[N as usize]<32) as i32;
    }

    pub fn encode_symbols(&mut self, in_bits: &Vec<u8>, out_bits: &mut Vec<u8>, p: i32) {
        let mut cdf: Vec<i32> = vec![ ((p << 15) + 256 - p) >> 8, 1<<15, 0 ];
        for bit in in_bits { self.encode_symbol(*bit as i32, out_bits, &mut cdf); }
    }

    pub fn encode_symbolp(&mut self, bit: u8, out_bits: &mut Vec<u8>, p: i32) {
        let mut cdf: Vec<i32> = vec![ ((p << 15) + 256 - p) >> 8, 1<<15, 0 ];
        self.encode_symbol(bit as i32, out_bits, &mut cdf);
    }

    // TODO skip cdf update
    pub fn encode_literal(&mut self, out_bits: &mut Vec<u8>, literal: u64, n: u8) {
        for i in (0..n).rev() {
            self.encode_symbolp(((literal>>i)&1) as u8, out_bits, 128);
        }
    }

    pub fn decode_literal(&mut self, in_bits: &mut VecDeque<u8>, n: u8) -> u64 {
        let mut x: u64 = 0;
        for _ in 0..n { x = (x<<1) + self.decode_symbolp(in_bits, 128); }
        x
    }

    pub fn encode_uvlc(&mut self, out_bits: &mut Vec<u8>, val: u32) {
        let v: u64 = (val + 1) as u64;
        let n = msb64(v);
        for _ in 0..n { out_bits.push(0) }
        out_bits.push(1);
        for i in (0..n).rev() { out_bits.push(((v>>i)&1) as u8); }
    }

    pub fn decode_uvlc(&mut self, in_bits: &mut VecDeque<u8>) -> u32 {
        let mut lz = 0;
        loop {
            if let Some(bit) = in_bits.pop_front() {
                let done: bool = bit>0;
                if done { break; }
                else { lz += 1; }
            }
        }
        if lz>=32 { return ((1u64<<32) - 1) as u32; }
        let mut v: u64 = 0;
        for _ in 0..lz {
            if let Some(bit) = in_bits.pop_front() { v = (v<<1) + bit as u64; }
            else { assert!(false); }
        }
        (v + (1<<lz) - 1) as u32
    }

    pub fn encode_le(&mut self, out_bits: &mut Vec<u8>, val: u64, n: u8) {
        for i in 0..n {
            let byte: u8 = ((val>>(i*8)) & 0xff) as u8;
            for j in (0..8).rev() {
                out_bits.push(((byte>>j) & 1) as u8);
            }
        }
    }

    pub fn decode_le(&mut self, in_bits: &mut VecDeque<u8>, n: u8) -> u64 {
        let mut t: u64 = 0;
        for i in 0..n {
            let mut byte: u8 = 0;
            for _ in 0..8 {
                if let Some(bit) = in_bits.pop_front() { byte = (byte<<1) + bit; }
                else { assert!(false); }
            }
            t += (byte as u64) << (i*8);
        }
        t
    }

    pub fn encode_leb128(&mut self, out_bits: &mut Vec<u8>, val: u64) {
        assert!(val<(1<<56));
        let mut v = val;
        while {
            out_bits.push((v>127) as u8);
            for i in (0..7).rev() {
                out_bits.push(((v>>i)&1) as u8);
            }
            v >>= 7;
            v>0
        } {}
    }

    pub fn decode_leb128(&mut self, in_bits: &mut VecDeque<u8>) -> u64 {
        let mut v: u64 = 0;
        for i in 0..8 {
            let mut leb: u64 = 0;
            for _ in 0..8 {
                if let Some(bit) = in_bits.pop_front() { leb = (leb<<1) + bit as u64; }
                else { assert!(false); }
            }
            v |= (leb & 0x7f) << (i*7);
            self.leb_128_bytes += 1;
            if leb & 0x80 == 0 { break; }
        }
        v
    }

    pub fn encode_su(&mut self, out_bits: &mut Vec<u8>, val: i64, n: u8) {
        let signed: bool = val < 0;
        out_bits.push(signed as u8);
        if signed {
            let v: u64 = (val + 2*(1<<(n-1))) as u64;
            for i in (0..(n-1)).rev() { out_bits.push(((v>>i)&1) as u8); }
        }
        else {
            let v: u64 = val.abs() as u64;
            for i in (0..(n-1)).rev() { out_bits.push(((v>>i)&1) as u8); }
        }
    }

    pub fn decode_su(&mut self, in_bits: &mut VecDeque<u8>, n: u8) -> i64 {
        let mut v: u64 = 0;
        for _ in 0..n {
            if let Some(bit) = in_bits.pop_front() { v = (v<<1) + bit as u64; }
            else { assert!(false); }
        }
        let sign_mask: u64 = 1 << (n-1);
        let ret = if v&sign_mask > 0 { v as i64 - 2*(sign_mask as i64) } else { v as i64 };
        ret
    }

    pub fn encode_ns(&mut self, out_bits: &mut Vec<u8>, val: u64, n: u64) {
        let w: u64 = FloorLog2!(n) + 1;
        let m: u64 = (1<<w) - n;
        if val<m {
            for i in (0..(w-1)).rev() {
                out_bits.push(((val>>i)&1) as u8);
            }
        }
        else {
            let v: u64 = val + m;
            for i in (1..w).rev() {
                out_bits.push(((v>>i)&1) as u8);
            }
            out_bits.push((v&1) as u8);
        }
    }

    pub fn decode_ns(&mut self, in_bits: &mut VecDeque<u8>, n: u64) -> u64 {
        let w: u64 = FloorLog2!(n) + 1;
        let m: u64 = (1<<w) - n;
        let mut v: u64 = 0;
        for _ in 0..(w-1) {
            if let Some(bit) = in_bits.pop_front() { v = (v<<1) + bit as u64; }
            else { assert!(false); }
        }
        if v<m { return v };
        if let Some(extra_bit) = in_bits.pop_front() { return (v<<1) - m + extra_bit as u64; }
        else { assert!(false); }
        0
    }

    pub fn exit_encoder(&mut self, out_bits: &mut Vec<u8>) {
        if self.bool_value_last_zero {
            out_bits.push(0);
            self.bit_count += 1;
        }
        for _ in 0..self.bool_value_zero_count {
            out_bits.push(1);
            self.bit_count += 1;
        }
        for i in 0..16 {
            out_bits.push(((self.bool_value >> (15 - i)) & 1) as u8);
        }

        let padding = (8 - self.bit_count % 8) % 8;
        for _ in 0..(padding as i32) {
            out_bits.push(0);
        }
    }
/*
    FIXME no Default_Intrabc_Cdf definition
    pub fn encode_intrabc(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF) {
        let cdf = &mut tile_cdf.use_intrabc;
    }
*/
    pub fn encode_intra_frame_y_mode(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, y_mode: YMode, above_mode: YMode, left_mode: YMode) {
        /*
        abovemode = Intra_Mode_Context[ AvailU ? YModes[ MiRow - 1 ][ MiCol ] : DC_PRED ]
        leftmode = Intra_Mode_Context[ AvailL ? YModes[ MiRow ][ MiCol - 1] : DC_PRED ]*/
        let am = intra_mode_context[above_mode as usize];
        let lm = intra_mode_context[left_mode as usize];
        let cdf = &mut tile_cdf.intra_frame_y_mode[am][lm];
        self.encode_symbol(y_mode as i32, out_bits, cdf);
    }

    pub fn encode_y_mode(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, y_mode: YMode, mi_size: BlockSize) {
        let ctx = size_group[mi_size as usize];
        let cdf = &mut tile_cdf.y_mode[ctx];
        self.encode_symbol(y_mode as i32, out_bits, cdf);
    }

    pub fn encode_uv_mode(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, uv_mode: YMode, mi_size: BlockSize, lossless: bool, subsampling_x: bool, subsampling_y: bool) {
        let cdf: &mut [i32] = if lossless && get_plane_residual_size(mi_size, 1, subsampling_x, subsampling_y)==BLOCK_4X4 {
            &mut tile_cdf.uv_mode_cfl_allowed[uv_mode as usize]
        } else if !lossless && Max!(block_width[mi_size as usize], block_height[mi_size as usize])<=32 {
            &mut tile_cdf.uv_mode_cfl_allowed[uv_mode as usize]
        } else {
            &mut tile_cdf.uv_mode_cfl_not_allowed[uv_mode as usize]
        };
        self.encode_symbol(uv_mode as i32, out_bits, cdf);
    }

    pub fn encode_angle_delta_y(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, y_mode: YMode, angle_delta_y: usize) {
        let cdf = &mut tile_cdf.angle_delta[y_mode as usize - V_PRED as usize];
        self.encode_symbol(angle_delta_y as i32, out_bits, cdf);
    }

    pub fn encode_angle_delta_uv(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, uv_mode: YMode, angle_delta_uv: usize) {
        let cdf = &mut tile_cdf.angle_delta[uv_mode as usize - V_PRED as usize];
        self.encode_symbol(angle_delta_uv as i32, out_bits, cdf);
    }

    pub fn encode_partition(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, bsize: BlockSize, above: bool, left: bool, partition: Partition) {
        let bsl = mi_width_log2[bsize as usize];
        let ctx = (left as usize)*2 + above as usize;
        let cdf: &mut [i32] = if bsl==1 {
            &mut tile_cdf.partition_w8[ctx]
        } else if bsl==2 {
            &mut tile_cdf.partition_w16[ctx]
        } else if bsl==3 {
            &mut tile_cdf.partition_w32[ctx]
        } else if bsl==4 {
            &mut tile_cdf.partition_w64[ctx]
        } else {
            &mut tile_cdf.partition_w128[ctx]
        };
        self.encode_symbol(partition as i32, out_bits, cdf);
    }

    pub fn encode_split_or_horz(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, bsize: BlockSize, above: bool, left: bool, split_or_horz: i32) {
        let bsl = mi_width_log2[bsize as usize];
        let ctx = (left as usize)*2 + above as usize;
        let partition_cdf: &mut [i32] = if bsl==1 {
            &mut tile_cdf.partition_w8[ctx]
        } else if bsl==2 {
            &mut tile_cdf.partition_w16[ctx]
        } else if bsl==3 {
            &mut tile_cdf.partition_w32[ctx]
        } else if bsl==4 {
            &mut tile_cdf.partition_w64[ctx]
        } else {
            &mut tile_cdf.partition_w128[ctx]
        };
        let mut psum = partition_cdf[PARTITION_VERT as usize] - partition_cdf[PARTITION_VERT as usize - 1] +
                   partition_cdf[PARTITION_VERT as usize] - partition_cdf[PARTITION_VERT as usize - 1] +
                   partition_cdf[PARTITION_VERT as usize] - partition_cdf[PARTITION_VERT as usize - 1] +
                   partition_cdf[PARTITION_VERT as usize] - partition_cdf[PARTITION_VERT as usize - 1] +
                   partition_cdf[PARTITION_VERT as usize] - partition_cdf[PARTITION_VERT as usize - 1];
        if bsize != BLOCK_128X128 {
            psum += partition_cdf[PARTITION_VERT_4 as usize] - partition_cdf[PARTITION_VERT_4 as usize - 1];
        }
        let mut cdf = [ 1<<15 - psum, 1<<15, 0 ];
        self.encode_symbol(split_or_horz, out_bits, &mut cdf);
    }

    pub fn encode_split_or_vert(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, bsize: BlockSize, above: bool, left: bool, split_or_vert: i32) {
        let bsl = mi_width_log2[bsize as usize];
        let ctx = (left as usize)*2 + above as usize;
        let partition_cdf: &mut [i32] = if bsl==1 {
            &mut tile_cdf.partition_w8[ctx]
        } else if bsl==2 {
            &mut tile_cdf.partition_w16[ctx]
        } else if bsl==3 {
            &mut tile_cdf.partition_w32[ctx]
        } else if bsl==4 {
            &mut tile_cdf.partition_w64[ctx]
        } else {
            &mut tile_cdf.partition_w128[ctx]
        };
        let mut psum = partition_cdf[PARTITION_VERT as usize] - partition_cdf[PARTITION_VERT as usize - 1] +
                   partition_cdf[PARTITION_SPLIT as usize] - partition_cdf[PARTITION_SPLIT as usize - 1] +
                   partition_cdf[PARTITION_HORZ_A as usize] - partition_cdf[PARTITION_HORZ_A as usize - 1] +
                   partition_cdf[PARTITION_VERT_A as usize] - partition_cdf[PARTITION_VERT_A as usize - 1] +
                   partition_cdf[PARTITION_VERT_B as usize] - partition_cdf[PARTITION_VERT_B as usize - 1];
        if bsize != BLOCK_128X128 {
            psum += partition_cdf[PARTITION_VERT_4 as usize] - partition_cdf[PARTITION_VERT_4 as usize - 1];
        }
        let mut cdf = [ 1<<15 - psum, 1<<15, 0 ];
        self.encode_symbol(split_or_vert, out_bits, &mut cdf);
    }

    pub fn encode_tx_depth(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, mi_row: usize, mi_col: usize, max_rect_tx_size: TxSize, max_tx_depth: usize, frame: &mut Frame, tx_depth: usize) {
        let max_tx_width = tx_width[max_rect_tx_size as usize];
        let max_tx_height = tx_height[max_rect_tx_size as usize];
        let above_w = match frame.get_above_mi(mi_row, mi_col) {
            Some(mi) => if mi.is_inter {
                block_width[mi.mi_size as usize]
            } else {
                tx_width[mi.tx_size as usize]
            },
            None => 0,
        };
        let left_h = match frame.get_left_mi(mi_row, mi_col) {
            Some(mi) => if mi.is_inter {
                block_width[mi.mi_size as usize]
            } else {
                tx_height[mi.tx_size as usize]
            },
            None => 0,
        };
        let ctx = (above_w>=max_tx_width) as usize + (left_h>=max_tx_height) as usize;
        let mut cdf: &mut [i32] = if max_tx_depth==4 {
            &mut tile_cdf.tx_64x64[ctx]
        } else if max_tx_depth==3 {
            &mut tile_cdf.tx_32x32[ctx]
        } else if max_tx_depth==2 {
            &mut tile_cdf.tx_16x16[ctx]
        } else {
            &mut tile_cdf.tx_8x8[ctx]
        };
        self.encode_symbol(tx_depth as i32, out_bits, &mut cdf);
    }

    pub fn encode_txfm_split(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, row: usize, col: usize, mi_row: usize, mi_col: usize, frame: &mut Frame, tx_size: TxSize, mi_size: BlockSize, tx_fm_split: usize) {
        let above = if row==mi_row {
            match frame.get_above_mi(row, col) {
                Some(mi) => if mi.skip && mi.is_inter {
                    block_width[mi.mi_size as usize]
                } else { tx_width[mi.tx_size as usize] },
                None => 64,
            }
        } else { tx_width[frame.get_above_mi(row, col).unwrap().tx_size as usize] };
        let left = if col==mi_col {
            match frame.get_left_mi(row, col) {
                Some(mi) => if mi.skip && mi.is_inter {
                    block_width[mi.mi_size as usize]
                } else { tx_width[mi.tx_size as usize] },
                None => 64,
            }
        } else { tx_width[frame.get_above_mi(row, col).unwrap().tx_size as usize] };
        let size = Min!(64, Max!(block_width[mi_size as usize], block_height[mi_size as usize]));
        let max_tx_size = find_tx_size(size, size);
        let tx_sz_sqr_up = tx_size_sqr_up[tx_size as usize];
        let ctx = (tx_sz_sqr_up!=max_tx_size && max_tx_size>TX_8X8) as usize * 3 + (TX_SIZES-1-max_tx_size as usize) * 6 + above + left;
        let mut cdf: &mut [i32] = &mut tile_cdf.txfm_split[ctx];
        self.encode_symbol(tx_fm_split as i32, out_bits, &mut cdf);
    }

    pub fn encode_segment_id(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, segment_id: usize, prev_u: isize, prev_l: isize, prev_ul: isize) {
        let ctx = if prev_ul<0 || prev_u<0 || prev_l<0 {
            0
        } else if prev_ul==prev_u && prev_ul==prev_l {
            2
        } else if prev_ul==prev_u || prev_ul==prev_l || prev_u==prev_l {
            1
        } else { 0 };
        let mut cdf: &mut [i32] = &mut tile_cdf.segment_id[ctx];
        self.encode_symbol(segment_id as i32, out_bits, &mut cdf);
    }

    pub fn encode_seg_id_predicted(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, seg_id_predicted: usize, mi_row: usize, mi_col: usize, frame: &mut Frame) {
        let ctx = match frame.get_left_mi(mi_row, mi_col) {
            Some(mi) => mi.seg_id_predicted,
            None => 0,
        } + match frame.get_above_mi(mi_row, mi_col) {
            Some(mi) => mi.seg_id_predicted,
            None => 0,
        };
        let mut cdf: &mut [i32] = &mut tile_cdf.segment_id_predicted[ctx];
        self.encode_symbol(seg_id_predicted as i32, out_bits, &mut cdf);
    }

    pub fn encode_new_mv(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, new_mv: usize, new_mv_context: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.new_mv[new_mv_context];
        self.encode_symbol(new_mv as i32, out_bits, &mut cdf);
    }

    pub fn encode_zero_mv(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, zero_mv: usize, zero_mv_context: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.zero_mv[zero_mv_context];
        self.encode_symbol(zero_mv as i32, out_bits, &mut cdf);
    }

    pub fn encode_ref_mv(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, ref_mv: usize, ref_mv_context: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.ref_mv[ref_mv_context];
        self.encode_symbol(ref_mv as i32, out_bits, &mut cdf);
    }

    pub fn encode_drl_mode(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, drl_mode: usize, idx: usize, drl_ctx_stack: &[usize]) {
        let mut cdf: &mut [i32] = &mut tile_cdf.drl_mode[drl_ctx_stack[idx]];
        self.encode_symbol(drl_mode as i32, out_bits, &mut cdf);
    }

    pub fn encode_is_inter(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, is_inter: bool, mi_row: usize, mi_col: usize, frame: &mut Frame) {
        let ctx = match (frame.get_left_mi(mi_row, mi_col), frame.get_above_mi(mi_row, mi_col)) {
            (Some(mi_l), Some(mi_u)) => if mi_l.intra&&mi_u.intra { 3 } else { (mi_l.intra||mi_u.intra) as usize },
            (Some(mi_l), _) => 2*(mi_l.intra as usize),
            (_, Some(mi_u)) => 2*(mi_u.intra as usize),
            _ => 0,
        };
        let mut cdf: &mut [i32] = &mut tile_cdf.is_inter[ctx];
        self.encode_symbol(is_inter as i32, out_bits, &mut cdf);
    }

    pub fn encode_use_filter_intra(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, use_filter_intra: usize, mi_size: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.filter_intra[mi_size];
        self.encode_symbol(use_filter_intra as i32, out_bits, &mut cdf);
    }

    pub fn encode_filter_intra_mode(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, filter_intra_mode: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.filter_intra_mode;
        self.encode_symbol(filter_intra_mode as i32, out_bits, &mut cdf);
    }

    pub fn encode_comp_mode(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, comp_mode: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let ctx = match (frame.get_above_mi(mi_row, mi_col), frame.get_left_mi(mi_row, mi_col)) {
            (Some(mi_u), Some(mi_l)) => if mi_u.single&&mi_l.single {
                (check_backward(mi_u.ref_frame[0]) ^ check_backward(mi_l.ref_frame[0])) as usize
            } else if mi_u.single {
                2 + (check_backward(mi_u.ref_frame[0]) || mi_u.intra) as usize
            } else if mi_l.single {
                2 + (check_backward(mi_l.ref_frame[0]) || mi_l.intra) as usize
            } else { 4 },
            (Some(mi_u), _) => if mi_u.single { check_backward(mi_u.ref_frame[0]) as usize } else { 3 },
            (_, Some(mi_l)) => if mi_l.single { check_backward(mi_l.ref_frame[0]) as usize } else { 3 },
            _ => 1,
        };
        let mut cdf: &mut [i32] = &mut tile_cdf.comp_mode[ctx];
        self.encode_symbol(comp_mode as i32, out_bits, &mut cdf);
    }

    pub fn encode_skip_mode(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, skip_mode: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let ctx = match frame.get_above_mi(mi_row, mi_col) {
            Some(mi) => mi.skip_mode,
            None => 0,   
        } + match frame.get_left_mi(mi_row, mi_col) {
            Some(mi) => mi.skip_mode,
            None => 0,   
        };
        let mut cdf: &mut [i32] = &mut tile_cdf.skip_mode[ctx];
        self.encode_symbol(skip_mode as i32, out_bits, &mut cdf);
    }

    pub fn encode_skip(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, skip: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let ctx = match frame.get_above_mi(mi_row, mi_col) {
            Some(mi) => mi.skip as usize,
            None => 0,   
        } + match frame.get_left_mi(mi_row, mi_col) {
            Some(mi) => mi.skip as usize,
            None => 0,   
        };
        let mut cdf: &mut [i32] = &mut tile_cdf.skip[ctx];
        self.encode_symbol(skip as i32, out_bits, &mut cdf);
    }

    pub fn encode_comp_ref(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, comp_ref: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let last12_count = count_refs(LAST_FRAME, frame, mi_row, mi_col) + count_refs(LAST2_FRAME, frame, mi_row, mi_col);
        let last3_gold_count = count_refs(LAST3_FRAME, frame, mi_row, mi_col) + count_refs(GOLDEN_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(last12_count, last3_gold_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.comp_ref[ctx][0];
        self.encode_symbol(comp_ref as i32, out_bits, &mut cdf);
    }

    pub fn encode_comp_ref_p1(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, comp_ref_p1: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let last_count = count_refs(LAST_FRAME, frame, mi_row, mi_col);
        let last2_count = count_refs(LAST2_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(last_count, last2_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.comp_ref[ctx][1];
        self.encode_symbol(comp_ref_p1 as i32, out_bits, &mut cdf);
    }

    pub fn encode_comp_ref_p2(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, comp_ref_p2: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let last3_count = count_refs(LAST3_FRAME, frame, mi_row, mi_col);
        let gold_count = count_refs(GOLDEN_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(last3_count, gold_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.comp_ref[ctx][2];
        self.encode_symbol(comp_ref_p2 as i32, out_bits, &mut cdf);
    }

    pub fn encode_comp_bwdref(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, comp_bwdref: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let brfarf2_count = count_refs(BWDREF_FRAME, frame, mi_row, mi_col) + count_refs(ALTREF2_FRAME, frame, mi_row, mi_col);
        let arf_count = count_refs(ALTREF_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(brfarf2_count, arf_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.comp_bwd_ref[ctx][0];
        self.encode_symbol(comp_bwdref as i32, out_bits, &mut cdf);
    }

    pub fn encode_comp_bwdref_p1(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, comp_bwdref_p1: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let brf_count = count_refs(BWDREF_FRAME, frame, mi_row, mi_col);
        let arf2_count = count_refs(ALTREF2_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(brf_count, arf2_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.comp_bwd_ref[ctx][1];
        self.encode_symbol(comp_bwdref_p1 as i32, out_bits, &mut cdf);
    }

    pub fn encode_single_ref_p1(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, single_ref_p1: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let fwd_count = count_refs(LAST_FRAME, frame, mi_row, mi_col) +
                        count_refs(LAST2_FRAME, frame, mi_row, mi_col) +
                        count_refs(LAST3_FRAME, frame, mi_row, mi_col) +
                        count_refs(GOLDEN_FRAME, frame, mi_row, mi_col);
        let bwd_count = count_refs(BWDREF_FRAME, frame, mi_row, mi_col) +
                        count_refs(ALTREF2_FRAME, frame, mi_row, mi_col) +
                        count_refs(ALTREF_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(fwd_count, bwd_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.single_ref[ctx][0];
        self.encode_symbol(single_ref_p1 as i32, out_bits, &mut cdf);
    }

    pub fn encode_single_ref_p2(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, single_ref_p2: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let brfarf2_count = count_refs(BWDREF_FRAME, frame, mi_row, mi_col) + count_refs(ALTREF2_FRAME, frame, mi_row, mi_col);
        let arf_count = count_refs(ALTREF_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(brfarf2_count, arf_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.single_ref[ctx][1];
        self.encode_symbol(single_ref_p2 as i32, out_bits, &mut cdf);
    }

    pub fn encode_single_ref_p3(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, single_ref_p3: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let last12_count = count_refs(LAST_FRAME, frame, mi_row, mi_col) + count_refs(LAST2_FRAME, frame, mi_row, mi_col);
        let last3_gold_count = count_refs(LAST3_FRAME, frame, mi_row, mi_col) + count_refs(GOLDEN_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(last12_count, last3_gold_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.single_ref[ctx][2];
        self.encode_symbol(single_ref_p3 as i32, out_bits, &mut cdf);
    }

    pub fn encode_single_ref_p4(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, single_ref_p4: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let last_count = count_refs(LAST_FRAME, frame, mi_row, mi_col);
        let last2_count = count_refs(LAST2_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(last_count, last2_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.single_ref[ctx][3];
        self.encode_symbol(single_ref_p4 as i32, out_bits, &mut cdf);
    }

    pub fn encode_single_ref_p5(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, single_ref_p5: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let last3_count = count_refs(LAST3_FRAME, frame, mi_row, mi_col);
        let gold_count = count_refs(GOLDEN_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(last3_count, gold_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.single_ref[ctx][4];
        self.encode_symbol(single_ref_p5 as i32, out_bits, &mut cdf);
    }

    pub fn encode_single_ref_p6(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, single_ref_p6: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let brf_count = count_refs(BWDREF_FRAME, frame, mi_row, mi_col);
        let arf2_count = count_refs(ALTREF2_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(brf_count, arf2_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.single_ref[ctx][5];
        self.encode_symbol(single_ref_p6 as i32, out_bits, &mut cdf);
    }

    pub fn encode_compound_mode(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, compound_mode: usize, ref_mv_context: usize, new_mv_context: usize) {
        const compound_mode_ctx_map: [[usize; COMP_NEWMV_CTXS]; 3] = [
            [ 0, 1, 1, 1, 1 ],
            [ 1, 2, 3, 4, 4 ],
            [ 4, 4, 5, 6, 7 ]
        ];
        let ctx = compound_mode_ctx_map[ref_mv_context>>1][Min!(new_mv_context, COMP_NEWMV_CTXS-1)];
        let mut cdf: &mut [i32] = &mut tile_cdf.compound_mode[ctx];
        self.encode_symbol(compound_mode as i32, out_bits, &mut cdf);
    }

    pub fn encode_interp_filter(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, interp_filter: InterpFilter, frame: &mut Frame, mi_row: usize, mi_col: usize, dir: usize) {
        let ref_frame = frame.get_mi(mi_row, mi_col).unwrap().ref_frame;
        let mut ctx = ((dir&1)*2 + (ref_frame[1]>INTRA_FRAME) as usize) * 4;
        let left_type = if let Some(mi) = frame.get_left_mi(mi_row, mi_col) {
            if mi.ref_frame[0]==ref_frame[0] || mi.ref_frame[1]==ref_frame[1] {
                mi.interp_filter[dir]
            } else { BILINEAR }
        } else { BILINEAR };
        let above_type = if let Some(mi) = frame.get_above_mi(mi_row, mi_col) {
            if mi.ref_frame[0]==ref_frame[0] || mi.ref_frame[1]==ref_frame[1] {
                mi.interp_filter[dir]
            } else { BILINEAR }
        } else { BILINEAR };
        ctx += if left_type==above_type {
            left_type as usize
        } else if left_type==BILINEAR {
            above_type as usize
        } else if above_type==BILINEAR {
            left_type as usize
        } else { BILINEAR as usize };
        let mut cdf: &mut [i32] = &mut tile_cdf.interp_filter[ctx];
        self.encode_symbol(interp_filter as i32, out_bits, &mut cdf);
    }

    pub fn encode_motion_mode(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, motion_mode: usize, mi_size: BlockSize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.motion_mode[mi_size as usize];
        self.encode_symbol(motion_mode as i32, out_bits, &mut cdf);
    }

    pub fn encode_mv_joint(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, mv_joint: MvJoint, mv_ctx: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.mv_joint[mv_ctx as usize];
        self.encode_symbol(mv_joint as i32, out_bits, &mut cdf);
    }

    pub fn encode_mv_sign(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, mv_sign: MvJoint, mv_ctx: usize, comp: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.mv_sign[mv_ctx as usize][comp];
        self.encode_symbol(mv_sign as i32, out_bits, &mut cdf);
    }

    pub fn encode_mv_class(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, mv_class: usize, mv_ctx: usize, comp: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.mv_class[mv_ctx as usize][comp];
        self.encode_symbol(mv_class as i32, out_bits, &mut cdf);
    }

    pub fn encode_mv_class0_bit(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, mv_class0_bit: usize, mv_ctx: usize, comp: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.mv_class0_bit[mv_ctx as usize][comp];
        self.encode_symbol(mv_class0_bit as i32, out_bits, &mut cdf);
    }

    pub fn encode_mv_class0_fr(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, mv_class0_fr: usize, mv_class0_bit: usize, mv_ctx: usize, comp: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.mv_class0_fr[mv_ctx as usize][comp][mv_class0_bit];
        self.encode_symbol(mv_class0_fr as i32, out_bits, &mut cdf);
    }

    pub fn encode_mv_class0_hp(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, mv_class0_hp: usize, mv_ctx: usize, comp: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.mv_class0_hp[mv_ctx as usize][comp];
        self.encode_symbol(mv_class0_hp as i32, out_bits, &mut cdf);
    }

    pub fn encode_mv_fr(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, mv_fr: usize, mv_ctx: usize, comp: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.mv_fr[mv_ctx as usize][comp];
        self.encode_symbol(mv_fr as i32, out_bits, &mut cdf);
    }

    pub fn encode_mv_hp(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, mv_hp: usize, mv_ctx: usize, comp: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.mv_hp[mv_ctx as usize][comp];
        self.encode_symbol(mv_hp as i32, out_bits, &mut cdf);
    }

    pub fn encode_mv_bit(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, mv_bit: usize, i: usize, mv_ctx: usize, comp: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.mv_bit[mv_ctx as usize][comp][i];
        self.encode_symbol(mv_bit as i32, out_bits, &mut cdf);
    }

    pub fn encode_all_zero(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, all_zero: usize, plane: usize, tx_size: TxSize, tx_size_ctx: usize, frame: &mut Frame, mi_size: BlockSize, mi_row: usize, mi_col: usize, w4: usize, h4: usize, x4: usize, y4: usize) {
        let mut max_x4 = frame.mi_cols();
        let mut max_y4 = frame.mi_rows();
        if plane>0 {
            max_x4 >>= frame.subsampling_x as usize;
            max_y4 >>= frame.subsampling_y as usize;
        }
        let w = tx_width[tx_size as usize];
        let h = tx_height[tx_size as usize];
        let bsize = get_plane_residual_size(mi_size, plane, frame.subsampling_x, frame.subsampling_y);
        let bw = block_width[bsize as usize];
        let bh = block_height[bsize as usize];
        let mut ctx;
        if plane==0 {
            let mut top = 0;
            let mut left = 0;
            for k in 0..w4 {
                if x4+k < max_x4 {
                    top = Max!(top, match frame.get_above_mi(mi_row, mi_col) {
                        Some(mi) => mi.level_context[plane],
                        None => 0,
                    });
                }
            }
            for k in 0..h4 {
                if y4+k < max_y4 {
                    left = Max!(left, match frame.get_left_mi(mi_row, mi_col) {
                        Some(mi) => mi.level_context[plane],
                        None => 0,
                    });
                }
            }
            top = Min!(top, 255);
            left = Min!(left, 255);
            ctx = if bw==w && bh==h {
                0
            } else if top==0 && left==0 {
                1
            } else if top==0 || left==0 {
                2 + (Max!(top, left)>3) as usize
            } else if Max!(top, left) <= 3 {
                4
            } else if Min!(top, left) <= 3 {
                5
            } else { 6 };
        }
        else {
            let mut above = 0;
            let mut left = 0;
            for i in 0..w4 {
                if x4+i < max_x4 {
                    above |= match frame.get_above_mi(mi_row, mi_col) {
                        Some(mi) => mi.level_context[plane],
                        None => 0,
                    } | match frame.get_above_mi(mi_row, mi_col) {
                        Some(mi) => mi.dc_context[plane],
                        None => 0,
                    };
                }
            }
            for i in 0..h4 {
                if y4+i < max_y4 {
                    left |= match frame.get_left_mi(mi_row, mi_col) {
                        Some(mi) => mi.level_context[plane],
                        None => 0,
                    } | match frame.get_left_mi(mi_row, mi_col) {
                        Some(mi) => mi.dc_context[plane],
                        None => 0,
                    };
                }
            }
            ctx = (above!=0) as usize + (left!=0) as usize;
            ctx += 7;
            if bw*bh > w*h { ctx += 3; }
        }
        let mut cdf: &mut [i32] = &mut tile_cdf.txb_skip[0/*FIXME*/][tx_size_ctx][ctx];
        self.encode_symbol(all_zero as i32, out_bits, &mut cdf);
    }

    pub fn encode_eob_pt_16(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, eob_pt_16: usize, frame: &mut Frame, mi_row: usize, mi_col: usize, uv_mode: YMode, lossless: bool, is_inter: bool, reduced_tx_set: bool, plane: usize, tx_size: TxSize, x4: usize, y4: usize, ptype: usize) {
        let tx_type = compute_tx_type(plane, tx_size, x4, y4, frame, lossless, is_inter, reduced_tx_set, uv_mode, mi_row, mi_col);
        let ctx = if get_tx_class(tx_type)==TX_CLASS_2D { 0 } else { 1 };
        let mut cdf: &mut [i32] = &mut tile_cdf.eob_pt_16[0/*FIXME*/][ptype][ctx];
        self.encode_symbol(eob_pt_16 as i32, out_bits, &mut cdf);
    }

    pub fn encode_eob_pt_32(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, eob_pt_32: usize, frame: &mut Frame, mi_row: usize, mi_col: usize, uv_mode: YMode, lossless: bool, is_inter: bool, reduced_tx_set: bool, plane: usize, tx_size: TxSize, x4: usize, y4: usize, ptype: usize) {
        let tx_type = compute_tx_type(plane, tx_size, x4, y4, frame, lossless, is_inter, reduced_tx_set, uv_mode, mi_row, mi_col);
        let ctx = if get_tx_class(tx_type)==TX_CLASS_2D { 0 } else { 1 };
        let mut cdf: &mut [i32] = &mut tile_cdf.eob_pt_32[0/*FIXME*/][ptype][ctx];
        self.encode_symbol(eob_pt_32 as i32, out_bits, &mut cdf);
    }

    pub fn encode_eob_pt_64(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, eob_pt_64: usize, frame: &mut Frame, mi_row: usize, mi_col: usize, uv_mode: YMode, lossless: bool, is_inter: bool, reduced_tx_set: bool, plane: usize, tx_size: TxSize, x4: usize, y4: usize, ptype: usize) {
        let tx_type = compute_tx_type(plane, tx_size, x4, y4, frame, lossless, is_inter, reduced_tx_set, uv_mode, mi_row, mi_col);
        let ctx = if get_tx_class(tx_type)==TX_CLASS_2D { 0 } else { 1 };
        let mut cdf: &mut [i32] = &mut tile_cdf.eob_pt_64[0/*FIXME*/][ptype][ctx];
        self.encode_symbol(eob_pt_64 as i32, out_bits, &mut cdf);
    }

    pub fn encode_eob_pt_128(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, eob_pt_128: usize, frame: &mut Frame, mi_row: usize, mi_col: usize, uv_mode: YMode, lossless: bool, is_inter: bool, reduced_tx_set: bool, plane: usize, tx_size: TxSize, x4: usize, y4: usize, ptype: usize) {
        let tx_type = compute_tx_type(plane, tx_size, x4, y4, frame, lossless, is_inter, reduced_tx_set, uv_mode, mi_row, mi_col);
        let ctx = if get_tx_class(tx_type)==TX_CLASS_2D { 0 } else { 1 };
        let mut cdf: &mut [i32] = &mut tile_cdf.eob_pt_128[0/*FIXME*/][ptype][ctx];
        self.encode_symbol(eob_pt_128 as i32, out_bits, &mut cdf);
    }

    pub fn encode_eob_pt_256(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, eob_pt_256: usize, frame: &mut Frame, mi_row: usize, mi_col: usize, uv_mode: YMode, lossless: bool, is_inter: bool, reduced_tx_set: bool, plane: usize, tx_size: TxSize, x4: usize, y4: usize, ptype: usize) {
        let tx_type = compute_tx_type(plane, tx_size, x4, y4, frame, lossless, is_inter, reduced_tx_set, uv_mode, mi_row, mi_col);
        let ctx = if get_tx_class(tx_type)==TX_CLASS_2D { 0 } else { 1 };
        let mut cdf: &mut [i32] = &mut tile_cdf.eob_pt_256[0/*FIXME*/][ptype][ctx];
        self.encode_symbol(eob_pt_256 as i32, out_bits, &mut cdf);
    }

    pub fn encode_eob_pt_512(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, eob_pt_512: usize, frame: &mut Frame, mi_row: usize, mi_col: usize, uv_mode: YMode, lossless: bool, is_inter: bool, reduced_tx_set: bool, plane: usize, tx_size: TxSize, x4: usize, y4: usize, ptype: usize) {
        let tx_type = compute_tx_type(plane, tx_size, x4, y4, frame, lossless, is_inter, reduced_tx_set, uv_mode, mi_row, mi_col);
        let ctx = if get_tx_class(tx_type)==TX_CLASS_2D { 0 } else { 1 };
        let mut cdf: &mut [i32] = &mut tile_cdf.eob_pt_512[0/*FIXME*/][ptype][ctx];
        self.encode_symbol(eob_pt_512 as i32, out_bits, &mut cdf);
    }

    pub fn encode_eob_pt_1024(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, eob_pt_1024: usize, frame: &mut Frame, mi_row: usize, mi_col: usize, uv_mode: YMode, lossless: bool, is_inter: bool, reduced_tx_set: bool, plane: usize, tx_size: TxSize, x4: usize, y4: usize, ptype: usize) {
        let tx_type = compute_tx_type(plane, tx_size, x4, y4, frame, lossless, is_inter, reduced_tx_set, uv_mode, mi_row, mi_col);
        let ctx = if get_tx_class(tx_type)==TX_CLASS_2D { 0 } else { 1 };
        let mut cdf: &mut [i32] = &mut tile_cdf.eob_pt_1024[0/*FIXME*/][ptype][ctx];
        self.encode_symbol(eob_pt_1024 as i32, out_bits, &mut cdf);
    }

    pub fn encode_eob_extra(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, eob_extra: usize, tx_size_ctx: usize, ptype: usize, eob_pt: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.eob_extra[0/*FIXME*/][tx_size_ctx][ptype][eob_pt];
        self.encode_symbol(eob_extra as i32, out_bits, &mut cdf);
    }

    pub fn encode_coeff_base(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, coeff_base: usize, tx_size: TxSize, tx_size_ctx: usize, ptype: usize, plane: usize, x4: usize, y4: usize, scan: &[usize], c: usize, frame: &mut Frame, lossless: bool, is_inter: bool, reduced_tx_set: bool, uv_mode: YMode, mi_row: usize, mi_col: usize, quant: &[isize]) {
        let ctx = get_coeff_base_ctx(tx_size, plane, x4, y4, scan[c], c, false, frame, lossless, is_inter, reduced_tx_set, uv_mode, mi_row, mi_col, quant);
        let mut cdf: &mut [i32] = &mut tile_cdf.coeff_base[0/*FIXME*/][tx_size_ctx][ptype][ctx];
        self.encode_symbol(coeff_base as i32, out_bits, &mut cdf);
    }

    pub fn encode_coeff_base_eob(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, coeff_base_eob: usize, tx_size: TxSize, tx_size_ctx: usize, ptype: usize, plane: usize, x4: usize, y4: usize, scan: &[usize], c: usize, frame: &mut Frame, lossless: bool, is_inter: bool, reduced_tx_set: bool, uv_mode: YMode, mi_row: usize, mi_col: usize, quant: &[isize]) {
        let ctx = get_coeff_base_ctx(tx_size, plane, x4, y4, scan[c], c, true, frame, lossless, is_inter, reduced_tx_set, uv_mode, mi_row, mi_col, quant) - SIG_COEF_CONTEXTS + SIG_COEF_CONTEXTS_EOB;
        let mut cdf: &mut [i32] = &mut tile_cdf.coeff_base_eob[0/*FIXME*/][tx_size_ctx][ptype][ctx];
        self.encode_symbol(coeff_base_eob as i32, out_bits, &mut cdf);
    }

    pub fn encode_dc_sign(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, ptype: usize, frame: &mut Frame, plane: usize, mi_row: usize, mi_col: usize, w4: usize, h4: usize, x4: usize, y4: usize) {
        let mut max_x4 = frame.mi_cols();
        let mut max_y4 = frame.mi_rows();
        if plane>0 {
            max_x4 >>= frame.subsampling_x as usize;
            max_y4 >>= frame.subsampling_y as usize;
        }
        let mut dc_sign: isize = 0;
        for k in 0..w4 {
            if x4+k < max_x4 {
                let sign = match frame.get_above_mi(mi_row, mi_col) {
                    Some(mi) => mi.dc_context[plane],
                    None => 0,
                };
                if sign==1 {
                    dc_sign -= 1;
                } else if sign==2 {
                    dc_sign += 1;
                }
            }
        }
        for k in 0..h4 {
            if y4+k < max_y4 {
                let sign = match frame.get_left_mi(mi_row, mi_col) {
                    Some(mi) => mi.dc_context[plane],
                    None => 0,
                };
                if sign==1 {
                    dc_sign -= 1;
                } else if sign==2 {
                    dc_sign += 1;
                }
            }
        }
        let ctx = if dc_sign<0 { 1 } else if dc_sign>0 {
            2
        } else { 0 };
        let mut cdf: &mut [i32] = &mut tile_cdf.dc_sign[0/*FIXME*/][ptype][ctx];
        self.encode_symbol(dc_sign as i32, out_bits, &mut cdf);
    }

    pub fn encode_coeff_br(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, coeff_br: usize, frame: &mut Frame, tx_size_ctx: usize, lossless: bool, is_inter: bool, reduced_tx_set: bool, uv_mode: YMode, pos: usize, ptype: usize, tx_size: TxSize, x4: usize, y4: usize, plane: usize, quant: &[isize], mi_row: usize, mi_col: usize) {
        let adj_tx_size = adjusted_tx_size[tx_size as usize];
        let bwl = tx_width_log2[adj_tx_size as usize];
        let txw = tx_width[adj_tx_size as usize];
        let txh = tx_height[adj_tx_size as usize];
        let row = pos >> bwl;
        let col = pos - (row << bwl);
        let mut mag = 0;
        let tx_type = compute_tx_type(plane, tx_size, x4, y4, frame, lossless, is_inter, reduced_tx_set, uv_mode, mi_row, mi_col);
        let tx_class = get_tx_class(tx_type);
        for idx in 0..3 {
            let ref_row = row + mag_ref_offset_with_tx_class[tx_class][idx][0];
            let ref_col = col + mag_ref_offset_with_tx_class[tx_class][idx][1];
            if /*ref_row>=0 && ref_col>=0 && FIXME*/ref_row<txh && ref_col<(1<<bwl) {
                mag += Min!(quant[ref_row*txw + ref_col] as usize, COEFF_BASE_RANGE + NUM_BASE_LEVELS + 1);
            }
        }
        mag = Min!((mag+1) >> 1, 6);
        let ctx = if pos==0 {
            mag
        } else if tx_class==0 {
            if row<2 && col<2 {
                mag + 7
            } else { mag + 14 }
        } else {
            if tx_class==1 {
                if col==0 {
                    mag + 7
                } else {
                    mag + 14
                }
            } else {
                if row==0 {
                    mag + 7
                } else {
                    mag + 14
                }
            }
        };
        let mut cdf: &mut [i32] = &mut tile_cdf.coeff_br[0/*FIXME*/][Min!(tx_size_ctx, TX_32X32 as usize)][ptype][ctx];
        self.encode_symbol(coeff_br as i32, out_bits, &mut cdf);
    }

    pub fn encode_has_palette_y(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, has_palette_y: usize, frame: &mut Frame, mi_row: usize, mi_col: usize, bsize_ctx: usize) {
        let mut ctx = 0;
        ctx += match frame.get_above_mi(mi_row, mi_col) {
            Some(mi) => if mi.palette_size[0]>0 { 1 } else { 0 },
            None => 0,
        } + match frame.get_left_mi(mi_row, mi_col) {
            Some(mi) => if mi.palette_size[0]>0 { 1 } else { 0 },
            None => 0,
        };
        let mut cdf: &mut [i32] = &mut tile_cdf.palette_y_mode[bsize_ctx][ctx];
        self.encode_symbol(has_palette_y as i32, out_bits, &mut cdf);
    }

    pub fn encode_has_palette_uv(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, has_palette_uv: usize, palette_size_y: usize) {
        let ctx = if palette_size_y>0 { 1 } else { 0 };
        let mut cdf: &mut [i32] = &mut tile_cdf.palette_uv_mode[ctx];
        self.encode_symbol(has_palette_uv as i32, out_bits, &mut cdf);
    }

    pub fn encode_palette_size_y_minus_2(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, palette_size_y_minus_2: usize, bsize_ctx: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.palette_y_size[bsize_ctx];
        self.encode_symbol(palette_size_y_minus_2 as i32, out_bits, &mut cdf);
    }

    pub fn encode_palette_size_uv_minus_2(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, palette_size_uv_minus_2: usize, bsize_ctx: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.palette_uv_size[bsize_ctx];
        self.encode_symbol(palette_size_uv_minus_2 as i32, out_bits, &mut cdf);
    }

    pub fn encode_palette_color_idx_y(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, palette_color_idx_y: usize, color_context_hash: usize, palette_size_y: usize) {
        let ctx = palette_color_context[color_context_hash];
        assert!(ctx>=0);
        let mut cdf: &mut [i32] = match palette_size_y {
            2 => &mut tile_cdf.palette_size_2_y_color[ctx as usize],
            3 => &mut tile_cdf.palette_size_3_y_color[ctx as usize],
            4 => &mut tile_cdf.palette_size_4_y_color[ctx as usize],
            5 => &mut tile_cdf.palette_size_5_y_color[ctx as usize],
            6 => &mut tile_cdf.palette_size_6_y_color[ctx as usize],
            7 => &mut tile_cdf.palette_size_7_y_color[ctx as usize],
            _ => &mut tile_cdf.palette_size_8_y_color[ctx as usize],
        };
        self.encode_symbol(palette_color_idx_y as i32, out_bits, &mut cdf);
    }

    pub fn encode_palette_color_idx_uv(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, palette_color_idx_uv: usize, color_context_hash: usize, palette_size_uv: usize) {
        let ctx = palette_color_context[color_context_hash];
        assert!(ctx>=0);
        let mut cdf: &mut [i32] = match palette_size_uv {
            2 => &mut tile_cdf.palette_size_2_uv_color[ctx as usize],
            3 => &mut tile_cdf.palette_size_3_uv_color[ctx as usize],
            4 => &mut tile_cdf.palette_size_4_uv_color[ctx as usize],
            5 => &mut tile_cdf.palette_size_5_uv_color[ctx as usize],
            6 => &mut tile_cdf.palette_size_6_uv_color[ctx as usize],
            7 => &mut tile_cdf.palette_size_7_uv_color[ctx as usize],
            _ => &mut tile_cdf.palette_size_8_uv_color[ctx as usize],
        };
        self.encode_symbol(palette_color_idx_uv as i32, out_bits, &mut cdf);
    }

    pub fn encode_delta_q_abs(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, delta_q_abs: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.delta_q;
        self.encode_symbol(delta_q_abs as i32, out_bits, &mut cdf);
    }

    pub fn encode_delta_lf_abs(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, delta_lf_abs: usize, delta_lf_multi: usize, i: usize) {
        let mut cdf: &mut [i32] = if delta_lf_multi==0 {
            &mut tile_cdf.delta_lf
        } else {
            &mut tile_cdf.delta_lf_multi[i]
        };
        self.encode_symbol(delta_lf_abs as i32, out_bits, &mut cdf);
    }

    pub fn encode_intra_tx_type(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, intra_tx_type: TxType, tx_set: usize, tx_size: TxSize, use_filter_intra: bool, filter_intra_mode: usize, y_mode: YMode) {
        let intra_dir = if use_filter_intra { filter_intra_mode_to_intra_dir[filter_intra_mode] as usize } else {
            y_mode as usize
        };
        let mut cdf: &mut [i32] = if tx_set==TX_SET_INTRA_1 {
            &mut tile_cdf.intra_tx_type_set1[tx_size_sqr[tx_size as usize] as usize][intra_dir]
        } else {
            &mut tile_cdf.intra_tx_type_set2[tx_size_sqr[tx_size as usize] as usize][intra_dir]
        };
        self.encode_symbol(intra_tx_type as i32, out_bits, &mut cdf);
    }

    pub fn encode_inter_tx_type(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, inter_tx_type: TxType, tx_size: TxSize, tx_set: usize) {
        let mut cdf: &mut [i32] = if tx_set==TX_SET_INTER_1 {
            &mut tile_cdf.inter_tx_type_set1[tx_size_sqr[tx_size as usize] as usize]
        } else if tx_set==TX_SET_INTER_2 {
            &mut tile_cdf.inter_tx_type_set2[tx_size_sqr[tx_size as usize] as usize]
        } else {
            &mut tile_cdf.inter_tx_type_set3[tx_size_sqr[tx_size as usize] as usize]
        };
        self.encode_symbol(inter_tx_type as i32, out_bits, &mut cdf);
    }

    pub fn encode_comp_ref_type(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, comp_ref_type: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let (avail_u, above0, above1, above_comp_inter, above_intra) = match frame.get_above_mi(mi_row, mi_col) {
            Some(mi) => (true, mi.ref_frame[0], mi.ref_frame[1], !mi.intra&&!mi.single, mi.intra),
            None => (false, INTRA_FRAME, INTRA_FRAME, false, false),
        };
        let (avail_l, left0, left1, left_comp_inter, left_intra) = match frame.get_left_mi(mi_row, mi_col) {
            Some(mi) => (true, mi.ref_frame[0], mi.ref_frame[1], !mi.intra&&!mi.single, mi.intra),
            None => (false, INTRA_FRAME, INTRA_FRAME, false, false),
        };
        let above_uni_comp = above_comp_inter && is_samedir_ref_pair(above0, above1);
        let left_uni_comp = left_comp_inter && is_samedir_ref_pair(left0, left1);
        let ctx = if avail_u && !above_intra && avail_l && !left_intra {
            let samedir = is_samedir_ref_pair(above0, left0) as usize;
            if !above_comp_inter && !left_comp_inter {
                1 + 2*samedir
            } else if !above_comp_inter {
                if !left_uni_comp { 1 } else { 3 + samedir }
            } else if !left_comp_inter {
                if !above_uni_comp { 1 } else { 3 + samedir }
            } else {
                if !above_uni_comp && !left_uni_comp {
                    0
                } else if !above_uni_comp || left_uni_comp {
                    2
                } else { 3 + ((above0==BWDREF_FRAME) as usize == (left0==BWDREF_FRAME) as usize) as usize }
            }
        } else if avail_u && avail_l {
            if above_comp_inter {
                1 + 2*above_uni_comp as usize
            } else if left_comp_inter {
                1 + 2*left_uni_comp as usize
            } else { 2 }
        } else if above_comp_inter {
            4 * above_uni_comp as usize
        } else if left_comp_inter {
            4 * left_uni_comp as usize
        } else { 2 };
        let mut cdf: &mut [i32] = &mut tile_cdf.comp_ref_type[ctx];
        self.encode_symbol(comp_ref_type as i32, out_bits, &mut cdf);
    }

    pub fn encode_uni_comp_ref(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, uni_comp_ref: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let fwd_count = count_refs(LAST_FRAME, frame, mi_row, mi_col) +
                        count_refs(LAST2_FRAME, frame, mi_row, mi_col) +
                        count_refs(LAST3_FRAME, frame, mi_row, mi_col) +
                        count_refs(GOLDEN_FRAME, frame, mi_row, mi_col);
        let bwd_count = count_refs(BWDREF_FRAME, frame, mi_row, mi_col) +
                        count_refs(ALTREF2_FRAME, frame, mi_row, mi_col) +
                        count_refs(ALTREF_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(fwd_count, bwd_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.uni_comp_ref[ctx][0];
        self.encode_symbol(uni_comp_ref as i32, out_bits, &mut cdf);
    }

    pub fn encode_uni_comp_ref_p1(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, uni_comp_ref_p1: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let last2_count = count_refs(LAST2_FRAME, frame, mi_row, mi_col);
        let last3_gold_count = count_refs(LAST3_FRAME, frame, mi_row, mi_col) + count_refs(GOLDEN_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(last2_count, last3_gold_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.uni_comp_ref[ctx][1];
        self.encode_symbol(uni_comp_ref_p1 as i32, out_bits, &mut cdf);
    }

    pub fn encode_uni_comp_ref_p2(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, uni_comp_ref_p2: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let last3_count = count_refs(LAST3_FRAME, frame, mi_row, mi_col);
        let gold_count = count_refs(GOLDEN_FRAME, frame, mi_row, mi_col);
        let ctx = ref_count_ctx(last3_count, gold_count);
        let mut cdf: &mut [i32] = &mut tile_cdf.uni_comp_ref[ctx][2];
        self.encode_symbol(uni_comp_ref_p2 as i32, out_bits, &mut cdf);
    }

    pub fn encode_comp_group_idx(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, comp_group_idx: usize, frame: &mut Frame, mi_row: usize, mi_col: usize) {
        let mut ctx = match frame.get_above_mi(mi_row, mi_col) {
            Some(mi) => if !mi.single { mi.comp_group_idx } else if mi.ref_frame[0]==ALTREF_FRAME { 3 } else { 0 },
            None => 0,
        } + match frame.get_left_mi(mi_row, mi_col) {
            Some(mi) => if !mi.single { mi.comp_group_idx } else if mi.ref_frame[0]==ALTREF_FRAME { 3 } else { 0 },
            None => 0,
        };
        ctx = Min!(5, ctx);
        let mut cdf: &mut [i32] = &mut tile_cdf.comp_group_idx[ctx];
        self.encode_symbol(comp_group_idx as i32, out_bits, &mut cdf);
    }

    pub fn encode_compound_idx(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, compound_idx: usize, order_hints: &[usize], order_hint: usize, ref_frame: &[RefFrame; 2], frame: &mut Frame, mi_row: usize, mi_col: usize, enable_order_hint: bool, order_hint_bits: usize) {
        let fwd = Abs!(get_relative_dist(order_hints[ref_frame[0] as usize], order_hint, enable_order_hint, order_hint_bits));
        let bwd = Abs!(get_relative_dist(order_hints[ref_frame[1] as usize], order_hint, enable_order_hint, order_hint_bits));
        let mut ctx = if fwd==bwd { 3 } else { 0 };
        ctx += match frame.get_above_mi(mi_row, mi_col) {
            Some(mi) => if !mi.single {
                mi.compound_idx
            } else if mi.ref_frame[0]==ALTREF_FRAME {
                1
            } else { 0 },
            None => 0,
        } + match frame.get_left_mi(mi_row, mi_col) {
            Some(mi) => if !mi.single {
                mi.compound_idx
            } else if mi.ref_frame[0]==ALTREF_FRAME {
                1
            } else { 0 },
            None => 0,
        };
        let mut cdf: &mut [i32] = &mut tile_cdf.compound_idx[ctx];
        self.encode_symbol(compound_idx as i32, out_bits, &mut cdf);
    }

    pub fn encode_compound_type(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, compound_type: usize, mi_size: BlockSize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.compound_type[mi_size as usize];
        self.encode_symbol(compound_type as i32, out_bits, &mut cdf);
    }

    pub fn encode_interintra(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, interintra: usize, mi_size: BlockSize) {
        let ctx = size_group[mi_size as usize];
        let mut cdf: &mut [i32] = &mut tile_cdf.inter_intra[ctx];
        self.encode_symbol(interintra as i32, out_bits, &mut cdf);
    }

    pub fn encode_interintra_mode(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, interintra_mode: usize, mi_size: BlockSize) {
        let ctx = size_group[mi_size as usize];
        let mut cdf: &mut [i32] = &mut tile_cdf.inter_intra_mode[ctx];
        self.encode_symbol(interintra_mode as i32, out_bits, &mut cdf);
    }

    pub fn encode_wedge_index(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, wedge_index: usize, mi_size: BlockSize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.wedge_index[mi_size as usize];
        self.encode_symbol(wedge_index as i32, out_bits, &mut cdf);
    }

    pub fn encode_wedge_interintra(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, wedge_interintra: usize, mi_size: BlockSize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.wedge_inter_intra[mi_size as usize];
        self.encode_symbol(wedge_interintra as i32, out_bits, &mut cdf);
    }

    pub fn encode_use_obmc(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, use_obmc: usize, mi_size: BlockSize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.use_obmc[mi_size as usize];
        self.encode_symbol(use_obmc as i32, out_bits, &mut cdf);
    }

    pub fn encode_cfl_alpha_signs(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, cfl_alpha_signs: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.cfl_sign;
        self.encode_symbol(cfl_alpha_signs as i32, out_bits, &mut cdf);
    }

    pub fn encode_cfl_alpha_u(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, cfl_alpha_u: usize, cfl_alpha_signs: usize) {
        let ctx = if cfl_alpha_signs>=2 { cfl_alpha_signs-1 } else { assert!(false); 0 };
        let mut cdf: &mut [i32] = &mut tile_cdf.cfl_alpha[ctx];
        self.encode_symbol(cfl_alpha_u as i32, out_bits, &mut cdf);
    }

    pub fn encode_cfl_alpha_v(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, cfl_alpha_v: usize, cfl_alpha_signs: usize) {
        let ctx = if cfl_alpha_signs%3<2 { (cfl_alpha_signs%3)*3+cfl_alpha_signs/3-1 } else { assert!(false); 0 };
        let mut cdf: &mut [i32] = &mut tile_cdf.cfl_alpha[ctx];
        self.encode_symbol(cfl_alpha_v as i32, out_bits, &mut cdf);
    }

    pub fn encode_use_wiener(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, use_wiener: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.use_wiener;
        self.encode_symbol(use_wiener as i32, out_bits, &mut cdf);
    }

    pub fn encode_use_sgrproj(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, use_sgrproj: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.use_sgrproj;
        self.encode_symbol(use_sgrproj as i32, out_bits, &mut cdf);
    }

    pub fn encode_restoration_type(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF, restoration_type: usize) {
        let mut cdf: &mut [i32] = &mut tile_cdf.restoration_type;
        self.encode_symbol(restoration_type as i32, out_bits, &mut cdf);
    }

    //pub fn encode_hoge(&mut self, out_bits: &mut Vec<u8>, tile_cdf: &mut CDF) {
    //    let mut cdf: &mut [i32] = &mut tile_cdf.;
    //    self.encode_symbol( as i32, out_bits, &mut cdf);
    //}
    /*
    void init_cdf();
    void load_cdf(int ctx);
    void save_cdf(int ctx);
    void average_cdf();
    */
}