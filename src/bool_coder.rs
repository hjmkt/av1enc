#![allow(dead_code)]
#![allow(non_snake_case)]

use std::collections::VecDeque;
use util::*;

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

    fn push_bit(out_bits: &mut Vec<u8>, bit: u8) {
        out_bits.push(bit);
    }

    fn push_bits_with_size(out_bits: &mut Vec<u8>, bits: u32, size: u8) {
        for i in (0..size).rev() { out_bits.push(((bits>>i)&1) as u8); }
    }

    fn push_bits(out_bits: &mut Vec<u8>, bits: &Vec<u8>) {
        for bit in bits { out_bits.push(*bit); }
    }

    fn push_bits_align(out_bits: &mut Vec<u8>) {
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

    fn decode_symbolp(&mut self, in_bits: &mut VecDeque<u8>, p: i32) -> u64 {
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
    void init_cdf();
    void load_cdf(int ctx);
    void save_cdf(int ctx);
    void average_cdf();
    void encode_intra_frame_y_mode(std::vector<int> &out_bits, IntraMode mode, IntraMode above_mode, IntraMode left_mode);
    void encode_y_mode(std::vector<int> &out_bits, int mode, int mi_size);
    void encode_uv_mode(std::vector<int> &out_bits, int mode, int mi_size);
    void encode_angle_delta_y(std::vector<int> &out_bits, int angle_delta_y, int mi_size);
    void encode_angle_delta_uv(std::vector<int> &out_bits, int angle_delta_uv, int mi_size);
    void encode_partition(std::vector<int> &out_bits, Partition partition, int b_size, Block* above_block, Block* left_block);
    void encode_split_or_horz(std::vector<int> &out_bits, int split, int b_size, Block* above_block, Block* left_block);
    void encode_split_or_vert(std::vector<int> &out_bits, int split, int b_size, Block* above_block, Block* left_block);
    void encode_tx_size(std::vector<int> &out_bits, TxSize tx_size, TxSize max_tx_size, Block* above_block, Block* left_block);
    void encode_txfm_split(std::vector<int> &out_bits, int split, int mi_size, TxSize tx_size, Block *above_block, Block *left_block);
    void encode_segment_id(std::vector<int> &out_bits, int segment_id);
    void encode_seg_id_predicted(std::vector<int> &out_bits, int seg_id_predicted, Block *above_block, Block *left_block);
    void encode_new_mv(std::vector<int> &out_bits, int new_mv, int new_mv_context);
    void encode_zero_mv(std::vector<int> &out_bits, int zero_mv, int zero_mv_context);
    void encode_ref_mv(std::vector<int> &out_bits, int ref_mv, int ref_mv_context);
    void encode_drl_mode(std::vector<int> &out_bits, int drl_mode, int drl_context);
    void encode_is_inter(std::vector<int> &out_bits, int is_inter, Block *above_block, Block *left_block);
    void encode_comp_mode(std::vector<int> &out_bits, int comp_mode, Block *above_block, Block *left_block);
    void encode_skip(std::vector<int> &out_bits, int skip, Block *above_block, Block *left_block);
    void encode_comp_ref(std::vector<int> &out_bits, int comp_ref, Block *above_block, Block *left_block);
    void encode_comp_ref_p1(std::vector<int> &out_bits, int comp_ref, int fwd_ref_sign_idx, Block *above_block, Block *left_block);
    void encode_comp_ref_p2(std::vector<int> &out_bits, int comp_ref, int fwd_ref_sign_idx, Block *above_block, Block *left_block);
    void encode_comp_bwdref(std::vector<int> &out_bits, int comp_bwdref, Block* block);
    void encode_comp_bwdref_p1(std::vector<int> &out_bits, int comp_bwdref, Block* block);
    void encode_single_ref_p1(std::vector<int> &out_bits, int single_ref, Block *above_block, Block *left_block);
    void encode_single_ref_p2(std::vector<int> &out_bits, int single_ref, Block *block);
    void encode_single_ref_p3(std::vector<int> &out_bits, int single_ref, Block *above_block, Block *left_block);
    void encode_single_ref_p4(std::vector<int> &out_bits, int single_ref, Block *above_block, Block *left_block);
    void encode_single_ref_p5(std::vector<int> &out_bits, int single_ref, Block *above_block, Block *left_block);
    void encode_single_ref_p6(std::vector<int> &out_bits, int single_ref, Block *block);
    void encode_compound_mode(std::vector<int> &out_bits, int compound_mode, int compound_context);
    void encode_interp_filter(std::vector<int> &out_bits, int dir, Block* block, Block *above_block, Block *left_block);
    void encode_motion_mode(std::vector<int> &out_bits, int motion_mode, int mi_size);
    void encode_mv_joint(std::vector<int> &out_bits, MvJoint mv_joint, int mv_ctx);
    void encode_mv_sign(std::vector<int> &out_bits, int mv_sign, int mv_ctx);
    void encode_mv_class(std::vector<int> &out_bits, MvClass mv_class, int mv_ctx, int comp);
    void encode_mv_class0_bit(std::vector<int> &out_bits, int mv_class0_bit, int mv_ctx, int comp);
    void encode_mv_class0_fr(std::vector<int> &out_bits, int mv_fr, int mv_ctx, int mv_class0_bit, int comp);
    void encode_mv_class0_hp(std::vector<int> &out_bits, int mv_class0_hp, int comp);
    void encode_mv_fr(std::vector<int> &out_bits, int mv_fr, int mv_ctx, int comp);
    void encode_mv_hp(std::vector<int> &out_bits, int mv_hp, int mv_ctx, int comp);
    void encode_mv_bit(std::vector<int> &out_bits, int mv_bit, int mv_ctx, int comp);
    void encode_tail_token(std::vector<int> &out_bits, Token token, Block* above_block, Block* left_block, int c, TxSize tx_size, int plane, bool is_inter, int band, int sx, int sy, int start_x, int start_y, int mi_cols, int mi_rows, int pos, TxType plane_tx_type, int* token_cache);
    void encode_dc_head_token(std::vector<int> &out_bits, Token token, Block* above_block, Block* left_block, int c, TxSize tx_size, int plane, bool is_inter, int band, int sx, int sy, int start_x, int start_y, int mi_cols, int mi_rows, int pos, TxType plane_tx_type, int* token_cache);
    void encode_ac_head_token(std::vector<int> &out_bits, Token token, Block* above_block, Block* left_block, int c, TxSize tx_size, int plane, bool is_inter, int band, int sx, int sy, int start_x, int start_y, int mi_cols, int mi_rows, int pos, TxType plane_tx_type, int* token_cache);
    void encode_has_palette_y(std::vector<int> &out_bits, int has_palette_y, MiSize mi_size, Block* above_block, Block* left_block);
    void encode_has_palette_uv(std::vector<int> &out_bits, int has_palette_uv, int palette_size_y, MiSize mi_size, Block* above_block, Block* left_block);
    void encode_palette_size_y_minus_2(std::vector<int> &out_bits, int palette_size_y_minus_2, MiSize mi_size);
    void encode_palette_size_uv_minus_2(std::vector<int> &out_bits, int palette_size_uv_minus_2, MiSize mi_size);
    void encode_palette_color_idx_y(std::vector<int> &out_bits, int palette_color_idx_y, int palette_size_y, int color_context_hash);
    void encode_palette_color_idx_uv(std::vector<int> &out_bits, int palette_color_idx_uv, int palette_size_uv, int color_context_hash);
    void encode_delta_q_abs(std::vector<int> &out_bits, int delta_q_abs);
    void encode_delta_lf_abs(std::vector<int> &out_bits, int delta_lf_abs, int delta_lf_multi);
    void encode_intra_tx_type(std::vector<int> &out_bits, int set, TxType tx_type, TxSize tx_size, IntraMode y_mode);
    void encode_inter_tx_type(std::vector<int> &out_bits, int set, TxType tx_type, TxSize tx_size);
    void encode_comp_ref_type(std::vector<int> &out_bits, CompRefType comp_ref_type, Block* above_block, Block* left_block);
    void encode_uni_comp_ref(std::vector<int> &out_bits, int uni_comp_ref, Block* block);
    void encode_uni_comp_ref_p1(std::vector<int> &out_bits, int uni_comp_ref, Block* block);
    void encode_uni_comp_ref_p2(std::vector<int> &out_bits, int uni_comp_ref, Block* block);
    void encode_compound_type(std::vector<int> &out_bits, CompoundType compound_type, MiSize mi_size);
    void encode_interintra(std::vector<int> &out_bits, int interintra, MiSize mi_size);
    void encode_interintra_mode(std::vector<int> &out_bits, InterIntraMode interintra_mode, MiSize mi_size);
    void encode_wedge_interintra(std::vector<int> &out_bits, int wedge_interintra, MiSize mi_size);
    void encode_use_obmc(std::vector<int> &out_bits, int use_obmc, MiSize mi_size);
    void encode_cfl_alpha_signs(std::vector<int> &out_bits, int cfl_alpha_signs);
    void encode_cfl_alpha_u(std::vector<int> &out_bits, int cfl_alpha_u, int sign_u, int sign_v);
    void encode_cfl_alpha_v(std::vector<int> &out_bits, int cfl_alpha_u, int sign_u, int sign_v);
    void encode_use_wiener(std::vector<int> &out_bits, int use_wiener);
    void encode_use_sgrproj(std::vector<int> &out_bits, int use_sgrproj);
    void encode_restoration_type(std::vector<int> &out_bits, int restoration_type);
    void encode_all_zero(std::vector<int> &out_bits, int all_zero, TxSize tx_size_ctx, int plane, TxSize tx_size, int x, int y, Block* block);
    void encode_is_eob(std::vector<int> &out_bits, int is_eob, int plane, TxSize tx_size, TxSize tx_size_ctx, int ptype, int eob_pt, int start_x, int start_y, Block* block);
    void encode_eob_extra(std::vector<int> &out_bits, int eob_extra, TxSize tx_size_ctx, int ptype, int eob_pt);
    void encode_coeff_base(std::vector<int> &out_bits, int coeff_base, int plane, int start_x, int start_y, const int* scan, int c, TxSize tx_size, TxSize tx_size_ctx, int ptype, int16_t* dequant, Block* block);
    void encode_coeff_base_eob(std::vector<int> &out_bits, int coeff_base_eob, int plane, int start_x, int start_y, const int* scan, int c, TxSize tx_size, TxSize tx_size_ctx, int ptype, int16_t* dequant, Block* block);
    void encode_dc_sign(std::vector<int> &out_bits, int dc_sign, int plane, int x, int y, TxSize tx_size, int ptype, Block* block);
    void encode_coeff_br(std::vector<int> &out_bits, int coeff_br, TxSize tx_size, TxSize tx_size_ctx, int ptype, int pos, int16_t* dequant);
    */
}