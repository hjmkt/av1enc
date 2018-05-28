extern crate num;

use constants::*;
use constants::TxSize::*;
use constants::RefFrame::*;
use self::num::FromPrimitive;
use frame::*;
use constants::TxType::*;

pub fn get_plane_residual_size(subsize: BlockSize, plane: usize, subsampling_x: bool, subsampling_y: bool) -> BlockSize {
    let subx = if plane>0 { subsampling_x as usize } else { 0 };
    let suby = if plane>0 { subsampling_y as usize } else { 0 };
    subsampled_size[subsize as usize][subx][suby]
}

pub fn find_tx_size(w: usize, h: usize) -> TxSize {
    for tx_size in 0..TX_SIZES_ALL {
        if tx_width[tx_size]==w && tx_height[tx_size]==h {
            return TxSize::from_usize(tx_size).unwrap();
        }
    }
    TX_INVALID
}

pub fn check_backward(ref_frame: RefFrame) -> bool {
    ref_frame>=BWDREF_FRAME && ref_frame<=ALTREF_FRAME
}

pub fn count_refs(ref_frame: RefFrame, frame: &mut Frame, mi_row: usize, mi_col: usize) -> usize {
    return match frame.get_above_mi(mi_row, mi_col) {
        Some(mi) => (mi.ref_frame[0]==ref_frame) as usize + (mi.ref_frame[1]==ref_frame) as usize,
        None => 0,
    } + match frame.get_left_mi(mi_row, mi_col) {
        Some(mi) => (mi.ref_frame[0]==ref_frame) as usize + (mi.ref_frame[1]==ref_frame) as usize,
        None => 0,
    };
}

pub fn ref_count_ctx(counts0: usize, counts1: usize) -> usize {
    if counts0<counts1 { 0 }
    else if counts0==counts1 { 1 }
    else { 2 }
}

pub fn get_tx_set(tx_size: TxSize, is_inter: bool, reduced_tx_set: bool) -> usize {
    let tx_sz_sqr = tx_size_sqr[tx_size as usize];
    let tx_sz_sqr_up = tx_size_sqr_up[tx_size as usize];
    if tx_sz_sqr_up > TX_32X32 {
        return TX_SET_DCTONLY;
    }
    if is_inter {
        if reduced_tx_set || tx_sz_sqr_up==TX_32X32 {
            return TX_SET_INTER_3;
        } else if tx_sz_sqr==TX_16X16 {
            return TX_SET_INTER_2;
        } else { return TX_SET_INTER_1; }
    } else {
        if tx_sz_sqr_up==TX_32X32 {
            return TX_SET_DCTONLY;
        } else if reduced_tx_set {
            return TX_SET_INTRA_2;
        } else if tx_sz_sqr==TX_16X16 {
            return TX_SET_INTRA_2;
        } else { return TX_SET_INTRA_1; }
    }
}

pub fn is_tx_type_in_set(tx_set: usize, tx_type: TxType, is_inter: bool) -> bool {
    return if is_inter { tx_type_in_set_inter[tx_set as usize][tx_type as usize] } else {
        tx_type_in_set_intra[tx_set as usize][tx_type as usize]
    };
}

pub fn compute_tx_type(plane: usize, tx_size: TxSize, block_x: usize, block_y: usize, frame: &mut Frame, lossless: bool, is_inter: bool, reduced_tx_set: bool, uv_mode: YMode, mi_row: usize, mi_col: usize) -> TxType {
    let tx_sz_sqr_up = tx_size_sqr_up[tx_size as usize];
    if lossless || tx_sz_sqr_up>TX_32X32 {
        return DCT_DCT;
    }
    let tx_set = get_tx_set(tx_size, is_inter, reduced_tx_set);
    if plane==0 {
        let tx_type = frame.get_mi(mi_row, mi_col).unwrap().tx_type;
        if !is_tx_type_in_set(tx_set, tx_type, is_inter) {
            return DCT_DCT;
        } else {
            return tx_type;
        }
    }

    if is_inter {
        let x4 = Max!(mi_col, block_x << frame.subsampling_x as usize);
        let y4 = Max!(mi_row, block_y << frame.subsampling_y as usize);
        let tx_type = frame.get_mi(y4, x4).unwrap().tx_type;
        if !is_tx_type_in_set(tx_set, tx_type, is_inter) {
            return DCT_DCT;
        } else {
            return tx_type;
        }
    }

    let tx_type = mode_to_txfm[uv_mode as usize];
    if !is_tx_type_in_set(tx_set, tx_type, is_inter) {
        return DCT_DCT;
    } else {
        return tx_type;
    }
}

pub fn get_tx_class(tx_type: TxType) -> usize {
    return  match tx_type {
        V_DCT | V_ADST | V_FLIPADST => TX_CLASS_VERT,
        H_DCT | H_ADST | H_FLIPADST => TX_CLASS_HORIZ,
        _ => TX_CLASS_2D,
    };
}

pub fn get_coeff_base_ctx(tx_size: TxSize, plane: usize, block_x: usize, block_y: usize, pos: usize, c: usize, is_eob: bool, frame: &mut Frame, lossless: bool, is_inter: bool, reduced_tx_set: bool, uv_mode: YMode, mi_row: usize, mi_col: usize, quant: &[isize]) -> usize {
    let adj_tx_size = adjusted_tx_size[tx_size as usize];
    let bwl = tx_width_log2[adj_tx_size as usize];
    let width = 1 << bwl;
    let height = tx_height[adj_tx_size as usize];
    let tx_type = compute_tx_type(plane, tx_size, block_x, block_y, frame, lossless, is_inter, reduced_tx_set, uv_mode, mi_row, mi_col);
    if is_eob {
        if c == 0 { return SIG_COEF_CONTEXTS - 4; }
        if c <= (height<<bwl) / 8 { return SIG_COEF_CONTEXTS - 3; }
        if c <= (height<<bwl) / 4 { return SIG_COEF_CONTEXTS - 2; }
        else { return SIG_COEF_CONTEXTS - 1; }
    }
    let tx_class = get_tx_class(tx_type);
    let row = pos >> bwl;
    let col = pos - (row << bwl);
    let mut mag = 0;
    for idx in 0..SIG_REF_DIFF_OFFSET_NUM {
        let ref_row = row + sig_ref_diff_offset[tx_class][idx][0];
        let ref_col = col + sig_ref_diff_offset[tx_class][idx][1];
        if /*ref_row>=0 && ref_col>=0 && FIXME*/ref_row<height && ref_col<width {
            mag += Min!(Abs!(quant[(ref_row<<bwl) + ref_col]) as usize, 3);
        }
    }
    let ctx: usize = Min!((mag+1) >> 1, 4);
    if tx_class == TX_CLASS_2D {
        if row==0 && col==0 {
            return 0;
        } else {
            return ctx + coeff_base_ctx_offset[tx_size as usize][Min!(row, 4)][Min!(col, 4)];
        }
    }
    let idx = if tx_class==TX_CLASS_VERT { row } else { col };
    return ctx + coeff_base_pos_ctx_offset[Min!(idx, 2)];
}

pub fn is_samedir_ref_pair(ref0: RefFrame, ref1: RefFrame) -> bool {
    return if ref0<=INTRA_FRAME || ref1<=INTRA_FRAME {
        false
    } else {
        (ref0>=BWDREF_FRAME) as usize == (ref1>=BWDREF_FRAME) as usize
    };
}

pub fn get_relative_dist(a: usize, b: usize, enable_order_hint: bool, order_hint_bits: usize) -> isize {
    if !enable_order_hint {
        return 0;
    }
    let mut diff: isize = a as isize - b as isize;
    let m = 1 << (order_hint_bits - 1);
    diff = (diff & (m-1)) - (diff & m);
    return diff;
}