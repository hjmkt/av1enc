use std::rc::Rc;
use std::cell::RefCell;
use binary_writer::*;
use bool_coder::*;
use encoder::*;
use frame::*;
use constants::*;
use constants::BlockSize::*;
use constants::FrameRestorationType::*;
use block::*;


pub struct TileEncoder<'a, 'b> where 'a: 'b {
    writer: &'b mut BinaryWriter<'a>,
    coder: &'b mut BoolCoder,
    ectx: Rc<RefCell<EncoderContext>>,
    pub mi_row_start: usize,
    pub mi_row_end: usize,
    pub mi_col_start: usize,
    pub mi_col_end: usize,
}

impl<'a, 'b> TileEncoder<'a, 'b> where 'a: 'b {
    pub fn new(ectx: &Rc<RefCell<EncoderContext>>, writer: &'b mut BinaryWriter<'a>, coder: &'b mut BoolCoder) -> TileEncoder<'a, 'b> { TileEncoder{
        writer: writer,
        coder: coder,
        ectx: Rc::clone(ectx),
        mi_row_start: 0,
        mi_row_end: 0,
        mi_col_start: 0,
        mi_col_end: 0,
    }}

    pub fn encode_tile(&mut self, frame: &mut Frame, tile_num: usize, out_bits: &mut Vec<u8>) -> () {
        let mut sb_size;
        {
            let mut ecx = self.ectx.borrow_mut();
            let tile_col = ecx.tile_num / frame.tile_cols();
            let tile_row = ecx.tile_num / frame.tile_rows();
            self.mi_row_start = ecx.mi_row_starts[tile_row];
            self.mi_row_end = ecx.mi_row_starts[tile_row+1];
            self.mi_col_start = ecx.mi_col_starts[tile_col];
            self.mi_col_end = ecx.mi_col_starts[tile_col+1];
            ecx.current_q_index = frame.quantization_params.base_q_idx as usize;

            // init symbol

            // clear above context

            for i in 0..FRAME_LF_COUNT {
                ecx.delta_lf[i] = 0;
            }
            for plane in 0..ecx.num_planes {
                for pass in 0..2 {
                    ecx.ref_sgr_xqd[plane][pass] = sgrproj_xqd_mid[pass];
                    for i in 0..WIENER_COEFFS {
                        ecx.ref_lr_wiener[plane][pass][i] = wiener_taps_mid[i];
                    }
                }
            }
            sb_size = if ecx.use_128x128_superblock { BLOCK_128X128 } else { BLOCK_64X64 };
        }
        let sb_size4 = num_4x4_blocks_wide[sb_size as usize];
        let mut r = self.mi_row_start;
        while r < self.mi_row_end {
            // clear left context

            let mut c = self.mi_col_start;
            while c < self.mi_col_end {
                {
                    let mut ecx = self.ectx.borrow_mut();
                    ecx.read_deltas = frame.delta_q_params.delta_q_present;
                }
                // clear cdef
                // clear block decoded flags (r, c, sb_size4)

                self.encode_lr(out_bits, frame, r, c, sb_size);

                // encode partition r, c, sb_size

                c += sb_size4;
            }
            r += sb_size4;
        }

        // exit symbol
    }

    pub fn encode_lr(&mut self, out_bits: &mut Vec<u8>, frame: &mut Frame, r: usize, c: usize, bsize: BlockSize) -> () {
        if frame.allow_intrabc {
            return;
        }
        let w = num_4x4_blocks_wide[bsize as usize];
        let h = num_4x4_blocks_high[bsize as usize];
        let mut ecx = self.ectx.borrow_mut();
        let sb_size = if ecx.use_128x128_superblock { BLOCK_128X128 } else { BLOCK_64X64 };
        let sb_size4 = num_4x4_blocks_wide[sb_size as usize];
        let sb = &mut frame.superblocks[r/sb_size4][c/sb_size4];
        for plane in 0..ecx.num_planes {
            if ecx.frame_restoration_type[plane] != RESTORE_NONE {
                let subx = if plane==0 { 0 } else { frame.subsampling_x as usize };
                let suby = if plane==0 { 0 } else { frame.subsampling_y as usize };
                let unit_size = ecx.loop_restoration_size[plane];

                let count_units_in_frame = |unit_size: usize, tile_size: usize| -> usize {
                    return Max!((tile_size+(unit_size>>1)) / unit_size, 1);
                };

                let unit_rows = count_units_in_frame(unit_size, Round2!(frame.frame_height, suby));
                let unit_cols = count_units_in_frame(unit_size, Round2!(ecx.upscaled_width, subx));
                let unit_row_start = (r * (MI_SIZE>>suby) + unit_size-1) / unit_size;
                let unit_row_end = Min!(unit_rows, ((r+h) * (MI_SIZE>>suby) + unit_size-1) / unit_size);

                let (numerator, denominator) = if frame.use_superres {
                    ((MI_SIZE>>subx)*ecx.superres_denom, unit_size*SUPERRES_NUM)
                } else {
                    (MI_SIZE>>subx, unit_size)
                };
                let unit_col_start = (c*numerator + denominator-1) / denominator;
                let unit_col_end = Min!(unit_cols, ((c+w)*numerator + denominator-1) / denominator);
                for unit_row in unit_row_start..unit_row_end {
                    for unit_col in unit_col_start..unit_col_end {
                        // encode lr unit
                        let restoration_type = match ecx.frame_restoration_type[plane] {
                            RESTORE_NONE => RESTORE_NONE,
                            RESTORE_WIENER => {
                                let use_wiener = sb.lr_units[plane][unit_row-unit_row_start][unit_col-unit_col_start].use_wiener;
                                self.coder.push_bit(out_bits, use_wiener as u8);
                                if use_wiener { RESTORE_WIENER }  else { RESTORE_NONE }
                            },
                            RESTORE_SGRPROG => {
                                let use_sgrproj = sb.lr_units[plane][unit_row-unit_row_start][unit_col-unit_col_start].use_sgrproj;
                                self.coder.push_bit(out_bits, use_sgrproj as u8);
                                if use_sgrproj { RESTORE_SGRPROJ }  else { RESTORE_NONE }
                            },
                            _ => sb.lr_units[plane][unit_row-unit_row_start][unit_col-unit_col_start].restoration_type,
                        };
                        frame.lr_type[plane][unit_row][unit_col] = restoration_type;

                        let inverse_recenter = |r: usize, v: usize| -> usize {
                            if v > 2*r {
                                v
                            } else if v&1 > 0 {
                                r - ((v+1)>>1)
                            } else {
                                r + (v>>1)
                            }
                        };
                        let mut encode_subexp_bool = |v: usize, num_syms: usize, k: usize, coder: &mut BoolCoder, out_bits: &mut Vec<u8>| -> () {
                            let mut i = 0;
                            let mut mk = 0;
                            loop {
                                let b2 = if i>0 { k+i-1 } else { k };
                                let a = 1 << b2;
                                if num_syms <= mk+3*a {
                                    let tv = v - mk;
                                    coder.encode_ns(out_bits, tv as u64, (num_syms-mk) as u64);
                                    break;
                                } else {
                                    let tv = v - mk;
                                    if tv < 1<<b2 {
                                        coder.encode_literal(out_bits, 0, 1); // subexp_more_bools
                                        coder.encode_literal(out_bits, tv as u64, b2 as u8);
                                        break;
                                    } else {
                                        coder.encode_literal(out_bits, 1, 1); // subexp_more_bools
                                        i += 1;
                                        mk += a;
                                    }
                                }
                            }
                        };
                        let mut encode_unsigned_subexp_with_ref_bool = |v: usize, mx: usize, k: usize, r: usize, coder: &mut BoolCoder, out_bits: &mut Vec<u8>| -> () {
                            if (r<<1) <= mx {
                                let tv = inverse_recenter(r, v);
                                encode_subexp_bool(tv, mx, k, coder, out_bits);
                            } else {
                                let mut tv = mx - 1 - v;
                                tv = inverse_recenter(mx-1-r, tv);
                                encode_subexp_bool(tv, mx-1-r, k, coder, out_bits);
                            }
                        };
                        let mut encode_signed_subexp_with_ref_bool = |v: isize, low: isize, high: isize, k: usize, r: isize, coder: &mut BoolCoder, out_bits: &mut Vec<u8>| -> () {
                            let x = v - low;
                            encode_unsigned_subexp_with_ref_bool(x as usize, (high-low) as usize, k, (r-low) as usize, coder, out_bits);
                        };
                        if restoration_type==RESTORE_WIENER {
                            for pass in 0..2 {
                                let mut first_coeff;
                                if plane>0 {
                                    first_coeff = 1;
                                    frame.lr_wiener[plane][unit_row][unit_col][pass][0] = 0;
                                } else {
                                    first_coeff = 0;
                                }
                                for j in first_coeff..3 {
                                    let min = wiener_taps_min[j];
                                    let max = wiener_taps_max[j];
                                    let k = wiener_taps_k[j];
                                    let v = sb.lr_units[plane][unit_row-unit_row_start][unit_col-unit_col_start].lr_wiener[pass][j];
                                    encode_signed_subexp_with_ref_bool(v, min, max+1, k as usize, ecx.ref_lr_wiener[plane][pass][j], self.coder, out_bits);
                                    frame.lr_wiener[plane][unit_row][unit_col][pass][j] = v;
                                    let mut ecx = self.ectx.borrow_mut();
                                    ecx.ref_lr_wiener[plane][pass][j] = v;
                                }
                            }
                        } else if restoration_type==RESTORE_SGRPROJ {
                            let lr_sgr_set = sb.lr_units[plane][unit_row-unit_row_start][unit_col-unit_col_start].lr_sgr_set;
                            self.coder.encode_literal(out_bits, lr_sgr_set as u64, SGRPROJ_PARAMS_BITS as u8);
                            frame.lr_sgr_set[plane][unit_row][unit_col] = lr_sgr_set;
                            for i in 0..2 {
                                let radius = sgr_params[lr_sgr_set][i*2];
                                let min = sgrproj_xqd_min[i];
                                let max = sgrproj_xqd_max[i];
                                let mut v;
                                if radius>0 {
                                    v = sb.lr_units[plane][unit_row-unit_row_start][unit_col-unit_col_start].lr_sgr_xqd[i];
                                    encode_signed_subexp_with_ref_bool(v, min, max+1, SGRPROJ_PRJ_SUBEXP_K as usize, ecx.ref_sgr_xqd[plane][i], self.coder, out_bits);
                                } else {
                                    v = 0;
                                    if i==1 {
                                        v = Clip3!(min, max, (1<<SGRPROJ_PARAMS_BITS) - ecx.ref_sgr_xqd[plane][0]);
                                    }
                                }
                                frame.lr_sgr_xqd[plane][unit_row][unit_col][i] = v;
                                ecx.ref_sgr_xqd[plane][i] = v;
                            }
                        }
                    }
                }
            }
        }
    }
}