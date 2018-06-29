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
use constants::Partition::*;
use cdf::*;
use util::*;
use constants::RefFrame::*;
extern crate num;
extern crate enum_primitive;
use self::num::FromPrimitive;
use constants::CompMode::*;
use constants::CompRefType::*;
use constants::YMode::*;
use common::*;
use constants::MvJoint::*;
use constants::MVClass::*;


//pub struct TileEncoder<'a, 'b> where 'a: 'b {
pub struct TileEncoder<'b> {
//    writer: &'b mut BinaryWriter<'a>,
    coder: &'b mut BoolCoder,
    ectx: Rc<RefCell<EncoderContext>>,
    pub mi_row_start: usize,
    pub mi_row_end: usize,
    pub mi_col_start: usize,
    pub mi_col_end: usize,
    pub avail_u: bool,
    pub avail_l: bool,
    pub cdf: CDF,
    pub mi_row: usize,
    pub mi_col: usize,
    pub mi_size: BlockSize,
    pub has_chroma: bool,
    pub avail_u_chroma: bool,
    pub avail_l_chroma: bool,
    pub y_mode: usize,
    pub lossless: bool,
    pub last_active_seg_id: usize,
    pub skip: bool,
    pub palette_size_y: usize,
    pub palette_size_uv: usize,
    pub ref_mv_context: usize,
    pub new_mv_context: usize,
    pub zero_mv_context: usize,
    pub ref_mv_idx: usize,
    pub num_mv_found: usize,
    pub drl_ctx_stack: [usize; 8], // FIXME
    pub pred_mv: [MotionVector; 2],
    pub ref_stack_mv: [[MotionVector; 2]; 8], // FIXME
    pub global_mvs: [MotionVector; 2],
    pub mv_ctx: usize,
    pub new_mv_count: usize,
    pub found_match: bool,
    pub weight_stack: [usize; 8], // FIXME
    pub close_matches: usize,
}

//impl<'a, 'b> TileEncoder<'a, 'b> where 'a: 'b {
impl<'b> TileEncoder<'b> {
//    pub fn new(ectx: &Rc<RefCell<EncoderContext>>, writer: &'b mut BinaryWriter<'a>, coder: &'b mut BoolCoder) -> TileEncoder<'a, 'b> { TileEncoder{
    pub fn new(ectx: &Rc<RefCell<EncoderContext>>, coder: &'b mut BoolCoder) -> TileEncoder<'b> { TileEncoder{
//        writer: writer,
        coder: coder,
        ectx: Rc::clone(ectx),
        mi_row_start: 0,
        mi_row_end: 0,
        mi_col_start: 0,
        mi_col_end: 0,
        avail_u: false,
        avail_l: false,
        cdf: default_cdf.clone(),
        mi_row: 0,
        mi_col: 0,
        mi_size: BLOCK_128X128,
        has_chroma: true,
        avail_u_chroma: false,
        avail_l_chroma: false,
        y_mode: 0,
        lossless: false,
        last_active_seg_id: 0,
        skip: false,
        palette_size_y: 0,
        palette_size_uv: 0,
        ref_mv_context: 0,
        new_mv_context: 0,
        zero_mv_context: 0,
        ref_mv_idx: 0,
        num_mv_found: 0,
        drl_ctx_stack: [0; 8],
        pred_mv: [MotionVector::new(0, 0); 2],
        ref_stack_mv: [[MotionVector::new(0, 0); 2]; 8],
        global_mvs: [MotionVector::new(0, 0); 2],
        mv_ctx: 0,
        new_mv_count: 0,
        found_match: false,
        weight_stack: [0; 8],
        close_matches: 0,
    }}

    pub fn is_inside(&self, r: isize, c: isize) -> bool {
        return c>=self.mi_col_start as isize && c<self.mi_col_end as isize && r>=self.mi_row_start as isize && r<self.mi_row_end as isize;
    }

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

                    // clear cdef r c
                    frame.cdef_idx[r][c] = -1;
                    if ecx.use_128x128_superblock {
                        let cdef_size4 = num_4x4_blocks_wide[BLOCK_64X64 as usize];
                        frame.cdef_idx[r][c+cdef_size4] = -1;
                        frame.cdef_idx[r+cdef_size4][c] = -1;
                        frame.cdef_idx[r+cdef_size4][c+cdef_size4] = -1;
                    }
                }

                // clear block decoded flags (r, c, sb_size4)

                self.encode_lr(out_bits, frame, r, c, sb_size);

                // encode partition r, c, sb_size
                let mut superblock = Superblock::new(sb_size, c*4, r*4);
                self.encode_superblock(&mut superblock, frame, out_bits);

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
                            RESTORE_SGRPROJ => {
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

    pub fn encode_superblock(&mut self, sb: &mut Superblock, frame: &mut Frame, out_bits: &mut Vec<u8>) {
        let r = sb.y/4;
        let c = sb.x/4;
        let mut part_tree = PartitionTree::new(sb.mi_size);
        self.encode_partition(sb, r, c, frame, out_bits, &mut part_tree);
        sb.part_tree = part_tree;
    }

    pub fn encode_partition(&mut self, sb: &mut Superblock, r: usize, c: usize, frame: &mut Frame, out_bits: &mut Vec<u8>, part_tree: &mut PartitionTree) {
        let mi_rows = frame.mi_rows();
        let mi_cols = frame.mi_cols();
        if r>=mi_rows || c>=mi_cols {
            return;
        }

        self.avail_u = self.is_inside(r as isize - 1, c as isize);
        self.avail_l = self.is_inside(r as isize, c as isize - 1);
        let num4x4 = num_4x4_blocks_wide[part_tree.bsize as usize];
        let half_block4x4 = num4x4 >> 1;
        let quarter_block4x4 = half_block4x4 >> 1;
        let has_rows = r+half_block4x4 < mi_rows;
        let has_cols = c+half_block4x4 < mi_cols;

        let mut best_part_bits: Vec<u8> = vec![];
        let mut best_part_tree = part_tree.clone();
        let mut best_cost = <u64>::max_value();
        let mut best_coder = self.coder.clone();
        let mut best_cdf = self.cdf.clone();
        for partition in [
            PARTITION_NONE, PARTITION_HORZ, PARTITION_VERT, PARTITION_SPLIT,
            PARTITION_HORZ_A, PARTITION_HORZ_B, PARTITION_VERT_A, PARTITION_VERT_B, PARTITION_HORZ_4, PARTITION_VERT_4
        ].iter() {
            let mut  part_bits: Vec<u8> = vec![];
            if part_tree.bsize < BLOCK_8X8 {
                if *partition != PARTITION_NONE {
                    continue;
                }
            } else if has_rows && has_cols {
                self.coder.encode_partition(&mut part_bits, &mut self.cdf, part_tree.bsize, self.avail_u, self.avail_l, *partition);
            } else if has_cols {
                match *partition {
                    PARTITION_SPLIT => {
                        self.coder.encode_split_or_horz(&mut part_bits, &mut self.cdf, part_tree.bsize, self.avail_u, self.avail_l, 1);
                    },
                    PARTITION_HORZ => {
                        self.coder.encode_split_or_horz(&mut part_bits, &mut self.cdf, part_tree.bsize, self.avail_u, self.avail_l, 0);
                    },
                    _ => { continue; }
                }
            } else if has_rows {
                match *partition {
                    PARTITION_SPLIT => {
                        self.coder.encode_split_or_vert(&mut part_bits, &mut self.cdf, part_tree.bsize, self.avail_u, self.avail_l, 1);
                    },
                    PARTITION_VERT => {
                        self.coder.encode_split_or_vert(&mut part_bits, &mut self.cdf, part_tree.bsize, self.avail_u, self.avail_l, 0);
                    },
                    _ => { continue; }
                }
            } else {
                if *partition != PARTITION_SPLIT {
                    continue;
                }
            }
            part_tree.part_type = *partition;
            let sub_size = partition_subsize[*partition as usize][part_tree.bsize as usize];
            let split_size = partition_subsize[PARTITION_SPLIT as usize][part_tree.bsize as usize];

            match *partition {
                PARTITION_NONE => {
                    self.encode_block(sb, r, c, frame, &mut part_bits, part_tree);
                },
                PARTITION_HORZ => {
                    part_tree.partitions = vec![ PartitionTree::new(sub_size) ];
                    self.encode_block(sb, r, c, frame, &mut part_bits, &mut part_tree.partitions[0]);
                    if has_rows {
                        part_tree.partitions.push(PartitionTree::new(sub_size));
                        self.encode_block(sb, r+half_block4x4, c, frame, &mut part_bits, &mut part_tree.partitions[1]);
                    }
                },
                PARTITION_VERT => {
                    part_tree.partitions = vec![ PartitionTree::new(sub_size) ];
                    self.encode_block(sb, r, c, frame, &mut part_bits, &mut part_tree.partitions[0]);
                    if has_cols {
                        part_tree.partitions.push(PartitionTree::new(sub_size));
                        self.encode_block(sb, r, c+half_block4x4, frame, &mut part_bits, &mut part_tree.partitions[1]);
                    }
                },
                PARTITION_SPLIT => {
                    part_tree.partitions = vec![
                        PartitionTree::new(sub_size),
                        PartitionTree::new(sub_size),
                        PartitionTree::new(sub_size),
                        PartitionTree::new(sub_size),
                    ];
                    self.encode_partition(sb, r, c, frame, &mut part_bits, &mut part_tree.partitions[0]);
                    self.encode_partition(sb, r, c+half_block4x4, frame, &mut part_bits, &mut part_tree.partitions[1]);
                    self.encode_partition(sb, r+half_block4x4, c, frame, &mut part_bits, &mut part_tree.partitions[2]);
                    self.encode_partition(sb, r+half_block4x4, c+half_block4x4, frame, &mut part_bits, &mut part_tree.partitions[3]);
                },
                PARTITION_HORZ_A => {
                    part_tree.partitions = vec![
                        PartitionTree::new(split_size),
                        PartitionTree::new(split_size),
                        PartitionTree::new(sub_size),
                    ];
                    self.encode_block(sb, r, c, frame, &mut part_bits, &mut part_tree.partitions[0]);
                    self.encode_block(sb, r, c+half_block4x4, frame, &mut part_bits, &mut part_tree.partitions[1]);
                    self.encode_block(sb, r+half_block4x4, c, frame, &mut part_bits, &mut part_tree.partitions[2]);
                },
                PARTITION_HORZ_B => {
                    part_tree.partitions = vec![
                        PartitionTree::new(sub_size),
                        PartitionTree::new(split_size),
                        PartitionTree::new(split_size),
                    ];
                    self.encode_block(sb, r, c, frame, &mut part_bits, &mut part_tree.partitions[0]);
                    self.encode_block(sb, r+half_block4x4, c, frame, &mut part_bits, &mut part_tree.partitions[1]);
                    self.encode_block(sb, r, c+half_block4x4, frame, &mut part_bits, &mut part_tree.partitions[2]);
                },
                PARTITION_VERT_A => {
                    part_tree.partitions = vec![
                        PartitionTree::new(split_size),
                        PartitionTree::new(split_size),
                        PartitionTree::new(sub_size),
                    ];
                    self.encode_block(sb, r, c, frame, &mut part_bits, &mut part_tree.partitions[0]);
                    self.encode_block(sb, r+half_block4x4, c, frame, &mut part_bits, &mut part_tree.partitions[1]);
                    self.encode_block(sb, r, c+half_block4x4, frame, &mut part_bits, &mut part_tree.partitions[2]);
                },
                PARTITION_VERT_B => {
                    part_tree.partitions = vec![
                        PartitionTree::new(sub_size),
                        PartitionTree::new(split_size),
                        PartitionTree::new(split_size),
                    ];
                    self.encode_block(sb, r, c, frame, &mut part_bits, &mut part_tree.partitions[0]);
                    self.encode_block(sb, r, c+half_block4x4, frame, &mut part_bits, &mut part_tree.partitions[1]);
                    self.encode_block(sb, r+half_block4x4, c, frame, &mut part_bits, &mut part_tree.partitions[2]);
                },
                PARTITION_HORZ_4 => {
                    part_tree.partitions = vec![
                        PartitionTree::new(sub_size),
                        PartitionTree::new(sub_size),
                        PartitionTree::new(sub_size),
                    ];
                    self.encode_block(sb, r+quarter_block4x4*0, c, frame, &mut part_bits, &mut part_tree.partitions[0]);
                    self.encode_block(sb, r+quarter_block4x4*1, c, frame, &mut part_bits, &mut part_tree.partitions[1]);
                    self.encode_block(sb, r+quarter_block4x4*2, c, frame, &mut part_bits, &mut part_tree.partitions[2]);
                    if r+quarter_block4x4*3 < mi_rows {
                        part_tree.partitions.push(PartitionTree::new(sub_size));
                        self.encode_block(sb, r+quarter_block4x4*3, c, frame, &mut part_bits, &mut part_tree.partitions[3]);
                    }
                },
                PARTITION_VERT_4 => {
                    part_tree.partitions = vec![
                        PartitionTree::new(sub_size),
                        PartitionTree::new(sub_size),
                        PartitionTree::new(sub_size),
                    ];
                    self.encode_block(sb, r, c+quarter_block4x4*0, frame, &mut part_bits, &mut part_tree.partitions[0]);
                    self.encode_block(sb, r, c+quarter_block4x4*1, frame, &mut part_bits, &mut part_tree.partitions[1]);
                    self.encode_block(sb, r, c+quarter_block4x4*2, frame, &mut part_bits, &mut part_tree.partitions[2]);
                    if c+quarter_block4x4*3 < mi_cols {
                        part_tree.partitions.push(PartitionTree::new(sub_size));
                        self.encode_block(sb, r+quarter_block4x4*3, c, frame, &mut part_bits, &mut part_tree.partitions[3]);
                    }
                },
            }

            let cost = 0; // TODO
            if cost < best_cost {
                best_cost = cost;
                best_part_tree = part_tree.clone();
                best_part_bits = part_bits;
                best_coder = self.coder.clone();
                best_cdf = self.cdf.clone();
            }
        }
        *self.coder = best_coder;
        self.cdf = best_cdf;
        *part_tree = best_part_tree;
        self.coder.push_bits(out_bits, &best_part_bits);
    }

    pub fn encode_block(&mut self, sb: &mut Superblock, r: usize, c: usize, frame: &mut Frame, out_bits: &mut Vec<u8>, part_tree: &mut PartitionTree) -> () {
        self.mi_row = r;
        self.mi_col = c;
        self.mi_size = part_tree.bsize;
        let bw4 = num_4x4_blocks_wide[part_tree.bsize as usize];
        let bh4 = num_4x4_blocks_high[part_tree.bsize as usize];
        if bh4==1 && frame.subsampling_y && (self.mi_row&1)==0 {
            self.has_chroma = false;
        } else if bw4==1 && frame.subsampling_x && (self.mi_col&1)==0 {
            self.has_chroma = false;
        } else {
            let mut ecx = self.ectx.borrow_mut();
            self.has_chroma = ecx.num_planes > 1;
        }
        self.avail_u = self.is_inside(r as isize - 1, c as isize);
        self.avail_l = self.is_inside(r as isize, c as isize - 1);
        self.avail_u_chroma = self.avail_u;
        self.avail_l_chroma = self.avail_l;
        if self.has_chroma {
            if frame.subsampling_y && bh4==1 {
                self.avail_u_chroma = self.is_inside(r as isize - 2, c as isize);
            }
            if frame.subsampling_x && bw4==1 {
                self.avail_l_chroma = self.is_inside(r as isize, c as isize - 2);
            }
        } else {
            self.avail_u_chroma = false;
            self.avail_l_chroma = false;
        }

        // mode info
        self.encode_mode_info(sb, frame, out_bits, part_tree);
    }

    pub fn encode_segment_id(&mut self, sb: &mut Superblock, frame: &mut Frame, out_bits: &mut Vec<u8>, part_tree: &mut PartitionTree) -> () {
        let prev_ul = if self.avail_u && self.avail_l {
            frame.mis[self.mi_row-1][self.mi_col-1].unwrap().segment_id as isize
        } else { -1 };
        let prev_u = if self.avail_u {
            frame.mis[self.mi_row-1][self.mi_col].unwrap().segment_id as isize
        } else { -1 };
        let prev_l = if self.avail_l {
            frame.mis[self.mi_row][self.mi_col-1].unwrap().segment_id as isize
        } else { -1 };
        let pred = if prev_u==-1 {
            if prev_l==-1 { 0 } else { prev_l }
        } else if prev_l==-1 {
            prev_u
        } else {
            if prev_ul==prev_u { prev_u } else { prev_l }
        };
        if self.skip {
            part_tree.segment_id = pred as usize;
        } else {
            let neg_interleave = |v: usize, rf: usize, max: usize| -> usize {
                if rf==0{
                    v
                } else if rf>=(max-1) {
                    max - 1 - v
                } else if 2*rf < max {
                    let d1 = ((v-rf)<<1) - 1;
                    let d2 = (rf-v)<<1;
                    if d1 <= 2*rf {
                        d1
                    } else if d2 <= 2*rf {
                        d2
                    } else {
                        v
                    }
                } else {
                    let d1 = ((v-rf)<<1) - 1;
                    let d2 = (rf-v)<<1;
                    if d1 <= 2*(max-rf-1) {
                        d1
                    } else if d2 <= 2*(max-rf-1) {
                        d2
                    } else {
                        rf-v-1
                    }
                }
            };
            let v = neg_interleave(part_tree.segment_id, pred as usize, self.last_active_seg_id+1);
            self.coder.encode_segment_id(out_bits, &mut self.cdf, v, prev_u as isize, prev_l as isize, prev_ul as isize);
        }
    }

    pub fn encode_mode_info(&mut self, sb: &mut Superblock, frame: &mut Frame, out_bits: &mut Vec<u8>, part_tree: &mut PartitionTree) -> () {
        let frame_is_intra = self.ectx.borrow_mut().frame_is_intra;
        if frame_is_intra {
            // intra frame mode info
            self.skip = false;
            let seg_id_pre_skip = self.ectx.borrow_mut().seg_id_pre_skip;
            if seg_id_pre_skip {
                // intra segment id
                if frame.segmentation_params.segmentation_enabled {
                    self.encode_segment_id(sb, frame, out_bits, part_tree);
                } else {
                    part_tree.segment_id = 0;
                }
                self.lossless = self.ectx.borrow_mut().lossless_array[part_tree.segment_id];
            }

            let skip_mode = 0;
            // encode skip
            let seg_feature_active_idx = |idx: usize, feature: usize, frame: &mut Frame| -> bool {
                return frame.segmentation_params.segmentation_enabled && frame.segmentation_params.feature_value[idx][feature].is_some()
            };
            if seg_id_pre_skip && seg_feature_active_idx(part_tree.segment_id, SEG_LVL_SKIP, frame) {
                part_tree.skip = true;
            } else {
                self.coder.encode_skip(out_bits, &mut self.cdf, part_tree.skip as usize, frame, self.mi_row, self.mi_col);
            }

            if !seg_id_pre_skip {
                // intra segment id
                if frame.segmentation_params.segmentation_enabled {
                    self.encode_segment_id(sb, frame, out_bits, part_tree);
                } else {
                    part_tree.segment_id = 0;
                }
                self.lossless = self.ectx.borrow_mut().lossless_array[part_tree.segment_id];
            }

            // encode cdef
            let coded_lossless = self.ectx.borrow_mut().coded_lossless;
            let enable_cdef = self.ectx.borrow_mut().enable_cdef;
            if !(self.skip || coded_lossless || !enable_cdef || frame.allow_intrabc) {
                let cdef_size4 = num_4x4_blocks_wide[BLOCK_64X64 as usize];
                let cdef_mask4 = !(cdef_size4 - 1);
                let r = self.mi_row & cdef_mask4;
                let c = self.mi_col & cdef_mask4;
                if frame.cdef_idx[r][c] == -1 {
                    self.coder.encode_literal(out_bits, part_tree.cdef_idx as u64, frame.cdef_params.cdef_bits as u8);
                    let w4 = num_4x4_blocks_wide[part_tree.bsize as usize];
                    let h4 = num_4x4_blocks_high[part_tree.bsize as usize];
                    {
                        let mut i = 0;
                        while i < r+h4 {
                            let mut j = 0;
                            while j < c+w4 {
                                frame.cdef_idx[i][j] = part_tree.cdef_idx;
                                j += cdef_size4;
                            }
                            i += cdef_size4;
                        }
                    }
                }
            }

            // encode delta_q index
            let use_128x128_superblock = self.ectx.borrow_mut().use_128x128_superblock;
            let sb_size = if use_128x128_superblock { BLOCK_128X128 } else { BLOCK_64X64 };
            let read_deltas = self.ectx.borrow_mut().read_deltas;
            if !(part_tree.bsize==sb_size && part_tree.skip) && read_deltas {
                let mut current_q_index = &mut self.ectx.borrow_mut().current_q_index;
                let delta_q_res = frame.delta_q_params.delta_q_res;
                let delta_q_index = (part_tree.q_index as isize - *current_q_index as isize) >> delta_q_res;
                let delta_q_abs = Abs!(delta_q_index) as usize;
                if delta_q_index != 0 {
                    if delta_q_abs < DELTA_Q_SMALL {
                        self.coder.encode_delta_q_abs(out_bits, &mut self.cdf, delta_q_abs);
                    } else {
                        self.coder.encode_delta_q_abs(out_bits, &mut self.cdf, DELTA_Q_SMALL);
                        let delta_q_rem_bits = msb16((delta_q_abs-1) as u16);
                        self.coder.encode_literal(out_bits, (delta_q_rem_bits-1) as u64, 3);
                        let delta_q_abs_bits = (delta_q_abs-1) - (1<<(delta_q_rem_bits-1));
                        self.coder.encode_literal(out_bits, delta_q_abs_bits as u64, delta_q_rem_bits);
                    }
                    let delta_q_sign_bit = if delta_q_index < 0 { 1 } else { 0 };
                    self.coder.encode_literal(out_bits, delta_q_sign_bit as u64, 1);
                } else {
                    self.coder.encode_delta_q_abs(out_bits, &mut self.cdf, delta_q_abs);
                }
                *current_q_index = part_tree.q_index;
            }

            // encode delta lf
            if !(part_tree.bsize==sb_size && part_tree.skip) && read_deltas && frame.delta_lf_params.delta_lf_present {
                let frame_lf_count = if frame.delta_lf_params.delta_lf_multi {
                    let num_planes = self.ectx.borrow_mut().num_planes;
                    if num_planes > 1 { FRAME_LF_COUNT } else { FRAME_LF_COUNT-2 }
                } else { 1 };
                for i in 0..frame_lf_count {
                    let current_delta_lf = &mut self.ectx.borrow_mut().delta_lf[i];
                    let delta_lf_res = frame.delta_lf_params.delta_lf_res;
                    let delta_lf = (part_tree.delta_lf[i]-*current_delta_lf) >> delta_lf_res;
                    if delta_lf != 0 {
                        let delta_lf_abs = Abs!(delta_lf) as usize;
                        if delta_lf_abs < DELTA_LF_SMALL {
                            self.coder.encode_delta_lf_abs(out_bits, &mut self.cdf, delta_lf_abs, frame.delta_lf_params.delta_lf_multi as usize, i);
                        } else {
                            self.coder.encode_delta_lf_abs(out_bits, &mut self.cdf, DELTA_LF_SMALL, frame.delta_lf_params.delta_lf_multi as usize, i);
                            let delta_lf_rem_bits = msb16((delta_lf_abs-1) as u16);
                            self.coder.encode_literal(out_bits, (delta_lf_rem_bits-1) as u64, 3);
                            let delta_lf_abs_bits = (delta_lf_abs-1) - (1<<(delta_lf_rem_bits-1));
                            self.coder.encode_literal(out_bits, (delta_lf_abs_bits-1) as u64, delta_lf_rem_bits);
                        }
                        let delta_lf_sign_bits = if delta_lf < 0 { 1 } else { 0 };
                        self.coder.encode_literal(out_bits, delta_lf_sign_bits as u64, 1);
                    } else {
                        self.coder.encode_delta_lf_abs(out_bits, &mut self.cdf, 0, frame.delta_lf_params.delta_lf_multi as usize, i);
                    }
                    *current_delta_lf = part_tree.delta_lf[i];
                }
            }

            self.ectx.borrow_mut().read_deltas = false;

            // encode is_inter
            if skip_mode > 0 {
                part_tree.is_inter = true;
            } else {
                let seg_feature_active_idx = |idx: usize, feature: usize, frame: &mut Frame| -> bool {
                    return frame.segmentation_params.segmentation_enabled && frame.segmentation_params.feature_value[idx][feature].is_some()
                };
                if seg_feature_active_idx(part_tree.segment_id, SEG_LVL_REF_FRAME, frame) {
                    part_tree.is_inter = frame.segmentation_params.feature_value[part_tree.segment_id][SEG_LVL_REF_FRAME].unwrap() != INTRA_FRAME as isize;
                } else if seg_feature_active_idx(part_tree.segment_id, SEG_LVL_GLOBALMV, frame) {
                    part_tree.is_inter = true;
                } else {
                    self.coder.encode_is_inter(out_bits, &mut self.cdf, part_tree.is_inter, self.mi_row, self.mi_col, frame);
                }
            }

            if part_tree.is_inter {
                // inter block mode info
                self.palette_size_y = 0;
                self.palette_size_uv = 0;

                // encode ref frames
                if skip_mode > 0 {
                    let ecx = self.ectx.borrow_mut();
                    part_tree.ref_frame[0] = ecx.skip_mode_frame[0];
                    part_tree.ref_frame[1] = ecx.skip_mode_frame[1];
                } else {
                    let seg_feature_active_idx = |idx: usize, feature: usize, frame: &mut Frame| -> bool {
                        return frame.segmentation_params.segmentation_enabled && frame.segmentation_params.feature_value[idx][feature].is_some()
                    };
                    if seg_feature_active_idx(part_tree.segment_id, SEG_LVL_REF_FRAME, frame) {
                        part_tree.ref_frame[0] = RefFrame::from_isize(frame.segmentation_params.feature_value[part_tree.segment_id][SEG_LVL_REF_FRAME].unwrap()).unwrap();
                        part_tree.ref_frame[1] = NONE;
                    } else if seg_feature_active_idx(part_tree.segment_id, SEG_LVL_SKIP, frame) ||
                              seg_feature_active_idx(part_tree.segment_id, SEG_LVL_GLOBALMV, frame){
                        part_tree.ref_frame[0] = LAST_FRAME;
                        part_tree.ref_frame[1] = NONE;
                    } else {
                        let bw4 = num_4x4_blocks_wide[part_tree.bsize as usize];
                        let bh4 = num_4x4_blocks_high[part_tree.bsize as usize];
                        if frame.reference_select && Min!(bw4, bh4) >= 2 {
                            self.coder.encode_comp_mode(out_bits, &mut self.cdf, part_tree.comp_mode as usize, frame, self.mi_row, self.mi_col);
                        } else {
                            part_tree.comp_mode = SINGLE_REFERENCE;
                        }
                        if part_tree.comp_mode==COMPOUND_REFERENCE {
                            self.coder.encode_comp_ref_type(out_bits, &mut self.cdf, part_tree.comp_ref_type as usize, frame, self.mi_row, self.mi_col);
                            if part_tree.comp_ref_type==UNIDIR_COMP_REFERENCE {
                                let uni_comp_ref = part_tree.ref_frame[0]==BWDREF_FRAME && part_tree.ref_frame[1]==ALTREF_FRAME;
                                self.coder.encode_uni_comp_ref(out_bits, &mut self.cdf, uni_comp_ref as usize, frame, self.mi_row, self.mi_col);
                                if !uni_comp_ref {
                                    let uni_comp_ref_p1 = part_tree.ref_frame[0]==LAST_FRAME && [GOLDEN_FRAME, LAST3_FRAME].contains(&part_tree.ref_frame[1]);
                                    self.coder.encode_uni_comp_ref_p1(out_bits, &mut self.cdf, uni_comp_ref_p1 as usize, frame, self.mi_row, self.mi_col);
                                }
                            } else {
                                let comp_ref = [GOLDEN_FRAME, LAST3_FRAME].contains(&part_tree.ref_frame[0]);
                                self.coder.encode_comp_ref(out_bits, &mut self.cdf, comp_ref as usize, frame, self.mi_row, self.mi_col);
                                if !comp_ref {
                                    let comp_ref_p1 = part_tree.ref_frame[0]==LAST2_FRAME;
                                    self.coder.encode_comp_ref_p1(out_bits, &mut self.cdf, comp_ref_p1 as usize, frame, self.mi_row, self.mi_col);
                                } else {
                                    let comp_ref_p2 = part_tree.ref_frame[0]==GOLDEN_FRAME;
                                    self.coder.encode_comp_ref_p2(out_bits, &mut self.cdf, comp_ref_p2 as usize, frame, self.mi_row, self.mi_col);
                                }

                                let comp_bwdref = part_tree.ref_frame[1]==ALTREF_FRAME;
                                self.coder.encode_comp_bwdref(out_bits, &mut self.cdf, comp_bwdref as usize, frame, self.mi_row, self.mi_col);
                                if !comp_ref {
                                    let comp_bwdref_p1 = part_tree.ref_frame[1]==ALTREF2_FRAME;
                                    self.coder.encode_comp_ref_p1(out_bits, &mut self.cdf, comp_bwdref_p1 as usize, frame, self.mi_row, self.mi_col);
                                }
                            }
                        } else {
                            let single_ref_p1 = [ALTREF2_FRAME, BWDREF_FRAME, ALTREF_FRAME].contains(&part_tree.ref_frame[0]);
                            self.coder.encode_single_ref_p1(out_bits, &mut self.cdf, single_ref_p1 as usize, frame, self.mi_row, self.mi_col);
                            if single_ref_p1 {
                                let single_ref_p2 = part_tree.ref_frame[0]==ALTREF_FRAME;
                                self.coder.encode_single_ref_p2(out_bits, &mut self.cdf, single_ref_p2 as usize, frame, self.mi_row, self.mi_col);
                                if !single_ref_p2 {
                                    let single_ref_p6 = part_tree.ref_frame[0]==ALTREF2_FRAME;
                                    self.coder.encode_single_ref_p6(out_bits, &mut self.cdf, single_ref_p6 as usize, frame, self.mi_row, self.mi_col);
                                }
                            } else {
                                let single_ref_p3 = [GOLDEN_FRAME, LAST3_FRAME].contains(&part_tree.ref_frame[0]);
                                self.coder.encode_single_ref_p3(out_bits, &mut self.cdf, single_ref_p3 as usize, frame, self.mi_row, self.mi_col);
                                if single_ref_p3 {
                                    let single_ref_p5 = part_tree.ref_frame[0]==GOLDEN_FRAME;
                                    self.coder.encode_single_ref_p5(out_bits, &mut self.cdf, single_ref_p5 as usize, frame, self.mi_row, self.mi_col);
                                } else {
                                    let single_ref_p4 = part_tree.ref_frame[0]==LAST2_FRAME;
                                    self.coder.encode_single_ref_p4(out_bits, &mut self.cdf, single_ref_p4 as usize, frame, self.mi_row, self.mi_col);
                                }
                            }
                            part_tree.ref_frame[1] = NONE;
                        }
                    }
                }

                let is_compound = part_tree.ref_frame[1] > INTRA_FRAME;

                self.find_mv_stack(is_compound, frame, part_tree);

                if skip_mode > 0 {
                    self.y_mode = NEAREST_NEARESTMV as usize;
                } else {
                    let seg_feature_active_idx = |idx: usize, feature: usize, frame: &mut Frame| -> bool {
                        return frame.segmentation_params.segmentation_enabled && frame.segmentation_params.feature_value[idx][feature].is_some()
                    };
                    if seg_feature_active_idx(part_tree.segment_id, SEG_LVL_SKIP, frame) ||
                       seg_feature_active_idx(part_tree.segment_id, SEG_LVL_GLOBALMV, frame){
                        self.y_mode = GLOBALMV as usize;
                    } else if is_compound {
                        let compound_mode = part_tree.compound_mode as usize - NEAREST_NEARESTMV as usize;
                        self.coder.encode_compound_mode(out_bits, &mut self.cdf, compound_mode, self.ref_mv_context, self.new_mv_context);
                    } else {
                        self.coder.encode_new_mv(out_bits, &mut self.cdf, part_tree.new_mv as usize, self.new_mv_context);
                        if !part_tree.new_mv {
                            self.y_mode = NEWMV as usize;
                        } else {
                            self.coder.encode_zero_mv(out_bits, &mut self.cdf, part_tree.zero_mv as usize, self.zero_mv_context);
                            if !part_tree.zero_mv {
                                self.y_mode = GLOBALMV as usize;
                            } else {
                                self.coder.encode_ref_mv(out_bits, &mut self.cdf, part_tree.ref_mv as usize, self.ref_mv_context);
                                self.y_mode = if !part_tree.ref_mv { NEARESTMV as usize } else { NEARMV as usize };
                            }
                        }
                    }
                    self.ref_mv_idx = 0;
                    let has_nearmv = |y_mode: usize| -> bool {
                        return y_mode==NEARMV as usize || y_mode==NEAR_NEARMV as usize ||
                               y_mode==NEAR_NEWMV as usize || y_mode==NEW_NEARMV as usize;
                    };
                    if self.y_mode==NEWMV as usize || self.y_mode==NEW_NEWMV as usize {
                        for idx in 0..2 {
                            if self.num_mv_found > idx+1 {
                                self.coder.encode_drl_mode(out_bits, &mut self.cdf, part_tree.drl_mode[idx], idx, &self.drl_ctx_stack);
                                if part_tree.drl_mode[idx]==0 {
                                    self.ref_mv_idx = idx;
                                    break;
                                }
                                self.ref_mv_idx = idx + 1;
                            }
                        }
                    } else if has_nearmv(self.y_mode) {
                        self.ref_mv_idx = 1;
                        for idx in 1..3 {
                            if self.num_mv_found > idx+1 {
                                self.coder.encode_drl_mode(out_bits, &mut self.cdf, part_tree.drl_mode[idx], idx, &self.drl_ctx_stack);
                                if part_tree.drl_mode[idx]==0 {
                                    self.ref_mv_idx = idx;
                                    break;
                                }
                                self.ref_mv_idx = idx + 1;
                            }
                        }
                    }

                    self.assign_mv(is_compound, frame, out_bits, part_tree);
                }
            } else {
                // intra block mode info

            }

        } else {
            // inter frame mode info

        }
    }

    pub fn find_mv_stack(&mut self, is_compound: bool, frame: &mut Frame, part_tree: &mut PartitionTree) {
        let bw4 = num_4x4_blocks_wide[part_tree.bsize as usize];
        let bh4 = num_4x4_blocks_high[part_tree.bsize as usize];
        self.num_mv_found = 0;
        self.new_mv_count = 0;

        {
            let mut setup_global_mv = |ref_list: usize| -> () {
                let rf = part_tree.ref_frame[ref_list];
                let ecx = self.ectx.borrow_mut();
                let typ = if rf!=INTRA_FRAME { ecx.gm_type[rf as usize] } else { 0 };
                let bw = block_width[part_tree.bsize as usize];
                let bh = block_height[part_tree.bsize as usize];
                let mut cand_mv = if rf==INTRA_FRAME || typ==IDENTITY {
                    MotionVector::new(0, 0)
                } else if typ==TRANSLATION {
                    MotionVector::new(
                        ecx.gm_params[rf as usize][0] >> (WARPEDMODEL_PREC_BITS-3),
                        ecx.gm_params[rf as usize][1] >> (WARPEDMODEL_PREC_BITS-3)
                    )
                } else {
                    let x = (self.mi_col * MI_SIZE + bw/2 - 1) as isize;
                    let y = (self.mi_row * MI_SIZE + bh/2 - 1) as isize;
                    let xc = (ecx.gm_params[rf as usize][2]-(1<<WARPEDMODEL_PREC_BITS)) * x +
                            ecx.gm_params[rf as usize][3] * y +
                            ecx.gm_params[rf as usize][0];
                    let yc = ecx.gm_params[rf as usize][4] * x +
                            (ecx.gm_params[rf as usize][5]-(1<<WARPEDMODEL_PREC_BITS)) * y +
                            ecx.gm_params[rf as usize][1];
                    if frame.allow_high_precision_mv {
                        MotionVector::new(
                            Round2Signed!(xc, WARPEDMODEL_PREC_BITS-3),
                            Round2Signed!(yc, WARPEDMODEL_PREC_BITS-3)
                        )
                    } else {
                        MotionVector::new(
                            Round2Signed!(xc, WARPEDMODEL_PREC_BITS-2) * 2,
                            Round2Signed!(yc, WARPEDMODEL_PREC_BITS-2) * 2
                        )
                    }
                };
                let lower_mv_precision = |v: isize| -> isize {
                    let mut r = v;
                    if frame.force_integer_mv>0 {
                        let a = Abs!(v);
                        let a_int = (a+3) >> 3;
                        if v>0 {
                            a_int << 3
                        } else {
                            -(a_int << 3)
                        }
                    } else {
                        if v&1 != 0 {
                            if v>0 {
                                v-1
                            } else {
                                v+1
                            }
                        } else { v }
                    }
                };
                self.global_mvs[ref_list] = MotionVector::new(
                    lower_mv_precision(cand_mv.x),
                    lower_mv_precision(cand_mv.y)
                )
            };
            setup_global_mv(0);
            if is_compound {
                setup_global_mv(1);
            }
        }

        self.found_match = false;
        {
            let mut scan_row = |delta_row: isize, is_compound: bool| -> () {
                let bw4 = num_4x4_blocks_wide[part_tree.bsize as usize];
                let end4 = Min!(Min!(bw4, frame.mi_cols()-self.mi_col), 16);
                let mut d_col = 0;
                let use_step16 = bw4 >= 16;
                let mut d_row = delta_row as isize;
                if Abs!(d_row as isize) > 1 {
                    d_row += self.mi_row as isize & 1;
                    d_col = 1 - (self.mi_col as isize)&1;
                }
                let mut i: isize = 0;
                while i < end4 as isize {
                    let mv_row = self.mi_row as isize + d_row;
                    let mv_col = self.mi_col as isize + d_col + i;
                    if !self.is_inside(mv_row, mv_col) {
                        break;
                    }
                    let mi = frame.get_mi(mv_row as usize, mv_col as usize).unwrap();
                    let mut len = Min!(bw4, num_4x4_blocks_wide[mi.mi_size as usize]);
                    if Abs!(d_row) > 1 {
                        len = Max!(2, len);
                    }
                    if use_step16 {
                        len = Max!(4, len);
                    }
                    let weight = len * 2;

                    self.add_ref_mv_candidate(mv_row as usize, mv_col as usize, is_compound, weight, &mi, frame, part_tree);

                    i += len as isize;
                }
            };
            scan_row(-1, is_compound);
        }
        let mut found_above_match = self.found_match;
        self.found_match = false;
        {
            let mut scan_col = |delta_col: isize, is_compound: bool| -> () {
                let bh4 = num_4x4_blocks_high[part_tree.bsize as usize];
                let end4 = Min!(Min!(bh4, frame.mi_rows()-self.mi_row), 16);
                let mut d_row = 0;
                let use_step16 = bh4 >= 16;
                let mut d_col = delta_col as isize;
                if Abs!(d_col as isize) > 1 {
                    d_col += self.mi_col as isize & 1;
                    d_row = 1 - (self.mi_row as isize)&1;
                }
                let mut i: isize = 0;
                while i < end4 as isize {
                    let mv_row = self.mi_row as isize + d_row + i;
                    let mv_col = self.mi_col as isize + d_col;
                    if !self.is_inside(mv_row, mv_col) {
                        break;
                    }
                    let mi = frame.get_mi(mv_row as usize, mv_col as usize).unwrap();
                    let mut len = Min!(bh4, num_4x4_blocks_high[mi.mi_size as usize]);
                    if Abs!(d_col) > 1 {
                        len = Max!(2, len);
                    }
                    if use_step16 {
                        len = Max!(4, len);
                    }
                    let weight = len * 2;

                    self.add_ref_mv_candidate(mv_row as usize, mv_col as usize, is_compound, weight, &mi, frame, part_tree);

                    i += len as isize;
                }
            };
            scan_col(-1, is_compound);
        }
        let mut found_left_match = self.found_match;
        self.found_match = false;
        if Max!(bw4, bh4) <= 16 {
            // scan point
            let delta_row: isize = -1;
            let delta_col: isize = bw4 as isize;
            let mv_row = self.mi_row as isize + delta_row;
            let mv_col = self.mi_col as isize + delta_col;
            let weight = 4;
            let mi = frame.get_mi(mv_row as usize, mv_col as usize).unwrap();
            if self.is_inside(mv_row, mv_col) && mi.ref_frame[0]!=NONE { // FIXME
                self.add_ref_mv_candidate(mv_row as usize, mv_col as usize, is_compound, weight, &mi, frame, part_tree);
            }
        }
        if self.found_match {
            found_above_match = true;
        }
        self.close_matches = found_above_match as usize + found_left_match as usize;
        let num_nearest = self.num_mv_found;
        let num_new = self.new_mv_count;
        if num_nearest > 0 {
            for idx in 0..num_nearest {
                self.weight_stack[idx] += REF_CAT_LEVEL;
            }
        }
        self.zero_mv_context = 0;
        {
            // temporal scan
            let bw4 = num_4x4_blocks_wide[part_tree.bsize as usize];
            let bh4 = num_4x4_blocks_high[part_tree.bsize as usize];
            let step_w4 = if bw4>=16 { 4 } else { 2 };
            let step_h4 = if bh4>=16 { 4 } else { 2 };
            let mut delta_row = 0;
            while delta_row < Min!(bh4, 16) {
                let mut delta_col = 0;
                while delta_col < Min!(bw4, 16) {
                    // temporal sample
                    let mv_row = (self.mi_row+delta_row) | 1;
                    let mv_col = (self.mi_col+delta_col) | 1;
                    if self.is_inside(mv_row as isize, mv_col as isize) {
                        let x8 = mv_col >> 1;
                        let y8 = mv_row >> 1;
                        if delta_row==0 && delta_col==0 {
                            self.zero_mv_context = 1;
                        }
                        let lower_mv_precision = |v: isize| -> isize {
                            let mut r = v;
                            if frame.force_integer_mv>0 {
                                let a = Abs!(v);
                                let a_int = (a+3) >> 3;
                                if v>0 {
                                    a_int << 3
                                } else {
                                    -(a_int << 3)
                                }
                            } else {
                                if v&1 != 0 {
                                    if v>0 {
                                        v-1
                                    } else {
                                        v+1
                                    }
                                } else { v }
                            }
                        };
                        if !is_compound {
                            let mut cand_mv = self.ectx.borrow_mut().motion_field_mvs[part_tree.ref_frame[0] as usize][y8][x8];
                            if cand_mv.x != -1 << 15 {
                                cand_mv = MotionVector::new(
                                    lower_mv_precision(cand_mv.x),
                                    lower_mv_precision(cand_mv.y)
                                );
                                if delta_row==0 && delta_col==0 {
                                    if Abs!(cand_mv.x-self.global_mvs[0].x) >= 16 || Abs!(cand_mv.y-self.global_mvs[0].y) >= 16 {
                                        self.zero_mv_context = 1;
                                    } else {
                                        self.zero_mv_context = 0;
                                    }
                                }
                                let mut idx = 0;
                                while idx < self.num_mv_found {
                                    if cand_mv.x==self.ref_stack_mv[idx][0].x && cand_mv.y==self.ref_stack_mv[idx][0].y {
                                        break;
                                    }
                                    idx += 1;
                                }
                                if idx < self.num_mv_found {
                                    self.weight_stack[idx] += 2;
                                } else if self.num_mv_found < MAX_REF_MV_STACK_SIZE {
                                    self.ref_stack_mv[self.num_mv_found][0] = cand_mv;
                                    self.weight_stack[self.num_mv_found] = 2;
                                    self.num_mv_found += 1;
                                }
                            }
                        } else {
                            let mut cand_mv0 = self.ectx.borrow_mut().motion_field_mvs[part_tree.ref_frame[0] as usize][y8][x8];
                            let mut cand_mv1 = self.ectx.borrow_mut().motion_field_mvs[part_tree.ref_frame[1] as usize][y8][x8];
                            if cand_mv0.x != -1<<15 && cand_mv1.x != -1<<15 {
                                cand_mv0 = MotionVector::new(
                                    lower_mv_precision(cand_mv0.x),
                                    lower_mv_precision(cand_mv0.y)
                                );
                                cand_mv1 = MotionVector::new(
                                    lower_mv_precision(cand_mv1.x),
                                    lower_mv_precision(cand_mv1.y)
                                );
                                if delta_row==0 && delta_col==0 {
                                    if Abs!(cand_mv0.x-self.global_mvs[0].x) >= 16 ||
                                       Abs!(cand_mv0.y-self.global_mvs[0].y) >= 16 ||
                                       Abs!(cand_mv1.x-self.global_mvs[1].x) >= 16 ||
                                       Abs!(cand_mv1.y-self.global_mvs[1].y) >= 16 {
                                        self.zero_mv_context = 1;
                                    } else {
                                        self.zero_mv_context = 0;
                                    }
                                }
                                let mut idx = 0;
                                while idx < self.num_mv_found {
                                    if cand_mv0.x==self.ref_stack_mv[idx][0].x &&
                                       cand_mv0.y==self.ref_stack_mv[idx][0].y &&
                                       cand_mv1.x==self.ref_stack_mv[idx][1].x &&
                                       cand_mv1.y==self.ref_stack_mv[idx][1].y {
                                        break;
                                    }
                                    idx += 1;
                                }
                                if idx < self.num_mv_found {
                                    self.weight_stack[idx] += 2;
                                } else if self.num_mv_found < MAX_REF_MV_STACK_SIZE {
                                    self.ref_stack_mv[self.num_mv_found][0] = cand_mv0;
                                    self.ref_stack_mv[self.num_mv_found][1] = cand_mv1;
                                    self.weight_stack[self.num_mv_found] = 2;
                                    self.num_mv_found += 1;
                                }
                            }
                        }
                    }

                    delta_col += step_w4
                }
                delta_row += step_h4;
            }
            let allow_extension = bh4 >= num_4x4_blocks_high[BLOCK_8X8 as usize] &&
                                  bh4 < num_4x4_blocks_high[BLOCK_64X64 as usize] &&
                                  bh4 >= num_4x4_blocks_wide[BLOCK_8X8 as usize] &&
                                  bh4 < num_4x4_blocks_wide[BLOCK_64X64 as usize];
            if allow_extension {
                let tpl_sample_pos: [[isize; 2]; 3] = [
                    [bh4 as isize, -2], [bh4 as isize, bw4 as isize], [bh4 as isize-2, bw4 as isize]
                ];
                for i in 0..3 {
                    let delta_row = tpl_sample_pos[i][0];
                    let delta_col = tpl_sample_pos[i][1];

                }
            }
        }
    }

    pub fn add_ref_mv_candidate(&mut self, mv_row: usize, mv_col: usize, is_compound: bool, weight: usize, mi: &ModeInfo, frame: &mut Frame, part_tree: &mut PartitionTree) {
        if mi.is_inter {
            if !is_compound {
                for cand_list in 0..2 {
                    if mi.ref_frame[cand_list]==part_tree.ref_frame[0] {
                        // search stack
                        let cand_mode = mi.y_mode;
                        let cand_size = mi.mi_size;
                        let large = Min!(block_width[cand_size as usize], block_height[cand_size as usize]) >= 8;
                        let mut cand_mv = if (cand_mode==GLOBALMV || cand_mode==GLOBAL_GLOBALMV) &&
                                            self.ectx.borrow_mut().gm_type[part_tree.ref_frame[0] as usize] > TRANSLATION &&
                                            large {
                            self.global_mvs[0]
                        } else {
                            mi.mvs[cand_list]
                        };
                        let lower_mv_precision = |v: isize| -> isize {
                            let mut r = v;
                            if frame.force_integer_mv>0 {
                                let a = Abs!(v);
                                let a_int = (a+3) >> 3;
                                if v>0 {
                                    a_int << 3
                                } else {
                                    -(a_int << 3)
                                }
                            } else {
                                if v&1 != 0 {
                                    if v>0 {
                                        v-1
                                    } else {
                                        v+1
                                    }
                                } else { v }
                            }
                        };
                        cand_mv = MotionVector::new(
                            lower_mv_precision(cand_mv.x),
                            lower_mv_precision(cand_mv.y)
                        );
                        let has_newmv = |y_mode: usize| -> bool {
                            return y_mode==NEWMV as usize || y_mode==NEW_NEWMV as usize ||
                                y_mode==NEAR_NEWMV as usize || y_mode==NEW_NEARMV as usize ||
                                y_mode==NEAREST_NEWMV as usize || y_mode==NEW_NEARESTMV as usize;
                        };
                        if has_newmv(cand_mode as usize) {
                            self.new_mv_count += 1;
                        }
                        self.found_match = true;
                        let mut found = false;
                        let mut idx = 0;
                        while idx < self.num_mv_found {
                            if cand_mv==self.ref_stack_mv[idx][0] {
                                found = true;
                                break;
                            }
                            idx += 1;
                        }
                        if found {
                            self.weight_stack[idx] += weight;
                        } else if self.num_mv_found < MAX_REF_MV_STACK_SIZE {
                            self.ref_stack_mv[self.num_mv_found][0] = cand_mv;
                            self.weight_stack[self.num_mv_found] = weight;
                            self.num_mv_found += 1;
                        }
                    }
                }
            } else {
                if mi.ref_frame[0]==part_tree.ref_frame[0] && mi.ref_frame[1]==part_tree.ref_frame[1] {
                    // compound search stack
                    let mut cand_mvs = mi.mvs;
                    let cand_mode = mi.y_mode;
                    let cand_size = mi.mi_size;
                    if cand_mode==GLOBAL_GLOBALMV {
                        for ref_list in 0..2 {
                            if self.ectx.borrow_mut().gm_type[part_tree.ref_frame[ref_list] as usize] > TRANSLATION {
                                cand_mvs[ref_list] = self.global_mvs[ref_list];
                            }
                        }
                    }
                    let lower_mv_precision = |v: isize| -> isize {
                        let mut r = v;
                        if frame.force_integer_mv>0 {
                            let a = Abs!(v);
                            let a_int = (a+3) >> 3;
                            if v>0 {
                                a_int << 3
                            } else {
                                -(a_int << 3)
                            }
                        } else {
                            if v&1 != 0 {
                                if v>0 {
                                    v-1
                                } else {
                                    v+1
                                }
                            } else { v }
                        }
                    };
                    for i in 0..2 {
                        cand_mvs[i] = MotionVector::new(
                            lower_mv_precision(cand_mvs[i].x),
                            lower_mv_precision(cand_mvs[i].y)
                        );
                    }
                    self.found_match = true;
                    let mut found = false;
                    let mut idx = 0;
                    while idx < self.num_mv_found {
                        if cand_mvs[0]==self.ref_stack_mv[idx][0] && cand_mvs[1]==self.ref_stack_mv[idx][1] {
                            found = true;
                            break;
                        }
                        idx += 1;
                    }
                    if found {
                        self.weight_stack[idx] += weight;
                    } else if self.num_mv_found < MAX_REF_MV_STACK_SIZE {
                        for i in 0..2 {
                            self.ref_stack_mv[self.num_mv_found][i] = cand_mvs[i];
                        }
                        self.weight_stack[self.num_mv_found] = weight;
                        self.num_mv_found += 1;
                    }
                    let has_newmv = |y_mode: usize| -> bool {
                        return y_mode==NEWMV as usize || y_mode==NEW_NEWMV as usize ||
                            y_mode==NEAR_NEWMV as usize || y_mode==NEW_NEARMV as usize ||
                            y_mode==NEAREST_NEWMV as usize || y_mode==NEW_NEARESTMV as usize;
                    };
                    if has_newmv(cand_mode as usize) {
                        self.new_mv_count += 1;
                    }
                }
            }
        }
    }

    pub fn assign_mv(&mut self, is_compound: bool, frame: &mut Frame, out_bits: &mut Vec<u8>, part_tree: &mut PartitionTree) {
        for i in 0..(1+is_compound as usize) {
            if part_tree.use_intrabc {
                part_tree.compound_mode = NEWMV;
            } else {
                let get_mode = |ref_list: usize| -> YMode {
                    let y_mode = self.y_mode;
                    return if ref_list==0 {
                        if y_mode < NEAREST_NEARESTMV as usize {
                            YMode::from_usize(y_mode).unwrap()
                        } else if y_mode==NEW_NEWMV as usize || y_mode==NEW_NEARESTMV as usize || y_mode==NEW_NEARMV as usize {
                            NEWMV
                        } else if y_mode==NEAREST_NEARESTMV as usize || y_mode==NEAREST_NEWMV as usize {
                            NEARESTMV
                        } else if y_mode==NEAR_NEARMV as usize || y_mode==NEAR_NEWMV as usize {
                            NEARMV
                        } else { GLOBALMV }
                    } else {
                        if y_mode==NEW_NEWMV as usize || y_mode==NEAREST_NEWMV as usize || y_mode==NEAR_NEWMV as usize {
                            NEWMV
                        } else if y_mode==NEAREST_NEARESTMV as usize || y_mode==NEW_NEARESTMV as usize {
                            NEARESTMV
                        } else if y_mode==NEAR_NEARMV as usize || y_mode==NEW_NEARMV as usize {
                            NEARMV
                        } else { GLOBALMV }
                    };
                };
                part_tree.compound_mode = get_mode(i);
            }
            if part_tree.use_intrabc {
                self.pred_mv[0] = self.ref_stack_mv[0][0];
                if self.pred_mv[0]==MotionVector::new(0, 0) {
                    self.pred_mv[0] = self.ref_stack_mv[1][0];
                }
                if self.pred_mv[0]==MotionVector::new(0, 0) {
                    let sb_size = if self.ectx.borrow_mut().use_128x128_superblock { BLOCK_128X128 } else { BLOCK_64X64 };
                    let sb_size4 = num_4x4_blocks_high[sb_size as usize];
                    if self.mi_row-sb_size4 < self.mi_row_start {
                        self.pred_mv[0] = MotionVector::new(
                            0,
                            -((sb_size4 * MI_SIZE + INTRABC_DELAY_PIXELS) as isize) * 8
                        );
                    } else {
                        self.pred_mv[0] = MotionVector::new(
                            -((sb_size4 * MI_SIZE * 8) as isize),
                            0
                        );
                    }
                }
            } else if part_tree.compound_mode==GLOBALMV {
                self.pred_mv[i] = self.global_mvs[i];
            } else {
                let pos = if part_tree.compound_mode==NEWMV && self.num_mv_found<=1 {
                    0
                } else if part_tree.compound_mode==NEARESTMV {
                    0
                } else { self.ref_mv_idx };
                self.pred_mv[i] = self.ref_stack_mv[pos][i];
            }
            if part_tree.compound_mode==NEWMV {
                // encode mv
                let mut diff_mv = part_tree.mvs[0] - self.pred_mv[0];
                if part_tree.use_intrabc {
                    self.mv_ctx = MV_INTRABC_CONTEXT;
                } else {
                    self.mv_ctx = 0;
                }
                let mv_joint = match diff_mv {
                    MotionVector {x: 0, y: 0} => MV_JOINT_ZERO,
                    MotionVector {x: 0, y: _} => MV_JOINT_HZVNZ,
                    MotionVector {x: _, y: 0} => MV_JOINT_HNZVZ,
                    _ => MV_JOINT_HNZVNZ,
                };
                self.coder.encode_mv_joint(out_bits, &mut self.cdf, mv_joint, self.mv_ctx);
                let encode_mv_component = |v: isize, comp: usize| -> () {
                    let mv_sign = if v<0 { 1 } else { 0 };
                    let mv_abs = Abs!(v) as usize;
                    let mv_class = MVClass::from_u8(msb16(mv_abs as u16)/8).unwrap();
                    self.coder.encode_mv_class(out_bits, &mut self.cdf, mv_class as usize, self.mv_ctx, comp);
                    if mv_class==MV_CLASS_0 {
                        let mv_class0_bit = (mv_abs >> 3) & 1;
                        self.coder.encode_mv_class0_bit(out_bits, &mut self.cdf, mv_class0_bit, self.mv_ctx, comp);
                        if frame.force_integer_mv==0 {
                            let mv_class0_fr = (mv_abs >> 1) & 3;
                            self.coder.encode_mv_class0_fr(out_bits, &mut self.cdf, mv_class0_fr, mv_class0_bit, self.mv_ctx, comp);
                        }
                        if frame.allow_high_precision_mv {
                            let mv_class0_hp = mv_abs & 1;
                            self.coder.encode_mv_class0_hp(out_bits, &mut self.cdf, mv_class0_hp, self.mv_ctx, comp);
                        }
                    } else {
                        let mag = mv_abs - (CLASS0_SIZE << (mv_class as usize + 2));
                        for i in 0..(mv_class as usize) {
                            let mv_bit = (mag>>(i+3)) & 1;
                            self.coder.encode_mv_bit(out_bits, &mut self.cdf, mv_bit, i, self.mv_ctx, comp);
                        }
                        if frame.force_integer_mv==0 {
                            let mv_fr = (mag >> 1) & 3;
                            self.coder.encode_mv_fr(out_bits, &mut self.cdf, mv_fr, self.mv_ctx, comp);
                        }
                        if frame.allow_high_precision_mv {
                            let mv_hp = mag & 1;
                            self.coder.encode_mv_hp(out_bits, &mut self.cdf, mv_hp, self.mv_ctx, comp);
                        }
                    }
                };
            } else {
                part_tree.mvs[i] = self.pred_mv[i];
            }
        }
    }
}