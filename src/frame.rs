#![allow(dead_code)]

use constants::*;
use util::*;
use sequence_header::DecoderModelInfo;
use constants::InterpFilter::*;


pub struct TemporalPointInfo {
    pub tu_presentation_delay: u32 // TODO
}

pub struct TileInfo {
    pub uniform_tile_spacing_flag: bool,
    // TODO
}

pub struct QuantizationParams {
    pub base_q_idx: u8,
    pub delta_q: Option<u8>,
    pub diff_uv_delta: bool,
    pub using_qmatrix: bool,
    pub qm_y: u8,
    pub qm_u: u8,
    pub qm_v: u8,
}

pub struct SegmentationParams {
    pub segmentation_enabled: bool,
    pub segmentation_update_map: bool,
    pub segmentation_temporal_update: bool,
    pub segmentation_update_data: bool,
    pub feature_value: [[Option<isize>; SEG_LVL_MAX as usize]; MAX_SEGMENTS as usize],
}

pub struct DeltaQParams {
    pub delta_q_present: bool,
    pub delta_q_res: u8
}

pub struct DeltaLFParams {
    pub delta_lf_present: bool,
    pub delta_lf_res: u8,
    pub delta_lf_multi: bool
}

pub struct LoopFilterParams {
    pub loop_filter_level: [u8; 4],
    pub loop_filter_sharpness: u8,
    pub loop_filter_delta_enabled: bool,
    pub loop_filter_delta_update: bool,
    pub update_ref_delta: [bool; TOTAL_REFS_PER_FRAME as usize],
    pub loop_filter_ref_deltas: [isize; TOTAL_REFS_PER_FRAME as usize],
    pub update_mode_delta: [bool; 2],
    pub loop_filter_mode_deltas: [u8; 2]
}

pub struct CDEFParams {
    pub cdef_damping: u8,
    pub cdef_bits: u8,
    pub cdef_y_pri_strength: Vec<u8>,
    pub cdef_y_sec_strength: Vec<u8>,
    pub cdef_uv_pri_strength: Vec<u8>,
    pub cdef_uv_sec_strength: Vec<u8>
}

pub struct LRParams {
    pub lr_type: u8,
    pub lr_unit_shift: usize,
    pub lr_unit_extra_shift: usize,
    pub lr_uv_shift: usize,
}

pub struct GlobalMotionParams {
    pub is_global: [bool; NUM_REF_FRAMES as usize],
    pub is_rot_zoom: [bool; NUM_REF_FRAMES as usize],
    pub is_translation: [bool; NUM_REF_FRAMES as usize],
}

pub struct FilmGrainParams {
    pub apply_grain: bool,
    pub grain_seed: u16,
    pub update_grain: bool,
    pub film_grain_params_ref_idx: usize,
    pub num_y_points: usize,
    pub point_y_value: Vec<u8>,
    pub point_y_scaling: Vec<u8>,
    pub chroma_scaling_from_luma: bool,
    pub num_cb_points: usize,
    pub num_cr_points: usize,
    pub point_cb_value: Vec<u8>,
    pub point_cb_scaling: Vec<u8>,
    pub point_cr_value: Vec<u8>,
    pub point_cr_scaling: Vec<u8>,
    pub grain_scaling: usize,
    pub ar_coeff_lag: usize,
    pub ar_coeffs_y: Vec<i16>,
    pub ar_coeffs_cb: Vec<i16>,
    pub ar_coeffs_cr: Vec<i16>,
    pub ar_coeff_shift: usize,
    pub grain_scale_shift: usize,
    pub cb_mult: usize,
    pub cb_luma_mult: usize,
    pub cb_offset: usize,
    pub cr_mult: usize,
    pub cr_luma_mult: usize,
    pub cr_offset: usize,
    pub overlap_flag: bool,
    pub clip_to_restricted_range: bool
}

impl FilmGrainParams {
    pub fn new() -> FilmGrainParams { FilmGrainParams {
        apply_grain: false,
        grain_seed: 0,
        update_grain: false,
        film_grain_params_ref_idx: 0,
        num_y_points: 0,
        point_y_value: vec![],
        point_y_scaling: vec![],
        chroma_scaling_from_luma: false,
        num_cb_points: 0,
        num_cr_points: 0,
        point_cb_value: vec![],
        point_cb_scaling: vec![],
        point_cr_value: vec![],
        point_cr_scaling: vec![],
        grain_scaling: 0,
        ar_coeff_lag: 0,
        ar_coeffs_y: vec![],
        ar_coeffs_cb: vec![],
        ar_coeffs_cr: vec![],
        ar_coeff_shift: 0,
        grain_scale_shift: 0,
        cb_mult: 0,
        cb_luma_mult: 0,
        cb_offset: 0,
        cr_mult: 0,
        cr_luma_mult: 0,
        cr_offset: 0,
        overlap_flag: false,
        clip_to_restricted_range: false,
    }}

    pub fn clone(&self) -> FilmGrainParams { FilmGrainParams {
        apply_grain: self.apply_grain,
        grain_seed: self.grain_seed,
        update_grain: self.update_grain,
        film_grain_params_ref_idx: self.film_grain_params_ref_idx,
        num_y_points: self.num_y_points,
        point_y_value: self.point_y_value.clone(),
        point_y_scaling: self.point_y_scaling.clone(),
        chroma_scaling_from_luma: self.chroma_scaling_from_luma,
        num_cb_points: self.num_cb_points,
        num_cr_points: self.num_cr_points,
        point_cb_value: self.point_cb_value.clone(),
        point_cb_scaling: self.point_cb_scaling.clone(),
        point_cr_value: self.point_cr_value.clone(),
        point_cr_scaling: self.point_cr_scaling.clone(),
        grain_scaling: self.grain_scaling,
        ar_coeff_lag: self.ar_coeff_lag,
        ar_coeffs_y: self.ar_coeffs_y.clone(),
        ar_coeffs_cb: self.ar_coeffs_cb.clone(),
        ar_coeffs_cr: self.ar_coeffs_cr.clone(),
        ar_coeff_shift: self.ar_coeff_shift,
        grain_scale_shift: self.grain_scale_shift,
        cb_mult: self.cb_mult,
        cb_luma_mult: self.cb_luma_mult,
        cb_offset: self.cb_offset,
        cr_mult: self.cr_mult,
        cr_luma_mult: self.cr_luma_mult,
        cr_offset: self.cr_offset,
        overlap_flag: self.overlap_flag,
        clip_to_restricted_range: self.clip_to_restricted_range,
    }}
}

pub struct Frame {
    pub frame_width: usize,
    pub frame_height: usize,
    pub render_width: usize,
    pub render_height: usize,
    pub frame_type: FrameType,
    pub subsampling_x: bool,
    pub subsampling_y: bool,
    pub tile_col_starts: Vec<usize>,
    pub tile_row_starts: Vec<usize>,
    pub tile_groups: Vec<usize>,

    pub original_y: Vec<i16>,
    pub original_u: Vec<i16>,
    pub original_v: Vec<i16>,

    pub mis: Vec<Vec<Option<ModeInfo>>>,
    pub decoder_model_info: Option<DecoderModelInfo>,
    pub temporal_point_info: Option<TemporalPointInfo>,
    pub display_frame_id: Option<u32>, // TODO
    pub error_resilient_mode: bool,
    pub allow_screen_content_tools: usize,
    pub force_integer_mv: usize,
    pub allow_high_precision_mv: bool,
    pub use_ref_frame_mvs: bool,
    pub order_hint: u8, // TODO
    pub primary_ref_frame: u8,
    pub allow_intrabc: bool,
    pub ref_order_hint: [usize; NUM_REF_FRAMES as usize], // TODO
    pub use_superres: bool,
    pub frame_refs_short_signaling: bool,
    pub last_frame_idx: u8,
    pub gold_frame_idx: u8,
    pub ref_frame_idx: [u8; REFS_PER_FRAME as usize],
    pub delta_frame_id: [usize; REFS_PER_FRAME as usize],
    pub is_filter_switchable: bool,
    pub interp_filter: InterpFilter,
    pub is_motion_mode_switchable: bool,
    pub tile_info: TileInfo,
    pub quantization_params: QuantizationParams,
    pub segmentation_params: SegmentationParams,
    pub delta_q_params: DeltaQParams,
    pub delta_lf_params: DeltaLFParams,
    pub loop_filter_params: LoopFilterParams,
    pub cdef_params: CDEFParams,
    pub lr_params: LRParams,
    pub tx_mode_select: bool,
    pub reference_select: bool,
    pub skip_mode_present: bool,
    pub allow_warped_motion: bool,
    pub reduced_tx_set: bool,
    pub global_motion_params: GlobalMotionParams,
}

impl Frame {
    pub fn new(frame_type: FrameType, width: usize, height: usize) -> Frame {
        Frame {
            frame_width: width,
            frame_height: height,
            render_width: width,
            render_height: height,
            frame_type: frame_type,
            subsampling_x: true,
            subsampling_y: true,
            tile_col_starts: vec![0],
            tile_row_starts: vec![0],
            tile_groups: vec![0],

            original_y: vec![0; width as usize * height as usize],
            original_u: vec![0; width as usize * height as usize / 4],
            original_v: vec![0; width as usize * height as usize / 4],

            mis: vec![vec![None; width as usize / 4]; height as usize / 4],
            decoder_model_info: None,
            temporal_point_info: None,
            display_frame_id: None, // TODO
            error_resilient_mode: false,
            allow_screen_content_tools: 0,
            force_integer_mv: 0,
            allow_high_precision_mv: true,
            use_ref_frame_mvs: true,
            order_hint: 0, // TODO
            primary_ref_frame: 0,
            allow_intrabc: false,
            ref_order_hint: [0; NUM_REF_FRAMES as usize], // TODO
            use_superres: false,
            frame_refs_short_signaling: false,
            last_frame_idx: 0,
            gold_frame_idx: 0,
            ref_frame_idx: [0; REFS_PER_FRAME as usize],
            delta_frame_id: [0; REFS_PER_FRAME as usize],
            is_filter_switchable: false,
            interp_filter: SWITCHABLE,
            is_motion_mode_switchable: true,
            tile_info: TileInfo {
                uniform_tile_spacing_flag: false,
            },
            quantization_params: QuantizationParams {
                base_q_idx: 0,
                delta_q: None,
                diff_uv_delta: false,
                using_qmatrix: false,
                qm_y: 0,
                qm_u: 0,
                qm_v: 0,
            },
            segmentation_params: SegmentationParams {
                segmentation_enabled: false,
                segmentation_update_map: false,
                segmentation_temporal_update: false,
                segmentation_update_data: false,
                feature_value: [[None; SEG_LVL_MAX as usize]; MAX_SEGMENTS as usize],
            },
            delta_q_params: DeltaQParams {
                delta_q_present: false,
                delta_q_res: 0
            },
            delta_lf_params: DeltaLFParams {
                delta_lf_present: false,
                delta_lf_res: 0,
                delta_lf_multi: false
            },
            loop_filter_params: LoopFilterParams {
                loop_filter_level: [0; 4],
                loop_filter_sharpness: 0,
                loop_filter_delta_enabled: false,
                loop_filter_delta_update: false,
                update_ref_delta: [false; TOTAL_REFS_PER_FRAME as usize],
                loop_filter_ref_deltas: [0; TOTAL_REFS_PER_FRAME as usize],
                update_mode_delta: [false; 2],
                loop_filter_mode_deltas: [0; 2]
            },
            cdef_params: CDEFParams {
                cdef_damping: 0,
                cdef_bits: 0,
                cdef_y_pri_strength: vec![],
                cdef_y_sec_strength: vec![],
                cdef_uv_pri_strength: vec![],
                cdef_uv_sec_strength: vec![]
            },
            lr_params: LRParams {
                lr_type: 0,
                lr_unit_shift: 0,
                lr_unit_extra_shift: 0,
                lr_uv_shift: 0,
            },
            tx_mode_select: true,
            reference_select: true,
            skip_mode_present: true,
            allow_warped_motion: false,
            reduced_tx_set: false,
            global_motion_params: GlobalMotionParams {
                is_global: [false; NUM_REF_FRAMES as usize],
                is_rot_zoom: [false; NUM_REF_FRAMES as usize],
                is_translation: [false; NUM_REF_FRAMES as usize],
            },
        }
    }

    pub fn get_mi(&mut self, mi_row: usize, mi_col: usize) -> Option<ModeInfo> {
        self.mis[mi_row][mi_col]
    }

    pub fn get_above_mi(&mut self, mi_row: usize, mi_col: usize) -> Option<ModeInfo> {
        if mi_row==0 { None }
        else { self.mis[mi_row-1][mi_col] }
    }

    pub fn get_left_mi(&mut self, mi_row: usize, mi_col: usize) -> Option<ModeInfo> {
        if mi_col==0 { None }
        else { self.mis[mi_row][mi_col-1] }
    }

    pub fn mi_rows(&self) -> usize {
        2 * ((self.frame_height+7)>>3)
    }

    pub fn mi_cols(&self) -> usize {
        2 * ((self.frame_width+7)>>3)
    }

    pub fn tile_cols(&self) -> usize {
        self.tile_col_starts.len()
    }

    pub fn tile_rows(&self) -> usize {
        self.tile_row_starts.len()
    }

    pub fn tile_cols_log2(&self) -> usize {
        msb16(self.tile_col_starts.len() as u16) as usize
    }

    pub fn tile_rows_log2(&self) -> usize {
        msb16(self.tile_row_starts.len() as u16) as usize
    }
}