#![allow(dead_code)]

use constants::*;
use frame_header::*;
//use block::*;

pub struct Frame<'a> {
    width: usize,
    height: usize,
    frame_type: FrameType,
    pub subsampling_x: bool,
    pub subsampling_y: bool,
    showable_frame: bool,
    error_resilient_mode: bool,
    allow_screen_content_tools: bool,
    force_integer_mv: bool,
    id: u32,
    order_hint: u8,
    use_superres: bool,
    allow_intrabc: bool,
    pub ref_frames: [Option<&'a Frame<'a>>; NUM_REF_FRAMES as usize],
    allow_high_precision_mv: bool,
    is_motion_mode_switchable: bool,
    use_ref_frame_mvs: bool,
    tile_info: TileInfo,
    quantization_params: QuantizationParams,
    segmentation_params: SegmentationParams,
    delta_q_params: DeltaQParams,
    delta_lf_params: DeltaLFParams,
    loop_filter_params: LoopFilterParams,
    cdef_params: CDEFParams,
    lr_params: LRParams,
    tx_mode_select: bool,
    reference_select: bool,
    skip_mode_present: bool,
    allow_warped_motion: bool,
    reduced_tx_set: bool,
    global_motion_params: GlobalMotionParams,
    film_frain_params: FilmGrainParams,

    pub original_y: Vec<i16>,
    pub original_u: Vec<i16>,
    pub original_v: Vec<i16>,

    pub mis: Vec<Vec<Option<ModeInfo>>>,
}

impl<'a> Frame<'a> {
    pub fn new(frame_type: FrameType, width: usize, height: usize) -> Frame<'a> {
        Frame {
            width: width,
            height: height,
            frame_type: frame_type,
            subsampling_x: true,
            subsampling_y: true,
            showable_frame: true,
            error_resilient_mode: false,
            allow_screen_content_tools: false,
            force_integer_mv: false,
            id: 0,
            order_hint: 0,
            use_superres: false,
            allow_intrabc: false,
            ref_frames: [None; NUM_REF_FRAMES as usize],
            allow_high_precision_mv: true,
            is_motion_mode_switchable: true,
            use_ref_frame_mvs: true,
            tile_info: TileInfo { uniform_tile_spacing_flag: true },
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
                lr_unit_shift: false,
                lr_unit_extra_shift: false,
                lr_uv_shift: false,
            },
            tx_mode_select: false,
            reference_select: false,
            skip_mode_present: false,
            allow_warped_motion: false,
            reduced_tx_set: false,
            global_motion_params: GlobalMotionParams {
                is_global: [false; NUM_REF_FRAMES as usize],
                is_rot_zoom: [false; NUM_REF_FRAMES as usize],
                is_translation: [false; NUM_REF_FRAMES as usize],
            },
            film_frain_params: FilmGrainParams {
                apply_grain: false,
                grain_seed: 0,
                update_grain: false,
                film_grain_params_ref_idx: 0,
                num_y_points: 0,
                point_y_value: vec![],
                point_y_scaling: vec![],
                chroma_scaling_from_luma: false,
                num_cb_points: 0,
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
                clip_to_restricted_range: false
            },

            original_y: vec![0; width as usize * height as usize],
            original_u: vec![0; width as usize * height as usize / 4],
            original_v: vec![0; width as usize * height as usize / 4],

            mis: vec![vec![None; width as usize / 4]; height as usize / 4],
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
        self.height / 4
    }

    pub fn mi_cols(&self) -> usize {
        self.width / 4
    }
}