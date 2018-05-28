#![allow(dead_code)]

use sequence_header::DecoderModelInfo;
use constants::FrameType;
use constants::*;

struct TemporalPointInfo {
    tu_presentation_delay: u32 // TODO
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
    pub feature_value: [[Option<bool>; SEG_LVL_MAX as usize]; MAX_SEGMENTS as usize],
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
    pub loop_filter_ref_deltas: [u8; TOTAL_REFS_PER_FRAME as usize],
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
    pub lr_unit_shift: bool,
    pub lr_unit_extra_shift: bool,
    pub lr_uv_shift: bool,
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
    pub film_grain_params_ref_idx: u8,
    pub num_y_points: u8,
    pub point_y_value: Vec<u8>,
    pub point_y_scaling: Vec<u8>,
    pub chroma_scaling_from_luma: bool,
    pub num_cb_points: u8,
    pub point_cb_value: Vec<u8>,
    pub point_cb_scaling: Vec<u8>,
    pub point_cr_value: Vec<u8>,
    pub point_cr_scaling: Vec<u8>,
    pub grain_scaling: u8,
    pub ar_coeff_lag: u8,
    pub ar_coeffs_y: Vec<i16>,
    pub ar_coeffs_cb: Vec<i16>,
    pub ar_coeffs_cr: Vec<i16>,
    pub ar_coeff_shift: u8,
    pub grain_scale_shift: u8,
    pub cb_mult: u8,
    pub cb_luma_mult: u8,
    pub cb_offset: u16,
    pub cr_mult: u8,
    pub cr_luma_mult: u8,
    pub cr_offset: u16,
    pub overlap_flag: bool,
    pub clip_to_restricted_range: bool
}

pub struct FrameHeader {
    show_existing_frame: bool,
    frame_to_show_map_idx: Option<u8>,
    display_frame_id: Option<u32>, // TODO
    temporal_point_info: Option<TemporalPointInfo>,
    frame_type: FrameType,
    show_frame: bool,
    showable_frame: bool,
    error_resilient_mode: bool,
    disable_cdf_update: bool,
    allow_screen_content_tools: bool,
    force_integer_mv: bool,
    current_frame_id: u32, // TODO
    frame_size_override_flag: bool,
    order_hint: u8, // TODO
    primary_ref_frame: u8,
    decoder_model_info: Option<DecoderModelInfo>,
    refresh_frame_flags: u8,
    ref_order_hint: [u32; NUM_REF_FRAMES as usize], // TODO
    frame_width: u8,
    frame_height: u8,
    use_superres: bool,
    coded_denom: u8,
    render_and_frame_size_different: bool,
    render_width: u8,
    render_height: u8,
    allow_intrabc: bool,
    frame_refs_short_signaling: bool,
    last_frame_idx: u8,
    gold_frame_idx: u8,
    ref_frame_idx: [u8; REFS_PER_FRAME as usize],
    delta_frame_id: [u8; REFS_PER_FRAME as usize],
    allow_high_precision_mv: bool,
    is_motion_mode_switchable: bool,
    use_ref_frame_mvs: bool,
    disable_frame_end_update_cdf: bool,
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
}