#![allow(dead_code)]

use constants::*;
use constants::ChromaSamplePosition::*;

pub struct OperatingPoint {
    pub idc: u16,
    pub seq_level_idx: u8,
    pub seq_tier: Option<bool>,
}

pub struct FrameID {
    pub delta_frame_id_length: u8,
    pub additional_frame_id_length: u8,
}

#[derive(PartialEq, Copy, Clone)]
pub struct ColorDescription {
    pub color_primaries: ColorPrimaries,
    pub transfer_characteristics: TransferCharacteristics,
    pub matrix_coefficients: MatrixCoefficients
}

pub struct ColorConfig {
    pub high_bitdepth: bool,
    pub twelve_bit: bool,
    pub mono_chrome: bool,
    pub color_description: Option<ColorDescription>,
    pub color_range: bool,
    pub subsampling_x: bool,
    pub subsampling_y: bool,
    pub chroma_sample_position: ChromaSamplePosition,
    pub separate_uv_delta_q: bool,
}

pub struct TimingInfo {
    pub num_units_in_display_tick: u32,
    pub time_scale: u32,
    pub equal_picture_interval: bool,
    pub num_ticks_per_picture: u32
}

pub struct DecoderModelInfo {
    pub bitrate_scale: u8,
    pub buffer_size_scale: u8,
    pub encoder_decoder_buffer_delay_length: u8,
    pub num_units_in_decoding_tick: u32,
    pub buffer_removal_delay_length: u8,
    pub frame_presentation_delay_length: u8
}

pub struct OperatingParametersInfo {
    pub bitrate: u32,
    pub buffer_size: u32,
    pub cbr_flag: bool,
    pub decoder_buffer_delay: u32,
    pub encoder_buffer_delay: u32,
    pub low_delay_mode_flag: bool
}

pub struct OperatingPointsDecoderModel {
    pub idc: Vec<u16>,
    pub initial_display_delay: Option<u8>,
    pub operating_parameters_info: Option<OperatingParametersInfo>
}

pub struct SequenceHeader {
    pub profile: u8,
    pub still_picture: bool,
    pub reduced_still_picture_header: bool,
    pub operating_points: Vec<OperatingPoint>,
    pub frame_width: usize,
    pub frame_height: usize,
    pub max_frame_width: usize,
    pub max_frame_height: usize,
    pub frame_id: Option<FrameID>,
    pub use_128x128_superblock: bool,
    pub enable_filter_intra: bool,
    pub enable_intra_edge_filter: bool,
    pub enable_interintra_compound: bool,
    pub enable_masked_compound: bool,
    pub enable_warped_motion: bool,
    pub enable_dual_filter: bool,
    pub enable_order_hint: bool,
    pub enable_jnt_comp: bool,
    pub enable_ref_frame_mvs: bool,
    pub seq_choose_screen_content_tools: bool,
    pub seq_force_screen_content_tools: usize,
    pub seq_choose_integer_mv: bool,
    pub seq_force_integer_mv: usize,
    pub order_hint_bits: u8,
    pub enable_superres: bool,
    pub enable_cdef: bool,
    pub enable_restoration: bool,
    pub color_config: ColorConfig,
    pub timing_info: Option<TimingInfo>,
    pub decoder_model_info: Option<DecoderModelInfo>,
    pub operating_points_decoder_model: Vec<OperatingPointsDecoderModel>,
    pub film_grain_params_present: bool
}

impl SequenceHeader {
    pub fn new(frame_width: usize, frame_height: usize) -> SequenceHeader { SequenceHeader {
        profile: 0,
        still_picture: false,
        reduced_still_picture_header: false,
        operating_points: vec![],
        frame_width: frame_width,
        frame_height: frame_height,
        max_frame_width: frame_width,
        max_frame_height: frame_height,
        frame_id: None,
        use_128x128_superblock: true,
        enable_filter_intra: false,
        enable_intra_edge_filter: false,
        enable_interintra_compound: false,
        enable_masked_compound: false,
        enable_warped_motion: false,
        enable_dual_filter: false,
        enable_order_hint: false,
        enable_jnt_comp: false,
        enable_ref_frame_mvs: false,
        seq_choose_screen_content_tools: false,
        seq_force_screen_content_tools: 0,
        seq_choose_integer_mv: false,
        seq_force_integer_mv: 0,
        order_hint_bits: 0,
        enable_superres: false,
        enable_cdef: false,
        enable_restoration: false,
        color_config: ColorConfig {
            high_bitdepth: false,
            twelve_bit: false,
            mono_chrome: false,
            color_description: None,
            color_range: false,
            subsampling_x: true,
            subsampling_y: true,
            chroma_sample_position: CSP_COLOCATED,
            separate_uv_delta_q: false,
        },
        timing_info: None,
        decoder_model_info: None,
        operating_points_decoder_model: vec![],
        film_grain_params_present: false
    }}
}