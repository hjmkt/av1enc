#![allow(dead_code)]

use constants::*;

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

#[derive(PartialEq, Copy, Clone)]
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

#[derive(Copy, Clone)]
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
    pub buffer_removal_delay: usize,
    pub buffer_removal_delay_length: u8,
    pub buffer_removal_delay_present: bool,
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
    pub idc: u16,
    pub initial_display_delay: Option<u8>,
    pub operating_parameters_info: Option<OperatingParametersInfo>
}

pub struct SequenceHeader {
    pub profile: u8,
    pub still_picture: bool,
    pub operating_points: Vec<OperatingPoint>,
    pub frame_width: usize,
    pub frame_height: usize,
    pub max_frame_width: usize,
    pub max_frame_height: usize,
    pub frame_id: Option<FrameID>,
    pub enable_filter_intra: bool,
    pub enable_intra_edge_filter: bool,
    pub enable_interintra_compound: bool,
    pub enable_masked_compound: bool,
    pub enable_dual_filter: bool,
    pub enable_jnt_comp: bool,
    pub seq_choose_screen_content_tools: bool,
    pub seq_choose_integer_mv: bool,
    pub enable_cdef: bool,
    pub decoder_model_info: Option<DecoderModelInfo>,
}

impl SequenceHeader {
    pub fn new(frame_width: usize, frame_height: usize) -> SequenceHeader { SequenceHeader {
        profile: 0,
        still_picture: false,
        operating_points: vec![],
        frame_width: frame_width,
        frame_height: frame_height,
        max_frame_width: frame_width,
        max_frame_height: frame_height,
        frame_id: None,
        enable_filter_intra: false,
        enable_intra_edge_filter: false,
        enable_interintra_compound: false,
        enable_masked_compound: false,
        enable_dual_filter: false,
        enable_jnt_comp: false,
        seq_choose_screen_content_tools: false,
        seq_choose_integer_mv: false,
        enable_cdef: false,
        decoder_model_info: None,
    }}
}