use constants::{ColorPrimaries, TransferCharacteristics, MatrixCoefficients};

struct OperatingPoint {
    idc: u16,
    seq_level_idx: u8,
    seq_tier: Option<bool>,
}

struct FrameID {
    delta_frame_id_length: u8,
    additional_frame_id_length: u8,
}

struct ColorDescription {
    color_primaries: ColorPrimaries,
    transfer_characteristics: TransferCharacteristics,
    matrix_coefficients: MatrixCoefficients
}

struct ColorConfig {
    high_bitdepth: bool,
    twelve_bit: bool,
    mono_chrome: bool,
    color_description: Option<ColorDescription>,
    color_range: bool,
    subsampling_x: bool,
    subsampling_y: bool,
    chroma_sample_position: u8,
    separate_uv_delta_q: bool,
}

struct TimingInfo {
    num_units_in_display_tick: u32,
    time_scale: u32,
    equal_picture_interval: bool,
    num_ticks_per_picture: u32
}

pub struct DecoderModelInfo {
    bitrate_scale: u8,
    buffer_size_scale: u8,
    encoder_decoder_buffer_delay_length: u8,
    num_units_in_decoding_tick: u32,
    buffer_removal_delay_length: u8,
    frame_presentation_delay_length: u8
}

struct OperatingParametersInfo {
    bitrate: u32,
    buffer_size: u32,
    cbr_flag: bool,
    decoder_buffer_delay: u32,
    encoder_buffer_delay: u32,
    low_delay_mode_flag: bool
}

struct OperatingPointsDecoderModel {
    idc: Vec<u16>,
    initial_display_delay: Option<u8>,
    operating_parameters_info: Option<OperatingParametersInfo>
}

pub struct SequenceHeader {
    profile: u8,
    still_picture: bool,
    reduced_still_picture_header: bool,
    operating_points: Vec<OperatingPoint>,
    frame_width: u16,
    frame_height: u16,
    max_frame_width: u16,
    max_frame_height: u16,
    frame_id: Option<FrameID>,
    use_128x128_superblock: bool,
    enable_filter_intra: bool,
    enable_intra_edge_filter: bool,
    enable_interintra_compound: bool,
    enable_masked_compound: bool,
    enable_warped_motion: bool,
    enable_dual_filter: bool,
    enable_order_hint: bool,
    enable_jnt_comp: bool,
    enable_ref_frame_mvs: bool,
    seq_choose_screen_content_tools: bool,
    seq_force_screen_content_tools: bool,
    seq_choose_integer_mv: bool,
    seq_force_integer_mv: bool,
    order_hint_bits: u8,
    enable_superres: bool,
    enable_cdef: bool,
    enable_restoration: bool,
    color_config: ColorConfig,
    timing_info: Option<TimingInfo>,
    decoder_model_info: Option<DecoderModelInfo>,
    operating_points_decoder_model: Vec<OperatingPointsDecoderModel>,
    film_grain_params_present: bool
}