use sequence_header::*;
use constants::*;
use constants::FrameType::*;
use frame::*;
use constants::ChromaSamplePosition::*;
use constants::FrameRestorationType::*;
use constants::TxMode::*;

pub struct EncoderContext {
    pub show_existing_frame: bool,
    pub show_frame: bool,
    pub showable_frame: bool,
    pub frame_to_show_map_idx: u8,
    pub frame_id_numbers_present_flag: bool,
    pub disable_cdf_update: bool,
    pub seen_frame_header: bool,
    pub additional_frame_id_length: usize,
    pub delta_frame_id_length: usize,
    pub reduced_still_picture_header: bool,
    pub frame_is_intra: bool,
    pub timing_info: Option<TimingInfo>,
    pub refresh_frame_flags: u8,
    pub ref_frame_type: [FrameType; NUM_REF_FRAMES],
    pub ref_valid: [bool; NUM_REF_FRAMES],
    pub ref_frame_id: [usize; NUM_REF_FRAMES],
    pub ref_order_hint: [usize; NUM_REF_FRAMES],
    pub gm_type: [usize; NUM_REF_FRAMES],
    pub gm_params: [[isize; 6]; NUM_REF_FRAMES],
    pub prev_gm_params: [[isize; 6]; NUM_REF_FRAMES],
    pub film_grain_params: Option<FilmGrainParams>,
    pub ref_film_grain_params: [FilmGrainParams; NUM_REF_FRAMES],
    pub order_hints: [usize; REFS_PER_FRAME],
    pub expected_frame_id: [usize; REFS_PER_FRAME],
    pub ref_frame_sign_bias: [usize; REFS_PER_FRAME],
    pub seq_force_screen_content_tools: usize,
    pub seq_force_integer_mv: usize,
    pub prev_frame_id: usize,
    pub current_frame_id: usize,
    pub operating_points_decoder_model: Vec<OperatingPointsDecoderModel>,
    pub temporal_id: usize,
    pub spatial_id: usize,
    pub enable_order_hint: bool,
    pub superres_denom: usize,
    pub upscaled_width: usize,
    pub enable_superres: bool,
    pub enable_ref_frame_mvs: bool,
    pub order_hint_bits: u8,
    pub use_128x128_superblock: bool,
    pub mi_col_starts: Vec<usize>,
    pub mi_row_starts: Vec<usize>,
    pub context_update_tile_id: usize,
    pub tile_size_bytes: usize,
    pub delta_q_y_dc: isize,
    pub delta_q_u_dc: isize,
    pub delta_q_u_ac: isize,
    pub delta_q_v_dc: isize,
    pub delta_q_v_ac: isize,
    pub num_planes: usize,
    pub color_config: ColorConfig,
    pub seg_id_pre_skip: bool,
    pub last_active_seg_id: usize,
    pub coded_lossless: bool,
    pub all_lossless: bool,
    pub lossless_array: [bool; MAX_SEGMENTS],
    pub seg_q_m_level: [[usize; MAX_SEGMENTS]; 3],
    pub current_q_index: usize,
    pub frame_restoration_type: [FrameRestorationType; 3],
    pub loop_restoration_size: [usize; 3],
    pub enable_restoration: bool,
    pub uses_lr: bool,
    pub uses_chroma_lr: bool,
    pub tx_mode: TxMode,
    pub skip_mode_allowed: bool,
    pub skip_mode_frame: [usize; 2],
    pub enable_warped_motion: bool,
    pub tg_start: usize,
    pub tg_end: usize,
    pub frame_size_override_flag: bool,
    pub disable_frame_end_update_cdf: bool,
}

impl EncoderContext {
    pub fn new() -> EncoderContext { EncoderContext {
        show_existing_frame: false,
        show_frame: true,
        showable_frame: true,
        frame_to_show_map_idx: 0,
        frame_id_numbers_present_flag: false,
        disable_cdf_update: false,
        seen_frame_header: false,
        additional_frame_id_length: 0,
        delta_frame_id_length: 0,
        reduced_still_picture_header: false,
        frame_is_intra: false,
        timing_info: None,
        refresh_frame_flags: 0,
        ref_frame_type: [KEY_FRAME; NUM_REF_FRAMES],
        ref_frame_id: [0; NUM_REF_FRAMES],
        ref_valid: [false; NUM_REF_FRAMES],
        ref_order_hint: [0; NUM_REF_FRAMES],
        gm_type: [0; NUM_REF_FRAMES],
        gm_params: [[0; 6]; NUM_REF_FRAMES],
        prev_gm_params: [[0; 6]; NUM_REF_FRAMES],
        film_grain_params: None,
        ref_film_grain_params: [ // FIXME
            FilmGrainParams::new(), FilmGrainParams::new(), FilmGrainParams::new(), FilmGrainParams::new(), FilmGrainParams::new(), FilmGrainParams::new(), FilmGrainParams::new(), FilmGrainParams::new()
        ],
        order_hints: [0; REFS_PER_FRAME],
        expected_frame_id: [0; REFS_PER_FRAME],
        ref_frame_sign_bias: [0; REFS_PER_FRAME],
        seq_force_screen_content_tools: 0,
        seq_force_integer_mv: 0,
        prev_frame_id: 0,
        current_frame_id: 0,
        operating_points_decoder_model: vec![],
        temporal_id: 0,
        spatial_id: 0,
        enable_order_hint: false,
        superres_denom: 1,
        upscaled_width: 0,
        enable_superres: false,
        enable_ref_frame_mvs: false,
        order_hint_bits: 0,
        use_128x128_superblock: true,
        mi_col_starts: vec![],
        mi_row_starts: vec![],
        context_update_tile_id: 0,
        tile_size_bytes: 0,
        delta_q_y_dc: 0,
        delta_q_u_dc: 0,
        delta_q_u_ac: 0,
        delta_q_v_dc: 0,
        delta_q_v_ac: 0,
        num_planes: 3,
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
        seg_id_pre_skip: false,
        last_active_seg_id: 0,
        coded_lossless: false,
        all_lossless: false,
        lossless_array: [false; MAX_SEGMENTS],
        seg_q_m_level: [[0; MAX_SEGMENTS]; 3],
        current_q_index: 0,
        frame_restoration_type: [RESTORE_NONE; 3],
        enable_restoration: false,
        uses_lr: false,
        uses_chroma_lr: false,
        loop_restoration_size: [0; 3],
        tx_mode: ONLY_4X4,
        skip_mode_allowed: true,
        skip_mode_frame: [0; 2],
        enable_warped_motion: false,
        tg_start: 0,
        tg_end: 0,
        frame_size_override_flag: false,
        disable_frame_end_update_cdf: false,
    }}
}