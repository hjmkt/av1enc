#![allow(dead_code)]

use constants::*;
use frame_header::*;

pub struct Frame<'a> {
    width: u16,
    height: u16,
    frame_type: FrameType,
    showable_frame: bool,
    error_resilient_mode: bool,
    allow_screen_content_tools: bool,
    force_integer_mv: bool,
    id: u32,
    order_hint: u8,
    use_superres: bool,
    allow_intrabc: bool,
    ref_frames: [&'a Frame<'a>; NUM_REF_FRAMES as usize],
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
}