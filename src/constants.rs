#![allow(non_camel_case_types)]
#![allow(dead_code)]

pub enum ColorPrimaries {
    CP_BT_709 = 1,
    CP_UNSPECIFIED = 2,
    CP_BT_470_M = 4,
    CP_BT_470_B_G = 5,
    CP_BT_601 = 6,
    CP_SMPTE_240 = 7,
    CP_GENERIC_FILM = 8,
    CP_BT_2020 = 9,
    CP_XYZ = 10,
    CP_SMPTE_431 = 11,
    CP_SMPTE_432 = 12,
    CP_EBU_3213 = 22,
}

pub enum TransferCharacteristics {
    TC_RESERVED_0 = 0,
    TC_BT_709 = 1,
    TC_UNSPECIFIED = 2,
    TC_RESERVED_3 = 3,
    TC_BT_470_M = 4,
    TC_BT_470_B_G = 5,
    TC_BT_601 = 6,
    TC_SMPTE_240 = 7,
    TC_LINEAR = 8,
    TC_LOG_100 = 9,
    TC_LOG_100_SQRT10 = 10,
    TC_IEC_61966 = 11,
    TC_BT_1361 = 12,
    TC_SRGB = 13,
    TC_BT_2020_10_BIT = 14,
    TC_BT_2020_12_BIT = 15,
    TC_SMPTE_2084 = 16,
    TC_SMPTE_428 = 17,
    TC_HLG = 18,
}

pub enum MatrixCoefficients {
    MC_IDENTITY = 0,
    MC_BT_709 = 1,
    MC_UNSPECIFIED = 2,
    MC_RESERVED_3 = 3,
    MC_FCC = 4,
    MC_BT_470_B_G = 5,
    MC_BT_601 = 6,
    MC_SMPTE_240 = 7,
    MC_SMPTE_YCGCO = 8,
    MC_BT_2020_NCL = 9,
    MC_BT_2020_CL = 10,
    MC_SMPTE_2085 = 11,
    MC_CHROMAT_NCL = 12,
    MC_CHROMAT_CL = 13,
    MC_ICTCP = 14,
}

pub enum FrameType {
    KEY_FRAME = 0,
    INTER_FRAME = 1,
    INTRA_ONLY_FRAME = 2,
    SWITCH_FRAME = 3,
}

pub const REFS_PER_FRAME: u8 = 7;
pub const TOTAL_REFS_PER_FRAME: u8 = 8;
pub const BLOCK_SIZE_GROUPS: u8 = 4;
pub const BLOCK_SIZES: u8 = 22;
pub const BLOCK_INVALID: u8 = 22;
pub const MAX_SB_SIZE: u8 = 128;
pub const MI_SIZE: u8 = 4;
pub const MI_SIZE_LOG2: u8 = 2;
pub const MAX_TILE_WIDTH: u16 = 4096;
pub const MAX_TILE_AREA: u32 = 4096 * 2304;
pub const MAX_TILE_ROWS: u8 = 64;
pub const MAX_TILE_COLS: u8 = 64;
pub const INTRABC_DELAY_PIXELS: u16 = 256;
pub const INTRABC_DELAY_SB64: u8 = 4;
pub const NUM_REF_FRAMES: u8 = 8;
pub const IS_INTER_CONTEXTS: u8 = 4;
pub const REF_CONTEXTS: u8 = 3;
pub const MAX_SEGMENTS	: u8 = 8;
pub const SEGMENT_ID_CONTEXTS: u8 = 3;
pub const SEG_LVL_ALT_Q: u8 = 0;
pub const SEG_LVL_ALT_LF_Y_V: u8 = 1;
pub const SEG_LVL_REF_FRAME: u8 = 5;
pub const SEG_LVL_SKIP: u8 = 6;
pub const SEG_LVL_GLOBALMV: u8 = 7;
pub const SEG_LVL_MAX: u8 = 7;
pub const PLANE_TYPES: u8 = 2;
pub const TX_SIZE_CONTEXTS: u8 = 3;
pub const INTERP_FILTERS: u8 = 3;
pub const INTERP_FILTER_CONTEXTS: u8 = 16;
pub const SKIP_MODE_CONTEXTS: u8 = 3;
pub const SKIP_CONTEXTS: u8 = 3;
pub const PARTITION_CONTEXTS: u8 = 4;
pub const TX_SIZES: u8 = 5;
pub const TX_SIZES_ALL: u8 = 19;
pub const TX_MODES: u8 = 3;
pub const DCT_DCT: u8 = 0;
pub const ADST_DCT: u8 = 1;
pub const DCT_ADST: u8 = 2;
pub const ADST_ADST: u8 = 3;
pub const FLIPADST_DCT: u8 = 4;
pub const DCT_FLIPADST: u8 = 5;
pub const FLIPADST_FLIPADST: u8 = 6;
pub const ADST_FLIPADST: u8 = 7;
pub const FLIPADST_ADST: u8 = 8;
pub const IDTX: u8 = 9;
pub const V_DCT: u8 = 10;
pub const H_DCT: u8 = 11;
pub const V_ADST: u8 = 12;
pub const H_ADST: u8 = 13;
pub const V_FLIPADST: u8 = 14;
pub const H_FLIPADST: u8 = 15;
pub const TX_TYPES: u8 = 16;
pub const MB_MODE_COUNT: u8 = 17;
pub const INTRA_MODES: u8 = 13;
pub const UV_INTRA_MODES_CFL_NOT_ALLOWED: u8 = 13;
pub const UV_INTRA_MODES_CFL_ALLOWED: u8 = 14;
pub const COMPOUND_MODES: u8 = 8;
pub const COMPOUND_MODE_CONTEXTS: u8 = 8;
pub const COMP_NEWMV_CTXS: u8 = 5;
pub const NEW_MV_CONTEXTS: u8 = 6;
pub const ZERO_MV_CONTEXTS: u8 = 2;
pub const REF_MV_CONTEXTS: u8 = 6;
pub const DRL_MODE_CONTEXTS: u8 = 3;
pub const MV_CONTEXTS: u8 = 2;
pub const MV_INTRABC_CONTEXT: u8 = 1;
pub const MV_JOINTS: u8 = 4;
pub const MV_CLASSES: u8 = 11;
pub const CLASS0_SIZE: u8 = 2;
pub const MV_OFFSET_BITS: u8 = 10;
pub const MAX_LOOP_FILTER: u8 = 63;
pub const REF_SCALE_SHIFT: u8 = 14;
pub const SUBPEL_BITS: u8 = 4;
pub const SUBPEL_MASK: u8 = 15;
pub const SCALE_SUBPEL_BITS: u8 = 10;
pub const MV_BORDER: u8 = 128;
pub const PALETTE_COLOR_CONTEXTS: u8 = 5;
pub const PALETTE_MAX_COLOR_CONTEXT_HASH: u8 = 8;
pub const PALETTE_BLOCK_SIZE_CONTEXTS: u8 = 7;
pub const PALETTE_Y_MODE_CONTEXTS: u8 = 3;
pub const PALETTE_UV_MODE_CONTEXTS: u8 = 2;
pub const PALETTE_SIZES: u8 = 7;
pub const PALETTE_COLORS: u8 = 8;
pub const PALETTE_NUM_NEIGHBORS: u8 = 3;
pub const DELTA_Q_SMALL: u8 = 3;
pub const DELTA_LF_SMALL: u8 = 3;
pub const QM_TOTAL_SIZE: u16 = 3344;
pub const MAX_ANGLE_DELTA: u8 = 3;
pub const DIRECTIONAL_MODES: u8 = 8;
pub const ANGLE_STEP: u8 = 3;
pub const TX_SET_TYPES_INTRA: u8 = 3;
pub const TX_SET_TYPES_INTER: u8 = 4;
pub const WARPEDMODEL_PREC_BITS: u8 = 16;
pub const IDENTITY: u8 = 0;
pub const TRANSLATION: u8 = 1;
pub const ROTZOOM: u8 = 2;
pub const AFFINE: u8 = 3;
pub const GM_ABS_TRANS_BITS: u8 = 12;
pub const GM_ABS_TRANS_ONLY_BITS: u8 = 9;
pub const GM_ABS_ALPHA_BITS: u8 = 12;
pub const DIV_LUT_PREC_BITS: u8 = 14;
pub const DIV_LUT_BITS: u8 = 8;
pub const DIV_LUT_NUM: u16 = 257;
pub const MOTION_MODES: u8 = 3;
pub const SIMPLE: u8 = 0;
pub const OBMC: u8 = 1;
pub const LOCALWARP: u8 = 2;
pub const LEAST_SQUARES_SAMPLES_MAX: u8 = 8;
pub const LS_MV_MAX: u16 = 256;
pub const WARPEDMODEL_TRANS_CLAMP: u32 = 1<<23;
pub const WARPEDMODEL_NONDIAGAFFINE_CLAMP: u16 = 1<<13;
pub const WARPEDPIXEL_PREC_SHIFTS: u8 = 1<<6;
pub const WARPEDDIFF_PREC_BITS: u8 = 10;
pub const GM_ALPHA_PREC_BITS: u8 = 15;
pub const GM_TRANS_PREC_BITS: u8 = 6;
pub const GM_TRANS_ONLY_PREC_BITS: u8 = 3;
pub const INTERINTRA_MODES: u8 = 4;
pub const MASK_MASTER_SIZE: u8 = 64;
pub const SEGMENT_ID_PREDICTED_CONTEXTS: u8 = 3;
pub const FWD_REFS: u8 = 4;
pub const BWD_REFS: u8 = 3;
pub const SINGLE_REFS: u8 = 7;
pub const UNIDIR_COMP_REFS: u8 = 4;
pub const COMPOUND_TYPES: u8 = 2;
pub const CFL_JOINT_SIGNS: u8 = 8;
pub const CFL_ALPHABET_SIZE: u8 = 16;
pub const COMP_INTER_CONTEXTS: u8 = 5;
pub const COMP_REF_TYPE_CONTEXTS: u8 = 5;
pub const CFL_ALPHA_CONTEXTS: u8 = 6;
pub const INTRA_MODE_CONTEXTS: u8 = 5;
pub const COMP_GROUP_IDX_CONTEXTS: u8 = 6;
pub const COMPOUND_IDX_CONTEXTS: u8 = 6;
pub const INTRA_EDGE_KERNELS: u8 = 3;
pub const INTRA_EDGE_TAPS: u8 = 5;
pub const FRAME_LF_COUNT: u8 = 4;
pub const MAX_VARTX_DEPTH: u8 = 2;
pub const TXFM_PARTITION_CONTEXTS: u8 = 21;
pub const REF_CAT_LEVEL: u16 = 640;
pub const MAX_REF_MV_STACK_SIZE: u8 = 8;
pub const MFMV_STACK_SIZE: u8 = 3;
pub const MAX_TX_DEPTH: u8 = 2;
pub const WEDGE_TYPES: u8 = 16;
pub const FILTER_BITS: u8 = 7;
pub const WIENER_COEFFS: u8 = 3;
pub const SGRPROJ_PARAMS_BITS: u8 = 4;
pub const SGRPROJ_PRJ_SUBEXP_K: u8 = 4;
pub const SGRPROJ_PRJ_BITS: u8 = 7;
pub const SGRPROJ_RST_BITS: u8 = 4;
pub const SGRPROJ_MTABLE_BITS: u8 = 20;
pub const SGRPROJ_RECIP_BITS: u8 = 12;
pub const SGRPROJ_SGR_BITS: u8 = 8;
pub const EC_PROB_SHIFT: u8 = 6;
pub const EC_MIN_PROB: u8 = 4;
pub const SELECT_SCREEN_CONTENT_TOOLS: u8 = 2;
pub const SELECT_INTEGER_MV: u8 = 2;
pub const RESTORATION_TILESIZE_MAX: u16 = 256;
pub const MAX_FRAME_DISTANCE: u8 = 31;
pub const MAX_OFFSET_WIDTH: u8 = 8;
pub const MAX_OFFSET_HEIGHT: u8 = 0;
pub const WARP_PARAM_REDUCE_BITS: u8 = 6;
pub const NUM_BASE_LEVELS: u8 = 2;
pub const COEFF_BASE_RANGE: u8 = 12;
pub const BR_CDF_SIZE: u8 = 4;
pub const SIG_COEF_CONTEXTS_EOB: u8 = 4;
pub const SIG_COEF_CONTEXTS_2D: u8 = 26;
pub const SIG_COEF_CONTEXTS: u8 = 42;
pub const SIG_REF_DIFF_OFFSET_NUM: u8 = 5;
pub const SUPERRES_NUM: u8 = 8;
pub const SUPERRES_DENOM_MIN: u8 = 9;
pub const SUPERRES_DENOM_BITS: u8 = 3;
pub const SUPERRES_FILTER_BITS: u8 = 6;
pub const SUPERRES_FILTER_SHIFTS: u8 = 1<<SUPERRES_FILTER_BITS;
pub const SUPERRES_FILTER_TAPS: u8 = 8;
pub const SUPERRES_FILTER_OFFSET: u8 = 3;
pub const SUPERRES_SCALE_BITS: u8 = 14;
pub const SUPERRES_SCALE_MASK: u16 = (1<<14) - 1;
pub const SUPERRES_EXTRA_BITS: u8 = 8;
pub const TXB_SKIP_CONTEXTS: u8 = 13;
pub const EOB_COEF_CONTEXTS: u8 = 22;
pub const DC_SIGN_CONTEXTS: u8 = 3;
pub const LEVEL_CONTEXTS: u8 = 21;
pub const TX_CLASS_2D: u8 = 0;
pub const TX_CLASS_HORIZ: u8 = 1;
pub const TX_CLASS_VERT: u8 = 2;
pub const REFMVS_LIMIT: u16 = (1<<12) - 1;
pub const INTRA_FILTER_SCALE_BITS: u8 = 4;
pub const INTRA_FILTER_MODES: u8 = 5;
pub const COEFF_CDF_Q_CTXS: u8 = 4;
pub const PRIMARY_REF_NONE: u8 = 7;