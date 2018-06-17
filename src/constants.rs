#![allow(non_camel_case_types)]
#![allow(non_upper_case_globals)]
#![allow(dead_code)]
#![allow(unused_imports)]

extern crate num;
extern crate enum_primitive;
use self::num::FromPrimitive;
use self::enum_primitive::*;
use self::TxType::*;
use constants::BlockSize::*;

#[derive(PartialEq, Copy, Clone)]
pub enum ColorPrimaries {
    CP_BT_709 = 1, // BT.709
    CP_UNSPECIFIED = 2, // Unspecified
    CP_BT_470_M = 4, // BT.470 System M (historical)
    CP_BT_470_B_G = 5, // BT.470 System B, G (historical)
    CP_BT_601 = 6, // BT.601
    CP_SMPTE_240 = 7, // SMPTE 240
    CP_GENERIC_FILM = 8, // Generic film (color filters using illuminant C)
    CP_BT_2020 = 9, // BT.2020, BT.2100
    CP_XYZ = 10, // SMPTE 428 (CIE 1921 XYZ)
    CP_SMPTE_431 = 11, // SMPTE RP 431-2
    CP_SMPTE_432 = 12, // SMPTE EG 432-1
    CP_EBU_3213 = 22, // EBU Tech. 3213-E
}

#[derive(PartialEq, Copy, Clone)]
pub enum TransferCharacteristics {
    TC_RESERVED_0 = 0, // For future use
    TC_BT_709 = 1, // BT.709
    TC_UNSPECIFIED = 2, // Unspecified
    TC_RESERVED_3 = 3, // For future use
    TC_BT_470_M = 4, // BT.470 System M (historical)
    TC_BT_470_B_G = 5, // BT.470 System B, G (historical)
    TC_BT_601 = 6, // BT.601
    TC_SMPTE_240 = 7, // SMPTE 240 M
    TC_LINEAR = 8, // Linear
    TC_LOG_100 = 9, // Logarithmic (100 : 1 range)
    TC_LOG_100_SQRT10 = 10, // Logarithmic (100 * Sqrt(10) : 1 range)
    TC_IEC_61966 = 11, // IEC 61966-2-4
    TC_BT_1361 = 12, // BT.1361
    TC_SRGB = 13, // sRGB or sYCC
    TC_BT_2020_10_BIT = 14, // BT.2020 10-bit systems
    TC_BT_2020_12_BIT = 15, // BT.2020 12-bit systems
    TC_SMPTE_2084 = 16, // SMPTE ST 2084, ITU BT.2100 PQ
    TC_SMPTE_428 = 17, // SMPTE ST 428
    TC_HLG = 18, // BT.2100 HLG, ARIB STD-B67
}

#[derive(PartialEq, Copy, Clone)]
pub enum MatrixCoefficients {
    MC_IDENTITY = 0, // Identity matrix
    MC_BT_709 = 1, // BT.709
    MC_UNSPECIFIED = 2, // Unspecified
    MC_RESERVED_3 = 3, // For future use
    MC_FCC = 4, // US FCC 73.628
    MC_BT_470_B_G = 5, // BT.470 System B, G (historical)
    MC_BT_601 = 6, // BT.601
    MC_SMPTE_240 = 7, // SMPTE 240 M
    MC_SMPTE_YCGCO = 8, // YCgCo
    MC_BT_2020_NCL = 9, // BT.2020 non-constant luminance, BT.2100 YCbCr
    MC_BT_2020_CL = 10, // BT.2020 constant luminance
    MC_SMPTE_2085 = 11, // SMPTE ST 2085 YDzDx
    MC_CHROMAT_NCL = 12, // Chromaticity-derived non-constant luminance
    MC_CHROMAT_CL = 13, // Chromaticity-derived constant luminance
    MC_ICTCP = 14, // BT.2100 ICtCp
}

#[derive(PartialEq, Copy, Clone)]
pub enum ChromaSamplePosition {
    CSP_UNKNOWN = 0,
    CSP_VERTICAL = 1,
    CSP_COLOCATED = 2,
    CSP_RESERVED = 3,
}

#[derive(PartialEq, Copy, Clone)]
pub enum FrameType {
    KEY_FRAME = 0,
    INTER_FRAME = 1,
    INTRA_ONLY_FRAME = 2,
    SWITCH_FRAME = 3,
}

#[derive(PartialEq, Copy, Clone)]
pub enum FrameRestorationType {
    RESTORE_NONE = 0, // lr_type=0
    RESTORE_WIENER = 1, // lr_type=2
    RESTORE_SGRPROJ = 2, // lr_type=3
    RESTORE_SWITCHABLE = 3, // lr_type=1
}

#[derive(PartialEq, Copy, Clone)]
pub enum Partition {
    PARTITION_NONE,
    PARTITION_HORZ,
    PARTITION_VERT,
    PARTITION_SPLIT,
    PARTITION_HORZ_A,
    PARTITION_HORZ_B,
    PARTITION_VERT_A,
    PARTITION_VERT_B,
    PARTITION_HORZ_4,
    PARTITION_VERT_4,
}

#[derive(PartialEq, Copy, Clone, PartialOrd)]
pub enum BlockSize {
    BLOCK_4X4,
    BLOCK_4X8,
    BLOCK_8X4,
    BLOCK_8X8,
    BLOCK_8X16,
    BLOCK_16X8,
    BLOCK_16X16,
    BLOCK_16X32,
    BLOCK_32X16,
    BLOCK_32X32,
    BLOCK_32X64,
    BLOCK_64X32,
    BLOCK_64X64,
    BLOCK_64X128,
    BLOCK_128X64,
    BLOCK_128X128,
    BLOCK_4X16,
    BLOCK_16X4,
    BLOCK_8X32,
    BLOCK_32X8,
    BLOCK_16X64,
    BLOCK_64X16,
    BLOCK_INVALID,
}

#[derive(PartialEq, Copy, Clone)]
pub enum YMode {
    DC_PRED,
    V_PRED,
    H_PRED,
    D45_PRED,
    D135_PRED,
    D113_PRED,
    D157_PRED,
    D203_PRED,
    D67_PRED,
    SMOOTH_PRED,
    SMOOTH_V_PRED,
    SMOOTH_H_PRED,
    PAETH_PRED, // __IntraFrameYMode
    UV_CFL_PRED, // __UVMode
    NEARESTMV, // ^^CompounMode
    NEARMV,
    GLOBALMV,
    NEWMV,
    NEAREST_NEARESTMV,
    NEAR_NEARMV,
    NEAREST_NEWMV,
    NEW_NEARESTMV,
    NEAR_NEWMV,
    NEW_NEARMV,
    GLOBAL_GLOBALMV,
    NEW_NEWMV,
}

enum_from_primitive! {
#[derive(PartialEq, Copy, Clone, PartialOrd)]
pub enum TxSize {
    TX_4X4 = 0,
    TX_8X8 = 1,
    TX_16X16 = 2,
    TX_32X32 = 3,
    TX_64X64 = 4,
    TX_4X8 = 5,
    TX_8X4 = 6,
    TX_8X16 = 7,
    TX_16X8 = 8,
    TX_16X32 = 9,
    TX_32X16 = 10,
    TX_32X64 = 11,
    TX_64X32 = 12,
    TX_4X16 = 13,
    TX_16X4 = 14,
    TX_8X32 = 15,
    TX_32X8 = 16,
    TX_16X64 = 17,
    TX_64X16 = 18,
    TX_INVALID = 19,
}
}

#[derive(PartialEq, Copy, Clone)]
pub enum TxMode {
    ONLY_4X4,
    TX_MODE_LARGEST,
    TX_MODE_SELECT,
}

#[derive(PartialEq, Copy, Clone)]
pub enum MVClass {
    MV_CLASS_0,
    MV_CLASS_1,
    MV_CLASS_2,
    MV_CLASS_3,
    MV_CLASS_4,
    MV_CLASS_5,
    MV_CLASS_6,
    MV_CLASS_7,
    MV_CLASS_8,
    MV_CLASS_9,
    MV_CLASS_10,
}

enum SignUV {
    CFL_SIGN_ZERO,
    CFL_SIGN_NEG,
    CFL_SIGN_POS,
}

const cfl_alpha_signs: [(SignUV, SignUV); 8] = [
    (SignUV::CFL_SIGN_ZERO, SignUV::CFL_SIGN_NEG),
    (SignUV::CFL_SIGN_ZERO, SignUV::CFL_SIGN_POS),
    (SignUV::CFL_SIGN_NEG, SignUV::CFL_SIGN_ZERO),
    (SignUV::CFL_SIGN_NEG, SignUV::CFL_SIGN_NEG),
    (SignUV::CFL_SIGN_NEG, SignUV::CFL_SIGN_POS),
    (SignUV::CFL_SIGN_POS, SignUV::CFL_SIGN_ZERO),
    (SignUV::CFL_SIGN_POS, SignUV::CFL_SIGN_NEG),
    (SignUV::CFL_SIGN_POS, SignUV::CFL_SIGN_POS),
];

#[derive(PartialEq, Copy, Clone, PartialOrd)]
pub enum InterpFilter {
    EIGHTTAP,
    EIGHTTAP_SMOOTH,
    EIGHTTAP_SHARP,
    BILINEAR,
    SWITCHABLE,
}

enum FilterIntraMode {
    FILTER_DC_PRED,
    FILTER_V_PRED,
    FILTER_H_PRED,
    FILTER_D157_PRED,
    FILTER_PAETH_PRED,
}

enum CompMode {
    SINGLE_REFERENCE,
    COMPOUND_REFERENCE,
}

enum CompRefType {
    UNIDIR_COMP_REFERENCE, // Both reference frames from the same group
    BIDIR_COMP_REFERENCE, // One from Group 1 and one from Group 2
}

#[derive(PartialEq, Copy, Clone, PartialOrd)]
pub enum RefFrame {
    NONE = -1, // ref[1]=NONE block uses single prediction
    INTRA_FRAME = 0, // ref[1]=INTRA_FRAME block uses interintra prediction
    LAST_FRAME = 1,
    LAST2_FRAME = 2,
    LAST3_FRAME = 3,
    GOLDEN_FRAME = 4,
    BWDREF_FRAME = 5,
    ALTREF2_FRAME = 6,
    ALTREF_FRAME = 7,
}

enum MotionMode {
    SIMPLE,
    OBMC,
    LOCALWARP,
}

enum InterIntraMode {
    II_DC_PRED,
    II_V_PRED,
    II_H_PRED,
    II_SMOOTH_PRED,
}

enum CompoundType {
    COMPOUND_WEDGE,
    COMPOUND_DIFFWTD,
    COMPOUND_AVERAGE,
    COMPOUND_INTRA,
    COMPOUND_DISTANCE,
}

enum MaskType {
    UNIFORM_45,
    UNIFORM_45_INV,
}

pub enum MvJoint {
    MV_JOINT_ZERO, // changes row = No, changes col = No
    MV_JOINT_HNZVZ, // changes row = No, changes col = Yes
    MV_JOINT_HZVNZ, // changes row = Yes, changes col = No
    MV_JOINT_HNZVNZ, // changes row = Yes, changes col = Yes
}

pub const REFS_PER_FRAME: usize = 7; // Number of reference frames that can be used for inter prediction
pub const TOTAL_REFS_PER_FRAME: usize = 8; // Number of reference frame types (including intra type)
pub const BLOCK_SIZE_GROUPS: usize = 4; // Number of contexts when decoding y_mode
pub const BLOCK_SIZES: usize = 22; // Number of different block sizes used
//pub const BLOCK_INVALID: usize = 22; // Sentinel value to mark partition choices that are not allowed
pub const MAX_SB_SIZE: usize = 128; // Maximum size of a superblock in luma samples
pub const MI_SIZE: usize = 4; // Smallest size of a mode info block in luma samples
pub const MI_SIZE_LOG2: usize = 2; // Base 2 logarithm of smallest size of a mode info block
pub const MAX_TILE_WIDTH: usize = 4096; // Maximum width of a tile in units of luma samples
pub const MAX_TILE_AREA: usize = 4096 * 2304; // Maximum area of a tile in units of luma samples
pub const MAX_TILE_ROWS: usize = 64; // Maximum number of tile rows
pub const MAX_TILE_COLS: usize = 64; // Maximum number of tile columns
pub const INTRABC_DELAY_PIXELS: usize = 256; // Number of horizontal luma samples before intra block copy can be used
pub const INTRABC_DELAY_SB64: usize = 4; // Number of 64 by 64 blocks before intra block copy can be used
pub const NUM_REF_FRAMES: usize = 8; // Number of frames that can be stored for future reference
pub const IS_INTER_CONTEXTS: usize = 4; // Number of contexts for is_inter
pub const REF_CONTEXTS: usize = 3; // Number of contexts for single_ref, comp_ref, comp_bwdref, uni_comp_ref, uni_comp_ref_p1 and uni_comp_ref_p2
pub const MAX_SEGMENTS	: usize = 8; // Number of segments allowed in segmentation map
pub const SEGMENT_ID_CONTEXTS: usize = 3; // Number of contexts for segment_id
pub const SEG_LVL_ALT_Q: usize = 0; // Index for quantizer segment feature
pub const SEG_LVL_ALT_LF_Y_V: usize = 1; // Index for vertical luma loop filter segment feature
pub const SEG_LVL_REF_FRAME: usize = 5; // Index for reference frame segment feature
pub const SEG_LVL_SKIP: usize = 6; // Index for skip segment feature
pub const SEG_LVL_GLOBALMV: usize = 7; // Index for global mv feature
pub const SEG_LVL_MAX: usize = 7; // Number of segment features
pub const PLANE_TYPES: usize = 2; // Number of different plane types (luma or chroma)
pub const TX_SIZE_CONTEXTS: usize = 3; // Number of contexts for transform size
pub const INTERP_FILTERS: usize = 3; // Number of values for interp_filter
pub const INTERP_FILTER_CONTEXTS: usize = 16; // Number of contexts for interp_filter
pub const SKIP_MODE_CONTEXTS: usize = 3; // Number of contexts for decoding skip_mode
pub const SKIP_CONTEXTS: usize = 3; // Number of contexts for decoding skip
pub const PARTITION_CONTEXTS: usize = 4; // Number of contexts when decoding partition
pub const TX_SIZES: usize = 5; // Number of square transform sizes
pub const TX_SIZES_ALL: usize = 19; // Number of transform sizes (including non-square sizes)
pub const TX_MODES: usize = 3; // Number of values for tx_mode
pub const TX_TYPES: usize = 16; // Number of inverse transform types
pub const MB_MODE_COUNT: usize = 17; // Number of values for YMode
pub const INTRA_MODES: usize = 13; // Number of values for y_mode
pub const UV_INTRA_MODES_CFL_NOT_ALLOWED: usize = 13; // Number of values for uv_mode when chroma from luma is not allowed
pub const UV_INTRA_MODES_CFL_ALLOWED: usize = 14; // Number of values for uv_mode when chroma from luma is allowed
pub const COMPOUND_MODES: usize = 8; // Number of values for compound_mode
pub const COMPOUND_MODE_CONTEXTS: usize = 8; // Number of contexts for compound_mode
pub const COMP_NEWMV_CTXS: usize = 5; // Number of new mv values used when constructing context for compound_mode
pub const NEW_MV_CONTEXTS: usize = 6; // Number of contexts for new_mv
pub const ZERO_MV_CONTEXTS: usize = 2; // Number of contexts for zero_mv
pub const REF_MV_CONTEXTS: usize = 6; // Number of contexts for ref_mv
pub const DRL_MODE_CONTEXTS: usize = 3; // Number of contexts for drl_mode
pub const MV_CONTEXTS: usize = 2; // Number of contexts for decoding motion vectors including one for intra block copy
pub const MV_INTRABC_CONTEXT: usize = 1; // Motion vector context used for intra block copy
pub const MV_JOINTS: usize = 4; // Number of values for mv_joint
pub const MV_CLASSES: usize = 11; // Number of values for mv_class
pub const CLASS0_SIZE: usize = 2; // Number of values for mv_class0_bit
pub const MV_OFFSET_BITS: usize = 10; // Maximum number of bits for decoding motion vectors
pub const MAX_LOOP_FILTER: usize = 63; // Maximum value used for loop filtering
pub const REF_SCALE_SHIFT: usize = 14; // Number of bits of precision when scaling reference frames
pub const SUBPEL_BITS: usize = 4; // Number of bits of precision when choosing an inter prediction filter kernel
pub const SUBPEL_MASK: usize = 15; // ( 1 << SUBPEL_BITS ) - 1
pub const SCALE_SUBPEL_BITS: usize = 10; // Number of bits of precision when computing inter prediction locations
pub const MV_BORDER: usize = 128; // Value used when clipping motion vectors
pub const PALETTE_COLOR_CONTEXTS: usize = 5; // Number of values for color contexts
pub const PALETTE_MAX_COLOR_CONTEXT_HASH: usize = 8; // Number of mappings between color context hash and color context
pub const PALETTE_BLOCK_SIZE_CONTEXTS: usize = 7; // Number of values for palette block size
pub const PALETTE_Y_MODE_CONTEXTS: usize = 3; // Number of values for palette Y plane mode contexts
pub const PALETTE_UV_MODE_CONTEXTS: usize = 2; // Number of values for palette U and V plane mode contexts
pub const PALETTE_SIZES: usize = 7; // Number of values for palette_size
pub const PALETTE_COLORS: usize = 8; // Number of values for palette_color
pub const PALETTE_NUM_NEIGHBORS: usize = 3; // Number of neighbors considered within palette computation
pub const DELTA_Q_SMALL: usize = 3; // Value indicating alternative encoding of quantizer index delta values
pub const DELTA_LF_SMALL: usize = 3; // Value indicating alternative encoding of loop filter delta values
pub const QM_TOTAL_SIZE: usize = 3344; // Number of values in the quantizer matrix
pub const MAX_ANGLE_DELTA: usize = 3; // Maximum magnitude of AngleDeltaY and AngleDeltaUV
pub const DIRECTIONAL_MODES: usize = 8; // Number of directional intra modes
pub const ANGLE_STEP: usize = 3; // Number of degrees of step per unit increase in AngleDeltaY or AngleDeltaUV.
pub const TX_SET_TYPES_INTRA: usize = 3; // Number of intra transform set types
pub const TX_SET_TYPES_INTER: usize = 4; // Number of inter transform set types
pub const WARPEDMODEL_PREC_BITS: usize = 16; // Internal precision of warped motion models
pub const IDENTITY: usize = 0; // Warp model is just an identity transform
pub const TRANSLATION: usize = 1; // Warp model is a pure translation
pub const ROTZOOM: usize = 2; // Warp model is a rotation + symmetric zoom + translation
pub const AFFINE: usize = 3; // Warp model is a general affine transform
pub const GM_ABS_TRANS_BITS: usize = 12; // Number of bits encoded for translational components of global motion models, if part of a ROTZOOM or AFFINE model
pub const GM_ABS_TRANS_ONLY_BITS: usize = 9; // Number of bits encoded for translational components of global motion models, if part of a TRANSLATION model
pub const GM_ABS_ALPHA_BITS: usize = 12; // Number of bits encoded for non-translational components of global motion models
pub const DIV_LUT_PREC_BITS: usize = 14; // Number of fractional bits of entries in divisor lookup table
pub const DIV_LUT_BITS: usize = 8; // Number of fractional bits for lookup in divisor lookup table
pub const DIV_LUT_NUM: usize = 257; // Number of entries in divisor lookup table
pub const MOTION_MODES: usize = 3; // Number of values for motion modes
pub const SIMPLE: usize = 0; // Use translation or global motion compensation
pub const OBMC: usize = 1; // Use overlapped block motion compensation
pub const LOCALWARP: usize = 2; // Use local warp motion compensation
pub const LEAST_SQUARES_SAMPLES_MAX: usize = 8; // Largest number of samples used when computing a local warp
pub const LS_MV_MAX: usize = 256; // Largest motion vector difference to include in local warp computation
pub const WARPEDMODEL_TRANS_CLAMP: usize = 1<<23; // Clamping value used for translation components of warp
pub const WARPEDMODEL_NONDIAGAFFINE_CLAMP: usize = 1<<13; // Clamping value used for matrix components of warp
pub const WARPEDPIXEL_PREC_SHIFTS: usize = 1<<6; // Number of phases used in warped filtering
pub const WARPEDDIFF_PREC_BITS: usize = 10; // Number of extra bits of precision in warped filtering
pub const GM_ALPHA_PREC_BITS: usize = 15; // Number of fractional bits for sending non-translational warp model coefficients
pub const GM_TRANS_PREC_BITS: usize = 6; // Number of fractional bits for sending translational warp model coefficients
pub const GM_TRANS_ONLY_PREC_BITS: usize = 3; // Number of fractional bits used for pure translational warps
pub const INTERINTRA_MODES: usize = 4; // Number of inter intra modes
pub const MASK_MASTER_SIZE: usize = 64; // Size of MasterMask array
pub const SEGMENT_ID_PREDICTED_CONTEXTS: usize = 3; // Number of contexts for segment_id_predicted
pub const FWD_REFS: usize = 4; // Number of contexts for is_inter
pub const BWD_REFS: usize = 3; // Number of contexts for skip
pub const SINGLE_REFS: usize = 7; // 	Number of syntax elements for single reference frames
pub const UNIDIR_COMP_REFS: usize = 4; // Number of syntax elements for unidirectional compound reference frames
pub const COMPOUND_TYPES: usize = 2; // Number of values for compound_type
pub const CFL_JOINT_SIGNS: usize = 8; // Number of values for cfl_alpha_signs
pub const CFL_ALPHABET_SIZE: usize = 16; // Number of values for cfl_alpha_u and cfl_alpha_v
pub const COMP_INTER_CONTEXTS: usize = 5; // Number of contexts for comp_mode
pub const COMP_REF_TYPE_CONTEXTS: usize = 5; // Number of contexts for comp_ref_type
pub const CFL_ALPHA_CONTEXTS: usize = 6; // Number of contexts for cfl_alpha_u and cfl_alpha_v
pub const INTRA_MODE_CONTEXTS: usize = 5; // Number of contexts for intra_frame_y_mode
pub const COMP_GROUP_IDX_CONTEXTS: usize = 6; // Number of contexts for comp_group_idx
pub const COMPOUND_IDX_CONTEXTS: usize = 6; // Number of contexts for compound_idx
pub const INTRA_EDGE_KERNELS: usize = 3; // Number of filter kernels for the intra edge filter
pub const INTRA_EDGE_TAPS: usize = 5; // Number of kernel taps for the intra edge filter
pub const FRAME_LF_COUNT: usize = 4; // Number of loop filter strength values
pub const MAX_VARTX_DEPTH: usize = 2; // Maximum depth for variable transform trees
pub const TXFM_PARTITION_CONTEXTS: usize = 21; // Number of contexts for txfm_split
pub const REF_CAT_LEVEL: usize = 640; // Bonus weight for close motion vectors
pub const MAX_REF_MV_STACK_SIZE: usize = 8; // Maximum number of motion vectors in the stack
pub const MFMV_STACK_SIZE: usize = 3; // Stack size for motion field motion vectors
pub const MAX_TX_DEPTH: usize = 2; // Number of contexts for tx_depth when the maximum transform size is 8x8
pub const WEDGE_TYPES: usize = 16; // Number of directions for the wedge mask process
pub const FILTER_BITS: usize = 7; // Number of bits used in Wiener filter coefficients
pub const WIENER_COEFFS: usize = 3; // Number of Wiener filter coefficients to read
pub const SGRPROJ_PARAMS_BITS: usize = 4; // Number of bits needed to specify self guided filter set
pub const SGRPROJ_PRJ_SUBEXP_K: usize = 4; // Controls how self guided deltas are read
pub const SGRPROJ_PRJ_BITS: usize = 7; // Precision bits during self guided restoration
pub const SGRPROJ_RST_BITS: usize = 4; // Restoration precision bits generated higher than source before projection
pub const SGRPROJ_MTABLE_BITS: usize = 20; // Precision of mtable division table
pub const SGRPROJ_RECIP_BITS: usize = 12; // Precision of division by n table
pub const SGRPROJ_SGR_BITS: usize = 8; // Internal precision bits for core selfguided_restoration
pub const EC_PROB_SHIFT: usize = 6; // Number of bits to reduce CDF precision during arithmetic coding
pub const EC_MIN_PROB: usize = 4; // Minimum probability assigned to each symbol during arithmetic coding
pub const SELECT_SCREEN_CONTENT_TOOLS: usize = 2; // Value that indicates the allow_screen_content_tools syntax element is coded
pub const SELECT_INTEGER_MV: usize = 2; // Value that indicates the force_integer_mv syntax element is coded
pub const RESTORATION_TILESIZE_MAX: usize = 256; // Maximum size of a loop restoration tile
pub const MAX_FRAME_DISTANCE: usize = 31; // Maximum distance when computing weighted prediction
pub const MAX_OFFSET_WIDTH: usize = 8; // Maximum horizontal offset of a projected motion vector
pub const MAX_OFFSET_HEIGHT: usize = 0; // Maximum vertical offset of a projected motion vector
pub const WARP_PARAM_REDUCE_BITS: usize = 6; // Rounding bitwidth for the parameters to the shear process
pub const NUM_BASE_LEVELS: usize = 2; // Number of quantizer base levels
pub const COEFF_BASE_RANGE: usize = 12; // The quantizer range above NUM_BASE_LEVELS above which the Exp-Golomb coding process is activated
pub const BR_CDF_SIZE: usize = 4; // Number of contexts for coeff_br
pub const SIG_COEF_CONTEXTS_EOB: usize = 4; // Number of contexts for coeff_base_eob
pub const SIG_COEF_CONTEXTS_2D: usize = 26; // Context offset for coeff_base for horizontal-only or vertical-only transforms.
pub const SIG_COEF_CONTEXTS: usize = 42; // Number of contexts for coeff_base
pub const SIG_REF_DIFF_OFFSET_NUM: usize = 5; // Maximum number of context samples to be used in determining the context index for coeff_base and coeff_base_eob.
pub const SUPERRES_NUM: usize = 8; // Numerator for upscaling ratio
pub const SUPERRES_DENOM_MIN: usize = 9; // Smallest denominator for upscaling ratio
pub const SUPERRES_DENOM_BITS: usize = 3; // Number of bits sent to specify denominator of upscaling ratio
pub const SUPERRES_FILTER_BITS: usize = 6; // Number of bits of fractional precision for upscaling filter selection
pub const SUPERRES_FILTER_SHIFTS: usize = 1<<SUPERRES_FILTER_BITS; // Number of phases of upscaling filters
pub const SUPERRES_FILTER_TAPS: usize = 8; // Number of taps of upscaling filters
pub const SUPERRES_FILTER_OFFSET: usize = 3; // Sample offset for upscaling filters
pub const SUPERRES_SCALE_BITS: usize = 14; // Number of fractional bits for computing position in upscaling
pub const SUPERRES_SCALE_MASK: usize = (1<<14) - 1; // Mask for computing position in upscaling
pub const SUPERRES_EXTRA_BITS: usize = 8; // Difference in precision between SUPERRES_SCALE_BITS and SUPERRES_FILTER_BITS
pub const TXB_SKIP_CONTEXTS: usize = 13; // Number of contexts for all_zero
pub const EOB_COEF_CONTEXTS: usize = 22; // Number of contexts for eob_extra
pub const DC_SIGN_CONTEXTS: usize = 3; // Number of contexts for dc_sign
pub const LEVEL_CONTEXTS: usize = 21; // Number of contexts for coeff_br
pub const TX_CLASS_2D: usize = 0; // Transform class for transform types performing non-identity transforms in both directions
pub const TX_CLASS_HORIZ: usize = 1; // Transform class for transforms performing only a horizontal non-identity transform
pub const TX_CLASS_VERT: usize = 2; // Transform class for transforms performing only a vertical non-identity transform
pub const REFMVS_LIMIT: usize = (1<<12) - 1; // Largest reference MV component that can be saved
pub const INTRA_FILTER_SCALE_BITS: usize = 4; // Scaling shift for intra filtering process
pub const INTRA_FILTER_MODES: usize = 5; // Number of types of intra filtering
pub const COEFF_CDF_Q_CTXS: usize = 4; // Number of selectable context types for the coeff( ) syntax structure
pub const PRIMARY_REF_NONE: usize = 7; // Value of primary_ref_frame indicating that there is no primary reference frame

pub const coeff_base_ctx_offset: [[[usize; 5]; 5]; TX_SIZES_ALL] = [
    [
        [ 0, 1, 6, 6, 0 ],
        [ 1, 6, 6, 21, 0 ],
        [ 6, 6, 21, 21, 0 ],
        [ 6, 21, 21, 21, 0 ],
        [ 0, 0, 0, 0, 0 ]
    ],
    [
        [ 0, 1, 6, 6, 21 ],
        [ 1, 6, 6, 21, 21 ],
        [ 6, 6, 21, 21, 21 ],
        [ 6, 21, 21, 21, 21 ],
        [ 21, 21, 21, 21, 21 ]
    ],
    [
        [ 0, 1, 6, 6, 21 ],
        [ 1, 6, 6, 21, 21 ],
        [ 6, 6, 21, 21, 21 ],
        [ 6, 21, 21, 21, 21 ],
        [ 21, 21, 21, 21, 21 ]
    ],
    [
        [ 0, 1, 6, 6, 21 ],
        [ 1, 6, 6, 21, 21 ],
        [ 6, 6, 21, 21, 21 ],
        [ 6, 21, 21, 21, 21 ],
        [ 21, 21, 21, 21, 21 ]
    ],
    [
        [ 0, 1, 6, 6, 21 ],
        [ 1, 6, 6, 21, 21 ],
        [ 6, 6, 21, 21, 21 ],
        [ 6, 21, 21, 21, 21 ],
        [ 21, 21, 21, 21, 21 ]
    ],
    [
        [ 0, 11, 11, 11, 0 ],
        [ 11, 11, 11, 11, 0 ],
        [ 6, 6, 21, 21, 0 ],
        [ 6, 21, 21, 21, 0 ],
        [ 21, 21, 21, 21, 0 ]
    ],
    [
        [ 0, 16, 6, 6, 21 ],
        [ 16, 16, 6, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 0, 0, 0, 0, 0 ]
    ],
    [
        [ 0, 11, 11, 11, 11 ],
        [ 11, 11, 11, 11, 11 ],
        [ 6, 6, 21, 21, 21 ],
        [ 6, 21, 21, 21, 21 ],
        [ 21, 21, 21, 21, 21 ]
    ],
    [
        [ 0, 16, 6, 6, 21 ],
        [ 16, 16, 6, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ]
    ],
    [
        [ 0, 11, 11, 11, 11 ],
        [ 11, 11, 11, 11, 11 ],
        [ 6, 6, 21, 21, 21 ],
        [ 6, 21, 21, 21, 21 ],
        [ 21, 21, 21, 21, 21 ]
    ],
    [
        [ 0, 16, 6, 6, 21 ],
        [ 16, 16, 6, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ]
    ],
    [
        [ 0, 11, 11, 11, 11 ],
        [ 11, 11, 11, 11, 11 ],
        [ 6, 6, 21, 21, 21 ],
        [ 6, 21, 21, 21, 21 ],
        [ 21, 21, 21, 21, 21 ]
    ],
    [
        [ 0, 16, 6, 6, 21 ],
        [ 16, 16, 6, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ]
    ],
    [
        [ 0, 11, 11, 11, 0 ],
        [ 11, 11, 11, 11, 0 ],
        [ 6, 6, 21, 21, 0 ],
        [ 6, 21, 21, 21, 0 ],
        [ 21, 21, 21, 21, 0 ]
    ],
    [
        [ 0, 16, 6, 6, 21 ],
        [ 16, 16, 6, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 0, 0, 0, 0, 0 ]
    ],
    [
        [ 0, 11, 11, 11, 11 ],
        [ 11, 11, 11, 11, 11 ],
        [ 6, 6, 21, 21, 21 ],
        [ 6, 21, 21, 21, 21 ],
        [ 21, 21, 21, 21, 21 ]
    ],
    [
        [ 0, 16, 6, 6, 21 ],
        [ 16, 16, 6, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ]
    ],
    [
        [ 0, 11, 11, 11, 11 ],
        [ 11, 11, 11, 11, 11 ],
        [ 6, 6, 21, 21, 21 ],
        [ 6, 21, 21, 21, 21 ],
        [ 21, 21, 21, 21, 21 ]
    ],
    [
        [ 0, 16, 6, 6, 21 ],
        [ 16, 16, 6, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ],
        [ 16, 16, 21, 21, 21 ]
    ]
];

pub const coeff_base_pos_ctx_offset: [usize; 3] = [
    SIG_COEF_CONTEXTS_2D,
    SIG_COEF_CONTEXTS_2D + 5,
    SIG_COEF_CONTEXTS_2D + 10
];

const subpel_filters: [[[isize; 8]; 16]; 6] = [
    [[0, 0, 0, 128, 0, 0, 0, 0],
     [0, 2, -6, 126, 8, -2, 0, 0],
     [0, 2, -10, 122, 18, -4, 0, 0],
     [0, 2, -12, 116, 28, -8, 2, 0],
     [0, 2, -14, 110, 38, -10, 2, 0],
     [0, 2, -14, 102, 48, -12, 2, 0],
     [0, 2, -16, 94, 58, -12, 2, 0],
     [0, 2, -14, 84, 66, -12, 2, 0],
     [0, 2, -14, 76, 76, -14, 2, 0],
     [0, 2, -12, 66, 84, -14, 2, 0],
     [0, 2, -12, 58, 94, -16, 2, 0],
     [0, 2, -12, 48, 102, -14, 2, 0],
     [0, 2, -10, 38, 110, -14, 2, 0],
     [0, 2, -8, 28, 116, -12, 2, 0],
     [0, 0, -4, 18, 122, -10, 2, 0],
     [0, 0, -2, 8, 126, -6, 2, 0]],
    [[0, 0, 0, 128, 0, 0, 0, 0],
     [0, 2, 28, 62, 34, 2, 0, 0],
     [0, 0, 26, 62, 36, 4, 0, 0],
     [0, 0, 22, 62, 40, 4, 0, 0],
     [0, 0, 20, 60, 42, 6, 0, 0],
     [0, 0, 18, 58, 44, 8, 0, 0],
     [0, 0, 16, 56, 46, 10, 0, 0],
     [0, -2, 16, 54, 48, 12, 0, 0],
     [0, -2, 14, 52, 52, 14, -2, 0],
     [0, 0, 12, 48, 54, 16, -2, 0],
     [0, 0, 10, 46, 56, 16, 0, 0],
     [0, 0, 8, 44, 58, 18, 0, 0],
     [0, 0, 6, 42, 60, 20, 0, 0],
     [0, 0, 4, 40, 62, 22, 0, 0],
     [0, 0, 4, 36, 62, 26, 0, 0],
     [0, 0, 2, 34, 62, 28, 2, 0]],
    [[0, 0, 0, 128, 0, 0, 0, 0],
     [-2, 2, -6, 126, 8, -2, 2, 0],
     [-2, 6, -12, 124, 16, -6, 4, -2],
     [-2, 8, -18, 120, 26, -10, 6, -2],
     [-4, 10, -22, 116, 38, -14, 6, -2],
     [-4, 10, -22, 108, 48, -18, 8, -2],
     [-4, 10, -24, 100, 60, -20, 8, -2],
     [-4, 10, -24, 90, 70, -22, 10, -2],
     [-4, 12, -24, 80, 80, -24, 12, -4],
     [-2, 10, -22, 70, 90, -24, 10, -4],
     [-2, 8, -20, 60, 100, -24, 10, -4],
     [-2, 8, -18, 48, 108, -22, 10, -4],
     [-2, 6, -14, 38, 116, -22, 10, -4],
     [-2, 6, -10, 26, 120, -18, 8, -2],
     [-2, 4, -6, 16, 124, -12, 6, -2],
     [0, 2, -2, 8, 126, -6, 2, -2]],
    [[0, 0, 0, 128, 0, 0, 0, 0],
     [0, 0, 0, 120, 8, 0, 0, 0],
     [0, 0, 0, 112, 16, 0, 0, 0],
     [0, 0, 0, 104, 24, 0, 0, 0],
     [0, 0, 0, 96, 32, 0, 0, 0],
     [0, 0, 0, 88, 40, 0, 0, 0],
     [0, 0, 0, 80, 48, 0, 0, 0],
     [0, 0, 0, 72, 56, 0, 0, 0],
     [0, 0, 0, 64, 64, 0, 0, 0],
     [0, 0, 0, 56, 72, 0, 0, 0],
     [0, 0, 0, 48, 80, 0, 0, 0],
     [0, 0, 0, 40, 88, 0, 0, 0],
     [0, 0, 0, 32, 96, 0, 0, 0],
     [0, 0, 0, 24, 104, 0, 0, 0],
     [0, 0, 0, 16, 112, 0, 0, 0],
     [0, 0, 0, 8, 120, 0, 0, 0]],
    [[0, 0, 0, 128, 0, 0, 0, 0],
     [0, 0, -4, 126, 8, -2, 0, 0],
     [0, 0, -8, 122, 18, -4, 0, 0],
     [0, 0, -10, 116, 28, -6, 0, 0],
     [0, 0, -12, 110, 38, -8, 0, 0],
     [0, 0, -12, 102, 48, -10, 0, 0],
     [0, 0, -14, 94, 58, -10, 0, 0],
     [0, 0, -12, 84, 66, -10, 0, 0],
     [0, 0, -12, 76, 76, -12, 0, 0],
     [0, 0, -10, 66, 84, -12, 0, 0],
     [0, 0, -10, 58, 94, -14, 0, 0],
     [0, 0, -10, 48, 102, -12, 0, 0],
     [0, 0, -8, 38, 110, -12, 0, 0],
     [0, 0, -6, 28, 116, -10, 0, 0],
     [0, 0, -4, 18, 122, -8, 0, 0],
     [0, 0, -2, 8, 126, -4, 0, 0]],
    [[0, 0, 0, 128, 0, 0, 0, 0],
     [0, 0, 30, 62, 34, 2, 0, 0],
     [0, 0, 26, 62, 36, 4, 0, 0],
     [0, 0, 22, 62, 40, 4, 0, 0],
     [0, 0, 20, 60, 42, 6, 0, 0],
     [0, 0, 18, 58, 44, 8, 0, 0],
     [0, 0, 16, 56, 46, 10, 0, 0],
     [0, 0, 14, 54, 48, 12, 0, 0],
     [0, 0, 12, 52, 52, 12, 0, 0],
     [0, 0, 12, 48, 54, 14, 0, 0],
     [0, 0, 10, 46, 56, 16, 0, 0],
     [0, 0, 8, 44, 58, 18, 0, 0],
     [0, 0, 6, 42, 60, 20, 0, 0],
     [0, 0, 4, 40, 62, 22, 0, 0],
     [0, 0, 4, 36, 62, 26, 0, 0],
     [0, 0, 2, 34, 62, 30, 0, 0]]
];

const warped_filters: [[isize; 8]; WARPEDPIXEL_PREC_SHIFTS*3+1] = [
    [ 0,   0, 127,   1,   0, 0, 0, 0 ], [ 0, - 1, 127,   2,   0, 0, 0, 0 ],
    [ 1, - 3, 127,   4, - 1, 0, 0, 0 ], [ 1, - 4, 126,   6, - 2, 1, 0, 0 ],
    [ 1, - 5, 126,   8, - 3, 1, 0, 0 ], [ 1, - 6, 125,  11, - 4, 1, 0, 0 ],
    [ 1, - 7, 124,  13, - 4, 1, 0, 0 ], [ 2, - 8, 123,  15, - 5, 1, 0, 0 ],
    [ 2, - 9, 122,  18, - 6, 1, 0, 0 ], [ 2, -10, 121,  20, - 6, 1, 0, 0 ],
    [ 2, -11, 120,  22, - 7, 2, 0, 0 ], [ 2, -12, 119,  25, - 8, 2, 0, 0 ],
    [ 3, -13, 117,  27, - 8, 2, 0, 0 ], [ 3, -13, 116,  29, - 9, 2, 0, 0 ],
    [ 3, -14, 114,  32, -10, 3, 0, 0 ], [ 3, -15, 113,  35, -10, 2, 0, 0 ],
    [ 3, -15, 111,  37, -11, 3, 0, 0 ], [ 3, -16, 109,  40, -11, 3, 0, 0 ],
    [ 3, -16, 108,  42, -12, 3, 0, 0 ], [ 4, -17, 106,  45, -13, 3, 0, 0 ],
    [ 4, -17, 104,  47, -13, 3, 0, 0 ], [ 4, -17, 102,  50, -14, 3, 0, 0 ],
    [ 4, -17, 100,  52, -14, 3, 0, 0 ], [ 4, -18,  98,  55, -15, 4, 0, 0 ],
    [ 4, -18,  96,  58, -15, 3, 0, 0 ], [ 4, -18,  94,  60, -16, 4, 0, 0 ],
    [ 4, -18,  91,  63, -16, 4, 0, 0 ], [ 4, -18,  89,  65, -16, 4, 0, 0 ],
    [ 4, -18,  87,  68, -17, 4, 0, 0 ], [ 4, -18,  85,  70, -17, 4, 0, 0 ],
    [ 4, -18,  82,  73, -17, 4, 0, 0 ], [ 4, -18,  80,  75, -17, 4, 0, 0 ],
    [ 4, -18,  78,  78, -18, 4, 0, 0 ], [ 4, -17,  75,  80, -18, 4, 0, 0 ],
    [ 4, -17,  73,  82, -18, 4, 0, 0 ], [ 4, -17,  70,  85, -18, 4, 0, 0 ],
    [ 4, -17,  68,  87, -18, 4, 0, 0 ], [ 4, -16,  65,  89, -18, 4, 0, 0 ],
    [ 4, -16,  63,  91, -18, 4, 0, 0 ], [ 4, -16,  60,  94, -18, 4, 0, 0 ],
    [ 3, -15,  58,  96, -18, 4, 0, 0 ], [ 4, -15,  55,  98, -18, 4, 0, 0 ],
    [ 3, -14,  52, 100, -17, 4, 0, 0 ], [ 3, -14,  50, 102, -17, 4, 0, 0 ],
    [ 3, -13,  47, 104, -17, 4, 0, 0 ], [ 3, -13,  45, 106, -17, 4, 0, 0 ],
    [ 3, -12,  42, 108, -16, 3, 0, 0 ], [ 3, -11,  40, 109, -16, 3, 0, 0 ],
    [ 3, -11,  37, 111, -15, 3, 0, 0 ], [ 2, -10,  35, 113, -15, 3, 0, 0 ],
    [ 3, -10,  32, 114, -14, 3, 0, 0 ], [ 2, - 9,  29, 116, -13, 3, 0, 0 ],
    [ 2, - 8,  27, 117, -13, 3, 0, 0 ], [ 2, - 8,  25, 119, -12, 2, 0, 0 ],
    [ 2, - 7,  22, 120, -11, 2, 0, 0 ], [ 1, - 6,  20, 121, -10, 2, 0, 0 ],
    [ 1, - 6,  18, 122, - 9, 2, 0, 0 ], [ 1, - 5,  15, 123, - 8, 2, 0, 0 ],
    [ 1, - 4,  13, 124, - 7, 1, 0, 0 ], [ 1, - 4,  11, 125, - 6, 1, 0, 0 ],
    [ 1, - 3,   8, 126, - 5, 1, 0, 0 ], [ 1, - 2,   6, 126, - 4, 1, 0, 0 ],
    [ 0, - 1,   4, 127, - 3, 1, 0, 0 ], [ 0,   0,   2, 127, - 1, 0, 0, 0 ],

    [ 0,  0,   0, 127,   1,   0,  0,  0], [ 0,  0,  -1, 127,   2,   0,  0,  0],
    [ 0,  1,  -3, 127,   4,  -2,  1,  0], [ 0,  1,  -5, 127,   6,  -2,  1,  0],
    [ 0,  2,  -6, 126,   8,  -3,  1,  0], [-1,  2,  -7, 126,  11,  -4,  2, -1],
    [-1,  3,  -8, 125,  13,  -5,  2, -1], [-1,  3, -10, 124,  16,  -6,  3, -1],
    [-1,  4, -11, 123,  18,  -7,  3, -1], [-1,  4, -12, 122,  20,  -7,  3, -1],
    [-1,  4, -13, 121,  23,  -8,  3, -1], [-2,  5, -14, 120,  25,  -9,  4, -1],
    [-1,  5, -15, 119,  27, -10,  4, -1], [-1,  5, -16, 118,  30, -11,  4, -1],
    [-2,  6, -17, 116,  33, -12,  5, -1], [-2,  6, -17, 114,  35, -12,  5, -1],
    [-2,  6, -18, 113,  38, -13,  5, -1], [-2,  7, -19, 111,  41, -14,  6, -2],
    [-2,  7, -19, 110,  43, -15,  6, -2], [-2,  7, -20, 108,  46, -15,  6, -2],
    [-2,  7, -20, 106,  49, -16,  6, -2], [-2,  7, -21, 104,  51, -16,  7, -2],
    [-2,  7, -21, 102,  54, -17,  7, -2], [-2,  8, -21, 100,  56, -18,  7, -2],
    [-2,  8, -22,  98,  59, -18,  7, -2], [-2,  8, -22,  96,  62, -19,  7, -2],
    [-2,  8, -22,  94,  64, -19,  7, -2], [-2,  8, -22,  91,  67, -20,  8, -2],
    [-2,  8, -22,  89,  69, -20,  8, -2], [-2,  8, -22,  87,  72, -21,  8, -2],
    [-2,  8, -21,  84,  74, -21,  8, -2], [-2,  8, -22,  82,  77, -21,  8, -2],
    [-2,  8, -21,  79,  79, -21,  8, -2], [-2,  8, -21,  77,  82, -22,  8, -2],
    [-2,  8, -21,  74,  84, -21,  8, -2], [-2,  8, -21,  72,  87, -22,  8, -2],
    [-2,  8, -20,  69,  89, -22,  8, -2], [-2,  8, -20,  67,  91, -22,  8, -2],
    [-2,  7, -19,  64,  94, -22,  8, -2], [-2,  7, -19,  62,  96, -22,  8, -2],
    [-2,  7, -18,  59,  98, -22,  8, -2], [-2,  7, -18,  56, 100, -21,  8, -2],
    [-2,  7, -17,  54, 102, -21,  7, -2], [-2,  7, -16,  51, 104, -21,  7, -2],
    [-2,  6, -16,  49, 106, -20,  7, -2], [-2,  6, -15,  46, 108, -20,  7, -2],
    [-2,  6, -15,  43, 110, -19,  7, -2], [-2,  6, -14,  41, 111, -19,  7, -2],
    [-1,  5, -13,  38, 113, -18,  6, -2], [-1,  5, -12,  35, 114, -17,  6, -2],
    [-1,  5, -12,  33, 116, -17,  6, -2], [-1,  4, -11,  30, 118, -16,  5, -1],
    [-1,  4, -10,  27, 119, -15,  5, -1], [-1,  4,  -9,  25, 120, -14,  5, -2],
    [-1,  3,  -8,  23, 121, -13,  4, -1], [-1,  3,  -7,  20, 122, -12,  4, -1],
    [-1,  3,  -7,  18, 123, -11,  4, -1], [-1,  3,  -6,  16, 124, -10,  3, -1],
    [-1,  2,  -5,  13, 125,  -8,  3, -1], [-1,  2,  -4,  11, 126,  -7,  2, -1],
    [ 0,  1,  -3,   8, 126,  -6,  2,  0], [ 0,  1,  -2,   6, 127,  -5,  1,  0],
    [ 0,  1,  -2,   4, 127,  -3,  1,  0], [ 0,  0,   0,   2, 127,  -1,  0,  0],

    [ 0, 0, 0,   1, 127,   0,   0, 0 ], [ 0, 0, 0, - 1, 127,   2,   0, 0 ],
    [ 0, 0, 1, - 3, 127,   4, - 1, 0 ], [ 0, 0, 1, - 4, 126,   6, - 2, 1 ],
    [ 0, 0, 1, - 5, 126,   8, - 3, 1 ], [ 0, 0, 1, - 6, 125,  11, - 4, 1 ],
    [ 0, 0, 1, - 7, 124,  13, - 4, 1 ], [ 0, 0, 2, - 8, 123,  15, - 5, 1 ],
    [ 0, 0, 2, - 9, 122,  18, - 6, 1 ], [ 0, 0, 2, -10, 121,  20, - 6, 1 ],
    [ 0, 0, 2, -11, 120,  22, - 7, 2 ], [ 0, 0, 2, -12, 119,  25, - 8, 2 ],
    [ 0, 0, 3, -13, 117,  27, - 8, 2 ], [ 0, 0, 3, -13, 116,  29, - 9, 2 ],
    [ 0, 0, 3, -14, 114,  32, -10, 3 ], [ 0, 0, 3, -15, 113,  35, -10, 2 ],
    [ 0, 0, 3, -15, 111,  37, -11, 3 ], [ 0, 0, 3, -16, 109,  40, -11, 3 ],
    [ 0, 0, 3, -16, 108,  42, -12, 3 ], [ 0, 0, 4, -17, 106,  45, -13, 3 ],
    [ 0, 0, 4, -17, 104,  47, -13, 3 ], [ 0, 0, 4, -17, 102,  50, -14, 3 ],
    [ 0, 0, 4, -17, 100,  52, -14, 3 ], [ 0, 0, 4, -18,  98,  55, -15, 4 ],
    [ 0, 0, 4, -18,  96,  58, -15, 3 ], [ 0, 0, 4, -18,  94,  60, -16, 4 ],
    [ 0, 0, 4, -18,  91,  63, -16, 4 ], [ 0, 0, 4, -18,  89,  65, -16, 4 ],
    [ 0, 0, 4, -18,  87,  68, -17, 4 ], [ 0, 0, 4, -18,  85,  70, -17, 4 ],
    [ 0, 0, 4, -18,  82,  73, -17, 4 ], [ 0, 0, 4, -18,  80,  75, -17, 4 ],
    [ 0, 0, 4, -18,  78,  78, -18, 4 ], [ 0, 0, 4, -17,  75,  80, -18, 4 ],
    [ 0, 0, 4, -17,  73,  82, -18, 4 ], [ 0, 0, 4, -17,  70,  85, -18, 4 ],
    [ 0, 0, 4, -17,  68,  87, -18, 4 ], [ 0, 0, 4, -16,  65,  89, -18, 4 ],
    [ 0, 0, 4, -16,  63,  91, -18, 4 ], [ 0, 0, 4, -16,  60,  94, -18, 4 ],
    [ 0, 0, 3, -15,  58,  96, -18, 4 ], [ 0, 0, 4, -15,  55,  98, -18, 4 ],
    [ 0, 0, 3, -14,  52, 100, -17, 4 ], [ 0, 0, 3, -14,  50, 102, -17, 4 ],
    [ 0, 0, 3, -13,  47, 104, -17, 4 ], [ 0, 0, 3, -13,  45, 106, -17, 4 ],
    [ 0, 0, 3, -12,  42, 108, -16, 3 ], [ 0, 0, 3, -11,  40, 109, -16, 3 ],
    [ 0, 0, 3, -11,  37, 111, -15, 3 ], [ 0, 0, 2, -10,  35, 113, -15, 3 ],
    [ 0, 0, 3, -10,  32, 114, -14, 3 ], [ 0, 0, 2, - 9,  29, 116, -13, 3 ],
    [ 0, 0, 2, - 8,  27, 117, -13, 3 ], [ 0, 0, 2, - 8,  25, 119, -12, 2 ],
    [ 0, 0, 2, - 7,  22, 120, -11, 2 ], [ 0, 0, 1, - 6,  20, 121, -10, 2 ],
    [ 0, 0, 1, - 6,  18, 122, - 9, 2 ], [ 0, 0, 1, - 5,  15, 123, - 8, 2 ],
    [ 0, 0, 1, - 4,  13, 124, - 7, 1 ], [ 0, 0, 1, - 4,  11, 125, - 6, 1 ],
    [ 0, 0, 1, - 3,   8, 126, - 5, 1 ], [ 0, 0, 1, - 2,   6, 126, - 4, 1 ],
    [ 0, 0, 0, - 1,   4, 127, - 3, 1 ], [ 0, 0, 0,   0,   2, 127, - 1, 0 ],

    [ 0, 0, 0,   0,   2, 127, - 1, 0 ]
];

const div_lut: [usize; DIV_LUT_NUM] = [
    16384, 16320, 16257, 16194, 16132, 16070, 16009, 15948, 15888, 15828, 15768,
    15709, 15650, 15592, 15534, 15477, 15420, 15364, 15308, 15252, 15197, 15142,
    15087, 15033, 14980, 14926, 14873, 14821, 14769, 14717, 14665, 14614, 14564,
    14513, 14463, 14413, 14364, 14315, 14266, 14218, 14170, 14122, 14075, 14028,
    13981, 13935, 13888, 13843, 13797, 13752, 13707, 13662, 13618, 13574, 13530,
    13487, 13443, 13400, 13358, 13315, 13273, 13231, 13190, 13148, 13107, 13066,
    13026, 12985, 12945, 12906, 12866, 12827, 12788, 12749, 12710, 12672, 12633,
    12596, 12558, 12520, 12483, 12446, 12409, 12373, 12336, 12300, 12264, 12228,
    12193, 12157, 12122, 12087, 12053, 12018, 11984, 11950, 11916, 11882, 11848,
    11815, 11782, 11749, 11716, 11683, 11651, 11619, 11586, 11555, 11523, 11491,
    11460, 11429, 11398, 11367, 11336, 11305, 11275, 11245, 11215, 11185, 11155,
    11125, 11096, 11067, 11038, 11009, 10980, 10951, 10923, 10894, 10866, 10838,
    10810, 10782, 10755, 10727, 10700, 10673, 10645, 10618, 10592, 10565, 10538,
    10512, 10486, 10460, 10434, 10408, 10382, 10356, 10331, 10305, 10280, 10255,
    10230, 10205, 10180, 10156, 10131, 10107, 10082, 10058, 10034, 10010, 9986,
    9963,  9939,  9916,  9892,  9869,  9846,  9823,  9800,  9777,  9754,  9732,
    9709,  9687,  9664,  9642,  9620,  9598,  9576,  9554,  9533,  9511,  9489,
    9468,  9447,  9425,  9404,  9383,  9362,  9341,  9321,  9300,  9279,  9259,
    9239,  9218,  9198,  9178,  9158,  9138,  9118,  9098,  9079,  9059,  9039,
    9020,  9001,  8981,  8962,  8943,  8924,  8905,  8886,  8867,  8849,  8830,
    8812,  8793,  8775,  8756,  8738,  8720,  8702,  8684,  8666,  8648,  8630,
    8613,  8595,  8577,  8560,  8542,  8525,  8508,  8490,  8473,  8456,  8439,
    8422,  8405,  8389,  8372,  8355,  8339,  8322,  8306,  8289,  8273,  8257,
    8240,  8224,  8208,  8192
];

const obmc_mask_1: [usize; 1] = [ 64 ];
const obmc_mask_2: [usize; 2] = [ 45, 64 ];
const obmc_mask_4: [usize; 4] = [ 39, 50, 59, 64 ];
const obmc_mask_8: [usize; 8] = [ 36, 42, 48, 53, 57, 61, 64, 64 ];
const obmc_mask_16: [usize; 16] = [
    34, 37, 40, 43, 46, 49, 52, 54,
    56, 58, 60, 61, 64, 64, 64, 64
];
const obmc_mask_32: [usize; 32] = [
    33, 35, 36, 38, 40, 41, 43, 44,
    45, 47, 48, 50, 51, 52, 53, 55,
    56, 57, 58, 59, 60, 60, 61, 62,
    64, 64, 64, 64, 64, 64, 64, 64
];

const wedge_master_oblique_odd: [usize; MASK_MASTER_SIZE] = [
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  2,  6,  18,
    37, 53, 60, 63, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64
];

const wedge_master_oblique_even: [usize; MASK_MASTER_SIZE] = [
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  4,  11, 27,
    46, 58, 62, 63, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64
];

const wedge_master_vertical: [usize; MASK_MASTER_SIZE] = [
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  7,  21,
    43, 57, 62, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64
];

enum WedgeDirection {
    WEDGE_HORIZONTAL,
    WEDGE_VERTICAL,
    WEDGE_OBLIQUE27,
    WEDGE_OBLIQUE63,
    WEDGE_OBLIQUE117,
    WEDGE_OBLIQUE153,
}

const ii_weights_1d: [usize; MAX_SB_SIZE] = [
    60, 58, 56, 54, 52, 50, 48, 47, 45, 44, 42, 41, 39, 38, 37, 35, 34, 33, 32,
    31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 22, 21, 20, 19, 19, 18, 18, 17, 16,
    16, 15, 15, 14, 14, 13, 13, 12, 12, 12, 11, 11, 10, 10, 10,  9,  9,  9,  8,
    8,  8,  8,  7,  7,  7,  7,  6,  6,  6,  6,  6,  5,  5,  5,  5,  5,  4,  4,
    4,  4,  4,  4,  4,  4,  3,  3,  3,  3,  3,  3,  3,  3,  3,  2,  2,  2,  2,
    2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  1,  1,  1,  1,  1,  1,  1,  1,
    1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1
];

const quant_dist_weight: [[usize; 2]; 4] = [
    [ 2, 3 ], [ 2, 5 ], [ 2, 7 ], [ 1, MAX_FRAME_DISTANCE ]
];

const quant_dist_lookup: [[usize; 2]; 4] = [
    [ 9, 7 ], [ 11, 5 ], [ 12, 4 ], [ 13, 3 ]
];

const cos128_lookup: [usize; 65] = [
    4096, 4095, 4091, 4085, 4076, 4065, 4052, 4036,
    4017, 3996, 3973, 3948, 3920, 3889, 3857, 3822,
    3784, 3745, 3703, 3659, 3612, 3564, 3513, 3461,
    3406, 3349, 3290, 3229, 3166, 3102, 3035, 2967,
    2896, 2824, 2751, 2675, 2598, 2520, 2440, 2359,
    2276, 2191, 2106, 2019, 1931, 1842, 1751, 1660,
    1567, 1474, 1380, 1285, 1189, 1092, 995, 897,
    799, 700, 601, 501, 401, 301, 201, 101, 0
];

const transform_row_shift: [usize; TX_SIZES_ALL] = [
    0, 1, 2, 2, 2, 0, 0, 1, 1,
    1, 1, 1, 1, 1, 1, 2, 2, 2, 2
];

const cdef_uv_dir: [[[usize; 8]; 2]; 2] = [
    [
        [0, 1, 2, 3, 4, 5, 6, 7],
        [1, 2, 2, 2, 3, 4, 6, 0]
    ],
    [
        [7, 0, 2, 4, 5, 6, 6, 6],
        [0, 1, 2, 3, 4, 5, 6, 7]
    ]
];

const cdef_pri_taps: [[usize; 2]; 2] = [
    [ 4, 2 ], [ 3, 3 ]
];

const cdef_sec_taps: [[usize; 2]; 2] = [
    [ 2, 1 ], [ 2, 1 ]
];

const cdef_direction: [[[isize; 2]; 2]; 8] = [
    [ [ -1, 1 ], [ -2,  2 ] ],
    [ [  0, 1 ], [ -1,  2 ] ],
    [ [  0, 1 ], [  0,  2 ] ],
    [ [  0, 1 ], [  1,  2 ] ],
    [ [  1, 1 ], [  2,  2 ] ],
    [ [  1, 0 ], [  2,  1 ] ],
    [ [  1, 0 ], [  2,  0 ] ],
    [ [  1, 0 ], [  2, -1 ] ]
];

const upscale_filter: [[isize; SUPERRES_FILTER_TAPS]; SUPERRES_FILTER_SHIFTS] = [
    [ 0, 0, 0, 128, 0, 0, 0, 0 ],        [ 0, 0, -1, 128, 2, -1, 0, 0 ],
    [ 0, 1, -3, 127, 4, -2, 1, 0 ],      [ 0, 1, -4, 127, 6, -3, 1, 0 ],
    [ 0, 2, -6, 126, 8, -3, 1, 0 ],      [ 0, 2, -7, 125, 11, -4, 1, 0 ],
    [ -1, 2, -8, 125, 13, -5, 2, 0 ],    [ -1, 3, -9, 124, 15, -6, 2, 0 ],
    [ -1, 3, -10, 123, 18, -6, 2, -1 ],  [ -1, 3, -11, 122, 20, -7, 3, -1 ],
    [ -1, 4, -12, 121, 22, -8, 3, -1 ],  [ -1, 4, -13, 120, 25, -9, 3, -1 ],
    [ -1, 4, -14, 118, 28, -9, 3, -1 ],  [ -1, 4, -15, 117, 30, -10, 4, -1 ],
    [ -1, 5, -16, 116, 32, -11, 4, -1 ], [ -1, 5, -16, 114, 35, -12, 4, -1 ],
    [ -1, 5, -17, 112, 38, -12, 4, -1 ], [ -1, 5, -18, 111, 40, -13, 5, -1 ],
    [ -1, 5, -18, 109, 43, -14, 5, -1 ], [ -1, 6, -19, 107, 45, -14, 5, -1 ],
    [ -1, 6, -19, 105, 48, -15, 5, -1 ], [ -1, 6, -19, 103, 51, -16, 5, -1 ],
    [ -1, 6, -20, 101, 53, -16, 6, -1 ], [ -1, 6, -20, 99, 56, -17, 6, -1 ],
    [ -1, 6, -20, 97, 58, -17, 6, -1 ],  [ -1, 6, -20, 95, 61, -18, 6, -1 ],
    [ -2, 7, -20, 93, 64, -18, 6, -2 ],  [ -2, 7, -20, 91, 66, -19, 6, -1 ],
    [ -2, 7, -20, 88, 69, -19, 6, -1 ],  [ -2, 7, -20, 86, 71, -19, 6, -1 ],
    [ -2, 7, -20, 84, 74, -20, 7, -2 ],  [ -2, 7, -20, 81, 76, -20, 7, -1 ],
    [ -2, 7, -20, 79, 79, -20, 7, -2 ],  [ -1, 7, -20, 76, 81, -20, 7, -2 ],
    [ -2, 7, -20, 74, 84, -20, 7, -2 ],  [ -1, 6, -19, 71, 86, -20, 7, -2 ],
    [ -1, 6, -19, 69, 88, -20, 7, -2 ],  [ -1, 6, -19, 66, 91, -20, 7, -2 ],
    [ -2, 6, -18, 64, 93, -20, 7, -2 ],  [ -1, 6, -18, 61, 95, -20, 6, -1 ],
    [ -1, 6, -17, 58, 97, -20, 6, -1 ],  [ -1, 6, -17, 56, 99, -20, 6, -1 ],
    [ -1, 6, -16, 53, 101, -20, 6, -1 ], [ -1, 5, -16, 51, 103, -19, 6, -1 ],
    [ -1, 5, -15, 48, 105, -19, 6, -1 ], [ -1, 5, -14, 45, 107, -19, 6, -1 ],
    [ -1, 5, -14, 43, 109, -18, 5, -1 ], [ -1, 5, -13, 40, 111, -18, 5, -1 ],
    [ -1, 4, -12, 38, 112, -17, 5, -1 ], [ -1, 4, -12, 35, 114, -16, 5, -1 ],
    [ -1, 4, -11, 32, 116, -16, 5, -1 ], [ -1, 4, -10, 30, 117, -15, 4, -1 ],
    [ -1, 3, -9, 28, 118, -14, 4, -1 ],  [ -1, 3, -9, 25, 120, -13, 4, -1 ],
    [ -1, 3, -8, 22, 121, -12, 4, -1 ],  [ -1, 3, -7, 20, 122, -11, 3, -1 ],
    [ -1, 2, -6, 18, 123, -10, 3, -1 ],  [ 0, 2, -6, 15, 124, -9, 3, -1 ],
    [ 0, 2, -5, 13, 125, -8, 2, -1 ],    [ 0, 1, -4, 11, 125, -7, 2, 0 ],
    [ 0, 1, -3, 8, 126, -6, 2, 0 ],      [ 0, 1, -3, 6, 127, -4, 1, 0 ],
    [ 0, 1, -2, 4, 127, -3, 1, 0 ],      [ 0, 0, -1, 2, 128, -1, 0, 0 ],
];

pub const sgr_params: [[usize; 4]; 1<<SGRPROJ_PARAMS_BITS] = [
    [ 2, 12, 1, 4 ],  [ 2, 15, 1, 6 ],  [ 2, 18, 1, 8 ],  [ 2, 21, 1, 9 ],
    [ 2, 24, 1, 10 ], [ 2, 29, 1, 11 ], [ 2, 36, 1, 12 ], [ 2, 45, 1, 13 ],
    [ 2, 56, 1, 14 ], [ 2, 68, 1, 15 ], [ 0, 0, 1, 5 ],   [ 0, 0, 1, 8 ],
    [ 0, 0, 1, 11 ],  [ 0, 0, 1, 14 ],  [ 2, 30, 0, 0 ],  [ 2, 75, 0, 0 ]
];

pub const sgrproj_xqd_min: [isize; 2] = [
    -96, -32
];

pub const sgrproj_xqd_max: [isize; 2] = [
    31, 95
];

// FIXME
pub const segmentation_feature_bits: [usize; SEG_LVL_MAX+1] = [
    8 , 6 , 6 , 6 , 6 , 3 , 0 , 0
];

pub const segmentation_feature_signed: [bool; SEG_LVL_MAX+1] = [
    true , true , true , true , true , false , false , false
];

pub const segmentation_feature_max: [usize; SEG_LVL_MAX+1] = [
    255, MAX_LOOP_FILTER, MAX_LOOP_FILTER, MAX_LOOP_FILTER, MAX_LOOP_FILTER, 7, 0, 0
];

pub const remap_lr_type: [FrameRestorationType; 4] = [
    FrameRestorationType::RESTORE_NONE,
    FrameRestorationType::RESTORE_SWITCHABLE,
    FrameRestorationType::RESTORE_WIENER,
    FrameRestorationType::RESTORE_SGRPROJ,
];

pub const wiener_taps_mid: [isize; 3] = [ 3,  -7,  15 ];
pub const sgrproj_xqd_mid: [isize; 2] = [ -32, 31 ];

const max_tx_depth: [usize; BLOCK_SIZES] = [
    0, 1, 1, 1,
    2, 2, 2, 3,
    3, 3, 4, 4,
    4, 4, 4, 4,
    2, 2, 3, 3,
    4, 4
];

const tx_mode_to_biggest_tx_size: [TxSize; TX_MODES] = [
    TxSize::TX_4X4,
    TxSize::TX_64X64,
    TxSize::TX_64X64
];

pub const subsampled_size: [[[BlockSize; 2]; 2]; BLOCK_SIZES] = [
    [ [ BlockSize::BLOCK_4X4,    BlockSize::BLOCK_4X4],      [BlockSize::BLOCK_4X4,     BlockSize::BLOCK_4X4] ],
    [ [ BlockSize::BLOCK_4X8,    BlockSize::BLOCK_4X4],      [BlockSize::BLOCK_INVALID, BlockSize::BLOCK_4X4] ],
    [ [ BlockSize::BLOCK_8X4,    BlockSize::BLOCK_INVALID],  [BlockSize::BLOCK_4X4,     BlockSize::BLOCK_4X4] ],
    [ [ BlockSize::BLOCK_8X8,    BlockSize::BLOCK_8X4],      [BlockSize::BLOCK_4X8,     BlockSize::BLOCK_4X4] ],
    [ [BlockSize::BLOCK_8X16,    BlockSize::BLOCK_8X8],      [BlockSize::BLOCK_INVALID, BlockSize::BLOCK_4X8] ],
    [ [BlockSize::BLOCK_16X8,    BlockSize::BLOCK_INVALID],  [BlockSize::BLOCK_8X8,     BlockSize::BLOCK_8X4] ],
    [ [BlockSize::BLOCK_16X16,   BlockSize::BLOCK_16X8],     [BlockSize::BLOCK_8X16,    BlockSize::BLOCK_8X8] ],
    [ [BlockSize::BLOCK_16X32,   BlockSize::BLOCK_16X16],    [BlockSize::BLOCK_INVALID, BlockSize::BLOCK_8X16] ],
    [ [BlockSize::BLOCK_32X16,   BlockSize::BLOCK_INVALID],  [BlockSize::BLOCK_16X16,   BlockSize::BLOCK_16X8] ],
    [ [BlockSize::BLOCK_32X32,   BlockSize::BLOCK_32X16],    [BlockSize::BLOCK_16X32,   BlockSize::BLOCK_16X16] ],
    [ [BlockSize::BLOCK_32X64,   BlockSize::BLOCK_32X32],    [BlockSize::BLOCK_INVALID, BlockSize::BLOCK_16X32] ],
    [ [BlockSize::BLOCK_64X32,   BlockSize::BLOCK_INVALID],  [BlockSize::BLOCK_32X32,   BlockSize::BLOCK_32X16] ],
    [ [BlockSize::BLOCK_64X64,   BlockSize::BLOCK_64X32],    [BlockSize::BLOCK_32X64,   BlockSize::BLOCK_32X32] ],
    [ [BlockSize::BLOCK_64X128,  BlockSize::BLOCK_64X64],    [BlockSize::BLOCK_INVALID, BlockSize::BLOCK_32X64] ],
    [ [BlockSize::BLOCK_128X64,  BlockSize::BLOCK_INVALID],  [BlockSize::BLOCK_64X64,   BlockSize::BLOCK_64X32] ],
    [ [BlockSize::BLOCK_128X128, BlockSize::BLOCK_128X64],   [BlockSize::BLOCK_64X128,  BlockSize::BLOCK_64X64] ],
    [ [BlockSize::BLOCK_4X16,    BlockSize::BLOCK_4X8],      [BlockSize::BLOCK_INVALID, BlockSize::BLOCK_4X8] ],
    [ [BlockSize::BLOCK_16X4,    BlockSize::BLOCK_INVALID],  [BlockSize::BLOCK_8X4,     BlockSize::BLOCK_8X4] ],
    [ [BlockSize::BLOCK_8X32,    BlockSize::BLOCK_8X16],     [BlockSize::BLOCK_INVALID, BlockSize::BLOCK_4X16] ],
    [ [BlockSize::BLOCK_32X8,    BlockSize::BLOCK_INVALID],  [BlockSize::BLOCK_16X8,    BlockSize::BLOCK_16X4] ],
    [ [BlockSize::BLOCK_16X64,   BlockSize::BLOCK_16X32],    [BlockSize::BLOCK_INVALID, BlockSize::BLOCK_8X32] ],
    [ [BlockSize::BLOCK_64X16,   BlockSize::BLOCK_INVALID],  [BlockSize::BLOCK_32X16,   BlockSize::BLOCK_32X8] ],
];

pub const tx_type_in_set_intra: [[bool; TX_TYPES]; TX_SET_TYPES_INTRA] = [
    [
        true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    ],
    [
        true, true, true, true, false, false, false, false, false, true, true, true, false, false, false, false,
    ],
    [
        true, true, true, true, false, false, false, false, false, true, false, false, false, false, false, false,
    ]
];

pub const tx_type_in_set_inter: [[bool; TX_TYPES]; TX_SET_TYPES_INTER] = [
    [
        true, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    ],
    [
        true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,
    ],
    [
        true, true, true, true, true, true, true, true, true, true, true, true, false, false, false, false,
    ],
    [
        true, false, false, false, false, false, false, false, false, true, false, false, false, false, false, false,
    ]
];

pub const tx_type_intra_inv_set1: [TxType; 7] = [
    IDTX, DCT_DCT, V_DCT, H_DCT, ADST_ADST, ADST_DCT, DCT_ADST
];

pub const tx_type_intra_inv_set2: [TxType; 5] = [
    IDTX, DCT_DCT, ADST_ADST, ADST_DCT, DCT_ADST
];

pub const tx_type_inter_inv_set1: [TxType; 16] = [
    IDTX, V_DCT, H_DCT, V_ADST, H_ADST, V_FLIPADST, H_FLIPADST,
    DCT_DCT, ADST_DCT, DCT_ADST, FLIPADST_DCT, DCT_FLIPADST, ADST_ADST,
    FLIPADST_FLIPADST, ADST_FLIPADST, FLIPADST_ADST
];

pub const tx_type_inter_inv_set12: [TxType; 12] = [
    IDTX, V_DCT, H_DCT, DCT_DCT, ADST_DCT, DCT_ADST, FLIPADST_DCT,
    DCT_FLIPADST, ADST_ADST, FLIPADST_FLIPADST, ADST_FLIPADST,
    FLIPADST_ADST
];

pub const tx_type_inter_inv_set3: [TxType; 2] = [
    IDTX, DCT_DCT
];

pub const wiener_taps_min: [isize; 3] = [ -5, -23, -17 ];
pub const wiener_taps_max: [isize; 3] = [ 10, 8, 46 ];
pub const wiener_taps_k: [isize; 3] = [ 1, 2, 3 ];

pub const sgrproj_taps_min: [isize; 2] = [ -96, -32 ];
pub const sgrproj_taps_max: [isize; 2] = [ 31, 95 ];

pub const intra_mode_context: [usize; INTRA_MODES] = [
    0, 1, 2, 3, 4, 4, 4, 4, 3, 0, 1, 2, 0
];

pub const size_group: [usize; BLOCK_SIZES] = [
    0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3,
    3, 3, 3, 3, 0, 0, 1, 1, 2, 2,
];

pub const block_width: [usize; BLOCK_SIZES] = [
    4, 4, 8, 8, 8, 16, 16, 16, 32, 32, 32, 64, 64, 64, 128, 128, 4, 16, 8, 32, 16, 64
];

pub const block_height: [usize; BLOCK_SIZES] = [
    4, 8, 4, 8, 16, 8, 16, 32, 16, 32, 64, 32, 64, 128, 64, 128, 16, 4, 32, 8, 64, 16
];

pub const mi_width_log2: [usize; BLOCK_SIZES] = [
    2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 7, 2, 4, 3, 5, 4, 6
];

pub const mi_height_log2: [usize; BLOCK_SIZES] = [
    2, 3, 2, 3, 4, 3, 4, 5, 4, 5, 6, 5, 6, 7, 6, 7, 4, 2, 5, 3, 6, 4
];

pub const num_4x4_blocks_wide: [usize; BLOCK_SIZES] = [
    1, 1, 2, 2, 2, 3, 3, 3, 8, 8, 8,
    16, 16, 16, 32, 32, 1, 4, 2, 8, 4, 16
];

pub const num_4x4_blocks_high: [usize; BLOCK_SIZES] = [
    1, 2, 1, 2, 4, 2, 4, 8, 4, 8, 16,
    8, 16, 32, 16, 32, 4, 1, 8, 2, 16, 4
];

pub const tx_width: [usize; TX_SIZES_ALL] = [
    4, 8, 16, 32, 64, 4, 8, 8, 16, 16, 32, 32, 64, 4, 16, 8, 32, 16, 64
];

pub const tx_width_log2: [usize; TX_SIZES_ALL] = [
    2, 3, 4, 5, 6, 2, 3, 3, 4, 4, 5, 5, 6, 2, 4, 3, 5, 4, 6
];

pub const tx_height: [usize; TX_SIZES_ALL] = [
    4, 8, 16, 32, 64, 8, 4, 16, 8, 32, 16, 64, 32, 16, 4, 32, 8, 64, 16
];

pub const tx_height_log2: [usize; TX_SIZES_ALL] = [
    2, 3, 4, 5, 6, 3, 2, 4, 3, 5, 4, 6, 5, 4, 2, 5, 3, 6, 4
];

#[derive(PartialEq, Copy, Clone)]
pub struct ModeInfo {
    pub mi_size: BlockSize,
    pub tx_size: TxSize,
    pub tx_type: TxType,
    pub is_inter: bool,
    pub skip: bool,
    pub intra: bool,
    pub single: bool,
    pub seg_id_predicted: usize,
    pub skip_mode: usize,
    pub ref_frame: [RefFrame; 2],
    pub interp_filter: [InterpFilter; 2],
    pub level_context: [usize; 3],
    pub dc_context: [usize; 3],
    pub palette_size: [usize; 3],
    pub comp_group_idx: usize,
    pub compound_idx: usize,
}

use self::TxSize::*;
pub const tx_size_sqr: [TxSize; TX_SIZES_ALL+1] = [
    TX_4X4,
    TX_8X8,
    TX_16X16,
    TX_32X32,
    TX_64X64,
    TX_4X4,
    TX_4X4,
    TX_8X8,
    TX_8X8,
    TX_16X16,
    TX_16X16,
    TX_32X32,
    TX_32X32,
    TX_4X4,
    TX_4X4,
    TX_8X8,
    TX_8X8,
    TX_16X16,
    TX_16X16,
    TX_INVALID,
];

pub const tx_size_sqr_up: [TxSize; TX_SIZES_ALL+1] = [
    TX_4X4,
    TX_8X8,
    TX_16X16,
    TX_32X32,
    TX_64X64,
    TX_8X8,
    TX_8X8,
    TX_16X16,
    TX_16X16,
    TX_32X32,
    TX_32X32,
    TX_64X64,
    TX_64X64,
    TX_16X16,
    TX_16X16,
    TX_32X32,
    TX_32X32,
    TX_64X64,
    TX_64X64,
    TX_INVALID,
];

pub const TX_SET_DCTONLY: usize = 0;
pub const TX_SET_INTRA_1: usize = 1;
pub const TX_SET_INTRA_2: usize = 2;
pub const TX_SET_INTER_1: usize = 0;
pub const TX_SET_INTER_2: usize = 1;
pub const TX_SET_INTER_3: usize = 2;

#[derive(PartialEq, Copy, Clone)]
pub enum TxType {
    DCT_DCT = 0, // Inverse transform rows with DCT and columns with DCT
    ADST_DCT = 1, // Inverse transform rows with DCT and columns with ADST
    DCT_ADST = 2, // Inverse transform rows with ADST and columns with DCT
    ADST_ADST = 3, // 	Inverse transform rows with ADST and columns with ADST
    FLIPADST_DCT = 4, // Inverse transform rows with FLIPADST and columns with DCT
    DCT_FLIPADST = 5, // Inverse transform rows with DCT and columns with FLIPADST
    FLIPADST_FLIPADST = 6, // Inverse transform rows with FLIPADST and columns with FLIPADST
    ADST_FLIPADST = 7, // Inverse transform rows with ADST and columns with FLIPADST
    FLIPADST_ADST = 8, // Inverse transform rows with FLIPADST and columns with ADST
    IDTX = 9, // Inverse transform rows with identity and columns with identity
    V_DCT = 10, // Inverse transform rows with identity and columns with DCT
    H_DCT = 11, // Inverse transform rows with DCT and columns with identity
    V_ADST = 12, // Inverse transform rows with identity and columns with ADST
    H_ADST = 13, // Inverse transform rows with ADST and columns with identity
    V_FLIPADST = 14, // Inverse transform rows with identity and columns with FLIPADST
    H_FLIPADST = 15, // Inverse transform rows with FLIPADST and columns with identity
}

// FIXME not in spec
pub const mode_to_txfm: [TxType; MB_MODE_COUNT] = [
    DCT_DCT, // DC
    ADST_DCT, // V
    DCT_ADST, // H
    DCT_DCT, // D45
    ADST_ADST, // D135
    ADST_DCT, // D117
    DCT_ADST, // D153
    DCT_ADST, // D207
    ADST_DCT, // D63
    ADST_ADST, // SMOOTH
    ADST_DCT,  // SMOOTH_V
    DCT_ADST,  // SMOOTH_H
    ADST_ADST, // TM
    DCT_DCT, // NEARESTMV
    DCT_DCT, // NEARMV
    DCT_DCT, // GLOBALMV
    DCT_DCT // NEWMV
];

pub const adjusted_tx_size: [TxSize; TX_SIZES_ALL] = [
    TX_4X4,
    TX_8X8,
    TX_16X16,
    TX_32X32,
    TX_32X32,
    TX_4X8,
    TX_8X4,
    TX_8X16,
    TX_16X8,
    TX_16X32,
    TX_32X16,
    TX_32X32,
    TX_32X32,
    TX_4X16,
    TX_16X4,
    TX_8X32,
    TX_32X8,
    TX_16X32,
    TX_32X16,
];

pub const sig_ref_diff_offset: [[[usize; 2]; SIG_REF_DIFF_OFFSET_NUM]; 3] = [
    [
        [0, 1], [1, 0], [1, 1], [0, 2], [2, 0]
    ],
    [
        [0, 1], [1, 0], [0, 2], [0, 3], [0, 4]
    ],
    [
        [0, 1], [1, 0], [2, 0], [3, 0], [4, 0]
    ],
];

pub const mag_ref_offset_with_tx_class: [[[usize; 2]; 3]; 3] = [
    [[0, 1], [1, 0], [1, 1]],
    [[0, 1], [1, 0], [0, 2]],
    [[0, 1], [1, 0], [2, 0]],
];

pub const palette_color_context: [isize; PALETTE_MAX_COLOR_CONTEXT_HASH+1] = [
    -1, -1, 0, -1, -1, 4, 3, 2, 1
];

use self::YMode::*;
pub const filter_intra_mode_to_intra_dir: [YMode; INTRA_FILTER_MODES] = [
    DC_PRED, V_PRED, H_PRED, D157_PRED, DC_PRED
];

pub const partition_subsize: [[BlockSize; BLOCK_SIZES]; 10] = [
    [
                                    BLOCK_4X4,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X8,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X16,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X32,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X64,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X128,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID
    ], [
                                    BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X8,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X32,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID
    ], [
                                    BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X16,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X64,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID
    ], [
                                    BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X4,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X8,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X16,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X32,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X64,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID
    ], [
                                    BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X8,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X32,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID
    ], [
                                    BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X4,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X8,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X16,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X32,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_128X64,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID
    ], [
                                    BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X16,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X64,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID
    ], [
                                    BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X8,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X16,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X32,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X64,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X128,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID
    ], [
                                    BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X4,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_32X8,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_64X16,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID
    ], [
                                    BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_4X16,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_8X32,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_16X64,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID,
        BLOCK_INVALID, BLOCK_INVALID, BLOCK_INVALID
    ]
];