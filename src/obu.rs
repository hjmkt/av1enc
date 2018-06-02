#![allow(dead_code)]
#![allow(non_camel_case_types)]

use sequence_header::SequenceHeader;
use frame_header::FrameHeader;
use frame::Frame;

struct SpatialLayer {
    max_width: u16,
    max_height: u16,
    ref_id: Option<u8>,
}

struct TemporalGroup {
    temporal_id: u8,
    temporal_switching_up_point_flag: u8,
    spatial_switching_up_point_flag: u8,
    ref_pic_diff: Vec<u8>,
}

enum ScalabilityModeIdc {
    SCALABILITY_L1T2,
    SCALABILITY_L1T3,
    SCALABILITY_L2T1,
    SCALABILITY_L2T2,
    SCALABILITY_L2T3,
    SCALABILITY_S2T1,
    SCALABILITY_S2T2,
    SCALABILITY_S2T3,
    SCALABILITY_L2T1h,
    SCALABILITY_L2T2h,
    SCALABILITY_L2T3h,
    SCALABILITY_S2T1h,
    SCALABILITY_S2T2h,
    SCALABILITY_S2T3h,
    SCALABILITY_SS,
}

pub struct ScalabilityStructure {
    spatial_layers: Vec<SpatialLayer>,
    temporal_group: Vec<TemporalGroup>,
}

struct Timecode {
    counting_type: u8,
    seconds: Option<u8>,
    minutes: Option<u8>,
    hours: Option<u8>,
    time_offset_value: Vec<u8>,
}

pub enum MetaData {
    METADATA_TYPE_ITUT_T35 { country_code: u16, payload_bytes: Vec<u16> },
    METADATA_TYPE_HDR_CLL { max_cll: u16, max_fall: u16 },
    METADATA_TYPE_HDR_MDCV {
        primary_chromaticity_x: [u16; 3],
        primary_chromaticity_y: [u16; 3],
        white_point_chromaticity_x: u16,
        white_point_chromaticity_y: u16,
        luminance_max: u32,
        luminance_min: u32,
    },
    METADATA_TYPE_SCALABILITY {
        mode_idc: u8,
        structure: Option<ScalabilityStructure>,
    },
    METADATA_TYPE_TIMECODE {
        counting_type: u8,
        seconds: Option<u8>,
        minutes: Option<u8>,
        hours: Option<u8>,
        time_offset_value: Vec<u8>,
    },
}

pub enum OBU<'a> {
    OBU_SEQUENCE_HEADER(SequenceHeader), // 1
    OBU_TEMPORAL_DELIMITER, // 2
    OBU_FRAME_HEADER(FrameHeader), // 3
    OBU_REDUNDANT_FRAME_HEADER, // 7
    OBU_TILE_GROUP, // 4
    OBU_METADATA(MetaData), // 5
    OBU_FRAME { // 6
        frame: Frame<'a>
    },
    OBU_PADDING(u32), // 15
}