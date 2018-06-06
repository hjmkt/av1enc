#![allow(dead_code)]
#![allow(non_camel_case_types)]

use sequence_header::SequenceHeader;
use frame::Frame;

pub struct SpatialLayerDimension {
    pub max_width: u16,
    pub max_height: u16,
}

pub struct SpatialLayerDescription {
    pub ref_id: u8,
}

pub struct TemporalGroup {
    pub temporal_id: u8,
    pub temporal_switching_up_point_flag: u8,
    pub spatial_switching_up_point_flag: u8,
    pub ref_pic_diffs: Vec<u8>,
}

#[derive(PartialEq, Copy, Clone)]
pub enum ScalabilityModeIdc {
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
    pub spatial_layer_dimensions: Vec<SpatialLayerDimension>,
    pub spatial_layer_descriptions: Vec<SpatialLayerDescription>,
    pub temporal_groups: Vec<TemporalGroup>,
}

pub enum MetaData {
    METADATA_TYPE_ITUT_T35 { country_code: u16, payload_bytes: Vec<u16> }, // 4
    METADATA_TYPE_HDR_CLL { max_cll: u16, max_fall: u16 }, // 1
    METADATA_TYPE_HDR_MDCV { // 2
        primary_chromaticity_x: [u16; 3],
        primary_chromaticity_y: [u16; 3],
        white_point_chromaticity_x: u16,
        white_point_chromaticity_y: u16,
        luminance_max: u32,
        luminance_min: u32,
    },
    METADATA_TYPE_SCALABILITY { // 3
        mode_idc: ScalabilityModeIdc,
        structure: Option<ScalabilityStructure>,
    },
    METADATA_TYPE_TIMECODE { // 5
        counting_type: u8,
        seconds: Option<u8>,
        minutes: Option<u8>,
        hours: Option<u8>,
        time_offset_length: u8,
        time_offset_value: u32,
    },
}

pub enum OBU {
    OBU_SEQUENCE_HEADER(SequenceHeader), // 1
    OBU_TEMPORAL_DELIMITER, // 2
    OBU_FRAME_HEADER(Frame), // 3
    OBU_REDUNDANT_FRAME_HEADER, // 7
    OBU_TILE_GROUP(Frame), // 4
    OBU_METADATA(MetaData), // 5
    OBU_FRAME { // 6
        frame: Frame
    },
    OBU_PADDING(u32), // 15
}