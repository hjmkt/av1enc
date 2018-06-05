use obu::*;
use binary_writer::*;
use bool_coder::*;
use obu::OBU;
use obu::OBU::*;
use util::*;
use constants::*;
use sequence_header::*;
use constants::ColorPrimaries::*;
use constants::TransferCharacteristics::*;
use constants::MatrixCoefficients::*;
use constants::ChromaSamplePosition::*;
use obu::MetaData::*;
use obu::ScalabilityModeIdc::*;

pub struct OBUEncoder<'a, 'b> where 'a: 'b {
    writer: &'b mut BinaryWriter<'a>,
    coder: &'b mut BoolCoder,
    //operating_point_idc: usize, // FIXME
    order_hint_bits: usize,
    bit_depth: usize,
    num_planes: usize,
    seen_frame_header: bool,
}

impl<'a, 'b> OBUEncoder<'a, 'b> where 'a: 'b {
    pub fn new(writer: &'b mut BinaryWriter<'a>, coder: &'b mut BoolCoder) -> OBUEncoder<'a, 'b> {
        OBUEncoder {
            writer: writer,
            coder: coder,
            //operating_point_idc: 0,
            order_hint_bits: 0,
            bit_depth: 8,
            num_planes: 3,
            seen_frame_header: false,
        }
    }

    pub fn encode_obu_header(&mut self, obu: &mut OBU, mut out_bits: &mut Vec<u8>) {
        let obu_forbidden_bit = 0;
        self.coder.push_bit(&mut out_bits, obu_forbidden_bit);
        let obu_type = match obu {
            &mut OBU_SEQUENCE_HEADER(_) => 1,
            &mut OBU_TEMPORAL_DELIMITER => 2,
            &mut OBU_FRAME_HEADER(_) => 3,
            &mut OBU_REDUNDANT_FRAME_HEADER => 7,
            &mut OBU_TILE_GROUP => 4,
            &mut OBU_METADATA(_) => 5,
            &mut OBU_FRAME { frame: _ } => 6,
            &mut OBU_PADDING(_) => 15,
        };
        self.coder.push_bits_with_size(&mut out_bits, obu_type, 4);
        let obu_extension_flag = 0;
        self.coder.push_bit(&mut out_bits, obu_extension_flag);
        let obu_has_size_field = 0;
        self.coder.push_bit(&mut out_bits, obu_has_size_field);
        let obu_reserved_1bit = 0;
        self.coder.push_bit(&mut out_bits, obu_reserved_1bit);
        if obu_extension_flag==1 {
            let temporal_id = 0;
            self.coder.push_bits_with_size(&mut out_bits, temporal_id, 3);
            let spatial_id = 0;
            self.coder.push_bits_with_size(&mut out_bits, spatial_id, 2);
            let extension_header_reserved_3bits = 0;
            self.coder.push_bits_with_size(&mut out_bits, extension_header_reserved_3bits, 3);
        }
    }

    pub fn encode_obu_payload(&mut self, obu: &mut OBU, out_bits: &mut Vec<u8>) -> () {
        // TODO operating point

        match obu {
            &mut OBU_SEQUENCE_HEADER(ref mut sequence_header) => {
                self.coder.push_bits_with_size(out_bits, sequence_header.profile as u32, 3);
                self.coder.push_bit(out_bits, sequence_header.still_picture as u8);
                self.coder.push_bit(out_bits, sequence_header.reduced_still_picture_header as u8);
                if sequence_header.reduced_still_picture_header {
                    // TODO
                } else {
                    // FIXME
                    self.coder.push_bits_with_size(out_bits, sequence_header.operating_points.len() as u32, 5);
                    for op in sequence_header.operating_points.iter() {
                        self.coder.push_bits_with_size(out_bits, op.idc as u32, 12);
                        self.coder.push_bits_with_size(out_bits, op.seq_level_idx as u32, 5);
                        if op.seq_level_idx>7 {
                            self.coder.push_bit(out_bits, op.seq_tier.unwrap() as u8);
                        } else {
                            assert!(op.seq_tier==None);
                        }
                    }
                }
                let frame_width_bits_minus_1 = msb16(sequence_header.frame_width as u16) - 1;
                let frame_height_bits_minus_1 = msb16(sequence_header.frame_height as u16) - 1;
                self.coder.push_bits_with_size(out_bits, frame_width_bits_minus_1 as u32, 4);
                self.coder.push_bits_with_size(out_bits, frame_height_bits_minus_1 as u32, 4);
                let frame_width_bits_minus_1 = msb16(sequence_header.max_frame_width as u16) - 1;
                let frame_height_bits_minus_1 = msb16(sequence_header.max_frame_height as u16) - 1;
                let nw = frame_width_bits_minus_1 + 1;
                let nh = frame_height_bits_minus_1 + 1;
                let max_frame_width_minus_1 = sequence_header.max_frame_width - 1;
                let max_frame_height_minus_1 = sequence_header.max_frame_height - 1;
                self.coder.push_bits_with_size(out_bits, max_frame_width_minus_1 as u32, nw);
                self.coder.push_bits_with_size(out_bits, max_frame_height_minus_1 as u32, nh);
                if sequence_header.reduced_still_picture_header {
                    assert!(sequence_header.frame_id.is_none());
                } else {
                    self.coder.push_bit(out_bits, if let Some(_) = sequence_header.frame_id { 1 } else { 0 } as u8);
                }
                if let Some(ref fid) = sequence_header.frame_id {
                    self.coder.push_bits_with_size(out_bits, (fid.delta_frame_id_length-2) as u32, 4);
                    self.coder.push_bits_with_size(out_bits, (fid.additional_frame_id_length-1) as u32, 3);
                }
                self.coder.push_bit(out_bits, sequence_header.use_128x128_superblock as u8);
                self.coder.push_bit(out_bits, sequence_header.enable_filter_intra as u8);
                self.coder.push_bit(out_bits, sequence_header.enable_intra_edge_filter as u8);
                if sequence_header.reduced_still_picture_header {
                    sequence_header.enable_interintra_compound = false;
                    sequence_header.enable_masked_compound = false;
                    sequence_header.enable_warped_motion = false;
                    sequence_header.enable_dual_filter = false;
                    sequence_header.enable_order_hint = false;
                    sequence_header.enable_jnt_comp = false;
                    sequence_header.enable_ref_frame_mvs = false;
                    sequence_header.seq_force_screen_content_tools = SELECT_SCREEN_CONTENT_TOOLS;
                    sequence_header.seq_force_integer_mv = SELECT_INTEGER_MV;
                    self.order_hint_bits = 0;
                } else {
                    self.coder.push_bit(out_bits, sequence_header.enable_interintra_compound as u8);
                    self.coder.push_bit(out_bits, sequence_header.enable_masked_compound as u8);
                    self.coder.push_bit(out_bits, sequence_header.enable_warped_motion as u8);
                    self.coder.push_bit(out_bits, sequence_header.enable_dual_filter as u8);
                    self.coder.push_bit(out_bits, sequence_header.enable_order_hint as u8);
                    if sequence_header.enable_order_hint {
                        self.coder.push_bit(out_bits, sequence_header.enable_jnt_comp as u8);
                        self.coder.push_bit(out_bits, sequence_header.enable_ref_frame_mvs as u8);
                    } else {
                        sequence_header.enable_jnt_comp = false;
                        sequence_header.enable_ref_frame_mvs = false;
                    }
                    self.coder.push_bit(out_bits, sequence_header.seq_choose_screen_content_tools as u8);
                    if sequence_header.seq_choose_screen_content_tools {
                        sequence_header.seq_force_screen_content_tools = SELECT_SCREEN_CONTENT_TOOLS;
                    } else {
                        self.coder.push_bit(out_bits, sequence_header.seq_force_screen_content_tools as u8);
                    }
                    if sequence_header.seq_force_screen_content_tools>0 {
                        self.coder.push_bit(out_bits, sequence_header.seq_choose_integer_mv as u8);
                        if sequence_header.seq_choose_integer_mv {
                            sequence_header.seq_force_integer_mv = SELECT_INTEGER_MV;
                        } else {
                            self.coder.push_bit(out_bits, sequence_header.seq_force_integer_mv as u8);
                        }
                    } else {
                        sequence_header.seq_force_integer_mv = SELECT_INTEGER_MV;
                    }
                    if sequence_header.enable_order_hint {
                        self.coder.push_bits_with_size(out_bits, (self.order_hint_bits-1) as u32, 3);
                    } else {
                        self.order_hint_bits = 0;
                    }
                }
                self.coder.push_bit(out_bits, sequence_header.enable_superres as u8);
                self.coder.push_bit(out_bits, sequence_header.enable_cdef as u8);
                self.coder.push_bit(out_bits, sequence_header.enable_restoration as u8);
                // color config
                self.coder.push_bit(out_bits, sequence_header.color_config.high_bitdepth as u8);
                if sequence_header.profile==2 && sequence_header.color_config.high_bitdepth {
                    self.coder.push_bit(out_bits, sequence_header.color_config.twelve_bit as u8);
                    self.bit_depth = if sequence_header.color_config.twelve_bit { 12 } else { 10 };
                } else if sequence_header.profile<=2 {
                    self.bit_depth = if sequence_header.color_config.high_bitdepth { 10 } else { 8 };
                }
                if sequence_header.profile==1 {
                    sequence_header.color_config.mono_chrome = false;
                } else {
                    self.coder.push_bit(out_bits, sequence_header.color_config.mono_chrome as u8);
                }
                self.num_planes = if sequence_header.color_config.mono_chrome { 1 } else { 3 };
                self.coder.push_bit(out_bits, sequence_header.color_config.color_description.is_some() as u8);
                match sequence_header.color_config.color_description {
                    Some(ref mut cd) => {
                        self.coder.push_bits_with_size(out_bits, cd.color_primaries as u32, 8);
                        self.coder.push_bits_with_size(out_bits, cd.transfer_characteristics as u32, 8);
                        self.coder.push_bits_with_size(out_bits, cd.matrix_coefficients as u32, 8);
                    },
                    None => {
                        sequence_header.color_config.color_description = Some( ColorDescription {
                            color_primaries: CP_UNSPECIFIED,
                            transfer_characteristics: TC_UNSPECIFIED,
                            matrix_coefficients: MC_UNSPECIFIED
                        });
                    },
                };
                if sequence_header.color_config.mono_chrome {
                    self.coder.push_bit(out_bits, sequence_header.color_config.color_range as u8);
                    sequence_header.color_config.subsampling_x = true;
                    sequence_header.color_config.subsampling_y = true;
                    sequence_header.color_config.chroma_sample_position = CSP_UNKNOWN;
                    sequence_header.color_config.separate_uv_delta_q = false;
                } else if sequence_header.color_config.color_description.unwrap() == (ColorDescription {
                    color_primaries: CP_BT_709,
                    transfer_characteristics: TC_SRGB,
                    matrix_coefficients: MC_IDENTITY
                }) {
                    sequence_header.color_config.color_range = true;
                    sequence_header.color_config.subsampling_x = false;
                    sequence_header.color_config.subsampling_y = false;
                } else {
                    self.coder.push_bit(out_bits, sequence_header.color_config.color_range as u8);
                    match sequence_header.profile {
                        0 => {
                            sequence_header.color_config.subsampling_x = true;
                            sequence_header.color_config.subsampling_y = true;
                        },
                        1 => {
                            sequence_header.color_config.subsampling_x = false;
                            sequence_header.color_config.subsampling_y = false;
                        },
                        _ => {
                            if self.bit_depth==12 {
                                self.coder.push_bit(out_bits, sequence_header.color_config.subsampling_x as u8);
                                if sequence_header.color_config.subsampling_x {
                                    self.coder.push_bit(out_bits, sequence_header.color_config.subsampling_y as u8);
                                } else {
                                    sequence_header.color_config.subsampling_y = false;
                                }
                            } else {
                                sequence_header.color_config.subsampling_x = true;
                                sequence_header.color_config.subsampling_y = false;
                            }
                        }
                    };
                    if sequence_header.color_config.subsampling_x && sequence_header.color_config.subsampling_y {
                        self.coder.push_bits_with_size(out_bits, sequence_header.color_config.chroma_sample_position as u32, 2);
                    }
                }
                self.coder.push_bit(out_bits, sequence_header.color_config.separate_uv_delta_q as u8);

                if sequence_header.reduced_still_picture_header {
                    sequence_header.timing_info = None;
                } else {
                    self.coder.push_bit(out_bits, if let Some(_) = sequence_header.timing_info { 1 } else { 0 } as u8);
                }
                match sequence_header.timing_info {
                    Some(ref ti) => {
                        // timing info
                        self.coder.push_bits_with_size(out_bits, ti.num_units_in_display_tick as u32, 32);
                        self.coder.push_bits_with_size(out_bits, ti.time_scale as u32, 32);
                        self.coder.push_bit(out_bits, ti.equal_picture_interval as u8);
                        if ti.equal_picture_interval {
                            self.coder.encode_uvlc(out_bits, ti.num_ticks_per_picture-1);
                        }

                        self.coder.push_bit(out_bits, sequence_header.decoder_model_info.is_some() as u8);
                        if let Some(ref dmi) = sequence_header.decoder_model_info {
                            self.coder.push_bits_with_size(out_bits, dmi.bitrate_scale as u32, 4);
                            self.coder.push_bits_with_size(out_bits, dmi.buffer_size_scale as u32, 4);
                            self.coder.push_bits_with_size(out_bits, (dmi.encoder_decoder_buffer_delay_length-1) as u32, 5);
                            self.coder.push_bits_with_size(out_bits, dmi.num_units_in_decoding_tick as u32, 32);
                            self.coder.push_bits_with_size(out_bits, (dmi.buffer_removal_delay_length-1) as u32, 5);
                            self.coder.push_bits_with_size(out_bits, (dmi.frame_presentation_delay_length-1) as u32, 5);
                        }
                    },
                    None => {
                        sequence_header.decoder_model_info = None;
                    },
                };
                self.coder.push_bit(out_bits, if sequence_header.operating_points_decoder_model.len()>0 { 1 } else { 0 } as u8);
                if sequence_header.operating_points_decoder_model.len()>0 {
                    self.coder.push_bits_with_size(out_bits, (sequence_header.operating_points_decoder_model.len()-1) as u32, 5);
                    for op in sequence_header.operating_points_decoder_model.iter() {
                        self.coder.push_bits_with_size(out_bits, op.idc as u32, 12);
                        self.coder.push_bit(out_bits, op.initial_display_delay.is_some() as u8);
                        if let Some(idd) = op.initial_display_delay {
                            self.coder.push_bits_with_size(out_bits, (idd-1) as u32, 4);
                        }
                        if sequence_header.decoder_model_info.is_some() {
                            self.coder.push_bit(out_bits, op.operating_parameters_info.is_some() as u8);
                            if let Some(ref opi) = op.operating_parameters_info {
                                self.coder.encode_uvlc(out_bits, opi.bitrate-1);
                                self.coder.encode_uvlc(out_bits, opi.buffer_size-1);
                                self.coder.push_bit(out_bits, opi.cbr_flag as u8);
                                let n = if let Some(ref dmi) = sequence_header.decoder_model_info { dmi.encoder_decoder_buffer_delay_length } else { 0 };
                                self.coder.push_bits_with_size(out_bits, opi.decoder_buffer_delay as u32, n);
                                self.coder.push_bits_with_size(out_bits, opi.encoder_buffer_delay as u32, n);
                                self.coder.push_bit(out_bits, opi.low_delay_mode_flag as u8);
                            }
                        }
                    }
                }
                self.coder.push_bit(out_bits, sequence_header.film_grain_params_present as u8);
            },
            &mut OBU_TEMPORAL_DELIMITER => {
                self.seen_frame_header = false;
            },
            &mut OBU_FRAME_HEADER(ref frame_header) => {

            },
            &mut OBU_REDUNDANT_FRAME_HEADER => {

            },
            &mut OBU_TILE_GROUP => {

            },
            &mut OBU_METADATA(ref meta_data) => {
                match meta_data {
                    &METADATA_TYPE_ITUT_T35 { country_code: cc, payload_bytes: ref pb } => {
                        self.coder.encode_leb128(out_bits, 4);
                        self.coder.push_bits_with_size(out_bits, cc as u32, 8);
                        if cc==0xFF {
                            self.coder.push_bits_with_size(out_bits, 0 as u32, 8); // FIXME itu_t_t35_country_code_extension_byte
                        }
                        // FIXME itu_t_t35_payload_bytes not in spec
                    },
                    &METADATA_TYPE_HDR_CLL { max_cll: mc, max_fall: mf } => {
                        self.coder.encode_leb128(out_bits, 1);
                        self.coder.push_bits_with_size(out_bits, mc as u32, 16);
                        self.coder.push_bits_with_size(out_bits, mf as u32, 16);
                    },
                    &METADATA_TYPE_HDR_MDCV {
                        primary_chromaticity_x: pcx,
                        primary_chromaticity_y: pcy,
                        white_point_chromaticity_x: wpcx,
                        white_point_chromaticity_y: wpcy,
                        luminance_max: lmax,
                        luminance_min: lmin,
                    } => {
                        self.coder.encode_leb128(out_bits, 2);
                        for i in 0..3 {
                            self.coder.push_bits_with_size(out_bits, pcx[i] as u32, 16);
                            self.coder.push_bits_with_size(out_bits, pcy[i] as u32, 16);
                        }
                            self.coder.push_bits_with_size(out_bits, wpcx as u32, 16);
                            self.coder.push_bits_with_size(out_bits, wpcy as u32, 16);
                            self.coder.push_bits_with_size(out_bits, lmax as u32, 32);
                            self.coder.push_bits_with_size(out_bits, lmin as u32, 32);
                    },
                    &METADATA_TYPE_SCALABILITY { mode_idc: mi, structure: ref st } => {
                        self.coder.encode_leb128(out_bits, 3);
                        self.coder.push_bits_with_size(out_bits, mi as u32, 8);
                        if mi==SCALABILITY_SS {
                            // scalability structure
                            if let &Some(ref s) = st {
                                self.coder.push_bits_with_size(out_bits, (Max!(s.spatial_layer_dimensions.len(), s.spatial_layer_descriptions.len())-1) as u32, 2);
                                self.coder.push_bit(out_bits, (s.spatial_layer_dimensions.len()>0) as u8);
                                self.coder.push_bit(out_bits, (s.spatial_layer_descriptions.len()>0) as u8);
                                self.coder.push_bit(out_bits, (s.temporal_groups.len()>0) as u8);
                                let scalability_structure_reserved_3bits = 0;
                                self.coder.push_bits_with_size(out_bits, scalability_structure_reserved_3bits, 3);
                                for dm in s.spatial_layer_dimensions.iter() {
                                    self.coder.push_bits_with_size(out_bits, dm.max_width as u32, 16);
                                    self.coder.push_bits_with_size(out_bits, dm.max_height as u32, 16);
                                }
                                for ds in s.spatial_layer_descriptions.iter() {
                                    self.coder.push_bits_with_size(out_bits, ds.ref_id as u32, 8);
                                }
                                if s.temporal_groups.len()>0 {
                                    self.coder.push_bits_with_size(out_bits, s.temporal_groups.len() as u32, 8);
                                    for tg in s.temporal_groups.iter() {
                                        self.coder.push_bits_with_size(out_bits, tg.temporal_id as u32, 3);
                                        self.coder.push_bit(out_bits, tg.temporal_switching_up_point_flag);
                                        self.coder.push_bit(out_bits, tg.spatial_switching_up_point_flag);
                                        self.coder.push_bits_with_size(out_bits, tg.ref_pic_diffs.len() as u32, 3);
                                        for d in tg.ref_pic_diffs.iter() {
                                            self.coder.push_bits_with_size(out_bits, *d as u32, 8);
                                        }
                                    }
                                }
                            } else {
                                assert!(false);
                            }
                        }
                    },
                    &METADATA_TYPE_TIMECODE {
                        counting_type: ct,
                        seconds: seconds,
                        minutes: minutes,
                        hours: hours,
                        time_offset_length: tol,
                        time_offset_value: tov,
                    } => {
                        self.coder.encode_leb128(out_bits, 3);
                        self.coder.push_bits_with_size(out_bits, ct as u32, 5);
                        if let (Some(s), Some(m), Some(h)) = (seconds, minutes, hours) {
                            self.coder.push_bit(out_bits, 1);
                            self.coder.push_bits_with_size(out_bits, s as u32, 6);
                            self.coder.push_bits_with_size(out_bits, m as u32, 6);
                            self.coder.push_bits_with_size(out_bits, h as u32, 5);
                        } else {
                            self.coder.push_bit(out_bits, 0);
                            if let Some(s) = seconds {
                                self.coder.push_bit(out_bits, 1);
                                self.coder.push_bits_with_size(out_bits, s as u32, 6);
                                if let Some(m) = minutes {
                                    self.coder.push_bit(out_bits, 1);
                                    self.coder.push_bits_with_size(out_bits, m as u32, 6);
                                    if let Some(h) = hours {
                                        self.coder.push_bit(out_bits, 1);
                                        self.coder.push_bits_with_size(out_bits, h as u32, 5);
                                    } else {
                                        self.coder.push_bit(out_bits, 0);
                                    }
                                } else {
                                    self.coder.push_bit(out_bits, 0);
                                }
                            } else {
                                self.coder.push_bit(out_bits, 0);
                            }
                        }
                        self.coder.push_bits_with_size(out_bits, tol as u32, 5);
                        self.coder.push_bits_with_size(out_bits, tov as u32, tol);
                    },
                    _ => { assert!(false); }
                };
            },
            &mut OBU_FRAME { ref frame } => {

            },
            &mut OBU_PADDING(padding) => {
                for i in 0..padding {
                    self.coder.push_bits_with_size(out_bits, 0, 8);
                }
            },
        };
    }

    pub fn encode_obu(&mut self, obu: &mut OBU, out_bits: &mut Vec<u8>) -> () {
        self.encode_obu_header(obu, out_bits);
        let mut obu_payload_bits = vec![];
        self.encode_obu_payload(obu, &mut obu_payload_bits);
        /*
        if obu.obu_has_size_field {
            let obu_size = obu_payload_bits.len();
        }*/
        self.writer.write_bits(&obu_payload_bits);
        self.writer.byte_align();
    }

    pub fn encode_temporal_unit(&mut self, mut obu: &mut OBU) -> () {
        let mut obu_bits = vec![];
        self.encode_obu(&mut obu, &mut obu_bits);
        let obu_length = obu_bits.len();
        let mut frame_unit_bits = vec![];
        self.coder.encode_leb128(&mut frame_unit_bits, obu_length as u64);
        self.coder.push_bits(&mut frame_unit_bits, &obu_bits);
        let frame_unit_size = frame_unit_bits.len();
        let mut temporal_unit_bits = vec![];
        self.coder.encode_leb128(&mut temporal_unit_bits, frame_unit_size as u64);
        self.coder.push_bits(&mut temporal_unit_bits, &frame_unit_bits);
        let temporal_unit_size = temporal_unit_bits.len();
        let mut temporal_unit_size_bits = vec![];
        self.coder.encode_leb128(&mut temporal_unit_size_bits, temporal_unit_size as u64);

        self.writer.write_bits(&temporal_unit_size_bits);
        self.writer.write_bits(&temporal_unit_bits);
    }
}