use std::rc::Rc;
use std::cell::RefCell;
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
use encoder::*;
use constants::FrameType::*;
use constants::RefFrame::*;
use constants::InterpFilter::*;
use constants::FrameRestorationType::*;
use constants::TxMode::*;
use frame::*;
use tile_encoder::*;


pub struct OBUEncoder<'a, 'b> where 'a: 'b {
    writer: &'b mut BinaryWriter<'a>,
    coder: &'b mut BoolCoder,
    //operating_point_idc: usize, // FIXME
    order_hint_bits: usize,
    bit_depth: usize,
    num_planes: usize,
    tile_num: usize,
    id_len: usize,
    all_frames: usize,
    ectx: Rc<RefCell<EncoderContext>>,
}

impl<'a, 'b> OBUEncoder<'a, 'b> where 'a: 'b {
    pub fn new(ectx: &Rc<RefCell<EncoderContext>>, writer: &'b mut BinaryWriter<'a>, coder: &'b mut BoolCoder) -> OBUEncoder<'a, 'b> {
        OBUEncoder {
            writer: writer,
            coder: coder,
            //operating_point_idc: 0,
            order_hint_bits: 0,
            bit_depth: 8,
            num_planes: 3,
            tile_num: 1,
            id_len: 0,
            all_frames: 0,
            ectx: Rc::clone(ectx),
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
            &mut OBU_TILE_GROUP(_) => 4,
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

                let mut ecx = self.ectx.borrow_mut();

                self.coder.push_bits_with_size(out_bits, sequence_header.profile as u32, 3);
                self.coder.push_bit(out_bits, sequence_header.still_picture as u8);
                self.coder.push_bit(out_bits, ecx.reduced_still_picture_header as u8);
                if ecx.reduced_still_picture_header {
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
                let frame_width_bits_minus_1 = msb16((sequence_header.frame_width-1) as u16) - 1;
                let frame_height_bits_minus_1 = msb16((sequence_header.frame_height-1) as u16) - 1;
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
                if ecx.reduced_still_picture_header {
                    assert!(sequence_header.frame_id.is_none());
                } else {
                    self.coder.push_bit(out_bits, if let Some(_) = sequence_header.frame_id { 1 } else { 0 } as u8);
                }
                if let Some(ref fid) = sequence_header.frame_id {
                    self.coder.push_bits_with_size(out_bits, (fid.delta_frame_id_length-2) as u32, 4);
                    self.coder.push_bits_with_size(out_bits, (fid.additional_frame_id_length-1) as u32, 3);
                }
                self.coder.push_bit(out_bits, ecx.use_128x128_superblock as u8);
                self.coder.push_bit(out_bits, sequence_header.enable_filter_intra as u8);
                self.coder.push_bit(out_bits, sequence_header.enable_intra_edge_filter as u8);
                if ecx.reduced_still_picture_header {
                    sequence_header.enable_interintra_compound = false;
                    sequence_header.enable_masked_compound = false;
                    ecx.enable_warped_motion = false;
                    sequence_header.enable_dual_filter = false;
                    ecx.enable_order_hint = false;
                    sequence_header.enable_jnt_comp = false;
                    ecx.enable_ref_frame_mvs = false;
                    ecx.seq_force_screen_content_tools = SELECT_SCREEN_CONTENT_TOOLS;
                    ecx.seq_force_integer_mv = SELECT_INTEGER_MV;
                    self.order_hint_bits = 0;
                } else {
                    self.coder.push_bit(out_bits, sequence_header.enable_interintra_compound as u8);
                    self.coder.push_bit(out_bits, sequence_header.enable_masked_compound as u8);
                    self.coder.push_bit(out_bits, ecx.enable_warped_motion as u8);
                    self.coder.push_bit(out_bits, sequence_header.enable_dual_filter as u8);
                    self.coder.push_bit(out_bits, ecx.enable_order_hint as u8);
                    if ecx.enable_order_hint {
                        self.coder.push_bit(out_bits, sequence_header.enable_jnt_comp as u8);
                        self.coder.push_bit(out_bits, ecx.enable_ref_frame_mvs as u8);
                    } else {
                        sequence_header.enable_jnt_comp = false;
                        ecx.enable_ref_frame_mvs = false;
                    }
                    self.coder.push_bit(out_bits, sequence_header.seq_choose_screen_content_tools as u8);
                    if sequence_header.seq_choose_screen_content_tools {
                        ecx.seq_force_screen_content_tools = SELECT_SCREEN_CONTENT_TOOLS;
                    } else {
                        self.coder.push_bit(out_bits, ecx.seq_force_screen_content_tools as u8);
                    }
                    if ecx.seq_force_screen_content_tools>0 {
                        self.coder.push_bit(out_bits, sequence_header.seq_choose_integer_mv as u8);
                        if sequence_header.seq_choose_integer_mv {
                            ecx.seq_force_integer_mv = SELECT_INTEGER_MV;
                        } else {
                            self.coder.push_bit(out_bits, ecx.seq_force_integer_mv as u8);
                        }
                    } else {
                        ecx.seq_force_integer_mv = SELECT_INTEGER_MV;
                    }
                    if ecx.enable_order_hint {
                        self.coder.push_bits_with_size(out_bits, (self.order_hint_bits-1) as u32, 3);
                    } else {
                        self.order_hint_bits = 0;
                    }
                }
                self.coder.push_bit(out_bits, ecx.enable_superres as u8);
                self.coder.push_bit(out_bits, ecx.enable_cdef as u8);
                self.coder.push_bit(out_bits, ecx.enable_restoration as u8);
                // color config
                self.coder.push_bit(out_bits, ecx.color_config.high_bitdepth as u8);
                if sequence_header.profile==2 && ecx.color_config.high_bitdepth {
                    self.coder.push_bit(out_bits, ecx.color_config.twelve_bit as u8);
                    self.bit_depth = if ecx.color_config.twelve_bit { 12 } else { 10 };
                } else if sequence_header.profile<=2 {
                    self.bit_depth = if ecx.color_config.high_bitdepth { 10 } else { 8 };
                }
                if sequence_header.profile==1 {
                    ecx.color_config.mono_chrome = false;
                } else {
                    self.coder.push_bit(out_bits, ecx.color_config.mono_chrome as u8);
                }
                self.num_planes = if ecx.color_config.mono_chrome { 1 } else { 3 };
                self.coder.push_bit(out_bits, ecx.color_config.color_description.is_some() as u8);
                if let Some(ref mut cd) = ecx.color_config.color_description {
                    self.coder.push_bits_with_size(out_bits, cd.color_primaries as u32, 8);
                    self.coder.push_bits_with_size(out_bits, cd.transfer_characteristics as u32, 8);
                    self.coder.push_bits_with_size(out_bits, cd.matrix_coefficients as u32, 8);
                }
                if ecx.color_config.color_description.is_none() {
                    ecx.color_config.color_description = Some( ColorDescription {
                        color_primaries: CP_UNSPECIFIED,
                        transfer_characteristics: TC_UNSPECIFIED,
                        matrix_coefficients: MC_UNSPECIFIED
                    });
                }
                if ecx.color_config.mono_chrome {
                    self.coder.push_bit(out_bits, ecx.color_config.color_range as u8);
                    ecx.color_config.subsampling_x = true;
                    ecx.color_config.subsampling_y = true;
                    ecx.color_config.chroma_sample_position = CSP_UNKNOWN;
                    ecx.color_config.separate_uv_delta_q = false;
                } else if ecx.color_config.color_description.unwrap() == (ColorDescription {
                    color_primaries: CP_BT_709,
                    transfer_characteristics: TC_SRGB,
                    matrix_coefficients: MC_IDENTITY
                }) {
                    ecx.color_config.color_range = true;
                    ecx.color_config.subsampling_x = false;
                    ecx.color_config.subsampling_y = false;
                } else {
                    self.coder.push_bit(out_bits, ecx.color_config.color_range as u8);
                    match sequence_header.profile {
                        0 => {
                            ecx.color_config.subsampling_x = true;
                            ecx.color_config.subsampling_y = true;
                        },
                        1 => {
                            ecx.color_config.subsampling_x = false;
                            ecx.color_config.subsampling_y = false;
                        },
                        _ => {
                            if self.bit_depth==12 {
                                self.coder.push_bit(out_bits, ecx.color_config.subsampling_x as u8);
                                if ecx.color_config.subsampling_x {
                                    self.coder.push_bit(out_bits, ecx.color_config.subsampling_y as u8);
                                } else {
                                    ecx.color_config.subsampling_y = false;
                                }
                            } else {
                                ecx.color_config.subsampling_x = true;
                                ecx.color_config.subsampling_y = false;
                            }
                        }
                    };
                    if ecx.color_config.subsampling_x && ecx.color_config.subsampling_y {
                        self.coder.push_bits_with_size(out_bits, ecx.color_config.chroma_sample_position as u32, 2);
                    }
                }
                self.coder.push_bit(out_bits, ecx.color_config.separate_uv_delta_q as u8);

                if ecx.reduced_still_picture_header {
                    ecx.timing_info = None;
                } else {
                    self.coder.push_bit(out_bits, if let Some(_) = ecx.timing_info { 1 } else { 0 } as u8);
                }
                match ecx.timing_info {
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
                self.coder.push_bit(out_bits, if ecx.operating_points_decoder_model.len()>0 { 1 } else { 0 } as u8);
                if ecx.operating_points_decoder_model.len()>0 {
                    self.coder.push_bits_with_size(out_bits, (ecx.operating_points_decoder_model.len()-1) as u32, 5);
                    for op in ecx.operating_points_decoder_model.iter() {
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
                self.coder.push_bit(out_bits, ecx.film_grain_params.is_some() as u8);
            },
            &mut OBU_TEMPORAL_DELIMITER => {
                let seen_frame_header = &mut self.ectx.borrow_mut().seen_frame_header;
                *seen_frame_header = false;
            },
            &mut OBU_FRAME_HEADER(ref mut frame) => {
                let mut ecx = self.ectx.borrow_mut();
                if ecx.seen_frame_header {
                    // TODO frame header copy
                } else {
                    ecx.seen_frame_header = true;

                    // uncompressed header
                    if ecx.frame_id_numbers_present_flag {
                        self.id_len = (ecx.additional_frame_id_length-1) + (ecx.delta_frame_id_length-2)+3;
                    }
                    self.all_frames = (1<<NUM_REF_FRAMES) - 1;
                    if ecx.reduced_still_picture_header {
                        ecx.show_existing_frame = false;
                        frame.frame_type = KEY_FRAME;
                        ecx.frame_is_intra = true;
                        ecx.show_frame = true;
                        ecx.showable_frame = false;
                    } else {
                        self.coder.push_bit(out_bits, ecx.show_existing_frame as u8);
                        if ecx.show_existing_frame {
                            self.coder.push_bits_with_size(out_bits, ecx.frame_to_show_map_idx as u32, 3);
                            if let Some(ref dmi) = frame.decoder_model_info {
                                if !(if let Some(ti) = ecx.timing_info { ti.equal_picture_interval } else { false }) {
                                    // temporal point info
                                    let n = dmi.frame_presentation_delay_length;
                                    if let Some(ref tpi) = frame.temporal_point_info {
                                        self.coder.push_bits_with_size(out_bits, tpi.tu_presentation_delay as u32, n);
                                    } else {
                                        assert!(false);
                                    }
                                }
                            }
                            ecx.refresh_frame_flags = 0;
                            if ecx.frame_id_numbers_present_flag {
                                self.coder.push_bits_with_size(out_bits, frame.display_frame_id.unwrap() as u32, self.id_len as u8);
                            }
                            frame.frame_type = ecx.ref_frame_type[ecx.frame_to_show_map_idx as usize];
                            if frame.frame_type==KEY_FRAME {
                                ecx.refresh_frame_flags = self.all_frames as u8;
                            }
                            if let Some(ref fgp) = ecx.film_grain_params {
                                // load grain params
                                // TODO
                            }
                            return;
                        }
                        self.coder.push_bits_with_size(out_bits, frame.frame_type as u32, 2);
                        ecx.frame_is_intra = frame.frame_type==INTRA_ONLY_FRAME || frame.frame_type==KEY_FRAME;
                        self.coder.push_bit(out_bits, ecx.show_frame as u8);
                        if ecx.show_frame {
                            if let Some(ref dmi) = frame.decoder_model_info {
                                if !(if let Some(ti) = ecx.timing_info { ti.equal_picture_interval } else { false }) {
                                    // temporal point info
                                    let n = dmi.frame_presentation_delay_length;
                                    if let Some(ref tpi) = frame.temporal_point_info {
                                        self.coder.push_bits_with_size(out_bits, tpi.tu_presentation_delay as u32, n);
                                    } else {
                                        assert!(false);
                                    }
                                }
                            }
                        }
                        if ecx.show_frame {
                            ecx.showable_frame = false;
                        } else {
                            self.coder.push_bit(out_bits, ecx.showable_frame as u8);
                        }
                        if frame.frame_type==SWITCH_FRAME || (frame.frame_type==KEY_FRAME && ecx.show_frame) {
                            frame.error_resilient_mode = false;
                        } else {
                            self.coder.push_bit(out_bits, frame.error_resilient_mode as u8);
                        }
                    }
                    if frame.frame_type==KEY_FRAME && ecx.show_frame {
                        for i in 0..NUM_REF_FRAMES {
                            ecx.ref_valid[i] = false;
                            ecx.ref_order_hint[i] = 0;
                        }
                        for i in 0..REFS_PER_FRAME {
                            ecx.order_hints[LAST_FRAME as usize + i] = 0;
                        }
                    }
                    self.coder.push_bit(out_bits, ecx.disable_cdf_update as u8);
                    if ecx.seq_force_screen_content_tools==SELECT_SCREEN_CONTENT_TOOLS {
                        self.coder.push_bit(out_bits, frame.allow_screen_content_tools as u8);
                    } else {
                        frame.allow_screen_content_tools = ecx.seq_force_screen_content_tools;
                    }
                    if frame.allow_screen_content_tools>0 {
                        if ecx.seq_force_integer_mv==SELECT_INTEGER_MV {
                            self.coder.push_bit(out_bits, frame.force_integer_mv as u8);
                        } else {
                            frame.force_integer_mv = ecx.seq_force_integer_mv;
                        }
                    } else {
                        frame.force_integer_mv = 0;
                    }
                    if ecx.frame_is_intra {
                        frame.force_integer_mv = 1;
                    }
                    if ecx.frame_id_numbers_present_flag {
                        ecx.prev_frame_id = ecx.current_frame_id;
                        self.coder.push_bits_with_size(out_bits, ecx.current_frame_id as u32, self.id_len as u8);
                        // mark ref frames
                        let diff_len = ecx.delta_frame_id_length;
                        for i in 0..NUM_REF_FRAMES {
                            if ecx.current_frame_id > (1<<diff_len) {
                                if ecx.ref_frame_id[i] > ecx.current_frame_id || ecx.ref_frame_id[i] < (ecx.current_frame_id-(1<<diff_len)) {
                                    ecx.ref_valid[i] = false;
                                }
                            } else {
                                if ecx.ref_frame_id[i] > ecx.current_frame_id && ecx.ref_frame_id[i] < ((1<<self.id_len)+ecx.current_frame_id-(1<<diff_len)) {
                                    ecx.ref_valid[i] = false;
                                }
                            }
                        }
                    } else {
                        ecx.current_frame_id = 0;
                    }
                    if frame.frame_type==SWITCH_FRAME {
                        ecx.frame_size_override_flag = true;
                    } else if ecx.reduced_still_picture_header {
                        ecx.frame_size_override_flag = false;
                    } else {
                        self.coder.push_bit(out_bits, ecx.frame_size_override_flag as u8);
                    }
                    self.coder.push_bits_with_size(out_bits, frame.order_hint as u32, self.order_hint_bits as u8);
                    // OrderHint = order_hint
                    if ecx.frame_is_intra || frame.error_resilient_mode {
                        frame.primary_ref_frame = PRIMARY_REF_NONE as u8;
                    } else {
                        self.coder.push_bits_with_size(out_bits, frame.primary_ref_frame as u32, 3);
                    }
                    if let Some(ref dmi) = frame.decoder_model_info {
                        self.coder.push_bit(out_bits, dmi.buffer_removal_delay_present as u8);
                        if dmi.buffer_removal_delay_present {
                            for op in ecx.operating_points_decoder_model.iter() {
                                let op_pt_idc = op.idc;
                                let in_temporal_layer = (op_pt_idc>>ecx.temporal_id) & 1 > 0;
                                let in_spatial_layer = (op_pt_idc>>(ecx.spatial_id+8)) & 1 > 0;
                                if op_pt_idc==0 || (in_temporal_layer&&in_spatial_layer) {
                                    let n = dmi.buffer_removal_delay_length;
                                    self.coder.push_bits_with_size(out_bits, dmi.buffer_removal_delay as u32, n);
                                }
                            }
                        }
                    }
                    frame.allow_high_precision_mv = false;
                    frame.use_ref_frame_mvs = false;
                    frame.allow_intrabc = false;
                    if frame.frame_type==SWITCH_FRAME || (frame.frame_type==KEY_FRAME && ecx.show_frame) {
                        ecx.refresh_frame_flags = self.all_frames as u8;
                    } else {
                        self.coder.push_bits_with_size(out_bits, ecx.refresh_frame_flags as u32, 8);
                    }
                    if !ecx.frame_is_intra || ecx.refresh_frame_flags != self.all_frames as u8 {
                        if frame.error_resilient_mode && ecx.enable_order_hint {
                            for i in 0..NUM_REF_FRAMES {
                                if frame.ref_order_hint[i] != ecx.ref_order_hint[i] {
                                    ecx.ref_valid[i] = false;
                                }
                            }
                        }
                    }
                    if frame.frame_type==KEY_FRAME {
                        // frame size
                        if ecx.frame_size_override_flag {
                            let nw = msb16((frame.frame_width-1) as u16);
                            let nh = msb16((frame.frame_height-1) as u16);
                            self.coder.push_bits_with_size(out_bits, (frame.frame_width-1) as u32, nw as u8);
                            self.coder.push_bits_with_size(out_bits, (frame.frame_height-1) as u32, nh as u8);
                        } else {
                            // TODO
                        }
                        {
                            // superres params
                            if ecx.enable_superres {
                                self.coder.push_bit(out_bits, frame.use_superres as u8);
                            } else {
                                frame.use_superres = false;
                            }
                            if frame.use_superres {
                                let coded_denom = ecx.superres_denom - SUPERRES_DENOM_MIN;
                                self.coder.push_bits_with_size(out_bits, coded_denom as u32, SUPERRES_DENOM_BITS as u8);
                            } else {
                                ecx.superres_denom = SUPERRES_NUM;
                            }
                            ecx.upscaled_width = frame.frame_width;
                            frame.frame_width = (ecx.upscaled_width*SUPERRES_NUM + ecx.superres_denom/2) / ecx.superres_denom;
                        }

                        // render size
                        let render_and_frame_size_different = frame.frame_width!=frame.render_width || frame.frame_height!=frame.render_height;
                        self.coder.push_bit(out_bits, render_and_frame_size_different as u8);
                        if render_and_frame_size_different {
                            self.coder.push_bits_with_size(out_bits, (frame.render_width-1) as u32, 16);
                            self.coder.push_bits_with_size(out_bits, (frame.render_height-1) as u32, 16);
                        } else {
                            frame.render_width = ecx.upscaled_width;
                            frame.render_height = frame.frame_height;
                        }

                        if frame.allow_screen_content_tools>0 {
                            self.coder.push_bit(out_bits, frame.allow_intrabc as u8);
                        }
                    } else {
                        if frame.frame_type==INTRA_ONLY_FRAME {
                            // frame size
                            if ecx.frame_size_override_flag {
                                let nw = msb16((frame.frame_width-1) as u16);
                                let nh = msb16((frame.frame_height-1) as u16);
                                self.coder.push_bits_with_size(out_bits, (frame.frame_width-1) as u32, nw as u8);
                                self.coder.push_bits_with_size(out_bits, (frame.frame_height-1) as u32, nh as u8);
                            } else {
                                // TODO
                            }
                            {
                                // superres params
                                if ecx.enable_superres {
                                    self.coder.push_bit(out_bits, frame.use_superres as u8);
                                } else {
                                    frame.use_superres = false;
                                }
                                if frame.use_superres {
                                    let coded_denom = ecx.superres_denom - SUPERRES_DENOM_MIN;
                                    self.coder.push_bits_with_size(out_bits, coded_denom as u32, SUPERRES_DENOM_BITS as u8);
                                } else {
                                    ecx.superres_denom = SUPERRES_NUM;
                                }
                                ecx.upscaled_width = frame.frame_width;
                                frame.frame_width = (ecx.upscaled_width*SUPERRES_NUM + ecx.superres_denom/2) / ecx.superres_denom;
                            }

                            // render size
                            let render_and_frame_size_different = frame.frame_width!=frame.render_width || frame.frame_height!=frame.render_height;
                            self.coder.push_bit(out_bits, render_and_frame_size_different as u8);
                            if render_and_frame_size_different {
                                self.coder.push_bits_with_size(out_bits, (frame.render_width-1) as u32, 16);
                                self.coder.push_bits_with_size(out_bits, (frame.render_height-1) as u32, 16);
                            } else {
                                frame.render_width = ecx.upscaled_width;
                                frame.render_height = frame.frame_height;
                            }

                            if frame.allow_screen_content_tools>0 && ecx.upscaled_width==frame.frame_width {
                                self.coder.push_bit(out_bits, frame.allow_intrabc as u8);
                            }
                        } else {
                            if !ecx.enable_order_hint {
                                frame.frame_refs_short_signaling = false;
                            } else {
                                self.coder.push_bit(out_bits, frame.frame_refs_short_signaling as u8);
                                if frame.frame_refs_short_signaling {
                                    self.coder.push_bits_with_size(out_bits, frame.last_frame_idx as u32, 3);
                                    self.coder.push_bits_with_size(out_bits, frame.gold_frame_idx as u32, 3);
                                    // set frame refs
                                    // TODO
                                }
                            }
                            for i in 0..REFS_PER_FRAME {
                                if !frame.frame_refs_short_signaling {
                                    self.coder.push_bits_with_size(out_bits, frame.ref_frame_idx[i] as u32, 3);
                                }
                                if ecx.frame_id_numbers_present_flag {
                                    let n = ecx.delta_frame_id_length;
                                    self.coder.push_bits_with_size(out_bits, (frame.delta_frame_id[i]-1) as u32, n as u8);
                                    ecx.expected_frame_id[i] = (ecx.current_frame_id+(1<<self.id_len)-frame.delta_frame_id[i]) % (1<<self.id_len);
                                }
                            }
                            if ecx.frame_size_override_flag && !frame.error_resilient_mode {
                                // frame size with refs
                                let mut found_ref = false;
                                for i in 0..REFS_PER_FRAME {
                                    // FIXME
                                    found_ref = true;
                                    self.coder.push_bit(out_bits, found_ref as u8);
                                    break;
                                }
                                if !found_ref {
                                    // frame size
                                    if ecx.frame_size_override_flag {
                                        let nw = msb16((frame.frame_width-1) as u16);
                                        let nh = msb16((frame.frame_height-1) as u16);
                                        self.coder.push_bits_with_size(out_bits, (frame.frame_width-1) as u32, nw as u8);
                                        self.coder.push_bits_with_size(out_bits, (frame.frame_height-1) as u32, nh as u8);
                                    } else {
                                        // TODO
                                    }
                                    {
                                        // superres params
                                        if ecx.enable_superres {
                                            self.coder.push_bit(out_bits, frame.use_superres as u8);
                                        } else {
                                            frame.use_superres = false;
                                        }
                                        if frame.use_superres {
                                            let coded_denom = ecx.superres_denom - SUPERRES_DENOM_MIN;
                                            self.coder.push_bits_with_size(out_bits, coded_denom as u32, SUPERRES_DENOM_BITS as u8);
                                        } else {
                                            ecx.superres_denom = SUPERRES_NUM;
                                        }
                                        ecx.upscaled_width = frame.frame_width;
                                        frame.frame_width = (ecx.upscaled_width*SUPERRES_NUM + ecx.superres_denom/2) / ecx.superres_denom;
                                    }

                                    // render size
                                    let render_and_frame_size_different = frame.frame_width!=frame.render_width || frame.frame_height!=frame.render_height;
                                    self.coder.push_bit(out_bits, render_and_frame_size_different as u8);
                                    if render_and_frame_size_different {
                                        self.coder.push_bits_with_size(out_bits, (frame.render_width-1) as u32, 16);
                                        self.coder.push_bits_with_size(out_bits, (frame.render_height-1) as u32, 16);
                                    } else {
                                        frame.render_width = ecx.upscaled_width;
                                        frame.render_height = frame.frame_height;
                                    }
                                } else {
                                    // superres params
                                    if ecx.enable_superres {
                                        self.coder.push_bit(out_bits, frame.use_superres as u8);
                                    } else {
                                        frame.use_superres = false;
                                    }
                                    if frame.use_superres {
                                        let coded_denom = ecx.superres_denom - SUPERRES_DENOM_MIN;
                                        self.coder.push_bits_with_size(out_bits, coded_denom as u32, SUPERRES_DENOM_BITS as u8);
                                    } else {
                                        ecx.superres_denom = SUPERRES_NUM;
                                    }
                                    ecx.upscaled_width = frame.frame_width;
                                    frame.frame_width = (ecx.upscaled_width*SUPERRES_NUM + ecx.superres_denom/2) / ecx.superres_denom;
                                }
                            } else {
                                // frame size
                                if ecx.frame_size_override_flag {
                                    let nw = msb16((frame.frame_width-1) as u16);
                                    let nh = msb16((frame.frame_height-1) as u16);
                                    self.coder.push_bits_with_size(out_bits, (frame.frame_width-1) as u32, nw as u8);
                                    self.coder.push_bits_with_size(out_bits, (frame.frame_height-1) as u32, nh as u8);
                                } else {
                                    // TODO
                                }
                                {
                                    // superres params
                                    if ecx.enable_superres {
                                        self.coder.push_bit(out_bits, frame.use_superres as u8);
                                    } else {
                                        frame.use_superres = false;
                                    }
                                    if frame.use_superres {
                                        let coded_denom = ecx.superres_denom - SUPERRES_DENOM_MIN;
                                        self.coder.push_bits_with_size(out_bits, coded_denom as u32, SUPERRES_DENOM_BITS as u8);
                                    } else {
                                        ecx.superres_denom = SUPERRES_NUM;
                                    }
                                    ecx.upscaled_width = frame.frame_width;
                                    frame.frame_width = (ecx.upscaled_width*SUPERRES_NUM + ecx.superres_denom/2) / ecx.superres_denom;
                                }

                                // render size
                                let render_and_frame_size_different = frame.frame_width!=frame.render_width || frame.frame_height!=frame.render_height;
                                self.coder.push_bit(out_bits, render_and_frame_size_different as u8);
                                if render_and_frame_size_different {
                                    self.coder.push_bits_with_size(out_bits, (frame.render_width-1) as u32, 16);
                                    self.coder.push_bits_with_size(out_bits, (frame.render_height-1) as u32, 16);
                                } else {
                                    frame.render_width = ecx.upscaled_width;
                                    frame.render_height = frame.frame_height;
                                }
                            }
                            if frame.force_integer_mv>0 {
                                frame.allow_high_precision_mv = false;
                            } else {
                                self.coder.push_bit(out_bits, frame.allow_high_precision_mv as u8);
                            }
                            // interpolation filter
                            self.coder.push_bit(out_bits, frame.is_filter_switchable as u8);
                            if frame.is_filter_switchable {
                                frame.interp_filter = SWITCHABLE;
                            } else {
                                self.coder.push_bits_with_size(out_bits, frame.interp_filter as u32, 2);
                            }

                            self.coder.push_bit(out_bits, frame.is_motion_mode_switchable as u8);
                            if frame.error_resilient_mode || !ecx.enable_ref_frame_mvs {
                                frame.use_ref_frame_mvs = false;
                            } else {
                                self.coder.push_bit(out_bits, frame.use_ref_frame_mvs as u8);
                            }
                        }
                    }

                    if !ecx.frame_is_intra {
                        for i in 0..REFS_PER_FRAME {
                            let ref_frame = LAST_FRAME as usize + i;
                            let hint = ecx.ref_order_hint[frame.ref_frame_idx[i] as usize];
                            ecx.order_hints[ref_frame] = hint;
                            if !ecx.enable_order_hint {
                                ecx.ref_frame_sign_bias[ref_frame] = 0;
                            } else {
                                ecx.ref_frame_sign_bias[ref_frame] = {
                                    // get relative dist
                                    if !ecx.enable_order_hint {
                                        0
                                    } else {
                                        let mut diff = hint - frame.order_hint as usize;
                                        let m = 1 << (ecx.order_hint_bits-1);
                                        diff = (diff&(m-1)) - (diff&m);
                                        diff
                                    }
                                };
                            }
                        }
                    }
                    if ecx.reduced_still_picture_header || ecx.disable_cdf_update {
                        ecx.disable_frame_end_update_cdf = true;
                    } else {
                        self.coder.push_bit(out_bits, ecx.disable_frame_end_update_cdf as u8);
                    }
                    if frame.primary_ref_frame==PRIMARY_REF_NONE as u8 {
                        // TODO
                        // init non coeff cdfs
                        // setup past independence
                    } else {
                        // TODO
                        // load cdfs( ref_frame_idx[ primary_ref_frame ] )
                        // load previous
                    }
                    if frame.use_ref_frame_mvs {
                        // TODO
                        // motion field estimation
                    }
                    // tile info
                    let mi_cols = frame.mi_cols();
                    let mi_rows = frame.mi_rows();
                    let sb_cols = if ecx.use_128x128_superblock {
                        (mi_cols+31) >> 5
                    } else { (mi_cols+15) >> 4 };
                    let sb_rows = if ecx.use_128x128_superblock {
                        (mi_rows+31) >> 5
                    } else { (mi_rows+15) >> 4 };
                    let sb_shift = if ecx.use_128x128_superblock { 5 } else { 4 };
                    let sb_size = sb_shift + 2;
                    let max_tile_width_sb = MAX_TILE_WIDTH >> sb_size;
                    let max_tile_area_sb = MAX_TILE_AREA >> (2*sb_size);
                    let tile_log2 = |blk_size: usize, target: usize| -> usize {
                        let mut k = 0;
                        while (blk_size<<k) < target {
                            k += 1;
                        }
                        k
                    };
                    let min_log2_tile_cols = tile_log2(max_tile_width_sb, sb_cols);
                    let max_log2_tile_cols = tile_log2(1, Min!(sb_cols, MAX_TILE_COLS));
                    let max_log2_tile_rows = tile_log2(1, Min!(sb_rows, MAX_TILE_ROWS));
                    let min_log2_tiles = Max!(min_log2_tile_cols, tile_log2(max_tile_area_sb, sb_rows*sb_cols));
                    self.coder.push_bit(out_bits, frame.tile_info.uniform_tile_spacing_flag as u8);
                    if frame.tile_info.uniform_tile_spacing_flag {
                        for _ in 0..(frame.tile_cols_log2()-min_log2_tile_cols) {
                            self.coder.push_bit(out_bits, 1 as u8); // increment_tile_cols_log2
                        }
                        self.coder.push_bit(out_bits, 0 as u8); // increment_tile_cols_log2
                        let tile_width_sb = (sb_cols+(1<<frame.tile_cols_log2())-1) >> frame.tile_cols_log2();
                        let mut i = 0;
                        let mut start_sb = 0;
                        while start_sb<sb_cols {
                            ecx.mi_col_starts[i] = start_sb << sb_shift;
                            i += 1;
                            start_sb += tile_width_sb
                        }
                        ecx.mi_col_starts[i] = mi_cols;
                        //frame.tile_cols() = i;

                        let min_log2_tile_rows = Max!(min_log2_tiles-frame.tile_cols_log2(), 0);
                        let max_tile_height_sb = sb_rows >> min_log2_tile_rows;
                        for _ in 0..(frame.tile_rows_log2()-min_log2_tile_rows) {
                            self.coder.push_bit(out_bits, 1 as u8); // increment_tile_rows_log2
                        }
                        self.coder.push_bit(out_bits, 0 as u8); // increment_tile_rows_log2
                        let tile_height_sb = (sb_rows+(1<<frame.tile_rows_log2())-1) >> frame.tile_rows_log2();
                        let mut i = 0;
                        start_sb = 0;
                        while start_sb<sb_rows {
                            ecx.mi_row_starts[i] = start_sb << sb_shift;
                            i += 1;
                            start_sb += tile_height_sb;
                        }
                        ecx.mi_row_starts[i] = mi_rows;
                        //ecx.tile_rows = i;
                    } else {
                        // TODO
                    }
                    if frame.tile_cols_log2()>0 || frame.tile_rows_log2()>0 {
                        self.coder.push_bits_with_size(out_bits, ecx.context_update_tile_id as u32, (frame.tile_rows_log2()+frame.tile_cols_log2()) as u8);
                        self.coder.push_bits_with_size(out_bits, (ecx.tile_size_bytes-1) as u32, 2);
                    } else {
                        ecx.context_update_tile_id = 0;
                    }

                    // quantization params
                    self.coder.push_bits_with_size(out_bits, frame.quantization_params.base_q_idx as u32, 8);
                    let mut encode_delta_q = |delta_q: isize, coder: &mut BoolCoder, out_bits: &mut Vec<u8>| -> () {
                        if delta_q != 0 {
                            coder.push_bit(out_bits, 1); // delta_coded
                            coder.encode_su(out_bits, delta_q as i64, 1+6);
                        } else {
                            coder.push_bit(out_bits, 0); // delta_coded
                        }
                    };
                    encode_delta_q(ecx.delta_q_y_dc, self.coder, out_bits);
                    if ecx.num_planes>1 {
                        if ecx.color_config.separate_uv_delta_q {
                            self.coder.push_bit(out_bits, frame.quantization_params.diff_uv_delta as u8);
                        } else {
                            frame.quantization_params.diff_uv_delta = false;
                        }
                        encode_delta_q(ecx.delta_q_u_dc, self.coder, out_bits);
                        encode_delta_q(ecx.delta_q_u_ac, self.coder, out_bits);
                        if frame.quantization_params.diff_uv_delta {
                            encode_delta_q(ecx.delta_q_v_dc, self.coder, out_bits);
                            encode_delta_q(ecx.delta_q_v_ac, self.coder, out_bits);
                        } else {
                            ecx.delta_q_v_dc = ecx.delta_q_u_dc;
                            ecx.delta_q_v_ac = ecx.delta_q_u_ac;
                        }
                    } else {
                        ecx.delta_q_u_dc = 0;
                        ecx.delta_q_u_ac = 0;
                        ecx.delta_q_v_dc = 0;
                        ecx.delta_q_v_ac = 0;
                    }
                    self.coder.push_bit(out_bits, frame.quantization_params.using_qmatrix as u8);
                    if frame.quantization_params.using_qmatrix {
                        self.coder.push_bits_with_size(out_bits, frame.quantization_params.qm_y as u32, 4);
                        self.coder.push_bits_with_size(out_bits, frame.quantization_params.qm_u as u32, 4);
                        if !ecx.color_config.separate_uv_delta_q {
                            frame.quantization_params.qm_v = frame.quantization_params.qm_u;
                        } else {
                            self.coder.push_bits_with_size(out_bits, frame.quantization_params.qm_v as u32, 4);
                        }
                    }

                    // segmentation params
                    self.coder.push_bit(out_bits, frame.segmentation_params.segmentation_enabled as u8);
                    if frame.segmentation_params.segmentation_enabled {
                        if frame.primary_ref_frame==PRIMARY_REF_NONE as u8 {
                            frame.segmentation_params.segmentation_update_map = true;
                            frame.segmentation_params.segmentation_temporal_update = false;
                            frame.segmentation_params.segmentation_update_data = true;
                        } else {
                            self.coder.push_bit(out_bits, frame.segmentation_params.segmentation_update_map as u8);
                            if frame.segmentation_params.segmentation_update_map {
                                self.coder.push_bit(out_bits, frame.segmentation_params.segmentation_temporal_update as u8);
                            }
                            self.coder.push_bit(out_bits, frame.segmentation_params.segmentation_update_data as u8);
                        }
                        if frame.segmentation_params.segmentation_update_data {
                            for i in 0..MAX_SEGMENTS {
                                for j in 0..SEG_LVL_MAX {
                                    self.coder.push_bit(out_bits, frame.segmentation_params.feature_value[i][j].is_some() as u8);
                                    let mut clipped_value = 0;
                                    if let Some(fv) = frame.segmentation_params.feature_value[i][j] {
                                        let bits_to_read = segmentation_feature_bits[j];
                                        let limit = segmentation_feature_max[j] as isize;
                                        if segmentation_feature_signed[j] {
                                            self.coder.encode_su(out_bits, fv as i64, 1+bits_to_read as u8);
                                            clipped_value = Clip3!(-limit, limit, fv);
                                        } else {
                                            self.coder.encode_su(out_bits, fv as i64, bits_to_read as u8);
                                            clipped_value = Clip3!(0, limit, fv);
                                        }
                                    }
                                    frame.segmentation_params.feature_value[i][j] = Some(clipped_value);
                                }
                            }
                        }
                    } else {
                        for i in 0..MAX_SEGMENTS {
                            for j in 0..SEG_LVL_MAX {
                                frame.segmentation_params.feature_value[i][j] = None;
                            }
                        }
                    }
                    ecx.seg_id_pre_skip = false;
                    ecx.last_active_seg_id = 0;
                    for i in 0..MAX_SEGMENTS {
                        for j in 0..SEG_LVL_MAX {
                            if frame.segmentation_params.feature_value[i][j].is_some() {
                                ecx.last_active_seg_id = i;
                                if j >= SEG_LVL_REF_FRAME {
                                    ecx.seg_id_pre_skip = true;
                                }
                            }
                        }
                    }

                    // delta q params
                    if frame.quantization_params.base_q_idx>0 {
                        self.coder.push_bit(out_bits, frame.delta_q_params.delta_q_present as u8);
                    } else {
                        frame.delta_q_params.delta_q_present = false;
                    }
                    if frame.delta_q_params.delta_q_present {
                        self.coder.push_bit(out_bits, frame.delta_q_params.delta_q_res as u8);
                    } else {
                        frame.delta_q_params.delta_q_res = 0;
                    }

                    // delta lf params
                    if frame.delta_q_params.delta_q_present {
                        if !frame.allow_intrabc {
                            self.coder.push_bit(out_bits, frame.delta_lf_params.delta_lf_present as u8);
                        } else {
                            frame.delta_lf_params.delta_lf_present = false;
                        }
                        if frame.delta_lf_params.delta_lf_present {
                            self.coder.push_bits_with_size(out_bits, frame.delta_lf_params.delta_lf_res as u32, 2);
                            self.coder.push_bit(out_bits, frame.delta_lf_params.delta_lf_multi as u8);
                        } else {
                            frame.delta_lf_params.delta_lf_res = 0;
                            frame.delta_lf_params.delta_lf_multi = false;
                        }
                    }

                    if frame.primary_ref_frame==PRIMARY_REF_NONE as u8 {
                        // init coeff cdfs
                    } else {
                        // load previous segment ids
                    }
                    ecx.coded_lossless = true;
                    let seg_feature_active_idx = |idx: usize, feature: usize, frame: &mut Frame| -> bool {
                        return frame.segmentation_params.segmentation_enabled && frame.segmentation_params.feature_value[idx][feature].is_some()
                    };
                    let get_qindex = |ignore_delta_q: bool, segment_id: usize, ecx: &mut EncoderContext, frame: &mut Frame| -> usize {
                        if seg_feature_active_idx(segment_id, SEG_LVL_ALT_Q, frame) {
                            let data = frame.segmentation_params.feature_value[segment_id][SEG_LVL_ALT_Q].unwrap();
                            let mut qindex = frame.quantization_params.base_q_idx as isize + data;
                            if !ignore_delta_q && frame.delta_q_params.delta_q_present {
                                qindex = ecx.current_q_index as isize + data;
                            }
                            return Clip3!(0isize, 255isize, qindex) as usize;
                        } else if ignore_delta_q {
                            return ecx.current_q_index;
                        } else {
                            return frame.quantization_params.base_q_idx as usize;
                        }
                    };
                    for segment_id in 0..MAX_SEGMENTS {
                        let qindex = get_qindex(true, segment_id, &mut ecx, frame);
                        ecx.lossless_array[segment_id] = qindex==0 && ecx.delta_q_y_dc==0 &&
                                                         ecx.delta_q_u_ac==0 && ecx.delta_q_u_dc==0 &&
                                                         ecx.delta_q_v_ac==0 && ecx.delta_q_v_dc==0;
                        if !ecx.lossless_array[segment_id] {
                            ecx.coded_lossless = false;
                        }
                        if frame.quantization_params.using_qmatrix {
                            if ecx.lossless_array[segment_id] {
                                ecx.seg_q_m_level[0][segment_id] = 15;
                                ecx.seg_q_m_level[1][segment_id] = 15;
                                ecx.seg_q_m_level[2][segment_id] = 15;
                            } else {
                                ecx.seg_q_m_level[0][segment_id] = frame.quantization_params.qm_y as usize;
                                ecx.seg_q_m_level[1][segment_id] = frame.quantization_params.qm_u as usize;
                                ecx.seg_q_m_level[2][segment_id] = frame.quantization_params.qm_v as usize;
                            }
                        }
                    }
                    ecx.all_lossless = ecx.coded_lossless && frame.frame_width==ecx.upscaled_width;

                    // loop filter params
                    if ecx.coded_lossless || frame.allow_intrabc {
                        frame.loop_filter_params.loop_filter_level[0] = 0;
                        frame.loop_filter_params.loop_filter_level[1] = 0;
                        frame.loop_filter_params.loop_filter_ref_deltas[INTRA_FRAME as usize] = 1;
                        frame.loop_filter_params.loop_filter_ref_deltas[LAST_FRAME as usize] = 0;
                        frame.loop_filter_params.loop_filter_ref_deltas[LAST2_FRAME as usize] = 0;
                        frame.loop_filter_params.loop_filter_ref_deltas[LAST3_FRAME as usize] = 0;
                        frame.loop_filter_params.loop_filter_ref_deltas[BWDREF_FRAME as usize] = 0;
                        frame.loop_filter_params.loop_filter_ref_deltas[GOLDEN_FRAME as usize] = -1;
                        frame.loop_filter_params.loop_filter_ref_deltas[ALTREF_FRAME as usize] = -1;
                        frame.loop_filter_params.loop_filter_ref_deltas[ALTREF2_FRAME as usize] = -1;
                        for i in 0..2 {
                            frame.loop_filter_params.loop_filter_mode_deltas[i] = 0;
                        }
                    } else {
                        self.coder.push_bits_with_size(out_bits, frame.loop_filter_params.loop_filter_level[0] as u32, 6);
                        self.coder.push_bits_with_size(out_bits, frame.loop_filter_params.loop_filter_level[1] as u32, 6);
                        if ecx.num_planes>1 {
                            if frame.loop_filter_params.loop_filter_level[0]>0 || frame.loop_filter_params.loop_filter_level[1]>0 {
                                self.coder.push_bits_with_size(out_bits, frame.loop_filter_params.loop_filter_level[2] as u32, 6);
                                self.coder.push_bits_with_size(out_bits, frame.loop_filter_params.loop_filter_level[3] as u32, 6);
                            }
                        }
                        self.coder.push_bits_with_size(out_bits, frame.loop_filter_params.loop_filter_sharpness as u32, 3);
                        self.coder.push_bit(out_bits, frame.loop_filter_params.loop_filter_delta_enabled as u8);
                        if frame.loop_filter_params.loop_filter_delta_enabled {
                            self.coder.push_bit(out_bits, frame.loop_filter_params.loop_filter_delta_update as u8);
                            if frame.loop_filter_params.loop_filter_delta_update {
                                for i in 0..TOTAL_REFS_PER_FRAME {
                                    self.coder.push_bit(out_bits, frame.loop_filter_params.update_ref_delta[i] as u8);
                                    if frame.loop_filter_params.update_ref_delta[i] {
                                        self.coder.encode_su(out_bits, frame.loop_filter_params.loop_filter_ref_deltas[i] as i64, 1+6);
                                    }
                                }
                                for i in 0..2 {
                                    self.coder.push_bit(out_bits, frame.loop_filter_params.update_mode_delta[i] as u8);
                                    if frame.loop_filter_params.update_ref_delta[i] {
                                        self.coder.encode_su(out_bits, frame.loop_filter_params.loop_filter_mode_deltas[i] as i64, 1+6);
                                    }
                                }
                            }
                        }
                    }

                    // cdef params
                    if ecx.coded_lossless || frame.allow_intrabc {
                        frame.cdef_params.cdef_bits = 0;
                        frame.cdef_params.cdef_y_pri_strength[0] = 0;
                        frame.cdef_params.cdef_y_sec_strength[0] = 0;
                        frame.cdef_params.cdef_uv_pri_strength[0] = 0;
                        frame.cdef_params.cdef_uv_sec_strength[0] = 0;
                        frame.cdef_params.cdef_damping = 3;
                    } else {
                        self.coder.push_bits_with_size(out_bits, (frame.cdef_params.cdef_damping-3) as u32, 2);
                        self.coder.push_bits_with_size(out_bits, frame.cdef_params.cdef_bits as u32, 2);
                        for i in 0..(1<<frame.cdef_params.cdef_bits) {
                            self.coder.push_bits_with_size(out_bits, frame.cdef_params.cdef_y_pri_strength[i] as u32, 4);
                            self.coder.push_bits_with_size(out_bits, frame.cdef_params.cdef_y_sec_strength[i] as u32, 2);
                            if frame.cdef_params.cdef_y_sec_strength[i]==3 {
                                frame.cdef_params.cdef_y_sec_strength[i] += 1;
                            }
                            if ecx.num_planes>1 {
                                self.coder.push_bits_with_size(out_bits, frame.cdef_params.cdef_uv_pri_strength[i] as u32, 4);
                                self.coder.push_bits_with_size(out_bits, frame.cdef_params.cdef_uv_sec_strength[i] as u32, 2);
                                if frame.cdef_params.cdef_uv_sec_strength[i]==3 {
                                    frame.cdef_params.cdef_uv_sec_strength[i] += 1;
                                }
                            }
                        }
                    }

                    // lr params
                    if ecx.all_lossless || frame.allow_intrabc || ecx.enable_restoration {
                        ecx.frame_restoration_type[0] = RESTORE_NONE;
                        ecx.frame_restoration_type[1] = RESTORE_NONE;
                        ecx.frame_restoration_type[2] = RESTORE_NONE;
                        ecx.uses_lr = false;
                    } else {
                        ecx.uses_lr = false;
                        ecx.uses_chroma_lr = false;
                        for i in 0..ecx.num_planes {
                            self.coder.push_bits_with_size(out_bits, frame.lr_params.lr_type as u32, 2);
                            ecx.frame_restoration_type[i] = remap_lr_type[frame.lr_params.lr_type as usize];
                            if ecx.frame_restoration_type[i] != RESTORE_NONE {
                                ecx.uses_lr = true;
                                if i>0 {
                                    ecx.uses_chroma_lr = true;
                                }
                            }
                        }
                        if ecx.uses_lr {
                            if ecx.use_128x128_superblock {
                                self.coder.push_bit(out_bits, (frame.lr_params.lr_unit_shift-1) as u8);
                            } else {
                                self.coder.push_bit(out_bits, frame.lr_params.lr_unit_shift as u8);
                                if frame.lr_params.lr_unit_shift>=1 {
                                    frame.lr_params.lr_unit_extra_shift = frame.lr_params.lr_unit_shift - 1;
                                    self.coder.push_bit(out_bits, frame.lr_params.lr_unit_extra_shift as u8);
                                }
                            }
                            ecx.loop_restoration_size[0] = RESTORATION_TILESIZE_MAX >> (2-frame.lr_params.lr_unit_shift);
                            if ecx.color_config.subsampling_x && ecx.color_config.subsampling_y && ecx.uses_chroma_lr {
                                self.coder.push_bit(out_bits, frame.lr_params.lr_uv_shift as u8);
                            } else {
                                frame.lr_params.lr_uv_shift = 0;
                            }
                            ecx.loop_restoration_size[1] >>= frame.lr_params.lr_uv_shift;
                            ecx.loop_restoration_size[2] >>= frame.lr_params.lr_uv_shift;
                        }
                    }

                    // tx mode
                    if ecx.coded_lossless {
                        ecx.tx_mode = ONLY_4X4;
                    } else {
                        self.coder.push_bit(out_bits, frame.tx_mode_select as u8);
                        if frame.tx_mode_select {
                            ecx.tx_mode = TX_MODE_SELECT;
                        } else {
                            ecx.tx_mode = TX_MODE_LARGEST;
                        }
                    }

                    // frame reference mode
                    if ecx.frame_is_intra {
                        frame.reference_select = false;
                    } else {
                        self.coder.push_bit(out_bits, frame.reference_select as u8);
                    }

                    // skip mode params
                    if ecx.frame_is_intra || !frame.reference_select || !ecx.enable_order_hint {
                        ecx.skip_mode_allowed = false;
                    } else {
                        let mut forward_idx: isize = -1;
                        let mut backward_idx: isize = -1;
                        let mut forward_hint = 0; // FIXME
                        let mut backward_hint = 0; // FIXME
                        let get_relative_dist = |a: isize, b: isize, ecx: &EncoderContext| -> isize {
                            if !ecx.enable_order_hint {
                                0
                            } else {
                                let mut diff = a - b;
                                let m = 1 << (ecx.order_hint_bits-1);
                                diff = (diff&(m-1)) - (diff&m);
                                diff
                            }
                        };
                        for i in 0..REFS_PER_FRAME {
                            let ref_hint = ecx.ref_order_hint[frame.ref_frame_idx[i] as usize];
                            if get_relative_dist(ref_hint as isize, frame.order_hint as isize, &ecx) < 0 {
                                if forward_idx<0 || get_relative_dist(ref_hint as isize, forward_hint as isize, &ecx)>0 {
                                    forward_idx = i as isize;
                                    forward_hint = ref_hint;
                                }
                            } else if get_relative_dist(ref_hint as isize, frame.order_hint as isize, &ecx) > 0 {
                                if backward_idx<0 || get_relative_dist(ref_hint as isize, backward_hint as isize, &ecx) < 0 {
                                    backward_idx = i as isize;
                                    backward_hint = ref_hint;
                                }
                            }
                        }
                        if forward_idx < 0 {
                            ecx.skip_mode_allowed = false;
                        } else if backward_idx >= 0 {
                            ecx.skip_mode_allowed = true;
                            ecx.skip_mode_frame[0] = (LAST_FRAME as isize + Min!(forward_idx, backward_idx)) as usize;
                            ecx.skip_mode_frame[1] = (LAST_FRAME as isize + Max!(forward_idx, backward_idx)) as usize;
                        } else {
                            let mut second_forward_idx: isize = -1;
                            let mut second_forward_hint = 0;
                            for i in 0..REFS_PER_FRAME {
                                let ref_hint = ecx.ref_order_hint[frame.ref_frame_idx[i] as usize];
                                if get_relative_dist(ref_hint as isize, forward_hint as isize, &ecx) < 0 {
                                    if second_forward_idx<0 || get_relative_dist(ref_hint as isize, second_forward_hint as isize, &ecx) > 0 {
                                        second_forward_idx = i as isize;
                                        second_forward_hint = ref_hint;
                                    }
                                }
                            }
                            if second_forward_idx < 0 {
                                ecx.skip_mode_allowed = false;
                            } else {
                                ecx.skip_mode_allowed = true;
                                ecx.skip_mode_frame[0] = (LAST_FRAME as isize + Min!(forward_idx, second_forward_idx)) as usize;
                                ecx.skip_mode_frame[1] = (LAST_FRAME as isize + Max!(forward_idx, second_forward_idx)) as usize;
                            }
                        }
                    }
                    if ecx.skip_mode_allowed {
                        self.coder.push_bit(out_bits, frame.skip_mode_present as u8);
                    } else {
                        frame.skip_mode_present = false;
                    }

                    if ecx.frame_is_intra || frame.error_resilient_mode || !ecx.enable_warped_motion {
                        frame.allow_warped_motion = false;
                    } else {
                        self.coder.push_bit(out_bits, frame.allow_warped_motion as u8);
                    }
                    self.coder.push_bit(out_bits, frame.reduced_tx_set as u8);

                    // global motion params
                    for rf in (LAST_FRAME as usize)..(ALTREF_FRAME as usize) {
                        ecx.gm_type[rf] = IDENTITY;
                        for i in 0..6 {
                            ecx.gm_params[rf][i] = if i%3==2 { 1<<WARPEDMODEL_PREC_BITS } else { 0 };
                        }
                    }
                    if !ecx.frame_is_intra {
                        for rf in (LAST_FRAME as usize)..(ALTREF_FRAME as usize) {
                            let mut typ;
                            self.coder.push_bit(out_bits, frame.global_motion_params.is_global[rf] as u8);
                            if frame.global_motion_params.is_global[rf] {
                                self.coder.push_bit(out_bits, frame.global_motion_params.is_rot_zoom[rf] as u8);
                                if frame.global_motion_params.is_rot_zoom[rf] {
                                    typ = ROTZOOM;
                                } else {
                                    self.coder.push_bit(out_bits, frame.global_motion_params.is_translation[rf] as u8);
                                    typ = if frame.global_motion_params.is_translation[rf] { TRANSLATION } else { AFFINE };
                                }
                            } else {
                                typ = IDENTITY;
                            }
                            ecx.gm_type[rf] = typ;

                            let mut encode_global_param = |typ: usize, rf: usize, idx: usize, ecx: &EncoderContext, coder: &mut BoolCoder| -> () {
                                let mut abs_bits = GM_ABS_ALPHA_BITS;
                                let mut prec_bits = GM_ALPHA_PREC_BITS;
                                if idx < 2 {
                                    if typ==TRANSLATION {
                                        abs_bits = GM_ABS_TRANS_ONLY_BITS - (!frame.allow_high_precision_mv) as usize;
                                        prec_bits = GM_TRANS_ONLY_PREC_BITS - (!frame.allow_high_precision_mv) as usize;
                                    } else {
                                        abs_bits = GM_ABS_TRANS_BITS;
                                        prec_bits = GM_TRANS_PREC_BITS;
                                    }
                                }
                                let prec_diff = WARPEDMODEL_PREC_BITS - prec_bits;
                                let round = if idx%3==2 { 1<<WARPEDMODEL_PREC_BITS } else { 0 };
                                let sub = if idx%3==2 { 1<<prec_bits } else { 0 };
                                let mx = 1 << abs_bits;
                                let r = (ecx.prev_gm_params[rf][idx]>>prec_diff) - sub;

                                let inverse_recenter = |r: usize, v: usize| -> usize {
                                    if v > 2*r {
                                        v
                                    } else if v&1 > 0 {
                                        r - ((v+1)>>1)
                                    } else {
                                        r + (v>>1)
                                    }
                                };
                                let mut encode_subexp = |v: usize, num_syms: usize, coder: &mut BoolCoder| -> () {
                                    let mut i = 0;
                                    let mut mk = 0;
                                    let mut k = 3;
                                    let mut tv = v;
                                    loop {
                                        let b2 = if i>0 { k+i-1 } else { k };
                                        let a = 1 << b2;
                                        if num_syms <= mk + 3*a {
                                            tv -= mk;
                                            coder.encode_ns(out_bits, tv as u64, (num_syms-mk) as u64);
                                            break;
                                        } else {
                                            if tv <= mk + (1<<b2) - 1 {
                                                coder.push_bit(out_bits, 0); // subexp_more_bits
                                                tv -= mk;
                                                coder.push_bits_with_size(out_bits, tv as u32, b2);
                                                break
                                            } else {
                                                coder.push_bit(out_bits, 1); // subexp_more_bits
                                                i += 1;
                                                mk += a;
                                            }
                                        }
                                    }
                                };
                                let mut encode_unsigned_subexp_with_ref = |v: usize, mx: usize, r: usize, coder: &mut BoolCoder| -> () {
                                    if (r<<1) <= mx {
                                        let tv = inverse_recenter(r, v);
                                        encode_subexp(tv, mx, coder);
                                    } else {
                                        let mut tv = mx - 1 - v;
                                        tv = inverse_recenter(mx-1-r, tv);
                                        encode_subexp(tv, mx-1-r, coder);
                                    }
                                };
                                let mut encode_signed_subexp_with_ref = |v: isize, low: isize, high: isize, r: isize, coder: &mut BoolCoder| -> () {
                                    let x = v - low;
                                    encode_unsigned_subexp_with_ref(x as usize, (high-low) as usize, (r-low) as usize, coder);
                                };
                                encode_signed_subexp_with_ref((ecx.gm_params[rf][idx]-round)>>prec_diff, -mx, mx+1, r, coder);
                            };

                            if typ >= ROTZOOM {
                                encode_global_param(typ, rf, 2, &ecx, self.coder);
                                encode_global_param(typ, rf, 3, &ecx, self.coder);
                                if typ==AFFINE {
                                encode_global_param(typ, rf, 4, &ecx, self.coder);
                                encode_global_param(typ, rf, 5, &ecx, self.coder);
                                } else {
                                    ecx.gm_params[rf][4] = -ecx.gm_params[rf][3];
                                    ecx.gm_params[rf][5] = -ecx.gm_params[rf][2];
                                }
                            }
                            if typ >= TRANSLATION {
                                encode_global_param(typ, rf, 0, &ecx, self.coder);
                                encode_global_param(typ, rf, 1, &ecx, self.coder);
                            }
                        }
                    }

                    // film grain params
                    // FIXME refactor
                    if ecx.film_grain_params.is_none() || (!ecx.show_frame && !ecx.showable_frame) {
                        ecx.film_grain_params = Some(FilmGrainParams::new());
                    } if ecx.film_grain_params.is_some() {
                        let tmp_fgp = ecx.ref_film_grain_params[ecx.film_grain_params.as_ref().unwrap().film_grain_params_ref_idx as usize].clone();
                        let tmp_cc = ecx.color_config.clone();
                        if let Some(ref mut fgp) = ecx.film_grain_params {
                            self.coder.push_bit(out_bits, fgp.apply_grain as u8);
                            if !fgp.apply_grain {
                                *fgp = FilmGrainParams::new();
                            } else {
                                self.coder.push_bits_with_size(out_bits, fgp.grain_seed as u32, 16);
                                if frame.frame_type==INTER_FRAME {
                                    self.coder.push_bit(out_bits, fgp.update_grain as u8);
                                } else {
                                    fgp.update_grain = true;
                                }
                                if !fgp.update_grain {
                                    self.coder.push_bits_with_size(out_bits, fgp.film_grain_params_ref_idx as u32, 3);
                                    let temp_grain_seed = fgp.grain_seed;
                                    *fgp =  tmp_fgp;
                                    fgp.grain_seed = temp_grain_seed;
                                } else {
                                    self.coder.push_bits_with_size(out_bits, fgp.num_y_points as u32, 4);
                                    for i in 0..fgp.num_y_points {
                                        self.coder.push_bits_with_size(out_bits, fgp.point_y_value[i] as u32, 8);
                                        self.coder.push_bits_with_size(out_bits, fgp.point_y_scaling[i] as u32, 8);
                                    }
                                    if tmp_cc.mono_chrome {
                                        fgp.chroma_scaling_from_luma = false;
                                    } else {
                                        self.coder.push_bit(out_bits, fgp.chroma_scaling_from_luma as u8);
                                    }
                                    if tmp_cc.mono_chrome || fgp.chroma_scaling_from_luma ||
                                       (tmp_cc.subsampling_x && tmp_cc.subsampling_y && fgp.num_y_points==0) {
                                        fgp.num_cb_points = 0;
                                        fgp.num_cr_points = 0;
                                    } else {
                                        self.coder.push_bits_with_size(out_bits, fgp.num_cb_points as u32, 4);
                                        for i in 0..fgp.num_cb_points {
                                            self.coder.push_bits_with_size(out_bits, fgp.point_cb_value[i] as u32, 8);
                                            self.coder.push_bits_with_size(out_bits, fgp.point_cb_scaling[i] as u32, 8);
                                        }
                                        self.coder.push_bits_with_size(out_bits, fgp.num_cr_points as u32, 4);
                                        for i in 0..fgp.num_cr_points {
                                            self.coder.push_bits_with_size(out_bits, fgp.point_cr_value[i] as u32, 8);
                                            self.coder.push_bits_with_size(out_bits, fgp.point_cr_scaling[i] as u32, 8);
                                        }
                                    }
                                    self.coder.push_bits_with_size(out_bits, (fgp.grain_scaling-8) as u32, 2);
                                    self.coder.push_bits_with_size(out_bits, fgp.ar_coeff_lag as u32, 2);
                                    let num_pos_luma = 2*fgp.ar_coeff_lag*(fgp.ar_coeff_lag+1);
                                    let mut num_pos_chroma;
                                    if fgp.num_y_points>0 {
                                        num_pos_chroma = num_pos_luma + 1;
                                        for i in 0..num_pos_luma {
                                            self.coder.push_bits_with_size(out_bits, (fgp.ar_coeffs_y[i]+128) as u32, 8);
                                        }
                                    } else {
                                        num_pos_chroma = num_pos_luma;
                                    }
                                    if fgp.chroma_scaling_from_luma || fgp.num_cb_points>0 {
                                        for i in 0..num_pos_chroma {
                                            self.coder.push_bits_with_size(out_bits, (fgp.ar_coeffs_cb[i]+128) as u32, 8);
                                        }
                                    }
                                    if fgp.chroma_scaling_from_luma || fgp.num_cr_points>0 {
                                        for i in 0..num_pos_chroma {
                                            self.coder.push_bits_with_size(out_bits, (fgp.ar_coeffs_cr[i]+128) as u32, 8);
                                        }
                                    }
                                    self.coder.push_bits_with_size(out_bits, (fgp.ar_coeff_shift-6) as u32, 2);
                                    self.coder.push_bits_with_size(out_bits, fgp.grain_scale_shift as u32, 2);
                                    if fgp.num_cb_points>0 {
                                        self.coder.push_bits_with_size(out_bits, fgp.cb_mult as u32, 8);
                                        self.coder.push_bits_with_size(out_bits, fgp.cb_luma_mult as u32, 8);
                                        self.coder.push_bits_with_size(out_bits, fgp.cb_offset as u32, 9);
                                    }
                                    if fgp.num_cr_points>0 {
                                        self.coder.push_bits_with_size(out_bits, fgp.cr_mult as u32, 8);
                                        self.coder.push_bits_with_size(out_bits, fgp.cr_luma_mult as u32, 8);
                                        self.coder.push_bits_with_size(out_bits, fgp.cr_offset as u32, 9);
                                    }
                                    self.coder.push_bit(out_bits, fgp.overlap_flag as u8);
                                    self.coder.push_bit(out_bits, fgp.clip_to_restricted_range as u8);
                                }
                            }
                        }
                    }


                    if ecx.show_existing_frame {
                        // frame wrapup
                        // TODO

                        ecx.seen_frame_header = false;
                    } else {
                        self.tile_num = 0;
                        ecx.seen_frame_header = true;
                    }
                }
            },
            &mut OBU_REDUNDANT_FRAME_HEADER => {

            },
            &mut OBU_TILE_GROUP(ref mut frame) => {
                let mut ecx = self.ectx.borrow_mut();
                let num_tiles = frame.tile_col_starts.len() * frame.tile_row_starts.len();
                if num_tiles>1 {
                    if frame.tile_groups.len()>0 {
                        self.coder.push_bit(out_bits, 1); // tile_start_and_end_present_flag
                    } else {
                        self.coder.push_bit(out_bits, 0); // tile_start_and_end_present_flag
                    }
                }
                if num_tiles==1 || frame.tile_groups.len()==0 {
                    ecx.tg_start = 0;
                    ecx.tg_end = num_tiles - 1;
                } else {
                    let tile_bits = frame.tile_cols_log2() + frame.tile_rows_log2();
                    self.coder.push_bits_with_size(out_bits, ecx.tg_start as u32, tile_bits as u8);
                    self.coder.push_bits_with_size(out_bits, ecx.tg_end as u32, tile_bits as u8);
                }
                self.coder.push_bits_align(out_bits);

                for n in ecx.tg_start..(ecx.tg_end+1) {
                    ecx.tile_num = n;
                    let last_tile = ecx.tile_num == ecx.tg_end;
                    let mut tile_bits: Vec<u8> = vec![];
                    {
                        //let mut tile_encoder = TileEncoder::new(&self.ectx.clone(), self.writer, self.coder);
                        let mut tile_encoder = TileEncoder::new(&self.ectx.clone(), self.coder);
                        tile_encoder.encode_tile(frame, ecx.tile_num, &mut tile_bits);
                    }
                    let tile_size_bytes = (tile_bits.len()+7) / 8;
                    if !last_tile {
                        self.coder.encode_le(out_bits, (tile_size_bytes-1) as u64, ecx.tile_size_bytes as u8);
                    }
                }

                if ecx.tg_end == frame.num_tiles()-1 {
                    if !ecx.disable_frame_end_update_cdf {
                        // frame_end_update_cdf
                        // TODO
                    }
                    // decode frame wrapup
                    // TODO
                    ecx.seen_frame_header = false;
                }
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
                        seconds: secs,
                        minutes: mins,
                        hours: hrs,
                        time_offset_length: tol,
                        time_offset_value: tov,
                    } => {
                        self.coder.encode_leb128(out_bits, 3);
                        self.coder.push_bits_with_size(out_bits, ct as u32, 5);
                        if let (Some(s), Some(m), Some(h)) = (secs, mins, hrs) {
                            self.coder.push_bit(out_bits, 1);
                            self.coder.push_bits_with_size(out_bits, s as u32, 6);
                            self.coder.push_bits_with_size(out_bits, m as u32, 6);
                            self.coder.push_bits_with_size(out_bits, h as u32, 5);
                        } else {
                            self.coder.push_bit(out_bits, 0);
                            if let Some(s) = secs {
                                self.coder.push_bit(out_bits, 1);
                                self.coder.push_bits_with_size(out_bits, s as u32, 6);
                                if let Some(m) = mins {
                                    self.coder.push_bit(out_bits, 1);
                                    self.coder.push_bits_with_size(out_bits, m as u32, 6);
                                    if let Some(h) = hrs {
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
                    //_ => { assert!(false); }
                };
            },
            &mut OBU_FRAME { ref mut frame } => {

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