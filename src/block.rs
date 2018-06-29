#![allow(dead_code)]

use constants::*;
use constants::Partition::*;
use constants::RefFrame::*;
use constants::CompMode::*;
use constants::CompRefType::*;
use constants::YMode::*;
use common::*;

pub struct LRUnit {
    pub use_wiener: bool,
    pub use_sgrproj: bool,
    pub restoration_type: FrameRestorationType,
    pub lr_sgr_set: usize,
    pub lr_wiener: [[isize; 3]; 2],
    pub lr_sgr_xqd: [isize; 2],
}

pub struct PartitionTree {
    pub part_type: Partition,
    pub partitions: Vec<PartitionTree>,
    pub bsize: BlockSize,
    pub segment_id: usize,
    pub skip: bool,
    pub cdef_idx: isize,
    pub q_index: usize,
    pub delta_lf: [isize; FRAME_LF_COUNT],
    pub is_inter: bool,
    pub ref_frame: [RefFrame; 2],
    pub comp_mode: CompMode,
    pub comp_ref_type: CompRefType,
    pub compound_mode: YMode,
    pub new_mv: bool,
    pub zero_mv: bool,
    pub ref_mv: bool,
    pub drl_mode: [usize; 8], // FIXME
    pub use_intrabc: bool,
    pub mvs: [MotionVector; 2],
}

impl PartitionTree {
    pub fn new(bsize: BlockSize) -> PartitionTree { PartitionTree {
        part_type: PARTITION_NONE,
        partitions: vec![],
        bsize: bsize,
        segment_id: 0,
        skip: false,
        cdef_idx: 0,
        q_index: 0,
        delta_lf: [0; FRAME_LF_COUNT],
        is_inter: false,
        ref_frame: [LAST_FRAME; 2],
        comp_mode: SINGLE_REFERENCE,
        comp_ref_type: UNIDIR_COMP_REFERENCE,
        compound_mode: NEAREST_NEARESTMV,
        new_mv: false,
        zero_mv: false,
        ref_mv: false,
        drl_mode: [0; 8],
        use_intrabc: false,
        mvs: [MotionVector::new(0, 0); 2],
    }}

    pub fn clone(&self) -> PartitionTree { PartitionTree {
        part_type: self.part_type,
        partitions: {
            let mut tmp = vec![];
            for p in self.partitions.iter() {
                tmp.push(p.clone());
            }
            tmp
        },
        bsize: self.bsize,
        segment_id: self.segment_id,
        skip: self.skip,
        cdef_idx: self.cdef_idx,
        q_index: self.q_index,
        delta_lf: self.delta_lf,
        is_inter: self.is_inter,
        ref_frame: self.ref_frame,
        comp_mode: self.comp_mode,
        comp_ref_type: self.comp_ref_type,
        compound_mode: self.compound_mode,
        new_mv: self.new_mv,
        zero_mv: self.zero_mv,
        ref_mv: self.ref_mv,
        drl_mode: self.drl_mode,
        use_intrabc: self.use_intrabc,
        mvs: self.mvs,
    }}
}

pub struct Superblock {
    pub mi_size: BlockSize,
    pub x: usize,
    pub y: usize,
    pub cdef_idx: usize,
    pub lr_units: Vec<Vec<Vec<LRUnit>>>,
    pub part_tree: PartitionTree,
}

impl Superblock {
    pub fn new(bsize: BlockSize, x: usize, y: usize) -> Superblock { Superblock {
        mi_size: bsize,
        x: x,
        y: y,
        cdef_idx: 0,
        lr_units: vec![],
        part_tree: PartitionTree::new(bsize),
    }}
}