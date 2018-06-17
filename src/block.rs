#![allow(dead_code)]

use constants::*;
use constants::Partition::*;

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
}

impl PartitionTree {
    pub fn new(bsize: BlockSize) -> PartitionTree { PartitionTree {
        part_type: PARTITION_NONE,
        partitions: vec![],
        bsize: bsize,
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