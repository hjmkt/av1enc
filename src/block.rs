#![allow(dead_code)]

use constants::*;
use constants::FrameRestorationType::*;

pub struct LRUnit {
    pub use_wiener: bool,
    pub use_sgrproj: bool,
    pub restoration_type: FrameRestorationType,
    pub lr_sgr_set: usize,
    pub lr_wiener: [[isize; 3]; 2],
    pub lr_sgr_xqd: [isize; 2],
}

pub struct Superblock {
    pub mi_size: BlockSize,
    pub x: usize,
    pub y: usize,
    pub cdef_idx: usize,
    pub lr_units: Vec<Vec<Vec<LRUnit>>>,
}