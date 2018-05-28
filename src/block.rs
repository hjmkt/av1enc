#![allow(dead_code)]

use constants::*;

pub struct Block {
    pub mi_size: BlockSize,
    pub is_inter: bool,
    pub skip: bool,
    pub tx_size: TxSize,
}