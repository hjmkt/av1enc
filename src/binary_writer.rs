#![allow(dead_code)]
#![allow(unused_imports)]

use std::fs::File;
use std::io::{self, BufWriter, Write};

pub struct BinaryWriter<'a> {
    output: Box<Write + 'a>,
    buf: u8,
    index: usize,
}

impl<'a> BinaryWriter<'a> {
    pub fn standard(stdout: &'a io::Stdout) -> BinaryWriter<'a> {
        BinaryWriter {
            output: Box::new(stdout.lock()),
            buf: 0,
            index: 0,
        }
    }

    pub fn file(path: &str) -> io::Result<BinaryWriter<'a>> {
        File::create(path).map(|file| BinaryWriter {
            output: Box::new(io::BufWriter::new(file)),
            buf: 0,
            index: 0,
        })
    }

    pub fn write_bit(&mut self, bit: u8) {
        self.buf = (self.buf<<1) + bit;
        self.index += 1;
        if self.index==8 {
            self.index = 0;
            let tmp = [self.buf];
            if let Err(e) = self.write(&tmp) {
                assert!(false, e);
            }
            self.buf = 0;
        }
    }

    pub fn write_bits(&mut self, bits: &Vec<u8>) {
        for bit in bits {
            self.write_bit(*bit);
        }
    }

    pub fn byte_align(&mut self) {
        let rem = if self.index>0 { 8 - self.index } else { 0 };
        for i in 0..rem {
            self.write_bit(0);
        }
    }
}

impl<'a> Write for BinaryWriter<'a> {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.output.write(buf)
    }

    fn flush(&mut self) -> io::Result<()> {
        self.output.flush()
    }
}