#![allow(dead_code)]
#![allow(unused_imports)]

use std::fs::File;
use std::io::{self, BufRead, Read};

pub struct BinaryReader<'a> {
    input: Box<BufRead + 'a>,
    buffer: u8
}

impl<'a> BinaryReader<'a> {
    pub fn standard(stdin: &'a io::Stdin) -> BinaryReader<'a> {
        BinaryReader {
            input: Box::new(stdin.lock()),
            buffer: 0
        }
    }

    pub fn file(path: &str) -> io::Result<BinaryReader<'a>> {
        File::open(path).map(|file| BinaryReader {
            input: Box::new(io::BufReader::new(file)),
            buffer: 0
        })
    }

    pub fn read_to_vec<T: From<u8>>(&mut self, v: &mut Vec<T>) -> io::Result<usize> {
        let len = v.len();
        let mut tmp: Vec<u8> = vec![0; len];
        match self.input.read(&mut tmp[..]) {
            Ok(s) => {
                for i in 0..len {
                    v[i] = T::from(tmp[i]);
                }
                Ok(s)
            },
            Err(e) => Err(e)
        }
    }
}

impl<'a> Read for BinaryReader<'a> {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        self.input.read(buf)
    }
}

impl<'a> BufRead for BinaryReader<'a> {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        self.input.fill_buf()
    }

    fn consume(&mut self, amt: usize) {
        self.input.consume(amt);
    }
}