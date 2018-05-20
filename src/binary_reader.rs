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