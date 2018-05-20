use std::fs::File;
use std::io::{self, BufWriter, Write};

pub struct BinaryWriter<'a> {
    output: Box<Write + 'a>,
}

impl<'a> BinaryWriter<'a> {
    pub fn standard(stdout: &'a io::Stdout) -> BinaryWriter<'a> {
        BinaryWriter {
            output: Box::new(stdout.lock()),
        }
    }

    pub fn file(path: &str) -> io::Result<BinaryWriter<'a>> {
        File::create(path).map(|file| BinaryWriter {
            output: Box::new(io::BufWriter::new(file)),
        })
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