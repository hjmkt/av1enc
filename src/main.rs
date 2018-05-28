#![allow(unused_imports)]
#![allow(unused_variables)]
#![allow(unused_mut)]

#[macro_use]
extern crate lazy_static;
extern crate colored;
mod option;
use option::*;
use colored::*;
use std::process;
mod binary_reader;
mod binary_writer;
use binary_reader::BinaryReader;
use binary_writer::BinaryWriter;
use std::io::{self, Read, Write};

mod constants;
mod cdf;
mod sequence_header;
mod obu;
mod frame_header;
mod util;
mod frame;
mod bool_coder;


fn main() {

    println!(
        "{} version {} Copyright (c) 2018-2018 {}",
        env!("CARGO_PKG_NAME").to_string(),
        env!("CARGO_PKG_VERSION").to_string(),
        env!("CARGO_PKG_AUTHORS").split(':').collect::<Vec<&str>>().join(", ").to_string()
    );

    let mut config: Config = Config::default();

    let args = std::env::args().skip(1).collect::<Vec<String>>();
    if args.len() == 0 { return; }

    CommandOption::init_config(&mut config);
    if let Err(e) = CommandOption::parse(args, &mut config) {
        println!("{}: failed to parse command-line arguments: {}", "error".red(), e);
        process::exit(0);
    }
    if let Err(e) = CommandOption::validate_config(&config) {
        println!("{}: invalid configuration: {}", "error".red(), e);
        process::exit(0);
    }

    // initialize binary reader
    let stdin = io::stdin();
    let mut reader = if config.input == Some("stdin".to_string()) {
        BinaryReader::standard(&stdin)
    }
    else {
        match BinaryReader::file(&config.input.unwrap()) {
            Ok(f) => f,
            Err(e) => {
                println!("{}: failed to open input file: {}", "error".red(), e);
                process::exit(0);
            }
        }
    };

    // initialize binary writer
    let stdout = io::stdout();
    let mut writer = if config.output == Some("stdout".to_string()) {
        BinaryWriter::standard(&stdout)
    }
    else {
        match BinaryWriter::file(&config.output.unwrap()) {
            Ok(f) => f,
            Err(e) => {
                println!("{}: failed to open output file: {}", "error".red(), e);
                process::exit(0);
            }
        }
    };

    
}
