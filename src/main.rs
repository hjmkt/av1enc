#[macro_use]
extern crate lazy_static;
extern crate colored;
mod option;
use option::*;
use colored::*;
use std::process;

fn main() {

    println!("{} version {} Copyright (c) 2018-2018 {}", env!("CARGO_PKG_NAME").to_string(), env!("CARGO_PKG_VERSION").to_string(), env!("CARGO_PKG_AUTHORS"));

    let mut config: Config = Config::default();

    let args = std::env::args().skip(1).collect::<Vec<String>>();
    if args.len() > 0 {
        CommandOption::init_config(&mut config);
        if let Err(e) = CommandOption::parse(args, &mut config) {
            println!("{}: failed to parse command-line arguments: {}", "error".red(), e);
            process::exit(0);
        }
        if let Err(e) = CommandOption::validate_config(&config) {
            println!("{}: invalid configuration: {}", "error".red(), e);
            process::exit(0);
        }
    }
}
