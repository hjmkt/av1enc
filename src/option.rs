extern crate num;

use self::num::rational::Ratio;
use std::str::FromStr;
use std::process;

pub struct Config {
    pub input: Option<String>,
    pub output: Option<String>,
    pub frame_width: Option<u16>,
    pub frame_height: Option<u16>,
    pub framerate: Option<Ratio<u16>>,
    pub frames: Option<u32>,
}

impl Default for Config {
    fn default() -> Config { Config {
        input: None,
        output: None,
        frame_width: None,
        frame_height: None,
        framerate: None,
        frames: None,
    }}
}

pub struct CommandOption {
    name: &'static str,
    short_name: Option<&'static str>,
    usage: &'static str,
    default: Option<&'static str>,
    parser: fn(&Vec<String>, &mut Config) -> Result<Vec<String>, String>,
    validator: fn(&Config) -> Result<i32, String>
}

impl<'a> CommandOption {
    fn parse_tuple<T: FromStr>(arg: &'a str, sep: char, len: usize, tuple: &mut [&mut Option<T>]) -> Result<i32, String> {
        let v : Vec<&str> = arg.split(sep).collect();
        if v.len() == len {
            for i in 0..v.len() {
                if let Ok(res) = T::from_str(v[i]) { *tuple[i] = Some(res); }
                else { return Err("wrong format".to_string()); }
            }
            Ok(0)
        }
        else{ Err("wrong format".to_string()) }
    }
}

macro_rules! parser_string {
    ($target: ident) => {
        | args: &Vec<String>, config: &mut Config | -> Result<Vec<String>, String> {
            if args.len() == 0 { return Err("".to_string()); }
            config.$target = Some(args[0].clone());
            Ok(args[1..].to_vec())
        }
    };
}

macro_rules! parser_num {
    ($type: ident, $target: ident) => {
        | args: &Vec<String>, config: &mut Config | -> Result<Vec<String>, String> {
            if args.len() == 0 { return Err("".to_string()); }
            if let Err(e) = CommandOption::parse_tuple::<$type>(&args[0], ' ', 1, &mut [&mut config.$target]) { Err(e) }
            else { Ok(args[1..].to_vec()) }
        }
    };
}

macro_rules! validator_skip {
    () => {
        | _config | -> Result<i32, String>  { Ok(0) }
    };
}

macro_rules! validator_default {
    ($target: ident) => {
        | config: &Config | -> Result<i32, String>  {
            if let Some(_) = config.$target { Ok(0) }
            else { Err(stringify!($target not specified.).to_string()) }
        }
    };
}

lazy_static! {
    static ref OPTIONS : Vec<CommandOption> = vec![
        CommandOption {
            name: "help",
            short_name: Some("h"),
            usage: "show help message",
            default: None,
            parser: | _args: &Vec<String>, _config: &mut Config | -> Result<Vec<String>, String> {
                CommandOption::show_help_and_exit();
                Ok(vec![])
            },
            validator: validator_skip!()
        },
        CommandOption {
            name: "input",
            short_name: Some("i"),
            usage: "input path or 'stdin' for standard input",
            default: None,
            parser: parser_string!(input),
            validator: validator_default!(input)
        },
        CommandOption {
            name: "output",
            short_name: Some("o"),
            usage: "output path or 'stdout' for standard output",
            default: None,
            parser: parser_string!(output),
            validator: validator_default!(output)
        },
        CommandOption {
            name: "frames",
            short_name: Some("fr"),
            usage: "number of maximum frames to encode",
            default: Some("4294967295"),
            parser: parser_num!(u32, frames),
            validator: validator_default!(frames)
        },
        CommandOption {
            name: "size",
            short_name: Some("s"),
            usage: "video resolution (WxH)",
            default: None,
            parser: | args: &Vec<String>, config: &mut Config | -> Result<Vec<String>, String> {
                if args.len() == 0 { return Err("".to_string()); }
                if let Err(e) = CommandOption::parse_tuple::<u16>(&args[0], 'x', 2, &mut [&mut config.frame_width, &mut config.frame_height]) { Err(e) }
                else { Ok(args[1..].to_vec()) }
            },
            validator: | config: &Config | -> Result<i32, String>  {
                if let (Some(w), Some(h)) = (config.frame_width, config.frame_height) {
                    if w>0 && w<8192 && h>0 && h<8192 { Ok(0) }
                    else { Err("invalid video resolution specified.".to_string()) }
                }
                else { Err("video resolution not specified.".to_string()) }
            }
        },
        CommandOption {
            name: "framerate",
            short_name: Some("r"),
            usage: "video frame rate (N/D, or N for N/1)",
            default: Some("30/1"),
            parser: | args: &Vec<String>, config: &mut Config | -> Result<Vec<String>, String> {
                if args.len() == 0 { return Err("".to_string()); }
                if let Ok(n) = u16::from_str(&args[0]) {
                    config.framerate = Some(Ratio::new(n, 1));
                    Ok(args[1..].to_vec())
                }
                else {
                    let (mut denom, mut numer) = (None, None);
                    if let Err(e) = CommandOption::parse_tuple::<u16>(&args[0], '/', 2, &mut [&mut numer, &mut denom]) { Err(e) }
                    else {
                        config.framerate = Some(Ratio::new(numer.unwrap(), denom.unwrap()));
                        Ok(vec![])
                    }
                }
            },
            validator: | config: &Config | -> Result<i32, String>  {
                if let Some(r) = config.framerate {
                    if r.numer() > &0 && r.denom() > &0 { Ok(0) }
                    else { Err("invalid video framerate specified.".to_string()) }
                }
                else { Err("video framerate not specified.".to_string()) }
            }
        }
    ];
}

impl CommandOption {
    fn set_option_by_name(name: &str, args: Vec<String>, config: &mut Config) -> Result<Vec<String>, String> {
        for opt in OPTIONS.iter() {
            if opt.name == name { return match (opt.parser)(&args, config) {
                Ok(v) => Ok(v),
                Err(e) => if args.len() > 0 {
                    Err(format!("invalid argument '{}' for option '--{}': {}", args[0], name, e))
                }
                else {
                    Err(format!("argument not specified for option '--{}': '{}'", name, e))
                }
            }}
        }
        Err(format!("option '--{}' not found.", name))
    }

    fn set_option_by_short_name(name: &str, args: Vec<String>, config: &mut Config) -> Result<Vec<String>, String> {
        for opt in OPTIONS.iter() {
            if opt.short_name == Some(name) { return match (opt.parser)(&args, config) {
                Ok(v) => Ok(v),
                Err(e) => if args.len() > 0 {
                    Err(format!("invalid argument '{}' for option '-{}': {}", args[0], name, e))
                }
                else {
                    Err(format!("argument not specified for option '-{}'", name))
                }
            }}
        }
        Err(format!("option '-{}' not found.", name))
    }

    pub fn init_config(config: &mut Config) {
        for opt in OPTIONS.iter() {
            if let Some(o) = opt.default {
                if let Err(e) = (opt.parser)(&[o.to_string()].to_vec(), config) { assert!(false, e) }
            }
        }
    }

    pub fn validate_config(config: &Config) -> Result<i32, String> {
        for opt in OPTIONS.iter() {
            if let Err(e) = (opt.validator)(config) { return Err(e) }
        }
        Ok(0)
    }

    fn show_help_and_exit() {
        for opt in OPTIONS.iter() {
            if let Some(sn) = opt.short_name { println!("--{} (-{}): {}", opt.name, sn, opt.usage) }
            else { println!("--{}: {}", opt.name, opt.usage); }
        }
        process::exit(0);
    }

    pub fn parse(args: Vec<String>, config: &mut Config) -> Result<Vec<String>, String> {
        if args.len() == 0 { Ok(vec![]) }
        else {
            let arg = &args[0];
            if &arg[0..1] == "-" {
                match {
                    if &arg[1..2] == "-" { CommandOption::set_option_by_name(&arg[2..], args[1..].to_vec(), config) }
                    else { CommandOption::set_option_by_short_name(&arg[1..], args[1..].to_vec(), config) }
                } {
                    Ok(v) => CommandOption::parse(v, config),
                    Err(e) => Err(e)
                }
            }
            else{ return Err("option name must start with '-'.".to_string()) }
        }
    }
}