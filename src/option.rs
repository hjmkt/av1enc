extern crate num;

use num::rational::Ratio;

struct Config {
    frame_width: i16,
    frame_height: i16,
    framerate: Ratio<i16>,
}

struct CommandOption {
    name: str,
    short_name: Option(str),
    usage: str,
    default: Option<str>,
    parser: fn(str, &mut Config) -> Result<i32, str>
}

impl CommandOption {
    fn parse_tuple<T: FromStr>(arg: str, sep: char, len: i32, tuple: &mut [T]) -> Result<i32, str> {
        let v = arg.split(sep).collect();
        if v.len() == len {
            for i in 0..(v.len()-1) {
                match T::from_str(v[i]) {
                    Ok(res) => tuple[i] = res,
                    Err(e) => return Err("Invalid argument")
                }
            }
            Ok(0)
        }
        else{
            Err("Invalid argument")
        }
    }
}

static options: Vec<CommandOption> = [
    CommandOption {
        name: "size",
        short_name: "s",
        usage: "video resolution (WxH)",
        default: None,
        parser: | arg: str, config: &mut Config | -> Result<i32, str> {
            CommandOption::parse_tuple::<u16>(arg, 'x', 2, [config.frame_width, config.frame_height])
        }
    },
    CommandOption {
        name: "framerate",
        short_name: "r",
        usage: "video frame rate (N/D, or N for N/1)",
        default: None,
        parser: | arg: str, config: &mut Config | -> Result<i32, str> {
            match u16::from_str(arg) {
                Ok(n) => {
                    config.framerate = Ratio(n, 1);
                    Ok(0)
                }
                Err(_) => {
                    CommandOption::parse_tuple::<u16>(arg, '/', 2, [config.frame_width, config.frame_height])
                }
            }
        }
    }
];