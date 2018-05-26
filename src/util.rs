#![allow(unused_macros)]

macro_rules! Abs {
    ($x: expr) => ( if x>=0 { $x } else { -$x } );
}

macro_rules! Clip3 {
    ($x: expr, $y: expr, $z: expr) => ( if $z<$x { $x } else if $z>$y { $y } else { $z } );
}

macro_rules! Clip1 {
    ($x: expr) => ( Clip3!(0, (1<<8)-1, $x) );
}

macro_rules! Min {
    ($x: expr, $y: expr) => ( if $x<$y { $x } else { $y } );
}

macro_rules! Max {
    ($x: expr, $y: expr) => ( if $x>$y { $x } else { $y } );
}

macro_rules! Round2 {
    ($x: expr, $n: expr) => ( if $n==0 { $x } else { ($x + (1 << ($n - 1))) >> $n } );
}

macro_rules! Round2Signed {
    ($x: expr, $n: expr) => ( if $x>=0 { Round2!($x, $n) } else { Round2!(-$x, $n) } );
}

macro_rules! FloorLog2 {
    ($x: expr) => (
        let mut s = 0;
        while $x != 0 {
            $x = ($x >> 1);
            s++;
        }
        s - 1
    );
}

macro_rules! CeilLog2 {
    ($x: expr) => (
        if (x < 2) { 0 }
        else {
            let mut i = 1;
            let mut p = 2;
            while p < $x {
                i++;
                p = p << 1;
            }
            i
        }
    );
}

pub fn msb16(n: u16) -> u8{
    static LUT: [u8; 16] = [0, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3];
    let mut bit = 0;
    let mut m0 = n >> 8;
    if m0>0 { bit += 8; }
    else { m0 = n; }
    let mut m1 = m0 >> 4;
    if m1>0 { bit += 4; }
    else { m1 = m0; }
    bit + LUT[m1 as usize]
}