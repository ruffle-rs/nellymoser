mod tables;

use bitstream_io::{read::BitRead, BitReader, LittleEndian};
use rustdct::DCTplanner;
use std::io::{Cursor, Read};
use tables::*;

const NELLY_BLOCK_LEN: usize = 64;
const NELLY_HEADER_BITS: usize = 116;
const NELLY_DETAIL_BITS: i32 = 198;
const NELLY_BUF_LEN: usize = 128;
const NELLY_FILL_LEN: usize = 124;
const NELLY_BIT_CAP: i16 = 6;
const NELLY_BASE_OFF: i32 = 4228;
const NELLY_BASE_SHIFT: i16 = 19;
const NELLY_SAMPLES: usize = NELLY_BUF_LEN * 2;

pub struct Decoder<R: Read> {
    reader: R,
    sample_rate: u32,
    state: [f32; NELLY_BUF_LEN],
    planner: DCTplanner<f32>,
    cur_frame: [f32; NELLY_SAMPLES],
    cur_sample: usize,
}

impl<R: Read> Decoder<R> {
    pub fn new(reader: R, sample_rate: u32) -> Self {
        Self {
            reader,
            sample_rate,
            state: [0.0; NELLY_BUF_LEN],
            planner: DCTplanner::new(),
            cur_frame: [0f32; NELLY_SAMPLES], // TODO: make uninitialized?
            cur_sample: 0,
        }
    }

    pub fn sample_rate(&self) -> u32 {
        self.sample_rate
    }

    fn next_frame(&mut self) -> Option<()> {
        let mut block = [0u8; NELLY_BLOCK_LEN];
        self.reader.read_exact(&mut block).ok()?;
        self.cur_frame = self.decode_block(&block);
        self.cur_sample = 0;
        Some(())
    }

    fn decode_block(&mut self, block: &[u8; NELLY_BLOCK_LEN]) -> [f32; NELLY_SAMPLES] {
        let mut buf = [0f32; NELLY_BUF_LEN];
        let mut pows = [0f32; NELLY_BUF_LEN];
        {
            let mut reader = BitReader::endian(Cursor::new(&block), LittleEndian);
            let mut val = NELLY_INIT_TABLE[reader.read::<u8>(6).unwrap() as usize] as f32;
            let mut ptr: usize = 0;
            for (i, x) in NELLY_BAND_SIZES_TABLE.iter().enumerate() {
                if i > 0 {
                    val += NELLY_DELTA_TABLE[reader.read::<u8>(5).unwrap() as usize] as f32;
                }

                let pval = (val / 2048.0).exp2();
                for _ in 0..*x {
                    buf[ptr] = val;
                    pows[ptr] = pval;
                    ptr += 1;
                }
            }
        }

        let bits = {
            let mut max = buf.iter().fold(0, |a, &b| a.max(b as i32));
            let mut shift = headroom(&mut max) - 16;

            let mut sbuf = [0i16; NELLY_BUF_LEN];
            for i in 0..NELLY_FILL_LEN {
                sbuf[i] = signed_shift(buf[i] as i32, shift as i32) as i16;
                sbuf[i] = ((3 * sbuf[i] as i32) >> 2) as i16;
            }

            let mut sum: i32 = sbuf.iter().map(|&s| s as i32).sum();

            shift += 11;
            let shift_saved = shift;
            sum -= NELLY_DETAIL_BITS << shift;
            shift += headroom(&mut sum);
            let mut small_off = (NELLY_BASE_OFF * (sum >> 16)) >> 15;
            shift = shift_saved - (NELLY_BASE_SHIFT + shift - 31);

            small_off = signed_shift(small_off, shift as i32);

            let mut bitsum = sum_bits(sbuf, shift_saved, small_off as i16);
            if bitsum != NELLY_DETAIL_BITS {
                let mut off = bitsum - NELLY_DETAIL_BITS;
                shift = 0;
                while off.abs() <= 16383 {
                    off *= 2;
                    shift += 1;
                }

                off = (off * NELLY_BASE_OFF) >> 15;
                shift = shift_saved - (NELLY_BASE_SHIFT + shift - 15);

                off = signed_shift(off, shift as i32);

                let mut last_off = small_off;
                let mut last_bitsum = bitsum;
                let mut last_j = 0;
                for j in 1..20 {
                    last_off = small_off;
                    small_off += off;
                    last_bitsum = bitsum;
                    last_j = j;

                    bitsum = sum_bits(sbuf, shift_saved, small_off as i16);

                    if (bitsum - NELLY_DETAIL_BITS) * (last_bitsum - NELLY_DETAIL_BITS) <= 0 {
                        break;
                    }
                }

                let mut big_off;
                let mut big_bitsum;
                let mut small_bitsum;
                if bitsum > NELLY_DETAIL_BITS {
                    big_off = small_off;
                    small_off = last_off;
                    big_bitsum = bitsum;
                    small_bitsum = last_bitsum;
                } else {
                    big_off = last_off;
                    big_bitsum = last_bitsum;
                    small_bitsum = bitsum;
                }

                while bitsum != NELLY_DETAIL_BITS && last_j <= 19 {
                    off = (big_off + small_off) >> 1;
                    bitsum = sum_bits(sbuf, shift_saved, off as i16);
                    if bitsum > NELLY_DETAIL_BITS {
                        big_off = off;
                        big_bitsum = bitsum;
                    } else {
                        small_off = off;
                        small_bitsum = bitsum;
                    }
                    last_j += 1;
                }

                if (big_bitsum - NELLY_DETAIL_BITS).abs()
                    >= (small_bitsum - NELLY_DETAIL_BITS).abs()
                {
                    bitsum = small_bitsum;
                } else {
                    small_off = big_off;
                    bitsum = big_bitsum;
                }
            }

            let mut bits = [0i32; NELLY_BUF_LEN];
            for i in 0..NELLY_FILL_LEN {
                let mut tmp = sbuf[i] as i32 - small_off;
                tmp = ((tmp >> (shift_saved - 1)) + 1) >> 1;
                bits[i] = if tmp < 0 {
                    0
                } else if tmp > NELLY_BIT_CAP as i32 {
                    NELLY_BIT_CAP as i32
                } else {
                    tmp
                };
            }

            if bitsum > NELLY_DETAIL_BITS {
                let mut i = 0;
                let mut tmp = 0;
                while tmp < NELLY_DETAIL_BITS {
                    tmp += bits[i];
                    i += 1;
                }

                bits[i - 1] -= tmp - NELLY_DETAIL_BITS;

                while i < NELLY_FILL_LEN {
                    bits[i] = 0;
                    i += 1;
                }
            }

            bits
        };

        let mut samples = [0f32; NELLY_SAMPLES];
        for i in 0..2 {
            let mut reader = BitReader::endian(Cursor::new(&block), LittleEndian);
            reader
                .skip(NELLY_HEADER_BITS as u32 + i * NELLY_DETAIL_BITS as u32)
                .unwrap();

            let input: Vec<f32> = (0..NELLY_BUF_LEN).map(|j| if j >= NELLY_FILL_LEN {
                0.0
            } else if bits[j] <= 0 {
                std::f32::consts::FRAC_1_SQRT_2
            } else {
                let v = reader.read::<u8>(bits[j] as u32).unwrap();
                NELLY_DEQUANTIZATION_TABLE[((1 << bits[j]) - 1 + v) as usize] as f32
            } * pows[j]).collect();

            let slice = &mut samples[i as usize * NELLY_BUF_LEN..(i as usize + 1) * NELLY_BUF_LEN];
            for (i, x) in self.state.iter_mut().enumerate() {
                slice[i] = *x;
                *x = 0.0;
            }

            let mdct = self.planner.plan_mdct(NELLY_BUF_LEN, window);
            mdct.process_imdct_split(&input, slice, &mut self.state);
        }

        samples
    }
}

#[inline]
fn signed_shift(i: i32, shift: i32) -> i32 {
    if shift > 0 {
        i << shift
    } else {
        i >> -shift
    }
}

fn sum_bits(buf: [i16; NELLY_BUF_LEN], shift: i16, off: i16) -> i32 {
    buf[0..NELLY_FILL_LEN].iter().fold(0i32, |ret, &i| {
        let b = i as i32 - off as i32;
        let b = ((b >> (shift - 1)) + 1) >> 1;
        ret + if b < 0 {
            0
        } else if b > NELLY_BIT_CAP as i32 {
            NELLY_BIT_CAP as i32
        } else {
            b
        }
    })
}

fn headroom(la: &mut i32) -> i16 {
    if *la == 0 {
        return 31;
    }

    let l = la.abs().leading_zeros() as i16 - 1;
    *la *= 1 << l;
    l
}

fn window(len: usize) -> Vec<f32> {
    (0..len)
        .map(|i| ((i as f32 + 0.5) / 128.0 * std::f32::consts::FRAC_PI_2).sin() / 8.0)
        .collect()
}

impl<R: AsRef<[u8]>> Decoder<Cursor<R>> {
    #[inline]
    pub fn reset(&mut self) {
        self.reader.set_position(0);
        self.state.iter_mut().for_each(|x| *x = 0.0);
    }
}

impl<R: Read> Iterator for Decoder<R> {
    type Item = f32;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cur_sample >= NELLY_SAMPLES {
            self.next_frame()?;
        }

        let sample = self.cur_frame[self.cur_sample];
        self.cur_sample += 1;
        Some(sample)
    }
}
