// Copyright 2020 Johannes Köster.
// Licensed under the GNU GPLv3 license (https://opensource.org/licenses/GPL-3.0)
// This file may not be copied, modified, or distributed
// except according to those terms.

use bio_types::genome;

/// Modeling of breakends.
/// 
/// # Examples
/// 
/// * `t[p[`: p=pos, t=ref_allele+seq, revcomp=false, side=RightOfPos
/// * `t]p]`: p=pos, t=ref_allele+seq, revcomp=true, side=
pub struct Breakend {
    pos: genome::Position,
    mate_pos: genome::Position,
    ref_allele: u8,
    seq: Vec<u8>,
    side_of_pos: Side,
    side_of_seq: 
    revcomp: bool,
}

pub enum Side {
    LeftOfPos,
    RightOfPos,
}