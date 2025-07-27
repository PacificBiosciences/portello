use rust_htslib::bam::record::Cigar;
use rust_vc_utils::cigar::{clean_up_cigar_edge_indels, compress_cigar, update_ref_and_read_pos};
use rust_vc_utils::{IntRange, get_indel_breakend_homology_info};

#[derive(Default)]
struct CigarBlockInfo {
    match_block_size: u32,

    shift_cigar: Vec<Cigar>,
}

impl CigarBlockInfo {
    fn add_match_len(&mut self, len: u32) {
        self.match_block_size += len;
    }

    fn add_non_match(&mut self, left_shift_len: u32, cigar_seg: Option<&Cigar>) {
        let actual_left_shift_len = std::cmp::min(self.match_block_size, left_shift_len);
        let shifted_match_block_size = self.match_block_size - actual_left_shift_len;
        if shifted_match_block_size > 0 {
            self.shift_cigar
                .push(Cigar::Match(shifted_match_block_size));
        }
        if let Some(x) = cigar_seg {
            self.shift_cigar.push(*x);
        }
        self.match_block_size = actual_left_shift_len;
    }
}

/// Does not preserve X/= match states in left-shifted output
///
/// Does not handle indel clusters like "10M2D2I10M"
///
pub fn left_shift_alignment(
    ref_pos: i64,
    cigar: &[Cigar],
    ref_seq: &[u8],
    read_seq: &[u8],
) -> (i64, Vec<Cigar>) {
    let ignore_hard_clip = false;

    let mut ref_head_pos = ref_pos;
    let mut read_head_pos = 0;

    let mut cigar_block_info = CigarBlockInfo::default();

    for c in cigar.iter() {
        use Cigar::*;

        match c {
            Match(len) | Equal(len) | Diff(len) => {
                cigar_block_info.add_match_len(*len);
            }
            Del(len) => {
                let ref_range = IntRange::from_pair(ref_head_pos, ref_head_pos + *len as i64);
                let read_range = IntRange::from_pair(read_head_pos as i64, read_head_pos as i64);
                let (range, _) =
                    get_indel_breakend_homology_info(ref_seq, &ref_range, read_seq, &read_range);
                let left_shift_len = std::cmp::max(0, -range.start) as u32;
                cigar_block_info.add_non_match(left_shift_len, Some(c));
            }
            Ins(len) => {
                let ref_range = IntRange::from_pair(ref_head_pos, ref_head_pos);
                let read_range =
                    IntRange::from_pair(read_head_pos as i64, read_head_pos as i64 + *len as i64);
                let (range, _) =
                    get_indel_breakend_homology_info(ref_seq, &ref_range, read_seq, &read_range);
                let left_shift_len = std::cmp::max(0, -range.start) as u32;
                cigar_block_info.add_non_match(left_shift_len, Some(c));
            }
            _ => {
                cigar_block_info.add_non_match(0, Some(c));
            }
        }
        update_ref_and_read_pos(c, &mut ref_head_pos, &mut read_head_pos, ignore_hard_clip);
    }
    cigar_block_info.add_non_match(0, None);

    let ref_pos_shift = clean_up_cigar_edge_indels(&mut cigar_block_info.shift_cigar);
    let shift_cigar = compress_cigar(&cigar_block_info.shift_cigar);
    (ref_pos + ref_pos_shift as i64, shift_cigar)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_left_shift_alignment() {
        use Cigar::*;

        {
            // Boring match
            let ref_pos = 2;
            let cigar = vec![Match(6)];
            let ref_seq = b"XXABCCDEXX";
            let read_seq = b"ABCCDE";

            let (shift_ref_pos, shift_cigar) =
                left_shift_alignment(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(shift_ref_pos, ref_pos);
            assert_eq!(shift_cigar, cigar);
        }

        {
            // Boring soft-clip
            let ref_pos = 4;
            let cigar = vec![SoftClip(2), Match(2), SoftClip(2)];
            let ref_seq = b"XXABCCDEXX";
            let read_seq = b"ABCCDE";

            let (shift_ref_pos, shift_cigar) =
                left_shift_alignment(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(shift_ref_pos, ref_pos);
            assert_eq!(shift_cigar, cigar);
        }

        {
            // Ins
            let ref_pos = 2;
            let cigar = vec![Match(3), Ins(1), Match(2)];
            let ref_seq = b"XXABCDEXX";
            let read_seq = b"ABCCDE";

            let (shift_ref_pos, shift_cigar) =
                left_shift_alignment(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(shift_ref_pos, ref_pos);
            assert_eq!(shift_cigar, vec![Match(2), Ins(1), Match(3)]);
        }

        {
            // Ins to edge
            let ref_pos = 4;
            let cigar = vec![Match(1), Ins(1), Match(2)];
            let ref_seq = b"XXABCDEXX";
            let read_seq = b"CCDE";

            let (shift_ref_pos, shift_cigar) =
                left_shift_alignment(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(shift_ref_pos, ref_pos);
            assert_eq!(shift_cigar, vec![SoftClip(1), Match(3)]);
        }

        {
            // Del
            let ref_pos = 2;
            let cigar = vec![Match(3), Del(1), Match(2)];
            let ref_seq = b"XXABCCDEXX";
            let read_seq = b"ABCDE";

            let (shift_ref_pos, shift_cigar) =
                left_shift_alignment(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(shift_ref_pos, ref_pos);
            assert_eq!(shift_cigar, vec![Match(2), Del(1), Match(3)]);
        }

        {
            // Del to edge
            let ref_pos = 4;
            let cigar = vec![Match(1), Del(1), Match(2)];
            let ref_seq = b"XXABCCDEXX";
            let read_seq = b"CDE";

            let (shift_ref_pos, shift_cigar) =
                left_shift_alignment(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(shift_ref_pos, 5);
            assert_eq!(shift_cigar, vec![Match(3)]);
        }

        {
            // Multi-indel
            let ref_pos = 2;
            let cigar = vec![Match(3), Ins(1), Match(2), Del(1), Match(1)];
            let ref_seq = b"XXABCDEEFXX";
            let read_seq = b"ABCCDEF";

            let (shift_ref_pos, shift_cigar) =
                left_shift_alignment(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(shift_ref_pos, ref_pos);
            assert_eq!(
                shift_cigar,
                vec![Match(2), Ins(1), Match(2), Del(1), Match(2)]
            );
        }
    }
}
