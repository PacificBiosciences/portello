use rust_htslib::bam::record::Cigar;
use rust_vc_utils::cigar::{clean_up_cigar_edge_indels, compress_cigar, update_ref_and_read_pos};

#[derive(Default)]
struct CigarBlockInfo {
    is_in_indel_block: bool,

    block_ref_start: i64,
    block_read_start: usize,

    block_del_size: u32,
    block_ins_size: u32,
}

impl CigarBlockInfo {
    fn _add_indel(&mut self, ref_pos: i64, read_pos: usize) {
        if !self.is_in_indel_block {
            self.is_in_indel_block = true;
            self.block_ref_start = ref_pos;
            self.block_read_start = read_pos;
        }
    }

    fn add_del(&mut self, len: u32, ref_pos: i64, read_pos: usize) {
        self._add_indel(ref_pos, read_pos);
        self.block_del_size += len;
    }

    fn add_ins(&mut self, len: u32, ref_pos: i64, read_pos: usize) {
        self._add_indel(ref_pos, read_pos);
        self.block_ins_size += len;
    }

    /// Return the simplified replacement for this cigar cluster
    fn end_indel(&mut self, ref_seq: &[u8], read_seq: &[u8]) -> Vec<Cigar> {
        let mut ret = Vec::new();
        if self.is_in_indel_block {
            self.is_in_indel_block = false;

            // Take care of simple indels
            match (self.block_del_size, self.block_ins_size) {
                (0, 0) => (),
                (0, len) => ret.push(Cigar::Ins(len)),
                (len, 0) => ret.push(Cigar::Del(len)),
                (1, 1) => {
                    // Don't even look at the sequence for this case, even if this is a SNP it changes edit distance 2 to 1:
                    ret.push(Cigar::Match(1))
                }
                (mut del_len, mut ins_len) => {
                    // The remaining cases benefit from looking at the sequence to find the best simplification:
                    let mut pre_match_len = 0;
                    let mut post_match_len = 0;

                    // First try to push as much insertion as possible onto the right-side match state:
                    while del_len > 0 && ins_len > 0 {
                        let last_del_base_ref_index =
                            self.block_ref_start as usize + del_len as usize - 1;
                        let last_del_base = ref_seq[last_del_base_ref_index];
                        let last_ins_base_read_index = self.block_read_start + ins_len as usize - 1;
                        let last_ins_base = read_seq[last_ins_base_read_index];
                        if last_del_base == last_ins_base {
                            del_len -= 1;
                            ins_len -= 1;
                            post_match_len += 1;
                        } else {
                            break;
                        }
                    }

                    // Check if more insertion can be pushed onto left-side match state:
                    while del_len > 0 && ins_len > 0 {
                        let first_del_base_ref_index =
                            self.block_ref_start as usize + pre_match_len as usize;
                        let first_del_base = ref_seq[first_del_base_ref_index];
                        let first_ins_base_read_index =
                            self.block_read_start + pre_match_len as usize;
                        let first_ins_base = read_seq[first_ins_base_read_index];
                        if first_del_base == first_ins_base {
                            del_len -= 1;
                            ins_len -= 1;
                            pre_match_len += 1;
                        } else {
                            break;
                        }
                    }

                    // Finally if we've simplified down to a SNP, then prefer 1 edit instead of 2
                    if del_len == 1 && ins_len == 1 {
                        del_len -= 1;
                        ins_len -= 1;
                        post_match_len += 1;
                    }

                    // Finished figuring out simplification, now return it:
                    let mut rpush = |c: Cigar| {
                        if !c.is_empty() {
                            ret.push(c);
                        }
                    };

                    rpush(Cigar::Match(pre_match_len));
                    rpush(Cigar::Ins(ins_len));
                    rpush(Cigar::Del(del_len));
                    rpush(Cigar::Match(post_match_len));
                }
            }
            self.block_ins_size = 0;
            self.block_del_size = 0;
        }
        ret
    }
}

/// For any complex cluster of multiple indel (ID) types in the cigar:
/// 1. See if any amount of insertion and deletion can be traded off to reduce total edit distance
/// 2. See if size of indel elements can be reduced.
/// 3. If simplification to one type is not possible, then standardize remaining cases to "nImD"
///
pub fn simplify_alignment_indels(
    ref_pos: i64,
    cigar: &[Cigar],
    ref_seq: &[u8],
    read_seq: &[u8],
) -> (i64, Vec<Cigar>) {
    let ignore_hard_clip = false;

    let mut ref_head_pos = ref_pos;
    let mut read_head_pos = 0;

    let mut cigar_block_info = CigarBlockInfo::default();

    let mut simple_cigar = Vec::new();

    for c in cigar.iter() {
        use Cigar::*;

        match c {
            Del(len) => {
                cigar_block_info.add_del(*len, ref_head_pos, read_head_pos);
            }
            Ins(len) => {
                cigar_block_info.add_ins(*len, ref_head_pos, read_head_pos);
            }
            _ => {
                simple_cigar.append(&mut cigar_block_info.end_indel(ref_seq, read_seq));
                simple_cigar.push(*c);
            }
        }
        update_ref_and_read_pos(c, &mut ref_head_pos, &mut read_head_pos, ignore_hard_clip);
    }
    simple_cigar.append(&mut cigar_block_info.end_indel(ref_seq, read_seq));

    let ref_pos_shift = clean_up_cigar_edge_indels(&mut simple_cigar);
    let simple_cigar = compress_cigar(&simple_cigar);
    (ref_pos + ref_pos_shift as i64, simple_cigar)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simplify_alignment_indels() {
        use Cigar::*;

        {
            // Boring match
            let ref_pos = 2;
            let cigar = vec![Match(6)];
            let ref_seq = b"XXABCCDEXX";
            let read_seq = b"ABCCDE";

            let (simple_ref_pos, simple_cigar) =
                simplify_alignment_indels(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(simple_ref_pos, ref_pos);
            assert_eq!(simple_cigar, cigar);
        }

        {
            // Boring Ins
            let ref_pos = 2;
            let cigar = vec![Match(2), Ins(1), Match(3)];
            let ref_seq = b"XXABCDEXX";
            let read_seq = b"ABCCDE";

            let (simple_ref_pos, simple_cigar) =
                simplify_alignment_indels(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(simple_ref_pos, ref_pos);
            assert_eq!(simple_cigar, cigar);
        }

        {
            // Boring Del
            let ref_pos = 2;
            let cigar = vec![Match(2), Del(1), Match(3)];
            let ref_seq = b"XXABCCDEXX";
            let read_seq = b"ABCDE";

            let (simple_ref_pos, simple_cigar) =
                simplify_alignment_indels(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(simple_ref_pos, ref_pos);
            assert_eq!(simple_cigar, vec![Match(2), Del(1), Match(3)]);
        }

        {
            // Boring InDel
            let ref_pos = 2;
            let cigar = vec![Match(2), Del(2), Ins(2), Match(3)];
            let ref_seq = b"XXABCCCDEXX";
            let read_seq = b"ABBBCDE";

            let (simple_ref_pos, simple_cigar) =
                simplify_alignment_indels(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(simple_ref_pos, ref_pos);
            assert_eq!(simple_cigar, vec![Match(2), Ins(2), Del(2), Match(3)]);
        }

        {
            // Simple Consolidation
            //
            // Test gives a choice of left or right side merge to assert that right-side is chosen first
            //
            let ref_pos = 2;
            let cigar = vec![Match(3), Ins(1), Del(2), Match(2)];
            let ref_seq = b"XXABCCCDEXX";
            let read_seq = b"ABCCDE";

            let (simple_ref_pos, simple_cigar) =
                simplify_alignment_indels(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(simple_ref_pos, ref_pos);
            assert_eq!(simple_cigar, vec![Match(3), Del(1), Match(3)]);
        }

        {
            // Left-side Consolidation
            //
            let ref_pos = 2;
            let cigar = vec![Match(3), Del(3), Ins(3), Match(1)];
            let ref_seq = b"XXABCCCDEXX";
            let read_seq = b"ABCCXXE";

            let (simple_ref_pos, simple_cigar) =
                simplify_alignment_indels(ref_pos, &cigar, ref_seq, read_seq);

            assert_eq!(simple_ref_pos, ref_pos);
            assert_eq!(simple_cigar, vec![Match(4), Ins(2), Del(2), Match(1)]);
        }
    }
}
