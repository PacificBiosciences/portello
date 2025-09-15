use std::cmp::Ordering;

use log::info;
use rust_htslib::bam::record::CigarString;
use rust_vc_utils::cigar::{
    clip_alignment_read_edges, get_cigar_read_offset, get_gap_compressed_identity_no_align_match,
    get_read_clip_positions,
};
use rust_vc_utils::{
    GenomeSegment, IntRange, SeqOrderSplitReadSegment, drop_true,
    get_read_segment_to_ref_pos_tree_map,
};

use super::{
    AllContigMappingInfo, ContigMappingInfo, ContigMappingSegmentInfo, print_split_read_summary,
};

fn get_seg_gap_compressed_identity(
    qname: &str,
    seg: &SeqOrderSplitReadSegment,
    seg_isec_seq_order_read_range: &IntRange,
) -> f64 {
    let ignore_hard_clip = false;

    let read_len = get_cigar_read_offset(&seg.cigar, ignore_hard_clip);

    // 1. Convert read intersect range from sequencing order to alignment order
    let seg_isec_read_range = if seg.is_fwd_strand {
        seg_isec_seq_order_read_range.clone()
    } else {
        seg_isec_seq_order_read_range.get_reverse_range(read_len as i64)
    };

    // 2. Get seg alignment trimmed to just the intersection range
    let (clipped_cigar, _) = clip_alignment_read_edges(
        &seg.cigar,
        seg_isec_read_range.start as usize,
        read_len - seg_isec_read_range.end as usize,
    );

    match get_gap_compressed_identity_no_align_match(&clipped_cigar) {
        Ok(x) => x,
        Err(msg) => {
            panic!(
                "Error generating gap-compressed identity for overlapping split read segment in assembly contig '{qname}': {msg}"
            );
        }
    }
}

/// Update the given split read segment to remove the intersection range from another split read
///
/// Return true if the split read should be eliminated entirely
fn clip_seg_isec_range(
    seg: &mut SeqOrderSplitReadSegment,
    seg_isec_seq_order_read_range: &IntRange,
) -> bool {
    let ignore_hard_clip = false;

    // Determine if we're clipping a prefix (suffix) in both sequence order and alignment order:
    let is_clipping_seq_order_prefix =
        { seg_isec_seq_order_read_range.start == seg.seq_order_read_start as i64 };
    let is_clipping_prefix = is_clipping_seq_order_prefix ^ (!seg.is_fwd_strand);

    // Convert read intersect range from sequencing order to alignment order:
    let read_len = get_cigar_read_offset(&seg.cigar, ignore_hard_clip);
    let mut seg_isec_read_range = if seg.is_fwd_strand {
        seg_isec_seq_order_read_range.clone()
    } else {
        seg_isec_seq_order_read_range.get_reverse_range(read_len as i64)
    };

    // Translate seg_isec_read_range into a cigar alignment clipping pattern:
    let (min_left_clip, min_right_clip) = if is_clipping_prefix {
        (seg_isec_read_range.end as usize, 0)
    } else {
        (0, read_len - seg_isec_read_range.start as usize)
    };
    let (shifted_cigar, ref_pos_shift) =
        clip_alignment_read_edges(&seg.cigar, min_left_clip, min_right_clip);
    seg.cigar = CigarString(shifted_cigar);
    seg.pos += ref_pos_shift;

    // Update seg_isec_read_range with the actual clip amount, which could exceed the min_left_clip/min_right_clip used above
    let (left_read_pos, right_read_pos, _) = get_read_clip_positions(&seg.cigar, ignore_hard_clip);

    // This indicates that the entire read is clipped out:
    if left_read_pos >= right_read_pos {
        return true;
    }

    if is_clipping_prefix {
        seg_isec_read_range.end = left_read_pos as i64;
    } else {
        seg_isec_read_range.start = right_read_pos as i64;
    };

    // Convert read intersect range from alignment order to sequencing order:
    let seg_isec_seq_order_read_range = if seg.is_fwd_strand {
        seg_isec_read_range.clone()
    } else {
        seg_isec_read_range.get_reverse_range(read_len as i64)
    };

    // Update seq_order_read_{start/end} based on seg_isec_seq_order_read_range
    if is_clipping_seq_order_prefix {
        seg.seq_order_read_start = seg_isec_seq_order_read_range.end as usize;
    } else {
        seg.seq_order_read_end = seg_isec_seq_order_read_range.start as usize;
    }
    false
}

/// Update the given split read segment info struct to remove the intersection range from another split read
///
/// Return true if the split read should be eliminated entirely
fn clip_seg_info_isec_range(
    seg_info: &mut ContigMappingSegmentInfo,
    seg_isec_seq_order_read_range: &IntRange,
) -> bool {
    let ignore_hard_clip = false;

    let eliminated = clip_seg_isec_range(
        &mut seg_info.seq_order_segment,
        seg_isec_seq_order_read_range,
    );

    if eliminated {
        true
    } else {
        let seg = &seg_info.seq_order_segment;
        seg_info.contig_to_ref_map =
            get_read_segment_to_ref_pos_tree_map(seg.pos, &seg.cigar, ignore_hard_clip);
        false
    }
}

/// Determine if two split read segments overlap, and if so, which segment should be clipped, and over what read-range
///
/// If an intersection is found, this returns a 2-tuple of:
/// 1. The intersection region of the two split alignments, expressed in contig fwd-strand coordinates
/// 2. Index of the split alignment to be clipped
///
fn get_seg_clip_info(
    contig_mapping_info: &ContigMappingInfo,
    split_index1: usize,
    split_index2: usize,
    debug: bool,
) -> Option<(IntRange, usize)> {
    let ordered_contig_segment_info = &contig_mapping_info.ordered_contig_segment_info;

    let seg1 = &ordered_contig_segment_info[split_index1].seq_order_segment;
    let seg2 = &ordered_contig_segment_info[split_index2].seq_order_segment;

    // determine if there's query overlap between the two split read alignments:
    if seg1.seq_order_read_end <= seg2.seq_order_read_start {
        return None;
    }

    let seg_isec_seq_order_read_range = IntRange::from_pair(
        seg2.seq_order_read_start as i64,
        seg1.seq_order_read_end as i64,
    );

    if debug {
        eprintln!(
            "Index {split_index1} and {split_index2} intersect over read_range {seg_isec_seq_order_read_range:?}!"
        );
    }

    let seg1_gci = get_seg_gap_compressed_identity(
        &contig_mapping_info.qname,
        seg1,
        &seg_isec_seq_order_read_range,
    );

    let seg2_gci = get_seg_gap_compressed_identity(
        &contig_mapping_info.qname,
        seg2,
        &seg_isec_seq_order_read_range,
    );

    let clip_seg1 = matches!(
        seg2_gci
            .partial_cmp(&seg1_gci)
            .unwrap()
            .then(seg2.mapq.cmp(&seg1.mapq)),
        Ordering::Greater
    );

    if debug {
        eprintln!("seg1 gci: {seg1_gci} mapq: {}", seg1.mapq);
        eprintln!("seg2 gci: {seg2_gci} mapq: {}", seg2.mapq);
        eprintln!("clip_seg1: {clip_seg1}");
    }

    let clip_split_index = if clip_seg1 {
        split_index1
    } else {
        split_index2
    };

    Some((seg_isec_seq_order_read_range, clip_split_index))
}

/// Clip repeated matches from contig split alignments
///
/// Repeated matches are regions where the same contig sequence range is mapped to multiple reference regions at a split
/// alignment junction.
///
/// This is a simple clipping routine which will choose to clip all of one side or the other, for each repeated match.
/// It could be improved in future to find some optimal clipping point along the repeated match range.
///
pub fn clip_repeated_contig_matches(result: &mut AllContigMappingInfo) {
    let debug = false;

    info!("Clipping repeated contig matches at split alignment segment boundaries");

    let mut segments_clipped = 0;

    for contig_mapping_info in result
        .iter_mut()
        .filter(|x| !x.ordered_contig_segment_info.is_empty())
    {
        if debug {
            print_split_read_summary(contig_mapping_info);
        }
        let split_read_count = contig_mapping_info.ordered_contig_segment_info.len();
        let mut split_read_eliminated = vec![false; split_read_count];
        for split_index1 in 0..split_read_count {
            for split_index2 in (split_index1 + 1)..split_read_count {
                if split_read_eliminated[split_index1] || split_read_eliminated[split_index2] {
                    continue;
                }

                let (seg_isec_seq_order_read_range, clip_split_index) =
                    match get_seg_clip_info(contig_mapping_info, split_index1, split_index2, debug)
                    {
                        Some(x) => x,
                        None => break,
                    };

                if debug {
                    let seg1 = &contig_mapping_info.ordered_contig_segment_info[split_index1]
                        .seq_order_segment;
                    let seg2 = &contig_mapping_info.ordered_contig_segment_info[split_index2]
                        .seq_order_segment;
                    let isec_segment = GenomeSegment {
                        chrom_index: seg1.chrom_index,
                        range: seg_isec_seq_order_read_range.clone(),
                    };
                    eprintln!("found clip for overlap segment: {isec_segment:?}");
                    eprintln!("seg1 index: {split_index1} info: {}", seg1.short_display());
                    eprintln!("seg2 index: {split_index2} info: {}", seg2.short_display());
                    eprintln!("clip_split_index: {clip_split_index}");
                }

                let eliminated = clip_seg_info_isec_range(
                    &mut contig_mapping_info.ordered_contig_segment_info[clip_split_index],
                    &seg_isec_seq_order_read_range,
                );

                if eliminated {
                    split_read_eliminated[clip_split_index] = true;
                }

                if debug {
                    let mut seg1_info = contig_mapping_info.ordered_contig_segment_info
                        [split_index1]
                        .seq_order_segment
                        .short_display();
                    let mut seg2_info = contig_mapping_info.ordered_contig_segment_info
                        [split_index2]
                        .seq_order_segment
                        .short_display();

                    if eliminated {
                        if clip_split_index == split_index1 {
                            seg1_info = "ELIMINATED".to_string();
                        } else {
                            seg2_info = "ELIMINATED".to_string();
                        }
                    }

                    eprintln!("After clipping:");
                    eprintln!("seg1 index: {split_index1} info: {}", &seg1_info);
                    eprintln!("seg2 index: {split_index2} info: {}", &seg2_info);
                    eprintln!("clip_split_index: {clip_split_index}");
                }

                segments_clipped += 1;
            }
        }

        // remove eliminated split reads:
        drop_true(
            &mut contig_mapping_info.ordered_contig_segment_info,
            &split_read_eliminated,
        );
    }

    info!("Clipped {segments_clipped} repeated contig match regions");
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::Cigar;

    #[test]
    fn test_clip_seg_isec_range() {
        use Cigar::*;

        // Simplest test scenario
        {
            let seg_isec_seq_order_read_range = IntRange::from_pair(10, 20);
            let mut seg = SeqOrderSplitReadSegment {
                seq_order_read_start: 0,
                seq_order_read_end: 20,
                chrom_index: 0,
                pos: 100,
                is_fwd_strand: true,
                cigar: CigarString(vec![Match(20)]),
                mapq: 60,
                from_primary_bam_record: false,
            };
            clip_seg_isec_range(&mut seg, &seg_isec_seq_order_read_range);

            assert_eq!(seg.pos, 100);
            assert_eq!(seg.cigar, CigarString(vec![Match(10), SoftClip(10)]));
            assert_eq!(seg.seq_order_read_start, 0);
            assert_eq!(seg.seq_order_read_end, 10);
        }

        // Simple scenario but reverse strand:
        {
            let seg_isec_seq_order_read_range = IntRange::from_pair(10, 20);
            let mut seg = SeqOrderSplitReadSegment {
                seq_order_read_start: 0,
                seq_order_read_end: 20,
                chrom_index: 0,
                pos: 100,
                is_fwd_strand: false,
                cigar: CigarString(vec![Match(20)]),
                mapq: 60,
                from_primary_bam_record: false,
            };
            clip_seg_isec_range(&mut seg, &seg_isec_seq_order_read_range);

            assert_eq!(seg.pos, 110);
            assert_eq!(seg.cigar, CigarString(vec![SoftClip(10), Match(10)]));
            assert_eq!(seg.seq_order_read_start, 0);
            assert_eq!(seg.seq_order_read_end, 10);
        }

        // Test scenario with edge insertion:
        {
            let seg_isec_seq_order_read_range = IntRange::from_pair(10, 20);
            let mut seg = SeqOrderSplitReadSegment {
                seq_order_read_start: 0,
                seq_order_read_end: 20,
                chrom_index: 0,
                pos: 100,
                is_fwd_strand: true,
                cigar: CigarString(vec![Match(5), Ins(10), Match(5)]),
                mapq: 60,
                from_primary_bam_record: false,
            };
            clip_seg_isec_range(&mut seg, &seg_isec_seq_order_read_range);

            assert_eq!(seg.pos, 100);
            assert_eq!(seg.cigar, CigarString(vec![Match(5), SoftClip(15)]));
            assert_eq!(seg.seq_order_read_start, 0);
            assert_eq!(seg.seq_order_read_end, 5);
        }

        // Test scenario with edge insertion & reverse aln:
        {
            let seg_isec_seq_order_read_range = IntRange::from_pair(10, 20);
            let mut seg = SeqOrderSplitReadSegment {
                seq_order_read_start: 0,
                seq_order_read_end: 20,
                chrom_index: 0,
                pos: 100,
                is_fwd_strand: false,
                cigar: CigarString(vec![Match(5), Ins(10), Match(5)]),
                mapq: 60,
                from_primary_bam_record: false,
            };
            clip_seg_isec_range(&mut seg, &seg_isec_seq_order_read_range);

            assert_eq!(seg.pos, 105);
            assert_eq!(seg.cigar, CigarString(vec![SoftClip(15), Match(5)]));
            assert_eq!(seg.seq_order_read_start, 0);
            assert_eq!(seg.seq_order_read_end, 5);
        }
    }
}
