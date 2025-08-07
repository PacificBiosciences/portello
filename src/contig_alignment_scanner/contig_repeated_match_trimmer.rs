use std::cmp::Ordering;

use rust_htslib::bam::record::CigarString;
use rust_vc_utils::IntRange;
use rust_vc_utils::SeqOrderSplitReadSegment;
use rust_vc_utils::cigar::{
    clip_alignment_read_edges, get_cigar_read_offset, get_gap_compressed_identity_no_align_match,
};

use super::{AllContigMappingInfo, ContigMappingInfo};

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
fn clip_seg_isec_range(
    seg: &mut SeqOrderSplitReadSegment,
    seg_isec_seq_order_read_range: &IntRange,
) {
    let ignore_hard_clip = false;

    // Determine if we're clipping a prefix (suffix) in both sequence order and alignment order:
    let is_clipping_seq_order_prefix =
        { seg_isec_seq_order_read_range.start == seg.seq_order_read_start as i64 };
    let is_clipping_prefix = is_clipping_seq_order_prefix ^ (!seg.is_fwd_strand);

    // Update seq_order_read_{start/end}
    if is_clipping_seq_order_prefix {
        seg.seq_order_read_start = seg_isec_seq_order_read_range.end as usize;
    } else {
        seg.seq_order_read_end = seg_isec_seq_order_read_range.start as usize;
    }

    // Convert read intersect range from sequencing order to alignment order
    let read_len = get_cigar_read_offset(&seg.cigar, ignore_hard_clip);
    let seg_isec_read_range = if seg.is_fwd_strand {
        seg_isec_seq_order_read_range.clone()
    } else {
        seg_isec_seq_order_read_range.get_reverse_range(read_len as i64)
    };

    // Get seg alignment trimmed to just the intersection range
    let (min_left_clip, min_right_clip) = if is_clipping_prefix {
        (seg_isec_read_range.end as usize, 0)
    } else {
        (0, read_len - seg_isec_read_range.start as usize)
    };
    let (shifted_cigar, ref_pos_shift) =
        clip_alignment_read_edges(&seg.cigar, min_left_clip, min_right_clip);
    seg.cigar = CigarString(shifted_cigar);
    seg.pos += ref_pos_shift;
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
    for contig_mapping_info in result.iter_mut() {
        let debug = false;
        if debug {
            eprintln!("contig split read summary:");
            for (i, x) in contig_mapping_info
                .ordered_contig_segment_info
                .iter()
                .enumerate()
            {
                eprintln!(
                    "split {i}: start: {} end: {}",
                    x.seq_order_segment.seq_order_read_start,
                    x.seq_order_segment.seq_order_read_end
                );
            }
        }
        let split_read_count = contig_mapping_info.ordered_contig_segment_info.len();
        for split_index1 in 0..split_read_count {
            for split_index2 in (split_index1 + 1)..split_read_count {
                let (seg_isec_seq_order_read_range, clip_split_index) =
                    match get_seg_clip_info(contig_mapping_info, split_index1, split_index2, debug)
                    {
                        Some(x) => x,
                        None => break,
                    };

                clip_seg_isec_range(
                    &mut contig_mapping_info.ordered_contig_segment_info[clip_split_index]
                        .seq_order_segment,
                    &seg_isec_seq_order_read_range,
                );
            }
        }
    }
}
