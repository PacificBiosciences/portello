use log::info;
use rust_htslib::bam::record::Cigar;
use rust_vc_utils::{
    SeqOrderSplitReadSegment,
    cigar::{get_cigar_ref_offset, strip_leading_clip, strip_trailing_clip},
    get_read_segment_to_ref_pos_tree_map,
};

use super::{AllContigMappingInfo, ContigMappingSegmentInfo, print_split_read_summary};

/// Get the reference distance gap between two adjacent split read segments
///
/// Assumes segments are on the same chromosome
///
fn get_seg_ref_gap(seg1: &SeqOrderSplitReadSegment, seg2: &SeqOrderSplitReadSegment) -> i64 {
    if seg1.is_fwd_strand {
        let seg1_read_end = seg1.pos + get_cigar_ref_offset(&seg1.cigar);
        seg2.pos - seg1_read_end
    } else {
        let seg2_read_end = seg2.pos + get_cigar_ref_offset(&seg2.cigar);
        seg1.pos - seg2_read_end
    }
}

/// Determine if these adjacent split read segments meet the co-linear join criteria:
///
fn are_segments_joinable(seg1: &SeqOrderSplitReadSegment, seg2: &SeqOrderSplitReadSegment) -> bool {
    // co-linearity and joinability criteria:
    //
    // 1. adjacent split reads
    // 2. Same chrom
    // 3. Same direction
    // 4. In forward direction: End of split read 1 preceeds start of split read 2, with a gap of less than max_segment_ref_gap bases
    //    In reverse direction: Start of split read 2 preceeds end of split read 1, with a gap of less than max_segment_ref_gap bases
    // 5. Same mapq (at least for first pass at this)
    //
    let max_segment_ref_gap = 1000;

    if seg1.chrom_index != seg2.chrom_index || seg1.is_fwd_strand != seg2.is_fwd_strand {
        false
    } else {
        let seg_ref_gap = get_seg_ref_gap(seg1, seg2);
        if seg_ref_gap < 0 || seg_ref_gap > max_segment_ref_gap {
            false
        } else {
            seg1.mapq == seg2.mapq
        }
    }
}

/// Join seg2 into seg1
///
/// seg1 must precede seg2 in sequencer read order, and must represent adjacent split read segments
///
/// seg1 and seg2 are assumed to not have any repeated match bases
///
fn join_segments(
    seg_info1: &mut ContigMappingSegmentInfo,
    mut seg_info2: ContigMappingSegmentInfo,
    debug: bool,
) {
    let contig_to_ref_map = {
        // first join the split read segment itself:
        let seg1 = &mut seg_info1.seq_order_segment;
        let seg2 = &mut seg_info2.seq_order_segment;

        let join_del_size = get_seg_ref_gap(seg1, seg2);
        assert!(join_del_size >= 0);
        let join_ins_size = {
            assert!(seg2.seq_order_read_start >= seg1.seq_order_read_end);
            seg2.seq_order_read_start - seg1.seq_order_read_end
        };

        if debug {
            eprintln!("join_segments ins/del size: {join_ins_size}/{join_del_size}");
        }

        /// Append Z-drop region and segment b cigar onto cigar a
        fn join_cigars(
            a: &mut Vec<Cigar>,
            b: &mut Vec<Cigar>,
            join_ins_size: usize,
            join_del_size: i64,
        ) {
            strip_trailing_clip(a);
            if join_ins_size > 0 {
                a.push(Cigar::Ins(join_ins_size as u32));
            }
            if join_del_size > 0 {
                a.push(Cigar::Del(join_del_size as u32));
            }
            strip_leading_clip(b);
            a.append(b);
        }

        if seg1.is_fwd_strand {
            join_cigars(
                &mut seg1.cigar.0,
                &mut seg2.cigar.0,
                join_ins_size,
                join_del_size,
            );
        } else {
            join_cigars(
                &mut seg2.cigar.0,
                &mut seg1.cigar.0,
                join_ins_size,
                join_del_size,
            );

            std::mem::swap(&mut seg1.cigar, &mut seg2.cigar);
            seg1.pos = seg2.pos;
        }

        seg1.seq_order_read_end = seg2.seq_order_read_end;

        get_read_segment_to_ref_pos_tree_map(seg1.pos, &seg1.cigar, false)
    };

    // Finally, update the contig ref map based on the updated cigar:
    seg_info1.contig_to_ref_map = contig_to_ref_map;
}

pub fn join_colinear_contig_segments(result: &mut AllContigMappingInfo) {
    let debug = false;

    info!("Joining colinear split alignment segments in each assembly contig");

    let mut segments_joined = 0;

    for contig_mapping_info in result
        .iter_mut()
        .filter(|x| !x.ordered_contig_segment_info.is_empty())
    {
        if debug {
            print_split_read_summary(contig_mapping_info);
        }

        let old_ordered_contig_segment_info =
            std::mem::take(&mut contig_mapping_info.ordered_contig_segment_info);
        let new_ordered_contig_segment_info = &mut contig_mapping_info.ordered_contig_segment_info;
        for segment in old_ordered_contig_segment_info.into_iter() {
            if new_ordered_contig_segment_info.is_empty() {
                new_ordered_contig_segment_info.push(segment);
                continue;
            }

            let last_segment = new_ordered_contig_segment_info.last_mut().unwrap();

            assert!(
                segment.seq_order_segment.seq_order_read_start
                    >= last_segment.seq_order_segment.seq_order_read_end,
                "Incomplete repeat trimming on qname: {} Segment1: {} Segment2: {}",
                contig_mapping_info.qname,
                last_segment.seq_order_segment.short_display(),
                segment.seq_order_segment.short_display()
            );

            if are_segments_joinable(&last_segment.seq_order_segment, &segment.seq_order_segment) {
                if debug {
                    eprintln!("Joinable segments found.");
                    eprintln!(
                        "Segment1: {}",
                        last_segment.seq_order_segment.short_display()
                    );
                    eprintln!("Segment2: {}", segment.seq_order_segment.short_display());
                }

                join_segments(last_segment, segment, debug);

                segments_joined += 1;

                if debug {
                    eprintln!(
                        "Segment1-after: {}",
                        last_segment.seq_order_segment.short_display()
                    );
                }
            } else {
                new_ordered_contig_segment_info.push(segment);
            }
        }
    }

    info!("Joined {segments_joined} colinear segments");
}
