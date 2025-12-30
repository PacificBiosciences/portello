use std::collections::BTreeSet;

use log::info;
use rust_vc_utils::cigar::get_cigar_ref_offset;
use rust_vc_utils::{IntRange, RegionMap};

use super::AllContigMappingInfo;

/// Track the span of a contig split-alignment segment over one chromosome
#[derive(Clone)]
struct ContigSegmentSpan {
    contig_index: usize,
    segment_index: usize,
    range: IntRange,
}

/// Add haplotype labels to each contig segment
///
/// These are used to set the HP tags of reads mapping to each contig segment
///
pub fn add_contig_segment_haplotype_labels(chrom_count: usize, result: &mut AllContigMappingInfo) {
    let debug = false;

    info!("Adding contig haplotype labels");

    // Step 1 - make a first pass through the contig segments to build a list of all contig segment reference intervals per chromosome:
    //
    let mut seg_spans = vec![Vec::new(); chrom_count];

    for (contig_index, contig_mapping_info) in result
        .iter()
        .enumerate()
        .filter(|(_, x)| !x.ordered_contig_segment_info.is_empty())
    {
        for (segment_index, seq_order_segment) in contig_mapping_info
            .ordered_contig_segment_info
            .iter()
            .map(|x| &x.seq_order_segment)
            .enumerate()
        {
            let start = seq_order_segment.pos;
            let end = start + get_cigar_ref_offset(&seq_order_segment.cigar);

            seg_spans[seq_order_segment.chrom_index].push(ContigSegmentSpan {
                contig_index,
                segment_index,
                range: IntRange { start, end },
            });
        }
    }

    // Step 2 - Update hap ids for each chromosome:
    //
    // 1. On each chromosome sort by segment size from largest to smallest
    // 2. Label HP in size order based on how many other HP ranges the segment overlaps
    //
    for mut chrom_seg_spans in seg_spans.into_iter() {
        if debug {
            eprintln!(
                "Starting new chrom_index. Span count: {}",
                chrom_seg_spans.len()
            );
        }
        chrom_seg_spans.sort_by_key(|x| std::cmp::Reverse(x.range.size()));

        let mut hap_id_map = RegionMap::default();

        for seg_span in chrom_seg_spans {
            let start = seg_span.range.start;
            let end = seg_span.range.end;

            // Allow for moderate overlaps between contig ends, if these exist we can still have the same HP tag
            //
            // This step isn't as important for scaffolded contigs, but has a large impact on the quality of phasing
            // annotation for unitigs, which tend to have more overlap.
            //
            let (adjusted_start, adjusted_end) = {
                let max_allowed_overlap_size = 24000;
                let size = end - start;
                let allowed_overlap_size = std::cmp::min(max_allowed_overlap_size, size / 3);
                (start + allowed_overlap_size, end - allowed_overlap_size)
            };

            let overlapping_hap_ids = hap_id_map
                .find_overlaps(adjusted_start, adjusted_end)
                .map(|x| *x.data())
                .collect::<BTreeSet<_>>();

            if debug {
                eprintln!(
                    "span: {:?} overlapping_hap_ids: {:?}",
                    &seg_span.range, &overlapping_hap_ids
                );
            }

            let candidate_hap_id = overlapping_hap_ids
                .iter()
                .enumerate()
                .find(|&(i, &x)| i != x as usize)
                .map(|(i, _)| i)
                .unwrap_or(overlapping_hap_ids.len());

            if debug {
                eprintln!("cand_hap_id: {candidate_hap_id}");
            }

            if candidate_hap_id > u8::MAX as usize {
                continue;
            }

            // Good candidate hap_id found, update label:
            hap_id_map.add_region_value(start, end, candidate_hap_id as u8);
            result[seg_span.contig_index].ordered_contig_segment_info[seg_span.segment_index]
                .hap_id = Some(candidate_hap_id as u8);
        }
    }
}
