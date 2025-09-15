use rust_vc_utils::{GenomeSegment, IntRange};

use super::AllContigMappingInfo;

/// Remove all split reads mapped to non-targeted regions
///
pub fn filter_non_targeted_segments(
    target_region: Option<&GenomeSegment>,
    result: &mut AllContigMappingInfo,
) {
    let target_region = match target_region {
        Some(x) => x,
        None => return,
    };

    for contig_mapping_info in result.iter_mut() {
        let mut old_ordered_contig_segment_info = Vec::new();
        std::mem::swap(
            &mut contig_mapping_info.ordered_contig_segment_info,
            &mut old_ordered_contig_segment_info,
        );

        for x in old_ordered_contig_segment_info.into_iter() {
            // Test if split read starts in the target region, and only keep those split reads
            //
            // Note that there's an upstream limitation to the targeted run where the supplementary reads for each
            // contigs will only be found if they start in the targeted region (rather than the more general criteria of
            // intersecting the region), so we need to replicate that limitation here.
            //
            let split_read_start_region = GenomeSegment {
                chrom_index: x.seq_order_segment.chrom_index,
                range: IntRange::from_int(x.seq_order_segment.pos),
            };
            if target_region.intersect(&split_read_start_region) {
                contig_mapping_info.ordered_contig_segment_info.push(x);
            }
        }
    }
}
