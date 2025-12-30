use std::collections::BTreeMap;

use rust_vc_utils::cigar::get_cigar_ref_offset;
use rust_vc_utils::depth_tracker::DepthTracker;
use rust_vc_utils::{ChromList, IntRange};

use super::{AllContigMappingInfo, clear_segment_phase_sets};
use crate::merged_range_set::MergedRangeSet;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq, PartialOrd, Ord)]
pub struct ContigSegmentInfo {
    pub contig_index: usize,
    pub segment_index: usize,
}

pub type ChromPhaseSetMap = BTreeMap<IntRange, Vec<ContigSegmentInfo>>;

/// Convert depth tracker into diploid ranges, and also force adjacent diploid ranges into separate groups
///
pub fn get_phased_diploid_ranges(
    depth_tracker: &DepthTracker,
    max_region_size: i64,
) -> MergedRangeSet {
    let mut mrs = MergedRangeSet::default();

    fn insert_range(mrs: &mut MergedRangeSet, start: i64, end: i64) {
        if start < end {
            mrs.insert(IntRange::from_pair(start, end));
        }
    }

    let mut total_depth = 0;
    let mut region_start_pos = None;
    for (pos, delta) in depth_tracker.depth_edges.iter() {
        total_depth += delta;
        let is_depth_condition = total_depth == 2;
        if region_start_pos.is_some() {
            let start_pos = region_start_pos.unwrap();
            let end_pos = std::cmp::min(*pos as i64, max_region_size);
            insert_range(&mut mrs, start_pos, end_pos);
            region_start_pos = None;
        }

        if is_depth_condition {
            region_start_pos = Some(*pos as i64);
        }
    }
    if let Some(region_start_pos) = region_start_pos {
        insert_range(&mut mrs, region_start_pos, max_region_size);
    }

    assert_eq!(total_depth, 0);

    mrs
}

fn merge_range_list(ranges: &[IntRange]) -> IntRange {
    assert!(!ranges.is_empty());

    let mut mr = ranges[0].clone();
    for r in ranges.iter().skip(1) {
        mr.merge(r);
    }
    mr
}

fn process_cluster(
    last_range: &IntRange,
    last_seg_list: &[ContigSegmentInfo],
    range_cluster: &mut Vec<IntRange>,
    remove_ranges: &mut Vec<IntRange>,
    add_entries: &mut BTreeMap<IntRange, Vec<ContigSegmentInfo>>,
) {
    if !range_cluster.is_empty() {
        range_cluster.push(last_range.clone());
        let merged_range = merge_range_list(range_cluster);

        // check that merged_size is some reasonable fraction of the components
        let component_range_size = range_cluster.iter().map(|x| x.size()).sum::<i64>();
        if component_range_size as f64 * 1.5 > merged_range.size() as f64 {
            remove_ranges.append(range_cluster);
            add_entries.insert(merged_range, last_seg_list.to_vec());
        } else {
            range_cluster.clear();
        }
    }
}

fn create_annotated_phase_set_ranges(
    chrom_list: &ChromList,
    all_contig_mapping_info: &AllContigMappingInfo,
) -> Vec<ChromPhaseSetMap> {
    // Step 1: Record contig depth across each chromosome
    //
    let chrom_count = chrom_list.data.len();
    let mut depth_tracker = vec![DepthTracker::default(); chrom_count];
    for contig_mapping_info in all_contig_mapping_info
        .iter()
        .filter(|x| !x.ordered_contig_segment_info.is_empty())
    {
        for seq_order_segment in contig_mapping_info
            .ordered_contig_segment_info
            .iter()
            .map(|x| &x.seq_order_segment)
        {
            let chrom_depth_tracker = &mut depth_tracker[seq_order_segment.chrom_index];
            chrom_depth_tracker
                .process_alignment(seq_order_segment.pos, &seq_order_segment.cigar.0);
        }
    }

    // Step 2: Convert contig depth into phase set ranges for each chromosome:
    //
    let mut phase_set_ranges = Vec::new();
    for (chrom_index, chrom_depth_tracker) in depth_tracker.into_iter().enumerate() {
        let chrom_info = &chrom_list.data[chrom_index];
        phase_set_ranges.push(get_phased_diploid_ranges(
            &chrom_depth_tracker,
            chrom_info.length as i64,
        ));
    }
    let phase_set_ranges = phase_set_ranges;

    // Step 3: Rescan over contigs to create annotated phase set ranges:
    //
    let mut annotated_phase_set_ranges = vec![ChromPhaseSetMap::default(); chrom_count];
    for (contig_index, contig_mapping_info) in all_contig_mapping_info
        .iter()
        .enumerate()
        .filter(|(_, x)| !x.ordered_contig_segment_info.is_empty())
    {
        for (segment_index, contig_mapping_segment_info) in contig_mapping_info
            .ordered_contig_segment_info
            .iter()
            .enumerate()
        {
            let seq_order_segment = &contig_mapping_segment_info.seq_order_segment;
            let chrom_phase_set_ranges = &phase_set_ranges[seq_order_segment.chrom_index];
            let chrom_annotated_phase_set_ranges =
                &mut annotated_phase_set_ranges[seq_order_segment.chrom_index];

            let start = seq_order_segment.pos;
            let end = start + get_cigar_ref_offset(&seq_order_segment.cigar);

            for phase_set_interval in chrom_phase_set_ranges
                .find_overlaps(&IntRange::from_pair(start, end))
                .iter()
            {
                let interval_value = chrom_annotated_phase_set_ranges
                    .entry(phase_set_interval.clone())
                    .or_default();

                interval_value.push(ContigSegmentInfo {
                    contig_index,
                    segment_index,
                });
            }
        }
    }

    annotated_phase_set_ranges
}

/// Refine annotated phase set ranges to remove small interupting haplotypes to otherwise contiguous diploid phase sets
fn refine_annotated_phase_set_ranges(annotated_phase_set_ranges: &mut [ChromPhaseSetMap]) {
    let debug = false;

    for chrom_apsr in annotated_phase_set_ranges.iter_mut() {
        let mut range_cluster = Vec::new();
        let mut remove_ranges = Vec::new();
        let mut add_entries = BTreeMap::new();

        let mut last_diploid_entry: Option<(&IntRange, &Vec<ContigSegmentInfo>)> = None;
        for (range, seg_list) in chrom_apsr.iter() {
            if seg_list.len() != 2 {
                continue;
            }

            if let Some((last_range, last_seg_list)) = last_diploid_entry {
                if debug {
                    eprintln!("seg_list: {seg_list:?}");
                    eprintln!("last_seg_list: {last_seg_list:?}");
                }
                if *seg_list == *last_seg_list {
                    range_cluster.push(last_range.clone());
                } else {
                    process_cluster(
                        last_range,
                        last_seg_list,
                        &mut range_cluster,
                        &mut remove_ranges,
                        &mut add_entries,
                    );
                }
            }
            last_diploid_entry = Some((range, seg_list));
        }

        if let Some((last_range, last_seg_list)) = last_diploid_entry {
            process_cluster(
                last_range,
                last_seg_list,
                &mut range_cluster,
                &mut remove_ranges,
                &mut add_entries,
            );
        }

        for x in remove_ranges {
            if debug {
                eprintln!("removing: {x:?}");
            }
            chrom_apsr.remove(&x);
        }

        for (k, v) in add_entries {
            if debug {
                eprintln!("adding: {k:?} {v:?}");
            }
            chrom_apsr.insert(k, v);
        }
    }
}

/// Update the contig segment phase_sets entry based on the most up to date annotated_phase_sets value
///
fn update_segment_phase_sets(
    annotated_phase_set_ranges: &[ChromPhaseSetMap],
    all_contig_mapping_info: &mut AllContigMappingInfo,
) {
    clear_segment_phase_sets(all_contig_mapping_info);

    for chrom_annotated_phase_set_ranges in annotated_phase_set_ranges.iter() {
        for (phase_set_range, phase_set_contig_segments) in chrom_annotated_phase_set_ranges.iter()
        {
            for contig_segment in phase_set_contig_segments {
                all_contig_mapping_info[contig_segment.contig_index].ordered_contig_segment_info
                    [contig_segment.segment_index]
                    .phase_sets
                    .add_region(phase_set_range.start, phase_set_range.end);
            }
        }
    }
}

/// Add initial phase set tags to contigs, based on regions with an expected contig depth
///
/// We separately return a structure containing phase set range information, which facilitates optional downstream
/// refinement of the phase set regions based on read evidence. This refinment step is intended for cases where
/// partially-phased contigs have been given as input.
///
pub fn add_contig_segment_phase_set_labels(
    chrom_list: &ChromList,
    all_contig_mapping_info: &mut AllContigMappingInfo,
) -> Vec<ChromPhaseSetMap> {
    let mut annotated_phase_set_ranges =
        create_annotated_phase_set_ranges(chrom_list, all_contig_mapping_info);

    refine_annotated_phase_set_ranges(&mut annotated_phase_set_ranges);

    update_segment_phase_sets(&annotated_phase_set_ranges, all_contig_mapping_info);

    annotated_phase_set_ranges
}
