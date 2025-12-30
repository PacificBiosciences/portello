use std::collections::BTreeMap;

use rust_htslib::bam::record::Cigar;

use crate::{
    IntRange, RegionMap,
    cigar::{get_alignment_ref_segments, get_cigar_ref_offset},
};

/// Track depth information from read alignments on one chromsome
///
/// A depth tracker is a data structure that enumerates read details to effeciently store an exact (base-resolution)
/// depth trace of the sample on one chromosome. It is mergable with other depth trackers (under the assumption that the
/// two merged trackers have been provided with independent read evidence).
///
/// The depth tracker can be projected into lower-resolution structures for different applications, such as depth bins
/// or regions defined on a depth predicate function.
///
#[derive(Clone, Default)]
pub struct DepthTracker {
    /// The smallest alignment CIGAR deletion size that will be accounted for while building the depth track
    /// (split reads are always accounted for).
    pub min_del_size: Option<u32>,

    /// Primary depth storage structure, key is chrom position, value is depth delta.
    ///
    pub depth_edges: BTreeMap<u64, i32>,
}

impl DepthTracker {
    pub fn new(min_del_size: Option<u32>) -> Self {
        Self {
            min_del_size,
            depth_edges: Default::default(),
        }
    }

    pub fn merge(&mut self, other: &Self) {
        assert_eq!(self.min_del_size, other.min_del_size);
        for (pos, delta) in other.depth_edges.iter() {
            *self.depth_edges.entry(*pos).or_insert(0) += delta;
        }
    }

    pub fn process_alignment(&mut self, ref_pos: i64, cigar: &[Cigar]) {
        let ref_segs = if let Some(min_del_size) = self.min_del_size {
            get_alignment_ref_segments(ref_pos, cigar, min_del_size)
        } else {
            let end = ref_pos + get_cigar_ref_offset(cigar);
            vec![IntRange::from_pair(ref_pos, end)]
        };

        for ref_seg in ref_segs {
            assert!(ref_seg.start <= ref_seg.end);

            let start = std::cmp::max(ref_seg.start, 0) as u64;
            let end = std::cmp::max(ref_seg.end, 0) as u64;

            *self.depth_edges.entry(start).or_default() += 1;
            *self.depth_edges.entry(end).or_default() -= 1;
        }
    }
}

/// Convert depth tracker into regions based on predicate function of depth
///
pub fn get_regions_from_depth_predicate<F: Fn(i32) -> bool>(
    depth_tracker: &DepthTracker,
    max_region_size: i64,
    depth_predicate: F,
) -> RegionMap {
    let mut region_map = RegionMap::default();

    let mut total_depth = 0;
    let mut region_start_pos = None;
    for (pos, delta) in depth_tracker.depth_edges.iter() {
        total_depth += delta;
        let is_depth_condition = depth_predicate(total_depth);
        if region_start_pos.is_some() {
            if !is_depth_condition {
                let start_pos = region_start_pos.unwrap();
                let end_pos = std::cmp::min(*pos as i64, max_region_size);
                if start_pos < end_pos {
                    region_map.add_region(start_pos, end_pos);
                }
                region_start_pos = None;
            }
        } else if is_depth_condition {
            region_start_pos = Some(*pos as i64);
        }
    }
    if let Some(region_start_pos) = region_start_pos
        && region_start_pos < max_region_size
    {
        region_map.add_region(region_start_pos, max_region_size);
    }

    assert_eq!(total_depth, 0);

    region_map
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_regions_from_depth_predicate() {
        let mut dt = DepthTracker::default();

        dt.process_alignment(10, &[Cigar::Match(100)]);
        dt.process_alignment(20, &[Cigar::Match(100)]);
        dt.process_alignment(30, &[Cigar::Match(100)]);

        let region_map = get_regions_from_depth_predicate(&dt, 200, |x| x >= 3);
        let val = region_map.find_overlaps(0, 200).next();
        assert!(val.is_some());
        let interval = val.unwrap().interval();
        assert_eq!(interval.start, 30);
        assert_eq!(interval.end, 110);
    }
}
