#![allow(unused)]

use std::collections::BTreeSet;
use std::collections::btree_set::Iter;
use std::ops::Bound::{Excluded, Unbounded};

use rust_vc_utils::IntRange;

/// A set of non-overlapping ranges over a single coordinate space
///
#[derive(Default)]
pub struct MergedRangeSet {
    ranges: BTreeSet<IntRange>,

    /// If true, merge all ranges which are adjacent
    ///
    /// Overlapping ranges are always merged
    merge_adjacent: bool,
}

impl MergedRangeSet {
    pub fn new(merge_adjacent: bool) -> Self {
        Self {
            ranges: Default::default(),
            merge_adjacent,
        }
    }

    /// Merge will join ranges which are overlapping (or adjacent if merge_adjacent is set)
    ///
    pub fn insert(&mut self, mut range: IntRange) {
        assert!(range.size() >= 1);

        let query_range = if self.merge_adjacent {
            // To support merging adjacent cases, we expand range by 1
            IntRange {
                start: range.start - 1,
                end: range.end + 1,
            }
        } else {
            range.clone()
        };
        for x in self.find_overlaps(&query_range) {
            if x.start < range.start {
                range.start = x.start;
            }
            if x.end > range.end {
                range.end = x.end;
            }

            self.ranges.remove(&x);
        }

        self.ranges.insert(range);
    }

    /// Remove will divide any existing ranges such that the input range is removed from the range set
    ///
    pub fn remove(&mut self, range: IntRange) {
        assert!(range.size() >= 1);

        for x in self.find_overlaps(&range) {
            if x.start < range.start {
                self.ranges
                    .insert(IntRange::from_pair(x.start, range.start));
            }
            if x.end > range.end {
                self.ranges.insert(IntRange::from_pair(range.end, x.end));
            }

            self.ranges.remove(&x);
        }
    }

    pub fn merge(&mut self, other: &Self) {
        assert_eq!(self.merge_adjacent, other.merge_adjacent);

        for x in other.iter() {
            self.insert(x.clone());
        }
    }

    pub fn iter(&self) -> Iter<'_, IntRange> {
        self.ranges.iter()
    }

    /// Find all members of the range set which overlap the query range
    ///
    /// This will not return adjacent ranges
    ///
    pub fn find_overlaps(&self, range: &IntRange) -> Vec<IntRange> {
        let mut ret = Vec::new();
        let mut iter = self.ranges.range((
            Unbounded,
            Excluded(&IntRange::from_pair(range.end, range.end)),
        ));

        while let Some(x) = iter.next_back() {
            if !range.intersect_range(x) {
                break;
            }
            ret.push(x.clone());
        }
        ret.reverse();
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_insert() {
        let mut mrs = MergedRangeSet::new(false);
        mrs.insert(IntRange::from_pair(10, 20));
        mrs.insert(IntRange::from_pair(20, 30));

        assert_eq!(
            mrs.iter().collect::<Vec<_>>(),
            vec![&IntRange::from_pair(10, 20), &IntRange::from_pair(20, 30)]
        );
    }

    #[test]
    fn test_insert_merge_adjacent() {
        let mut mrs = MergedRangeSet::new(true);
        mrs.insert(IntRange::from_pair(10, 20));
        mrs.insert(IntRange::from_pair(20, 30));

        assert_eq!(
            mrs.iter().collect::<Vec<_>>(),
            vec![&IntRange::from_pair(10, 30)]
        );

        let mut mrs = MergedRangeSet::new(true);
        mrs.insert(IntRange::from_pair(10, 20));
        mrs.insert(IntRange::from_pair(21, 30));

        assert_eq!(
            mrs.iter().collect::<Vec<_>>(),
            vec![&IntRange::from_pair(10, 20), &IntRange::from_pair(21, 30)]
        );

        let mut mrs = MergedRangeSet::new(true);
        mrs.insert(IntRange::from_pair(10, 20));
        mrs.insert(IntRange::from_pair(30, 40));
        mrs.insert(IntRange::from_pair(20, 30));

        assert_eq!(
            mrs.iter().collect::<Vec<_>>(),
            vec![&IntRange::from_pair(10, 40)]
        );
    }

    #[test]
    fn test_remove() {
        let mut mrs = MergedRangeSet::default();
        mrs.remove(IntRange::from_pair(10, 20));

        assert_eq!(
            mrs.iter().collect::<Vec<&IntRange>>(),
            Vec::<&IntRange>::new()
        );

        let mut mrs = MergedRangeSet::default();
        mrs.insert(IntRange::from_pair(10, 40));
        mrs.remove(IntRange::from_pair(20, 30));

        assert_eq!(
            mrs.iter().collect::<Vec<_>>(),
            vec![&IntRange::from_pair(10, 20), &IntRange::from_pair(30, 40)]
        );

        let mut mrs = MergedRangeSet::default();
        mrs.insert(IntRange::from_pair(10, 40));
        mrs.remove(IntRange::from_pair(30, 40));

        assert_eq!(
            mrs.iter().collect::<Vec<_>>(),
            vec![&IntRange::from_pair(10, 30)]
        );
    }

    #[test]
    fn test_find_overlaps() {
        let mut mrs = MergedRangeSet::default();

        let olap = mrs.find_overlaps(&IntRange::from_pair(10, 20));
        assert_eq!(olap, Vec::new());

        mrs.insert(IntRange::from_pair(20, 30));

        let olap = mrs.find_overlaps(&IntRange::from_pair(10, 20));
        assert_eq!(olap, Vec::new());

        let olap = mrs.find_overlaps(&IntRange::from_pair(30, 40));
        assert_eq!(olap, Vec::new());

        let olap = mrs.find_overlaps(&IntRange::from_pair(10, 21));
        assert_eq!(olap, vec![IntRange::from_pair(20, 30)]);

        let olap = mrs.find_overlaps(&IntRange::from_pair(29, 40));
        assert_eq!(olap, vec![IntRange::from_pair(20, 30)]);

        let olap = mrs.find_overlaps(&IntRange::from_pair(10, 40));
        assert_eq!(olap, vec![IntRange::from_pair(20, 30)]);
    }
}
