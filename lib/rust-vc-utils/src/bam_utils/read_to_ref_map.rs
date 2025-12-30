use std::collections::BTreeMap;

use rust_htslib::bam;
use rust_htslib::bam::record::Cigar;

use super::cigar::update_ref_and_read_pos;

/// Create an array of read length, mapping from the read to the reference position for the given read alignment
///
/// Allowing the ref_pos and cigar alignment input enables the processing of split alignment segments with this
/// function.
///
/// Read positions which are not mapped to a reference position are None (such as inside of an insertion)
///
/// Ref positions are 0-indexed
///
pub fn get_read_segment_to_ref_pos_map(
    seq_len: usize,
    mut ref_pos: i64,
    cigar: &[Cigar],
    ignore_hard_clip: bool,
) -> Vec<Option<i64>> {
    use rust_htslib::bam::record::Cigar::*;

    let mut read_to_ref = vec![None; seq_len];
    let mut read_pos = 0usize;

    for c in cigar.iter() {
        match c {
            Diff(len) | Equal(len) | Match(len) => {
                let len = *len as usize;
                for i in 0..len {
                    read_to_ref[read_pos + i] = Some(ref_pos + i as i64);
                }
            }
            _ => {}
        }
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, ignore_hard_clip);
    }
    read_to_ref
}

/// Create an array of read length, mapping from the read to the reference position for the read's primary alignment
/// segment
///
/// Read positions which are not mapped to a reference position are None (such as inside of an insertion)
///
/// Ref positions are 0-indexed
///
pub fn get_read_to_ref_pos_map(record: &bam::Record, ignore_hard_clip: bool) -> Vec<Option<i64>> {
    get_read_segment_to_ref_pos_map(
        record.seq_len(),
        record.pos(),
        &record.cigar(),
        ignore_hard_clip,
    )
}

/// Compressed representation of read to ref position mapping
///
#[derive(Clone, Default)]
pub struct ReadToRefTreeMap {
    // Key is the read_pos starting a block
    // Val is None for read not mapping to ref, or the ref_pos at the start of a match block
    //
    // This works because:
    // - a deletion is represented by adjacent keys where read_pos diff < ref_pos diff
    // - an insertion is represented by key with value of None, length is found by the read_pos diff to the adjacent key.
    map: BTreeMap<usize, Option<i64>>,
}

impl ReadToRefTreeMap {
    /// Given read_pos, return the corresponding reference position
    ///
    pub fn get_ref_pos(&self, read_pos: usize) -> Option<i64> {
        match self.map.range(..=read_pos).next_back() {
            Some((k, v)) => v.map(|ref_pos| ref_pos + (read_pos - k) as i64),
            None => None,
        }
    }

    pub fn get_ref_range(
        &self,
        read_start_pos: usize,
        read_end_pos: usize,
    ) -> std::collections::btree_map::Range<'_, usize, Option<i64>> {
        let read_start_block_pos = match self.map.range(..=read_start_pos).next_back() {
            Some((k, _)) => *k,
            None => read_start_pos,
        };

        self.map.range(read_start_block_pos..read_end_pos)
    }

    pub fn get_map(&self) -> &BTreeMap<usize, Option<i64>> {
        &self.map
    }

    /// Extended debug output
    pub fn dump(&self) {
        eprintln!("ReadToRefTreeMap");
        for (k, v) in self.map.iter() {
            eprintln!("k: {k} v: {v:?}");
        }
    }
}

/// Create a data structure enabling rapid lookup from the read to the reference position for the given read alignment
///
/// Allowing the ref_pos and cigar alignment input enables the processing of split alignment segments with this
/// function.
///
/// Read positions which are not mapped to a reference position are None (such as inside of an insertion)
///
/// Ref positions are 0-indexed
///
/// This function is designed to parallel the above `get_read_segment_to_ref_pos_map` function, but changes the return
/// value to the more efficient tree map
///
pub fn get_read_segment_to_ref_pos_tree_map(
    mut ref_pos: i64,
    cigar: &[Cigar],
    ignore_hard_clip: bool,
) -> ReadToRefTreeMap {
    use rust_htslib::bam::record::Cigar::*;

    let mut read_to_ref = ReadToRefTreeMap::default();
    let mut read_pos = 0usize;

    let mut update_map = |ref_pos: i64, read_pos: usize, match_len: &mut usize| {
        if *match_len > 0 {
            read_to_ref
                .map
                .insert(read_pos - *match_len, Some(ref_pos - *match_len as i64));
            read_to_ref.map.insert(read_pos, None);
            *match_len = 0;
        }
    };

    let mut match_len = 0;

    for c in cigar.iter() {
        match c {
            Diff(len) | Equal(len) | Match(len) => {
                match_len += *len as usize;
            }
            _ => {
                update_map(ref_pos, read_pos, &mut match_len);
            }
        }
        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, ignore_hard_clip);
    }
    update_map(ref_pos, read_pos, &mut match_len);

    read_to_ref
}

/// Compressed representation of ref to read position mapping
///
#[derive(Clone, Default)]
pub struct RefToReadTreeMap {
    // Key is the ref_pos starting a block
    // Val is None for ref not mapping to read (a deletion), or the read_pos at the start of a match block
    //
    map: BTreeMap<i64, Option<usize>>,
}

impl RefToReadTreeMap {
    /// Create tree map by inverting the read-to-ref tree map
    pub fn from_read_to_ref_tree_map(inmap: &ReadToRefTreeMap) -> Self {
        let mut map = BTreeMap::new();

        let mut is_last_entry_none = false;
        let mut last_block_start: Option<(usize, i64)> = None;
        for (read_pos, ref_pos) in inmap.map.iter() {
            if let Some(ref_pos) = ref_pos {
                // Determine if there's a read-to-ref deletion encoded between the last two blocks:
                if !is_last_entry_none && let Some((last_read_pos, last_ref_pos)) = last_block_start
                {
                    let read_block_size = *read_pos as i64 - last_read_pos as i64;
                    let next_ref_pos = last_ref_pos + read_block_size;
                    if next_ref_pos < *ref_pos {
                        map.insert(next_ref_pos, None);
                    }
                }
                map.insert(*ref_pos, Some(*read_pos));
                is_last_entry_none = false;
                last_block_start = Some((*read_pos, *ref_pos));
            } else {
                //
                if let Some((last_read_pos, last_ref_pos)) = last_block_start {
                    let read_block_size = *read_pos as i64 - last_read_pos as i64;
                    let next_ref_pos = last_ref_pos + read_block_size;
                    map.insert(next_ref_pos, None);
                    is_last_entry_none = true;
                }
            }
        }

        Self { map }
    }

    /// Given ref_pos, return the corresponding read position
    ///
    pub fn get_read_pos(&self, ref_pos: i64) -> Option<usize> {
        match self.map.range(..=ref_pos).next_back() {
            Some((k, v)) => v.map(|read_pos| read_pos + (ref_pos - k) as usize),
            None => None,
        }
    }

    pub fn get_read_range(
        &self,
        ref_start_pos: i64,
        ref_end_pos: i64,
    ) -> std::collections::btree_map::Range<'_, i64, Option<usize>> {
        let ref_start_block_pos = match self.map.range(..=ref_start_pos).next_back() {
            Some((k, _)) => *k,
            None => ref_start_pos,
        };

        self.map.range(ref_start_block_pos..ref_end_pos)
    }

    pub fn get_map(&self) -> &BTreeMap<i64, Option<usize>> {
        &self.map
    }

    /// Extended debug output
    pub fn dump(&self) {
        eprintln!("RefToReadTreeMap");
        for (k, v) in self.map.iter() {
            eprintln!("k: {k} v: {v:?}");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{Header, HeaderView, header};

    fn get_test_header() -> HeaderView {
        let mut _header = Header::new();
        _header.push_record(
            header::HeaderRecord::new(b"SQ")
                .push_tag(b"SN", "chr1")
                .push_tag(b"LN", 10000000),
        );
        HeaderView::from_header(&_header)
    }

    #[test]
    fn test_get_read_to_ref_pos_map() {
        let header = get_test_header();
        let sam_line = b"qname\t0\tchr1\t10\t60\t2H2M1I1M\t*\t0\t0\tACTG\tDDDD";
        let record = bam::Record::from_sam(&header, sam_line).unwrap();

        let read_to_ref_map = get_read_to_ref_pos_map(&record, true);
        assert_eq!(read_to_ref_map, vec![Some(9), Some(10), None, Some(11)]);
    }

    #[test]
    fn test_get_read_segment_to_ref_pos_tree_map() {
        let header = get_test_header();
        let sam_line = b"qname\t0\tchr1\t10\t60\t2H2M1I1M\t*\t0\t0\tACTG\tDDDD";
        let record = bam::Record::from_sam(&header, sam_line).unwrap();

        let read_to_ref_tmap =
            get_read_segment_to_ref_pos_tree_map(record.pos(), &record.cigar(), true);

        // Check that the tree map produces identical output to the uncompressed array output produced by
        // get_read_to_ref_pos_map:
        let rval = (0..4)
            .map(|x| read_to_ref_tmap.get_ref_pos(x))
            .collect::<Vec<_>>();
        assert_eq!(rval, vec![Some(9), Some(10), None, Some(11)]);

        // Check that the behavior off the end of the read makes sense:
        assert_eq!(read_to_ref_tmap.get_ref_pos(4), None);

        let rrange = read_to_ref_tmap.get_ref_range(0, 2).collect::<Vec<_>>();
        assert_eq!(rrange, vec![(&0, &Some(9))]);
    }

    #[test]
    fn test_ref_to_read_tree_map() {
        // First create test read_to_ref tree map, as before:
        let header = get_test_header();
        let sam_line = b"qname\t0\tchr1\t10\t60\t2H2M1I1M1D1M\t*\t0\t0\tACTGG\tDDDDD";
        let record = bam::Record::from_sam(&header, sam_line).unwrap();

        let read_to_ref_tmap =
            get_read_segment_to_ref_pos_tree_map(record.pos(), &record.cigar(), true);

        // Now convert this into a ref_to_read tree map:
        let ref_to_read_tmap = RefToReadTreeMap::from_read_to_ref_tree_map(&read_to_ref_tmap);

        assert_eq!(ref_to_read_tmap.get_read_pos(8), None);
        assert_eq!(ref_to_read_tmap.get_read_pos(9), Some(0));
        assert_eq!(ref_to_read_tmap.get_read_pos(11), Some(3));
        assert_eq!(ref_to_read_tmap.get_read_pos(12), None);
        assert_eq!(ref_to_read_tmap.get_read_pos(13), Some(4));
        assert_eq!(ref_to_read_tmap.get_read_pos(14), None);

        let rrange = ref_to_read_tmap.get_read_range(9, 12).collect::<Vec<_>>();
        assert_eq!(rrange, vec![(&9, &Some(0)), (&11, &Some(3))]);
    }

    #[test]
    fn test_ref_to_read_tree_map2() {
        // Snippet built from a failing real example:
        let mut map = BTreeMap::new();
        map.insert(80203, Some(61446));
        map.insert(81399, Some(62652));
        map.insert(81962, None);
        map.insert(82561, Some(63840));
        map.insert(84956, None);
        map.insert(84970, Some(66235));
        map.insert(91500, None);

        let read_to_ref_tmap = ReadToRefTreeMap { map };

        let ref_to_read_tmap = RefToReadTreeMap::from_read_to_ref_tree_map(&read_to_ref_tmap);

        assert_eq!(ref_to_read_tmap.get_read_pos(62653), Some(81400));
        assert_eq!(ref_to_read_tmap.get_read_pos(63553), None);
    }
}
