use rust_htslib::bam::record::{Cigar, CigarString};
use rust_vc_utils::ReadToRefTreeMap;
use rust_vc_utils::cigar::{
    clean_up_cigar_edge_indels, compress_cigar, is_alignment_match, update_ref_pos,
};

/// Given an ref1-ref2 adjacent block pair, determine how this needs to be translated into new ref2_cigar blocks
///
/// This is a private helper method for liftover_read_alignment
///
/// # Arguments
/// * `this_ref1_ref2_mapping_block` - The 'rightward' of the ref1-ref2 adjacent mapping block pair
///
/// * `last_ref1_ref2_mapping_block` - The 'leftward' of the ref1-ref2 adjacent mapping block pair
///
/// * `ref1_cigar_segment_end_pos` - The end position (in ref1 coordinates) of the ref1 cigar segment that we're
///   currently iterating ref1->ref2 blocks over
///
/// * `ref1_cigar_segment` - The ref1 cigar segment that we're currently iterating ref1->ref2 blocks over
///
/// * `block_ref1_pos` - Intermediate state variable that needs to be tracked over the ref1 cigar walk in the parent
///   liftover_read_alignment method. This marks the ref1 start position of the read-ref1 alignment within each
///   ref1->ref2 alignment block which we iterate through for each ref1 cigar segment.
///
/// * `ref2_start_pos` - Where the translated read-ref2 alignment starts in ref2 coordinates, this is will be the
///   persistent output fo the parent liftover_read_alignment method
///
/// * `ref2_end_pos` - Intermediate state variable that needs to be tracked over the ref1 cigar walk in the parent
///   liftover_read_alignment method
///
/// * `ref2_cigar` - The translated read-ref2 alignment, this is will be the persistent output fo the parent
///   liftover_read_alignment method
///
#[allow(clippy::too_many_arguments)]
fn update_ref2_cigar_segment(
    this_ref1_ref2_mapping_block: Option<(usize, Option<i64>)>,
    last_ref1_ref2_mapping_block: Option<(usize, Option<i64>)>,
    ref1_cigar_segment_end_pos: i64,
    ref1_cigar_segment: &Cigar,
    block_ref1_pos: &mut i64,
    ref2_start_pos: &mut Option<i64>,
    ref2_end_pos: &mut Option<i64>,
    ref2_cigar: &mut Vec<Cigar>,
) {
    let debug = false;
    if debug {
        eprintln!("starting update_ref2_cigar_segment");
        eprintln!(
            "this_ref1_ref2_mapping_block: {this_ref1_ref2_mapping_block:?} last_ref1_ref2_mapping_block: {last_ref1_ref2_mapping_block:?}"
        );
        eprintln!(
            "ref1_cigar_segment_end_pos: {ref1_cigar_segment_end_pos} ref1_cigar_segment: {ref1_cigar_segment}"
        );
        eprintln!(
            "START: block_ref1_pos: {block_ref1_pos} ref2_start_pos: {ref2_start_pos:?} ref2_end_pos: {ref2_end_pos:?} ref2_cigar: {}",
            CigarString(ref2_cigar.clone())
        );
    }

    // The end of the last block we're dealing with in ref1 coordinates is either the start of the next ref1-ref2 alignment block, or
    // the end of the current read-ref1 cigar segment, whichever comes first:
    let ref1_remapped_segment_end_pos =
        if let Some((this_block_ref1_start_pos, _)) = this_ref1_ref2_mapping_block {
            std::cmp::min(this_block_ref1_start_pos as i64, ref1_cigar_segment_end_pos)
        } else {
            ref1_cigar_segment_end_pos
        };

    if ref1_remapped_segment_end_pos > *block_ref1_pos {
        let remapped_segment_len = (ref1_remapped_segment_end_pos - *block_ref1_pos) as u32;
        let is_match_segment = is_alignment_match(ref1_cigar_segment);

        if debug {
            eprintln!(
                "ref1_remapped_segment_end_pos: {ref1_remapped_segment_end_pos} remapped_segment_len: {remapped_segment_len}"
            );
        }

        if let Some((last_block_ref1_start_pos, last_block_ref2_start_pos)) =
            last_ref1_ref2_mapping_block
        {
            match last_block_ref2_start_pos {
                Some(last_block_ref2_start_pos) => {
                    if is_match_segment && ref2_start_pos.is_none() {
                        let x = last_block_ref2_start_pos
                            + (*block_ref1_pos - last_block_ref1_start_pos as i64);
                        *ref2_start_pos = Some(x);
                    }

                    // Check for evidence of deletion since last match segment
                    if let Some(ref2_end_pos) = *ref2_end_pos {
                        let deletion_len = last_block_ref2_start_pos - ref2_end_pos;
                        if deletion_len > 0 && ref2_start_pos.is_some() {
                            ref2_cigar.push(Cigar::Del(deletion_len as u32));
                        }
                    }

                    let last_block_expected_ref2_len =
                        ref1_remapped_segment_end_pos - last_block_ref1_start_pos as i64;
                    *ref2_end_pos = Some(last_block_ref2_start_pos + last_block_expected_ref2_len);

                    if is_match_segment || ref2_start_pos.is_some() {
                        let new_segment = match ref1_cigar_segment {
                            Cigar::Del(_) => Cigar::Del(remapped_segment_len),
                            Cigar::RefSkip(_) => Cigar::RefSkip(remapped_segment_len),
                            _ => Cigar::Match(remapped_segment_len),
                        };
                        ref2_cigar.push(new_segment);
                    }
                }
                None => {
                    if is_match_segment {
                        ref2_cigar.push(Cigar::Ins(remapped_segment_len));
                    }
                }
            }
        } else {
            // If last_ref1_ref2_mapping_block is None, it can mean we're looking at a read-ref1 cigar segment prior to any
            // mapping to R2. In this case we convert read-consuming segments to clip:
            if is_match_segment {
                ref2_cigar.push(Cigar::SoftClip(remapped_segment_len));
            }
        }
        *block_ref1_pos = ref1_remapped_segment_end_pos;
    }

    if debug {
        eprintln!(
            "END: block_ref1_pos: {block_ref1_pos} ref2_start_pos: {ref2_start_pos:?} ref2_end_pos: {ref2_end_pos:?} ref2_cigar: {}",
            CigarString(ref2_cigar.clone())
        );
    }
}

/// Given an alignment from read to ref1, lift the alignment over to ref2
///
pub fn liftover_read_alignment(
    ref1_to_ref2_map: &ReadToRefTreeMap,
    mut ref1_cigar_segment_start_pos: i64,
    ref1_cigar: &[Cigar],
) -> Option<(i64, Vec<Cigar>)> {
    let mut ref2_start_pos = None;
    let mut ref2_end_pos = None;

    // walk through ref1_cigar, dynamically building ref2_cigar as we go
    //
    // The strategy involves creating lots of small cigar blocks, and compressing adjacent segments of the same type at the end.
    //
    // In the ref1 cigar walk, we divide handling between cigar segments that only consume read length without advancing the reference (like soft-clip and insert),
    // and all others.
    //
    // Segments which only consume read length can be directly transferred from the ref1 to the ref2 cigar.
    //
    let mut ref2_cigar = Vec::new();
    for ref1_cigar_segment in ref1_cigar.iter() {
        match ref1_cigar_segment {
            Cigar::Ins(_) | Cigar::SoftClip(_) | Cigar::HardClip(_) => {
                // Cigar segments that only consume read-length can be directly transferred from ref1 to ref2 cigar
                ref2_cigar.push(*ref1_cigar_segment);
            }
            Cigar::Diff(cigar_segment_len)
            | Cigar::Equal(cigar_segment_len)
            | Cigar::Match(cigar_segment_len)
            | Cigar::Del(cigar_segment_len)
            | Cigar::RefSkip(cigar_segment_len) => {
                // Cigar segments that consume reference-length all need to be reinterpreted from ref1 to ref coordinates

                // Iterate through "blocks" of the ref1-ref2 alignment that correspond to the read-ref1 alignment component
                // represented by the current cigar segment
                //
                let mut last_ref1_ref2_mapping_block: Option<(usize, Option<i64>)> = None;

                // ref1_cigar_segment_start_pos tracks the start of each ref1 cigar segment block as we walk through the cigar,
                // likewise block_ref1_pos marks the start of each ref1->ref2 block within the segment range as we iterate through
                // the blocks intersecting each cigar segment
                //
                let mut block_ref1_pos = ref1_cigar_segment_start_pos;
                let ref1_cigar_segment_end_pos =
                    ref1_cigar_segment_start_pos + *cigar_segment_len as i64;

                for this_ref1_ref2_mapping_block in ref1_to_ref2_map
                    .get_ref_range(
                        ref1_cigar_segment_start_pos as usize,
                        ref1_cigar_segment_end_pos as usize,
                    )
                    .map(|(&x, &y)| (x, y))
                {
                    update_ref2_cigar_segment(
                        Some(this_ref1_ref2_mapping_block),
                        last_ref1_ref2_mapping_block,
                        ref1_cigar_segment_end_pos,
                        ref1_cigar_segment,
                        &mut block_ref1_pos,
                        &mut ref2_start_pos,
                        &mut ref2_end_pos,
                        &mut ref2_cigar,
                    );

                    last_ref1_ref2_mapping_block = Some(this_ref1_ref2_mapping_block);
                }

                update_ref2_cigar_segment(
                    None,
                    last_ref1_ref2_mapping_block,
                    ref1_cigar_segment_end_pos,
                    ref1_cigar_segment,
                    &mut block_ref1_pos,
                    &mut ref2_start_pos,
                    &mut ref2_end_pos,
                    &mut ref2_cigar,
                );
            }
            Cigar::Pad(_) => (),
        };
        update_ref_pos(ref1_cigar_segment, &mut ref1_cigar_segment_start_pos);
    }

    ref2_start_pos.map(|x| {
        let shift = clean_up_cigar_edge_indels(&mut ref2_cigar);
        let ref2_cigar = compress_cigar(&ref2_cigar);
        (x + shift as i64, ref2_cigar)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_vc_utils::get_read_segment_to_ref_pos_tree_map;

    #[test]
    fn test_liftover_read_alignment_basic() {
        use Cigar::*;
        let read_to_ref1_cigar = vec![
            Match(10),
            Del(10),
            Match(10),
            Ins(10),
            Match(10),
            SoftClip(10),
        ];

        // Test case 1, ref1 doesn't map to ref2:
        let ref1_to_ref2_map = ReadToRefTreeMap::default();
        let result = liftover_read_alignment(&ref1_to_ref2_map, 10, &read_to_ref1_cigar);

        assert!(result.is_none());

        // Test case 2, ref1 is a simple exact match to ref2
        let ref1_to_ref2_cigar = vec![Match(100)];
        let ref1_to_ref2_map =
            get_read_segment_to_ref_pos_tree_map(1000, &ref1_to_ref2_cigar, false);
        let result = liftover_read_alignment(&ref1_to_ref2_map, 10, &read_to_ref1_cigar);

        let (ref2_pos, read_to_ref2_cigar) = result.unwrap();
        assert_eq!(ref2_pos, 1010);
        assert_eq!(read_to_ref2_cigar, read_to_ref1_cigar);

        // Test case 3, ref1 has dels compared to ref2
        let ref1_to_ref2_cigar = vec![
            Match(10),
            Del(1),
            Match(5),
            Del(1),
            Match(5),
            Del(1),
            Match(5),
            Del(1),
            Match(10),
            Del(1),
            Match(15),
            Del(1),
            Match(50),
        ];
        let ref1_to_ref2_map =
            get_read_segment_to_ref_pos_tree_map(1000, &ref1_to_ref2_cigar, false);
        let result = liftover_read_alignment(&ref1_to_ref2_map, 10, &read_to_ref1_cigar);

        let (ref2_pos, read_to_ref2_cigar) = result.unwrap();
        assert_eq!(ref2_pos, 1011);
        assert_eq!(
            read_to_ref2_cigar,
            vec![
                Match(5),
                Del(1),
                Match(5),
                Del(12),
                Match(5),
                Del(1),
                Match(5),
                Ins(10),
                Match(10),
                SoftClip(10),
            ]
        );

        // Test case 4, ref1 has ins compared to ref2
        let ref1_to_ref2_cigar = vec![
            Match(5),
            Ins(10),
            Match(10),
            Ins(5),
            Match(5),
            Ins(5),
            Match(3),
            Ins(5),
            Match(1),
            Ins(5),
            Match(46),
        ];
        let ref1_to_ref2_map =
            get_read_segment_to_ref_pos_tree_map(1000, &ref1_to_ref2_cigar, false);
        let result = liftover_read_alignment(&ref1_to_ref2_map, 10, &read_to_ref1_cigar);

        let (ref2_pos, read_to_ref2_cigar) = result.unwrap();
        assert_eq!(ref2_pos, 1005);
        assert_eq!(
            read_to_ref2_cigar,
            vec![
                SoftClip(5),
                Match(5),
                Del(5),
                Match(5),
                Ins(15),
                Match(3),
                Ins(5),
                Match(1),
                SoftClip(11),
            ]
        );
    }

    #[test]
    fn test_liftover_read_alignment_leading_clip() {
        use Cigar::*;
        let read_to_ref1_cigar = vec![
            Match(10),
            Del(10),
            Match(10),
            Ins(10),
            Match(10),
            SoftClip(10),
        ];

        // Test case 1, ref1 aligns with a leading clip against ref2
        let ref1_to_ref2_cigar = vec![SoftClip(30), Match(70)];
        let ref1_to_ref2_map =
            get_read_segment_to_ref_pos_tree_map(1000, &ref1_to_ref2_cigar, false);

        let result = liftover_read_alignment(&ref1_to_ref2_map, 0, &read_to_ref1_cigar);

        let (ref2_pos, read_to_ref2_cigar) = result.unwrap();
        assert_eq!(ref2_pos, 1000);
        assert_eq!(
            read_to_ref2_cigar,
            vec![SoftClip(30), Match(10), SoftClip(10),]
        );

        // Test case 2, ref1 aligns with a leading clip against ref2, but harder this time
        let ref1_to_ref2_cigar = vec![SoftClip(10), Match(90)];
        let ref1_to_ref2_map =
            get_read_segment_to_ref_pos_tree_map(1000, &ref1_to_ref2_cigar, false);

        let result = liftover_read_alignment(&ref1_to_ref2_map, 5, &read_to_ref1_cigar);

        let (ref2_pos, read_to_ref2_cigar) = result.unwrap();
        assert_eq!(ref2_pos, 1000);
        assert_eq!(
            read_to_ref2_cigar,
            vec![
                SoftClip(5),
                Match(5),
                Del(10),
                Match(10),
                Ins(10),
                Match(10),
                SoftClip(10),
            ]
        );

        // Test case 3, ref1 aligns with a leading clip against ref2, and a read-to-ref1 deletion spans the alignment start
        let read_to_ref1_cigar = vec![Match(10), Del(10), Match(10)];

        let ref1_to_ref2_cigar = vec![SoftClip(20), Match(90)];
        let ref1_to_ref2_map =
            get_read_segment_to_ref_pos_tree_map(1000, &ref1_to_ref2_cigar, false);

        let result = liftover_read_alignment(&ref1_to_ref2_map, 5, &read_to_ref1_cigar);

        let (ref2_pos, read_to_ref2_cigar) = result.unwrap();
        assert_eq!(ref2_pos, 1005);
        assert_eq!(read_to_ref2_cigar, vec![SoftClip(10), Match(10)]);
    }

    #[test]
    fn test_liftover_read_alignment_trailing_clip() {
        use Cigar::*;

        // Test case 1, ref1 aligns with a trailing clip from ref2
        let read_to_ref1_cigar = vec![Match(10), Del(10), Match(10)];

        let ref1_to_ref2_cigar = vec![Match(70), SoftClip(30)];
        let ref1_to_ref2_map =
            get_read_segment_to_ref_pos_tree_map(1000, &ref1_to_ref2_cigar, false);

        let result = liftover_read_alignment(&ref1_to_ref2_map, 45, &read_to_ref1_cigar);

        let (ref2_pos, read_to_ref2_cigar) = result.unwrap();
        assert_eq!(ref2_pos, 1045);
        assert_eq!(
            read_to_ref2_cigar,
            vec![Match(10), Del(10), Match(5), SoftClip(5),]
        );

        // Test case 2, ref1 aligns with a trailing clip from ref2,with insertion on boundary
        let read_to_ref1_cigar = vec![Match(10), Ins(10), Match(10)];

        let ref1_to_ref2_cigar = vec![Match(70), SoftClip(30)];
        let ref1_to_ref2_map =
            get_read_segment_to_ref_pos_tree_map(1000, &ref1_to_ref2_cigar, false);

        let result = liftover_read_alignment(&ref1_to_ref2_map, 60, &read_to_ref1_cigar);

        let (ref2_pos, read_to_ref2_cigar) = result.unwrap();
        assert_eq!(ref2_pos, 1060);
        assert_eq!(read_to_ref2_cigar, vec![Match(10), SoftClip(20),]);

        // Test case 3, ref1 aligns with a trailing clip from ref2, with a deletion spanning the clip boundary
        let read_to_ref1_cigar = vec![Match(10), Del(10), Match(10)];

        let ref1_to_ref2_cigar = vec![Match(70), SoftClip(30)];
        let ref1_to_ref2_map =
            get_read_segment_to_ref_pos_tree_map(1000, &ref1_to_ref2_cigar, false);
        /*
        let debug = true;
        if debug {
            eprintln!("TEST r1r2 map:");
            for x in ref1_to_ref2_map.get_map().iter().enumerate() {
                eprintln!("{x:?}");
            }
        }
        */

        let result = liftover_read_alignment(&ref1_to_ref2_map, 55, &read_to_ref1_cigar);

        let (ref2_pos, read_to_ref2_cigar) = result.unwrap();
        assert_eq!(ref2_pos, 1055);
        assert_eq!(read_to_ref2_cigar, vec![Match(10), SoftClip(10),]);
    }
}
