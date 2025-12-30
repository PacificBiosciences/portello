use std::collections::{BTreeMap, BTreeSet};
use std::ops::Deref;
use std::sync::mpsc::channel;

use camino::Utf8Path;
use log::info;
use rayon::ThreadPool;
use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::{self, Read};
use rust_vc_utils::cigar::{get_cigar_read_offset, update_ref_and_read_pos};
use rust_vc_utils::read_to_ref_map::{RefToReadTreeMap, get_read_segment_to_ref_pos_tree_map};
use rust_vc_utils::{ChromList, IntRange, SeqOrderSplitReadSegment, get_alignment_end, rev_comp};
use unwrap::unwrap;

use super::write_phased_het_vcf::{VcfHetSite, write_phased_het_vcf};
use crate::contig_alignment_scanner::{
    AllContigMappingInfo, ChromPhaseSetMap, ContigMappingSegmentInfo, ContigSegmentInfo,
    clear_segment_phase_sets,
};
use crate::region_map_anyvalue::RegionMapAnyValue;
use crate::worker_thread_data::BamReaderWorkerThreadDataSet;

/// Represents a variant from the reference
///
#[derive(Debug)]
struct VarSite {
    ref_range: IntRange,

    /// These are segment coordinates in seq order direction
    seg_seq_order_range: IntRange,

    /// This is the contig segment arranged back into reference alignment direction (similar to VCF ALT field for the variant)
    seg_ref_order_seq: Vec<u8>,
}

fn get_seg_ref_order_seq(contig_seq: &[u8], seg_range: &IntRange, is_fwd_strand: bool) -> Vec<u8> {
    let start = std::cmp::max(seg_range.start, 0) as usize;
    let end = std::cmp::max(seg_range.end, 0) as usize;
    let x = &contig_seq[start..end];
    if is_fwd_strand {
        x.to_vec()
    } else {
        rev_comp(x)
    }
}

/// Convert segment alignment into a vector of variant ranges
///
/// Output will always be sorted low->high on the reference genome coordinate
///
fn get_seg_variants(seg: &SeqOrderSplitReadSegment, contig_seq: &[u8]) -> Vec<VarSite> {
    use Cigar::*;

    let ignore_hard_clip = false;

    let contig_len = contig_seq.len();

    let mut var_sites = Vec::new();

    let mut ref_pos = seg.pos;
    let mut read_pos = 0;

    let mut indel_start_pos = None;
    let mut ins_len = 0;
    let mut del_len = 0;

    for c in seg.cigar.iter() {
        let mut process_indel = || {
            if let Some((indel_ref_pos, indel_read_pos)) = indel_start_pos {
                let ref_range = IntRange::from_pair(indel_ref_pos, indel_ref_pos + del_len);
                let aln_range =
                    IntRange::from_pair(indel_read_pos as i64, indel_read_pos as i64 + ins_len);
                let seg_range = translate_alignment_range_to_seq_order_read_range(
                    seg.is_fwd_strand,
                    contig_len,
                    aln_range,
                );
                let seg_ref_order_seq =
                    get_seg_ref_order_seq(contig_seq, &seg_range, seg.is_fwd_strand);

                var_sites.push(VarSite {
                    ref_range,
                    seg_seq_order_range: seg_range,
                    seg_ref_order_seq,
                });
            }
            indel_start_pos = None;
            ins_len = 0;
            del_len = 0;
        };

        match c {
            Del(len) => {
                if *len > 0 {
                    if indel_start_pos.is_none() {
                        indel_start_pos = Some((ref_pos, read_pos));
                    }
                    del_len += *len as i64;
                }
            }
            Ins(len) => {
                if *len > 0 {
                    if indel_start_pos.is_none() {
                        indel_start_pos = Some((ref_pos, read_pos));
                    }
                    ins_len += *len as i64;
                }
            }
            Diff(len) => {
                process_indel();

                let ref_range = IntRange::from_pair(ref_pos, ref_pos + *len as i64);
                let aln_range = IntRange::from_pair(read_pos as i64, read_pos as i64 + *len as i64);
                let seg_range = translate_alignment_range_to_seq_order_read_range(
                    seg.is_fwd_strand,
                    contig_len,
                    aln_range,
                );
                let seg_seq = get_seg_ref_order_seq(contig_seq, &seg_range, seg.is_fwd_strand);

                var_sites.push(VarSite {
                    ref_range,
                    seg_seq_order_range: seg_range,
                    seg_ref_order_seq: seg_seq,
                });
            }
            Match(_) => {
                panic!(
                    "Method assumes alignment CIGAR strings use seq match/mismatch (=/X) instead of alignment match (M)"
                );
            }
            _ => {
                process_indel();
            }
        }

        update_ref_and_read_pos(c, &mut ref_pos, &mut read_pos, ignore_hard_clip);
    }

    var_sites
}

#[derive(Clone, Debug, Default)]
struct TmpHetInfo {
    seg0_info: Option<(IntRange, Vec<u8>)>,
    seg1_info: Option<(IntRange, Vec<u8>)>,
}

type VarMap = BTreeMap<IntRange, TmpHetInfo>;

/// Identify and remove length 1 indels
///
/// Smaller indels can be more accurately called by a dedicated small variant caller like DeepVariant, but it's
/// challenging to phase directly off of these variants as detected by an assembler consensusing approach.
///
fn remove_len1_indels_from_var_map(var_map: &mut VarMap) {
    let mut indel_len1_list = Vec::new();
    for (ref_range, thi) in var_map.iter() {
        let rl = ref_range.size();
        match (&thi.seg0_info, &thi.seg1_info) {
            (Some((seg_range, _)), None) => {
                if (seg_range.size() - rl).abs() == 1 {
                    indel_len1_list.push(ref_range.clone());
                }
            }
            (None, Some((seg_range, _))) => {
                if (seg_range.size() - rl).abs() == 1 {
                    indel_len1_list.push(ref_range.clone());
                }
            }
            (Some((seg_range0, _)), Some((seg_range1, _))) => {
                if (seg_range0.size() - seg_range1.size()).abs() == 1 {
                    indel_len1_list.push(ref_range.clone());
                }
            }
            _ => (),
        }
    }

    for k in indel_len1_list.iter() {
        var_map.remove(k);
    }
}

/// Identify overlapping variants when the ref_range isn't identical
///
/// Overlapping ref_range can be transitive and potentially form some complex problems. To make progress without getting
/// stuck on corner cases, we only handle overlaps covering exactly two ref_range entries.
///
/// TODO: This function is unfinished WIP
///
#[allow(unused)]
fn merge_simple_variant_overlaps(ref_chrom_seq: &[u8], var_map: &mut VarMap) {
    #[derive(Clone)]
    struct OlapInfo {
        ref_range: IntRange,
        seg: Option<usize>,
        thi: TmpHetInfo,
    }

    let mut ref_range_clusters_to_merge = Vec::new();
    {
        let mut olap_end = 0;
        let mut olap_info_list = Vec::<OlapInfo>::new();
        for (ref_range, thi) in var_map.iter() {
            if ref_range.start > olap_end {
                // start a new cluster
                //
                // In addition to restricting cluster count to 2, we restrict it to two
                // simple cases where only one each of seg0 and seg1 are defined.
                if olap_info_list.len() == 2
                    && olap_info_list[0].seg.is_some()
                    && olap_info_list[1].seg.is_some()
                    && olap_info_list[0].seg != olap_info_list[1].seg
                {
                    ref_range_clusters_to_merge.push(olap_info_list.clone());
                }
                olap_info_list.clear();
            }
            olap_end = std::cmp::max(olap_end, ref_range.end);
            let seg = if thi.seg1_info.is_none() {
                Some(0)
            } else if thi.seg0_info.is_none() {
                Some(1)
            } else {
                None
            };
            olap_info_list.push(OlapInfo {
                ref_range: ref_range.clone(),
                seg,
                thi: thi.clone(),
            });
        }
    }

    let mut new_merged_var_map_entries = VarMap::default();
    for cluster in ref_range_clusters_to_merge.iter() {
        let mut merged_ref_range = cluster[0].ref_range.clone();
        for olap_info in cluster.iter().skip(1) {
            merged_ref_range.merge(&olap_info.ref_range);
        }

        let mut prefix_size = 0;
        let mut suffix_size = 0;
        for olap_info in cluster.iter() {
            prefix_size = std::cmp::max(
                prefix_size,
                olap_info.ref_range.start - merged_ref_range.start,
            );
            suffix_size =
                std::cmp::max(suffix_size, merged_ref_range.end - olap_info.ref_range.end);
        }

        let prefix_seq = &ref_chrom_seq
            [merged_ref_range.start as usize..(merged_ref_range.start + prefix_size) as usize];
        let suffix_seq = &ref_chrom_seq
            [(merged_ref_range.end - prefix_size) as usize..merged_ref_range.end as usize];

        /*
        let (mut merged_seg0_range, mut merged_seg1_range) = if cluster[0].seg == Some(0) {
            (cluster[0].thi.seg0_info.unwrap().0.clone(),cluster[1].thi.seg1_info.unwrap().0.clone())
        } else {
            (cluster[1].thi.seg0_info.unwrap().0.clone(),cluster[0].thi.seg1_info.unwrap().0.clone())
        };

        let value = TmpHetInfo {
            seg0_info: X,
            seg1_info: Y,
        };

        new_merged_var_map_entries.insert(merged_ref_range, value);
        */
    }

    // Remove old entries which have been merged:
    for cluster in ref_range_clusters_to_merge {
        for olap_info in cluster {
            var_map.remove(&olap_info.ref_range);
        }
    }

    // Add new merged entries
    for (k, v) in new_merged_var_map_entries {
        var_map.insert(k, v);
    }
}

fn get_het_var_map(seg0_vars: &[VarSite], seg1_vars: &[VarSite], debug: bool) -> VarMap {
    let mut var_map = VarMap::default();

    // First join variants using ref-range key:
    for x in seg0_vars.iter() {
        let entry = var_map.entry(x.ref_range.clone()).or_default();
        entry.seg0_info = Some((x.seg_seq_order_range.clone(), x.seg_ref_order_seq.clone()));
    }

    for x in seg1_vars.iter() {
        let entry = var_map.entry(x.ref_range.clone()).or_default();
        entry.seg1_info = Some((x.seg_seq_order_range.clone(), x.seg_ref_order_seq.clone()));
    }

    // Identify and remove (probable) homozygous variants from the set
    //
    let hom_list = {
        let mut x = Vec::new();
        for (ref_range, thi) in var_map.iter() {
            if let (Some((_, s0)), Some((_, s1))) = (&thi.seg0_info, &thi.seg1_info)
                && s0 == s1
            {
                x.push(ref_range.clone());
            }
        }
        x
    };

    if debug {
        eprintln!("var_map_len_before: {}", var_map.len());
        eprintln!("hom_list_len: {}", hom_list.len());
    }

    for k in hom_list.iter() {
        var_map.remove(k);
    }

    if debug {
        eprintln!("var_map_len_after: {}", var_map.len());

        eprintln!("----var_map after head:");
        for x in var_map.iter().take(20) {
            eprintln!("{x:?}");
        }
    }

    // Identify and remove length 1 indels
    remove_len1_indels_from_var_map(&mut var_map);

    // Identify overlapping variants when the ref_range isn't identical
    //merge_simple_variant_overlaps(ref_chrom_seq, &mut var_map);

    var_map
}

/// All data associated with each heterozygous variant site within a phase set block
///
/// The range for either segment are contig coordinates in the contig forward strand orientation
///
#[derive(Clone, Debug)]
struct HetSite {
    seg0_seq_order_range: IntRange,
    seg1_seq_order_range: IntRange,

    ref_range: IntRange,

    /// Haplotype sequence in the segment sequency order direction
    seg0_seq_order_seq: Vec<u8>,

    /// Haplotype sequence in the segment sequency order direction
    seg1_seq_order_seq: Vec<u8>,

    /// Haplotype sequence in the reference alignment direction. None if this haplotype is ref. Only used for phased het
    /// vcf output.
    seg0_ref_order_seq: Option<Vec<u8>>,

    /// Haplotype sequence in the reference alignment direction. None if this haplotype is ref. Only used for phased het
    /// vcf output.
    seg1_ref_order_seq: Option<Vec<u8>>,
}

impl HetSite {
    fn get_seg_sor(&self, is_seg0: bool) -> &IntRange {
        if is_seg0 {
            &self.seg0_seq_order_range
        } else {
            &self.seg1_seq_order_range
        }
    }

    fn get_seg_sos(&self, is_seg0: bool) -> &[u8] {
        if is_seg0 {
            &self.seg0_seq_order_seq
        } else {
            &self.seg1_seq_order_seq
        }
    }
}

/// Convert from alignment order range to sequence order read range
///
fn translate_alignment_range_to_seq_order_read_range(
    is_fwd_strand: bool,
    contig_len: usize,
    mut aln_range: IntRange,
) -> IntRange {
    if !is_fwd_strand {
        aln_range.reverse(contig_len as i64);
    }
    aln_range
}

/// Process a het site found on two contig alignment segments to the reference
///
/// This version handles the case where segA has a variant, and segB does not
///
#[allow(clippy::too_many_arguments)]
fn process_simple_het_site(
    seg_a_is_fwd_strand: bool,
    seg_b_info: &ContigMappingSegmentInfo,
    contig_b_seq: &[u8],
    seg_a_seq_order_range: &IntRange,
    seg_a_ref_order_seq: &[u8],
    ref_range: &IntRange,
    seg_a_is_seg0: bool,
    debug: bool,
) -> Option<HetSite> {
    if debug {
        eprintln!(
            "process_simple_het_site seg_a_seq_order_range {seg_a_seq_order_range:?} ref_range {ref_range:?}"
        );
    }

    let seg_a_seq_order_seq = if seg_a_is_fwd_strand {
        seg_a_ref_order_seq.to_vec()
    } else {
        rev_comp(seg_a_ref_order_seq)
    };

    if debug {
        eprintln!("seg_b c2r dump");
        seg_b_info.contig_to_ref_map.dump();
        eprintln!("seg_b r2c dump");
        seg_b_info.ref_to_contig_map.dump();
    }

    // Use ref_to_contig_map to translate ref_range back to the reference haplotype coordinate in the contig aligned in seg1;
    let seg_b_start = seg_b_info.ref_to_contig_map.get_read_pos(ref_range.start);
    let seg_b_end = seg_b_info.ref_to_contig_map.get_read_pos(ref_range.end);

    if debug {
        eprintln!("seg_b_start {seg_b_start:?} seg_b_end {seg_b_end:?}");
    }

    if let (Some(seg_b_start), Some(seg_b_end)) = (seg_b_start, seg_b_end) {
        let seg_b_range = IntRange::from_pair(seg_b_start as i64, seg_b_end as i64);
        // Translate segB_range into contig seq order coordinates:
        let seg_b_seq_order_range = translate_alignment_range_to_seq_order_read_range(
            seg_b_info.seq_order_segment.is_fwd_strand,
            contig_b_seq.len(),
            seg_b_range,
        );

        if debug {
            eprintln!("seg_b_seq_order_range {seg_b_seq_order_range:?}");
        }

        let seg_b_seq_order_seq = contig_b_seq
            [seg_b_seq_order_range.start as usize..seg_b_seq_order_range.end as usize]
            .to_vec();

        let (
            seg0_seq_order_range,
            seg0_seq_order_seq,
            seg0_ref_order_seq,
            seg1_seq_order_range,
            seg1_seq_order_seq,
            seg1_ref_order_seq,
        ) = if seg_a_is_seg0 {
            (
                seg_a_seq_order_range.clone(),
                seg_a_seq_order_seq,
                Some(seg_a_ref_order_seq.to_vec()),
                seg_b_seq_order_range,
                seg_b_seq_order_seq,
                None,
            )
        } else {
            (
                seg_b_seq_order_range,
                seg_b_seq_order_seq,
                None,
                seg_a_seq_order_range.clone(),
                seg_a_seq_order_seq,
                Some(seg_a_ref_order_seq.to_vec()),
            )
        };

        Some(HetSite {
            seg0_seq_order_range,
            seg1_seq_order_range,
            ref_range: ref_range.clone(),
            seg0_seq_order_seq,
            seg1_seq_order_seq,
            seg0_ref_order_seq,
            seg1_ref_order_seq,
        })
    } else {
        // TODO - add more sophisticated handling for Nones here later:
        None
    }
}

/// Given two overlapping contig split alignment segments from the same phase set, find their heterozygous variants
///
/// Het variant locations are reported back in the fwd-strand coordinates of each contig
///
fn get_segment_pair_het_sites(
    seg0_info: &ContigMappingSegmentInfo,
    seg1_info: &ContigMappingSegmentInfo,
    contig0_seq: &[u8],
    contig1_seq: &[u8],
    debug: bool,
) -> Vec<HetSite> {
    let seg0 = &seg0_info.seq_order_segment;
    let seg1 = &seg1_info.seq_order_segment;

    assert_eq!(seg0.chrom_index, seg1.chrom_index);

    if debug {
        eprintln!("seg0_info: hap_id {:?}", seg0_info.hap_id);
        eprintln!("seg1_info: hap_id {:?}", seg1_info.hap_id);
        eprintln!("seg0: {}", seg0.short_display());
        eprintln!("seg1: {}", seg1.short_display());

        eprintln!(
            "seg0_cigar_len {}",
            get_cigar_read_offset(&seg0.cigar, false)
        );
        eprintln!(
            "seg1_cigar_len {}",
            get_cigar_read_offset(&seg1.cigar, false)
        );
        eprintln!("contig0_seq_len: {}", contig0_seq.len());
        eprintln!("contig1_seq_len: {}", contig1_seq.len());
    }

    // Step 1. Get variant map from each contig alignment segment
    //
    let seg0_vars = get_seg_variants(seg0, contig0_seq);
    let seg1_vars = get_seg_variants(seg1, contig1_seq);

    if debug {
        eprintln!("seg0_var_count: {}", seg0_vars.len());
        eprintln!("seg1_var_count: {}", seg1_vars.len());

        let head_count = 20;
        eprintln!("---- seg0_vars head:");
        for x in seg0_vars.iter().take(head_count) {
            eprintln!("{x:?}");
        }
        eprintln!("---- seg1_vars head:");
        for x in seg1_vars.iter().take(head_count) {
            eprintln!("{x:?}");
        }
    }

    // Step 2. Join variants on common reference genome coordinates, and remove homozygous variants.
    //
    let var_map = get_het_var_map(&seg0_vars, &seg1_vars, debug);

    // Step 3. Identify and consolidate heterozygous variants
    //
    let mut het_sites = Vec::new();

    for (ref_range, thi) in var_map.iter() {
        assert!(thi.seg0_info.is_some() || thi.seg1_info.is_some());

        let het_site = match (&thi.seg0_info, &thi.seg1_info) {
            (Some((seg0_range, seg0_ref_order_seq)), None) => process_simple_het_site(
                seg0.is_fwd_strand,
                seg1_info,
                contig1_seq,
                seg0_range,
                seg0_ref_order_seq,
                ref_range,
                true,
                false,
            ),
            (None, Some((seg1_range, seg1_ref_order_seq))) => process_simple_het_site(
                seg1.is_fwd_strand,
                seg0_info,
                contig0_seq,
                seg1_range,
                seg1_ref_order_seq,
                ref_range,
                false,
                false,
            ),
            (Some((seg0_range, seg0_seq)), Some((seg1_range, seg1_seq))) => {
                let seg0_seq_order_seq = if seg0.is_fwd_strand {
                    seg0_seq.clone()
                } else {
                    rev_comp(seg0_seq)
                };
                let seg1_seq_order_seq = if seg1.is_fwd_strand {
                    seg1_seq.clone()
                } else {
                    rev_comp(seg1_seq)
                };
                Some(HetSite {
                    seg0_seq_order_range: seg0_range.clone(),
                    seg1_seq_order_range: seg1_range.clone(),
                    ref_range: ref_range.clone(),
                    seg0_seq_order_seq,
                    seg1_seq_order_seq,
                    seg0_ref_order_seq: Some(seg0_seq.clone()),
                    seg1_ref_order_seq: Some(seg1_seq.clone()),
                })
            }
            _ => None,
        };

        if let Some(het_site) = het_site {
            het_sites.push(het_site);
        }
    }

    het_sites
}

#[derive(Clone, Default)]
struct ReadSpanInfo {
    /// qname is only defined for debugging
    qname: Option<String>,

    matched_alleles: BTreeSet<IntRange>,
}

type ContigSegmentReadSpanInfo = RegionMapAnyValue<ReadSpanInfo>;

/// Returns informations about how reads span and support different het variant sites over the contig alignment segment
///
/// # Arguments
/// * `include_qnames` - include qnames in return data structure for downstream debug printing. The qnames are not
///   needed for normal phasing logic.
///
fn get_contig_read_span_info(
    bam_reader: &mut bam::IndexedReader,
    het_sites: &[HetSite],
    is_seg0: bool,
    contig_index: usize,
    contig_seg: &SeqOrderSplitReadSegment,
    include_qnames: bool,
) -> ContigSegmentReadSpanInfo {
    let ignore_hard_clip = false;

    let seg_eval_sites = get_seg_eval_sites(het_sites, is_seg0);

    let mut spans = ContigSegmentReadSpanInfo::default();

    bam_reader
        .fetch(bam::FetchDefinition::Region(
            contig_index as i32,
            contig_seg.seq_order_read_start as i64,
            contig_seg.seq_order_read_end as i64,
        ))
        .unwrap();

    let mut record = bam::Record::new();
    while let Some(r) = bam_reader.read(&mut record) {
        unwrap!(r, "Failed to parse alignment record");

        let contig_start = record.pos();
        let contig_end = get_alignment_end(&record);

        let read_seq = record.seq().as_bytes();
        let read_len = read_seq.len();

        let read_to_contig_map =
            get_read_segment_to_ref_pos_tree_map(contig_start, &record.cigar(), ignore_hard_clip);
        let contig_to_read_map = RefToReadTreeMap::from_read_to_ref_tree_map(&read_to_contig_map);

        let debug_read = false;
        let mut value = ReadSpanInfo::default();
        if include_qnames {
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            value.qname = Some(qname);
        }
        for x in seg_eval_sites.find_overlaps(contig_start, contig_end) {
            let data = x.data();
            let contig_variant_interval = &data.seg_seq_order_range;

            // Get the read sequence across contig interval
            let read_variant_start = contig_to_read_map.get_read_pos(contig_variant_interval.start);
            let read_variant_end = contig_to_read_map.get_read_pos(contig_variant_interval.end);

            if debug_read {
                eprintln!(
                    "contig_variant_interval: {contig_variant_interval:?} read_variant_start {read_variant_start:?} read_variant_end {read_variant_end:?}"
                );
            }

            if let (Some(read_start), Some(read_end)) = (read_variant_start, read_variant_end)
                && read_start <= read_len
                && read_end <= read_len
            {
                let allele_seq = read_seq[read_start..read_end].to_vec();

                if debug_read {
                    let a = std::str::from_utf8(&allele_seq).unwrap().to_string();
                    let d = std::str::from_utf8(&data.seg_seq_order_seq)
                        .unwrap()
                        .to_string();
                    eprintln!("allele_seq {a} data_seq {d}");
                }

                if allele_seq == data.seg_seq_order_seq {
                    value.matched_alleles.insert(data.ref_range.clone());
                }
            }
        }

        spans.add_region_value(contig_start, contig_end, value);
    }

    spans
}

fn is_het_pair_phased_by_seg_reads(
    contig_segment_read_span_info: &ContigSegmentReadSpanInfo,
    contig_seg: &SeqOrderSplitReadSegment,
    last_het_site: &HetSite,
    het_site: &HetSite,
    is_seg0: bool,
) -> bool {
    // Require that read not only touch each variant, but extend into it by a few bases:
    let min_overlap_dist = 5;

    // define the minimum interval a read must span to consider this het pair phased:
    let (max_start, min_end) = if contig_seg.is_fwd_strand {
        (
            last_het_site.get_seg_sor(is_seg0).start - min_overlap_dist,
            het_site.get_seg_sor(is_seg0).end + min_overlap_dist,
        )
    } else {
        (
            het_site.get_seg_sor(is_seg0).start - min_overlap_dist,
            last_het_site.get_seg_sor(is_seg0).end + min_overlap_dist,
        )
    };

    if max_start > min_end {
        let seg_index = if is_seg0 { 0 } else { 1 };
        panic!(
            "Unexpected het-site arrangement for seg{seg_index}. max_start {max_start} min_end {min_end}\n\
            last_het_site: {last_het_site:?}\n\
            het_site: {het_site:?}\n\
            contig_seg {}",
            contig_seg.short_display()
        );
    }

    let debug = false;
    for x in contig_segment_read_span_info.find_overlaps(max_start, min_end) {
        if x.interval().start <= max_start && x.interval().end >= min_end {
            // Now test that the read supports both of the expected het alleles:
            let data = x.data();
            if debug {
                eprintln!("Passing span found. qname {:?}", data.qname);
                eprintln!(
                    "last Ref-range {:?} found: {}",
                    last_het_site.ref_range,
                    data.matched_alleles.contains(&last_het_site.ref_range)
                );
                eprintln!(
                    "Ref-range {:?} found: {}",
                    het_site.ref_range,
                    data.matched_alleles.contains(&het_site.ref_range)
                );
            }
            if data.matched_alleles.contains(&het_site.ref_range)
                && data.matched_alleles.contains(&last_het_site.ref_range)
            {
                if debug {
                    eprintln!("Passing alleles found");
                }
                return true;
            }
        }
    }

    false
}

/// Return the intersection of 2 ranges, or None if the ranges don't interect at all
///
/// Allow for a zero-sized range to be contained inside of a longer range to still count as an intersection
///
fn intersect_ranges(a: &IntRange, b: &IntRange) -> Option<IntRange> {
    let start = std::cmp::max(a.start, b.start);
    let end = std::cmp::min(a.end, b.end);
    if end > start || (end == start && (a.intersect_pos(start) || b.intersect_pos(start))) {
        Some(IntRange::from_pair(start, end))
    } else {
        None
    }
}

fn append_alt_seq(ref_prefix: &[u8], seg_seq: Option<&[u8]>, alt_seqs: &mut Vec<Vec<u8>>) {
    if let Some(s) = seg_seq {
        if ref_prefix.is_empty() {
            alt_seqs.push(s.to_vec());
        } else {
            let mut alt_seq = ref_prefix.to_vec();
            alt_seq.extend(s);
            alt_seqs.push(alt_seq);
        }
    }
}

fn get_vcf_het_site(
    contig_segment0: &ContigMappingSegmentInfo,
    contig_segment1: &ContigMappingSegmentInfo,
    phase_set_range: &IntRange,
    ref_chrom_seq: &[u8],
    het_site: &HetSite,
) -> Option<VcfHetSite> {
    intersect_ranges(&het_site.ref_range, phase_set_range)?;

    let is_ref_prefix = het_site.ref_range.size() == 0
        || het_site.seg0_seq_order_range.size() == 0
        || het_site.seg1_seq_order_range.size() == 0;

    if is_ref_prefix && het_site.ref_range.start <= 0 {
        return None;
    }

    let chrom_index = contig_segment0.seq_order_segment.chrom_index;

    let pos = het_site.ref_range.start + if is_ref_prefix { -1 } else { 0 };

    let ref_seq = ref_chrom_seq[pos as usize..het_site.ref_range.end as usize].to_vec();
    let ref_prefix = if is_ref_prefix {
        &ref_chrom_seq[pos as usize..pos as usize + 1]
    } else {
        &[]
    };

    let mut alt_seqs = Vec::new();
    append_alt_seq(
        ref_prefix,
        het_site.seg0_ref_order_seq.as_deref(),
        &mut alt_seqs,
    );
    append_alt_seq(
        ref_prefix,
        het_site.seg1_ref_order_seq.as_deref(),
        &mut alt_seqs,
    );

    let is_alt0_hp0 = if het_site.seg0_ref_order_seq.is_some() {
        contig_segment0.hap_id.is_none_or(|hp| hp == 0)
    } else {
        contig_segment1.hap_id.is_none_or(|hp| hp == 0)
    };

    Some(VcfHetSite {
        chrom_index,
        pos,
        ref_seq,
        alt_seqs,
        is_alt0_hp0,
        phase_set: None,
    })
}

/// After het-sites have been found from contigs and read span info enumerated, this method puts the information
/// together to generate the final read-backed phase-sets.
///
/// * `het_sites`` - het variant list, sorted in ref-range order
/// * `contig_segment0_read_span_info` - Spanning information across het sites for reads mapped to contig_segment0
/// * `contig_segment1_read_span_info` - Spanning information across het sites for reads mapped to contig_segment1
/// * `contig_segment0` - Alignment details for contig_segment0
/// * `contig_segment1` - Alignment details for contig_segment1
///
#[allow(clippy::too_many_arguments)]
fn get_updated_phase_sets(
    het_sites: &[HetSite],
    contig_segment0_read_span_info: &ContigSegmentReadSpanInfo,
    contig_segment1_read_span_info: &ContigSegmentReadSpanInfo,
    contig_segment0: &ContigMappingSegmentInfo,
    contig_segment1: &ContigMappingSegmentInfo,
    phase_set_range: &IntRange,
    is_phased_het_vcf: bool,
    ref_chrom_seq: &[u8],
) -> (Vec<IntRange>, Vec<VcfHetSite>) {
    // Phasing method: walk through the het_sites list, define PS range as going from first to last phased het in
    // sequence. Currently configured to require the minimum supporting evidence of a single read in either haplotype to
    // consider adjacent hets as phased.

    let mut phase_sets = Vec::new();
    let mut vcf_het_sites = Vec::new();

    let mut phase_set_start = None;
    let mut last_het_site: Option<HetSite> = None;
    let mut phased_het_site_group = Vec::new();
    for het_site in het_sites.iter() {
        if is_phased_het_vcf {
            // Push in basic het site information before we know if the group is phased, or None if
            // we don't want to attempt a VCF record for this het for various reasons.
            //
            let vcf_het_site = get_vcf_het_site(
                contig_segment0,
                contig_segment1,
                phase_set_range,
                ref_chrom_seq,
                het_site,
            );
            phased_het_site_group.push(vcf_het_site);
        }

        if let Some(last_het_site) = last_het_site {
            // test for phase support on seg0 or seg1
            let het_pair_phased = is_het_pair_phased_by_seg_reads(
                contig_segment0_read_span_info,
                &contig_segment0.seq_order_segment,
                &last_het_site,
                het_site,
                true,
            ) || is_het_pair_phased_by_seg_reads(
                contig_segment1_read_span_info,
                &contig_segment1.seq_order_segment,
                &last_het_site,
                het_site,
                false,
            );

            if het_pair_phased {
                if phase_set_start.is_none() {
                    phase_set_start = Some(last_het_site.ref_range.start);
                }
            } else {
                if let Some(phase_set_start) = phase_set_start {
                    let new_ps = IntRange::from_pair(phase_set_start, last_het_site.ref_range.end);
                    if let Some(x) = intersect_ranges(&new_ps, phase_set_range) {
                        phase_sets.push(x);
                    }
                }

                if is_phased_het_vcf {
                    // Annotate PS tag value and transfer phased_het_site_group into vcf output buffer
                    assert!(!phased_het_site_group.is_empty());
                    let drain_end = phased_het_site_group.len() - 1;
                    for mut x in phased_het_site_group.drain(..drain_end).flatten() {
                        x.phase_set = phase_set_start;
                        vcf_het_sites.push(x);
                    }
                }

                phase_set_start = None;
            }
        }
        last_het_site = Some(het_site.clone());
    }

    // Handle the final phase block in case phasing extends into the last het-site:
    if let (Some(last_het_site), Some(phase_set_start)) = (last_het_site, phase_set_start) {
        let new_ps = IntRange::from_pair(phase_set_start, last_het_site.ref_range.end);
        if let Some(x) = intersect_ranges(&new_ps, phase_set_range) {
            phase_sets.push(x);
        }
    }

    if is_phased_het_vcf {
        // Annotate PS tag value and transfer phased_het_site_group into vcf output buffer
        for mut x in phased_het_site_group.drain(..).flatten() {
            x.phase_set = phase_set_start;
            vcf_het_sites.push(x);
        }
    }

    (phase_sets, vcf_het_sites)
}

struct SegmentPhaseSets {
    contig_info0: ContigSegmentInfo,
    contig_info1: ContigSegmentInfo,
    phase_sets: Vec<IntRange>,

    /// These will only be populated if the is_phased_het_vcf option is true
    vcf_het_sites: Vec<VcfHetSite>,
}

#[derive(Clone, Default)]
struct SegSpanInfo {
    seg_seq_order_range: IntRange,
    seg_seq_order_seq: Vec<u8>,
    ref_range: IntRange,
}

type SegEvalSites = RegionMapAnyValue<SegSpanInfo>;

fn get_seg_eval_sites(het_sites: &[HetSite], is_seg0: bool) -> SegEvalSites {
    let mut seg_eval_sites = RegionMapAnyValue::<SegSpanInfo>::default();
    for het_site in het_sites.iter() {
        let seg_range = het_site.get_seg_sor(is_seg0);
        let start = seg_range.start;
        let end = seg_range.end;

        // The interval tree data structure below won't retrieve empty intervals, so extend these by one (this won't
        // mess up the allele match because we copy in seg_seq_order_range for the client code to use)
        //
        let end = if start == end { end + 1 } else { end };
        let value = SegSpanInfo {
            seg_seq_order_range: het_site.get_seg_sor(is_seg0).clone(),
            seg_seq_order_seq: het_site.get_seg_sos(is_seg0).to_vec(),
            ref_range: het_site.ref_range.clone(),
        };
        seg_eval_sites.add_region_value(start, end, value);
    }
    seg_eval_sites
}

/// Phase a region with diploid assembly contig coverage
///
/// Assuming a partially-phased assembly, we're given a region that appears to be a single diploid phase-block by contig
/// structure and reference mapping pattern, but assumed to deliberatly contain switch errors via the
/// partially-phased/dual assembly output scheme.
///
/// This routine finds putative het variant sites directly from the two assembly contigs, then searches for read support
/// across each adjacent het variant pair to create read-backed phase set assignments for the segment.
///
/// # Arguments
/// * `phase_set_range` - Output phase sets should be bound by this input range. Most of the time the this range matches
///   the contig span and has no impact, but there are several more complicated cases where this bound is important to
///   prevent artifacts, such as a single diploid region spanned by two contigs, but interupted with a third haplotype
///   in the middle.
/// * `is_phased_het_vcf` - If true, return additional information which can be merged to write a het small variant VCF
///   file, which can be used to debug phasing accuracy.
///
#[allow(clippy::too_many_arguments)]
fn update_phase_set_for_diploid_segment(
    assembly_contig_list: &ChromList,
    worker_thread_dataset: BamReaderWorkerThreadDataSet,
    all_contig_mapping_info: &AllContigMappingInfo,
    contig_info0: ContigSegmentInfo,
    contig_info1: ContigSegmentInfo,
    phase_set_range: &IntRange,
    is_phased_het_vcf: bool,
    reference: &[Vec<u8>],
) -> SegmentPhaseSets {
    let debug = false;

    if debug {
        eprintln!(
            "contig0_name: {}",
            assembly_contig_list.data[contig_info0.contig_index].label
        );
        eprintln!(
            "contig1_name: {}",
            assembly_contig_list.data[contig_info1.contig_index].label
        );
    }

    let contig0 = &all_contig_mapping_info[contig_info0.contig_index];
    let contig1 = &all_contig_mapping_info[contig_info1.contig_index];
    let contig_segment0 = &contig0.ordered_contig_segment_info[contig_info0.segment_index];
    let contig_segment1 = &contig1.ordered_contig_segment_info[contig_info1.segment_index];

    let contig0_seq = unwrap!(
        contig0.contig_seq.as_ref(),
        "Expected contig_seq to be defined"
    );
    let contig1_seq = unwrap!(
        contig1.contig_seq.as_ref(),
        "Expected contig_seq to be defined"
    );

    let chrom_index = contig_segment0.seq_order_segment.chrom_index;
    let ref_chrom_seq = &reference[chrom_index];

    let het_sites = get_segment_pair_het_sites(
        contig_segment0,
        contig_segment1,
        contig0_seq,
        contig1_seq,
        debug,
    );

    if debug {
        eprintln!("het_sites count: {}", het_sites.len());
        for het_site in het_sites.iter().take(20) {
            eprintln!(
                "HetSite: ref_range: {:?} seg0_range: {:?} seg1_range: {:?}",
                &het_site.ref_range, &het_site.seg0_seq_order_range, &het_site.seg1_seq_order_range
            );
        }
    }

    let worker_id = rayon::current_thread_index().unwrap();
    let bam_reader = &mut worker_thread_dataset[worker_id].lock().unwrap().bam_reader;

    // Qnames are only used for certain debug operations
    let include_qnames_in_read_span_info = debug;

    let contig_segment0_read_span_info = get_contig_read_span_info(
        bam_reader,
        &het_sites,
        true,
        contig_info0.contig_index,
        &contig_segment0.seq_order_segment,
        include_qnames_in_read_span_info,
    );
    let contig_segment1_read_span_info = get_contig_read_span_info(
        bam_reader,
        &het_sites,
        false,
        contig_info1.contig_index,
        &contig_segment1.seq_order_segment,
        include_qnames_in_read_span_info,
    );

    let (phase_sets, vcf_het_sites) = get_updated_phase_sets(
        &het_sites,
        &contig_segment0_read_span_info,
        &contig_segment1_read_span_info,
        contig_segment0,
        contig_segment1,
        phase_set_range,
        is_phased_het_vcf,
        ref_chrom_seq,
    );

    SegmentPhaseSets {
        contig_info0,
        contig_info1,
        phase_sets,
        vcf_het_sites,
    }
}

/// Refine phase sets based on read support across heterozygous variants
///
/// Initial phase set ranges are based on the input assembly contigs representing one haplotype. This method assumes the
/// contigs result form a partial-phasing/dual assembly process and contain switch errors. The goal here is to run a
/// variation on read-backed phasing (RBP) to subdivide the contig into smaller phase-set blocks in diploid regions
/// based on support read evidence across heterozygous variants. A substantial difference from normal RBP is that the
/// phasing determinations are based on read alignments to the contigs and not directly to the reference.
///
/// Note that this isn't necessarily supposed to be better than traditional mapping-based RBP, but provides a
/// convenience to make portello realigned read output more directly ready for direct biological inference.
///
#[allow(clippy::too_many_arguments)]
pub fn segment_phase_set_blocks(
    assembly_contig_list: &ChromList,
    worker_thread_dataset: &BamReaderWorkerThreadDataSet,
    worker_pool: &ThreadPool,
    all_contig_mapping_info: &mut AllContigMappingInfo,
    annotated_phase_set_ranges: Vec<ChromPhaseSetMap>,
    phased_het_vcf_prefix: Option<&Utf8Path>,
    ref_filename: &Utf8Path,
    reference: &[Vec<u8>],
    ref_chrom_list: &ChromList,
    thread_count: usize,
) {
    let debug = false;

    info!("Scanning read to contig alignments to refine phase set region boundaries");

    clear_segment_phase_sets(all_contig_mapping_info);

    let is_phased_het_vcf = phased_het_vcf_prefix.is_some();

    let (tx, rx) = channel();

    {
        let annotated_phase_set_ranges = &annotated_phase_set_ranges;
        let all_contig_mapping_info = all_contig_mapping_info.deref();
        worker_pool.scope(move |scope| {
            // Iterate through existing phase blocks, walk through contigs cigar strings to find het variants and test for phasing support from the
            // reads aligned to those contig segments.
            //
            for (chrom_index, chrom_annotated_phase_set_ranges) in
            annotated_phase_set_ranges.iter().enumerate()
            {
                for (phase_set_range, phase_set_contig_segments) in chrom_annotated_phase_set_ranges.iter() {
                    if debug {
                        eprintln!("chrom_index: {chrom_index} phase_set_range: {phase_set_range:?} phase_set_contig_segments: {phase_set_contig_segments:?}");
                    }

                    if phase_set_contig_segments.len() == 2 {
                        let worker_thread_dataset = worker_thread_dataset.clone();
                        let tx = tx.clone();

                        scope.spawn(move |_| {
                            let segment_phase_sets = update_phase_set_for_diploid_segment(
                                assembly_contig_list,
                                worker_thread_dataset,
                                all_contig_mapping_info,
                                phase_set_contig_segments[0],
                                phase_set_contig_segments[1],
                                phase_set_range,
                                is_phased_het_vcf,
                                reference,
                            );

                            tx.send(segment_phase_sets).unwrap();
                        });
                    }
                }
            }
        });
    }

    // Copy all new phase sets into all_contig_mapping_info:
    let mut vcf_het_sites = Vec::new();
    for mut segment_phase_sets in rx {
        if debug {
            eprintln!(
                "Applying finished phase set to contig0: {:?} contig1: {:?}",
                segment_phase_sets.contig_info0, segment_phase_sets.contig_info1,
            );
            for x in segment_phase_sets.phase_sets.iter() {
                eprintln!("phase_set {x:?}");
            }
        }

        // The phase sets are merged in, rather than replaced, because the same contig alignment segment can show up in
        // more than one diploid phase block, paired against different partner contig alignment segments
        //
        for range in segment_phase_sets.phase_sets.iter() {
            all_contig_mapping_info[segment_phase_sets.contig_info0.contig_index]
                .ordered_contig_segment_info[segment_phase_sets.contig_info0.segment_index]
                .phase_sets
                .add_region(range.start, range.end);
            all_contig_mapping_info[segment_phase_sets.contig_info1.contig_index]
                .ordered_contig_segment_info[segment_phase_sets.contig_info1.segment_index]
                .phase_sets
                .add_region(range.start, range.end);
        }

        if is_phased_het_vcf {
            vcf_het_sites.append(&mut segment_phase_sets.vcf_het_sites);
        }
    }

    if let Some(phased_het_vcf_prefix) = phased_het_vcf_prefix {
        let mut phased_het_vcf = phased_het_vcf_prefix.to_string();
        phased_het_vcf.push_str(".vcf.gz");
        write_phased_het_vcf(
            ref_filename,
            ref_chrom_list,
            thread_count,
            &phased_het_vcf,
            vcf_het_sites,
        );
    }
}
