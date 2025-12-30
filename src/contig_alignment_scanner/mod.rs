mod contig_colinear_segment_joiner;
mod contig_repeated_match_trimmer;
mod contig_segment_hap_labeler;
mod contig_segment_phase_set_labeler;
mod non_targeted_segment_filter;

use std::collections::{BTreeMap, HashSet};
use std::sync::{Arc, Mutex};

use log::info;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::{self, Read};
use rust_vc_utils::cigar::get_read_clip_positions;
use rust_vc_utils::{
    ChromList, GenomeSegment, IntRange, ProgressReporter, ReadToRefTreeMap, RefToReadTreeMap,
    RegionMap, SeqOrderSplitReadSegment, get_read_segment_to_ref_pos_tree_map,
    get_region_segments_with_offset, get_seq_order_read_split_segments, rev_comp_in_place,
};
use unwrap::unwrap;

use self::contig_colinear_segment_joiner::join_colinear_contig_segments;
use self::contig_repeated_match_trimmer::clip_repeated_contig_matches;
use self::contig_segment_hap_labeler::add_contig_segment_haplotype_labels;
use self::contig_segment_phase_set_labeler::add_contig_segment_phase_set_labels;
pub use self::contig_segment_phase_set_labeler::{ChromPhaseSetMap, ContigSegmentInfo};
use self::non_targeted_segment_filter::filter_non_targeted_segments;
use crate::cli::{self, InputAssemblyMode};
use crate::worker_thread_data::{BamReaderWorkerThreadDataSet, get_bam_reader_worker_thread_data};

/// All data for a set of co-linear contig split segments
///
pub struct ContigMappingSegmentInfo {
    pub seq_order_segment: SeqOrderSplitReadSegment,

    /// Data structure derived from the contig split segment location and cigar string, used to quickly remap reads onto
    /// the contig split segment mapping to reference:
    ///
    pub contig_to_ref_map: ReadToRefTreeMap,

    /// Data structure derived from contig_to_ref_map, used to quickly map reference variant locations back to a contigs
    /// supporting reads to determine phasing support
    ///
    pub ref_to_contig_map: RefToReadTreeMap,

    /// The haplotype id for this contig segment, this will be used to set the bam aux HP tag if set
    pub hap_id: Option<u8>,

    /// Phase sets overlapping this contig segment, expressed as start-end ranges in coordinates of the reference
    /// chromosome this segment is mapped to.
    ///
    /// Multiple phase sets are allowed here because a dual-assembly contig that's partially-phased may need to be split
    /// into chunks.
    ///
    pub phase_sets: RegionMap,
}

/// All contig info that will ultimately be shared downstream outside of this scanner module. This includes all split
/// alignments, and the reverse contig sequence for any contigs with reverse mappings to the reference.
///
#[derive(Default)]
pub struct ContigMappingInfo {
    /// Contig name, only used to improve error messages
    pub qname: String,

    /// Primary contig mapping information structure
    pub ordered_contig_segment_info: Vec<ContigMappingSegmentInfo>,

    /// Contig sequence needs to be stored forgenerating read-backed phasing annotations given partially-phased assembly
    /// contig input
    pub contig_seq: Option<Vec<u8>>,

    /// Revcomp of the contig sequence, used to left-shift indel alignments on any reverse-mapped contigs
    pub rev_contig_seq: Option<Vec<u8>>,
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
struct SplitReadKey {
    pub chrom_index: usize,
    pub pos: i64,
    pub is_fwd_strand: bool,
    pub leading_cigar_clip: u32,
    pub trailing_cigar_clip: u32,
}

type SuppCigarMap = BTreeMap<SplitReadKey, (CigarString, ReadToRefTreeMap)>;

/// All data pertaining to a contig accumulated during the bam scan
///
#[derive(Default)]
struct ContigCell {
    pub contig_mapping_info: Option<ContigMappingInfo>,

    /// Accurate cigar strings taken directly from the split read alignments
    ///
    /// These are used to fill-in and replace the approximate cigar strings provided by minimap2 for split reads in the SA tag
    ///
    pub supp_cigars: SuppCigarMap,
}

type ParallelContigCell = Arc<Mutex<ContigCell>>;

/// Public format distributed outside of this module for read re-mapping
///
/// The outer vector is ordered by assembly contig index
///
pub type AllContigMappingInfo = Vec<ContigMappingInfo>;

/// Debug summary of all exported mapping info for a contig
pub fn print_split_read_summary(contig_mapping_info: &ContigMappingInfo) {
    eprintln!("contig '{}' split read summary:", contig_mapping_info.qname);
    for (i, x) in contig_mapping_info
        .ordered_contig_segment_info
        .iter()
        .enumerate()
    {
        eprintln!("split {i}: {}", x.seq_order_segment.short_display());
    }
}

/// Clear all existing phase sets from contig mapping info
///
pub fn clear_segment_phase_sets(all_contig_mapping_info: &mut AllContigMappingInfo) {
    for chrom_contig_mapping_info in all_contig_mapping_info.iter_mut() {
        for contig_mapping_segment_info in chrom_contig_mapping_info
            .ordered_contig_segment_info
            .iter_mut()
        {
            contig_mapping_segment_info.phase_sets.clear();
        }
    }
}

/// Use the primary read for each contig to create the core contig_mapping_info structure
///
fn add_primary_read(
    are_input_contigs_partially_phased: bool,
    ref_chrom_list: &ChromList,
    record: &bam::Record,
) -> ContigMappingInfo {
    let ordered_contig_segments = get_seq_order_read_split_segments(ref_chrom_list, record);

    let ordered_contig_segment_info = ordered_contig_segments
        .into_iter()
        .map(|seq_order_segment| {
            let contig_to_ref_map = if seq_order_segment.from_primary_bam_record {
                get_read_segment_to_ref_pos_tree_map(
                    seq_order_segment.pos,
                    &seq_order_segment.cigar,
                    false,
                )
            } else {
                Default::default()
            };

            ContigMappingSegmentInfo {
                seq_order_segment,
                contig_to_ref_map,
                ref_to_contig_map: Default::default(),
                hap_id: None,
                phase_sets: Default::default(),
            }
        })
        .collect::<Vec<_>>();

    let need_contig_seq = are_input_contigs_partially_phased;

    let contig_seq = if need_contig_seq {
        let mut seq = record.seq().as_bytes();
        if record.is_reverse() {
            rev_comp_in_place(&mut seq);
        }
        Some(seq)
    } else {
        None
    };

    let need_rev_contig_seq = ordered_contig_segment_info
        .iter()
        .any(|x| !x.seq_order_segment.is_fwd_strand);

    let rev_contig_seq = if need_rev_contig_seq {
        let mut rev_seq = record.seq().as_bytes();
        if !record.is_reverse() {
            rev_comp_in_place(&mut rev_seq);
        }
        Some(rev_seq)
    } else {
        None
    };

    let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
    ContigMappingInfo {
        qname,
        ordered_contig_segment_info,
        contig_seq,
        rev_contig_seq,
    }
}

fn add_split_read_cigar_to_supp_cigar_set(
    ref_chrom_list: &ChromList,
    parallel_contig_data: &[ParallelContigCell],
    record: &bam::Record,
    contig_id: usize,
) {
    let split_read_key = {
        let (leading_cigar_clip, trailing_cigar_clip) = {
            let ignore_hard_clip = false;
            let cigar_vec = &record.cigar().0;
            let (read_start, read_end, read_size) =
                get_read_clip_positions(cigar_vec, ignore_hard_clip);
            (read_start as u32, (read_size - read_end) as u32)
        };
        SplitReadKey {
            chrom_index: record.tid() as usize,
            pos: record.pos(),
            is_fwd_strand: !record.is_reverse(),
            leading_cigar_clip,
            trailing_cigar_clip,
        }
    };
    let read_to_ref_tree_map =
        get_read_segment_to_ref_pos_tree_map(record.pos(), &record.cigar(), false);
    let ret = parallel_contig_data[contig_id]
        .lock()
        .unwrap()
        .supp_cigars
        .insert(
            split_read_key,
            (record.cigar().take(), read_to_ref_tree_map),
        );

    if let Some((old_cigar, _)) = ret {
        let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
        eprintln!("Can't uniquely identify split read alignment info in contig '{qname}'");
        let ref_chrom_name = &ref_chrom_list.data[record.tid() as usize].label;
        eprintln!(
            "trying to insert split read at tid: {} chrom: {} pos: {} fwd_strand: {} cigar: {}",
            record.tid(),
            ref_chrom_name,
            record.pos(),
            !record.is_reverse(),
            record.cigar()
        );
        eprintln!("old insert cigar string: {old_cigar}");
        panic!("Can't uniquely identify split read alignment info in contig '{qname}'");
    }
}

#[allow(clippy::too_many_arguments)]
fn scan_chromosome_segment(
    scan_settings: &ScanSettings,
    bam_reader: &mut bam::IndexedReader,
    ref_chrom_list: &ChromList,
    assembly_contig_list: &ChromList,
    chrom_index: usize,
    begin: i64,
    end: i64,
    parallel_contig_data: &[ParallelContigCell],
    progress_reporter: &ProgressReporter,
) {
    bam_reader
        .fetch(bam::FetchDefinition::Region(chrom_index as i32, begin, end))
        .unwrap();

    let worker_region = IntRange::from_pair(begin, end);

    let mut record = bam::Record::new();
    while let Some(r) = bam_reader.read(&mut record) {
        unwrap!(r, "Failed to parse alignment record");

        // We don't expect to see unmapped reads here because of the above fetch call, but just in case, the expected
        // behavior is to skip them and infer their unmapped status from their absence in the contig map
        if record.is_unmapped() || record.is_secondary() {
            continue;
        }

        // Requiring the contig to start in the worker region prevents double entries
        //
        let contig_starts_in_worker_region = worker_region.intersect_pos(record.pos());
        if !contig_starts_in_worker_region {
            continue;
        }

        let contig_id =
            assembly_contig_list.label_to_index[std::str::from_utf8(record.qname()).unwrap()];

        if !record.is_supplementary() {
            let contig_mapping_info = add_primary_read(
                scan_settings.are_input_contigs_partially_phased,
                ref_chrom_list,
                &record,
            );
            parallel_contig_data[contig_id]
                .lock()
                .unwrap()
                .contig_mapping_info = Some(contig_mapping_info);
        } else {
            add_split_read_cigar_to_supp_cigar_set(
                ref_chrom_list,
                parallel_contig_data,
                &record,
                contig_id,
            );
        }
    }

    let kb_completed = (end - begin) as u64 / 1000;
    progress_reporter.inc(kb_completed);
}

#[allow(clippy::too_many_arguments)]
fn scan_chromosome_segments(
    worker_thread_dataset: BamReaderWorkerThreadDataSet,
    scan_settings: &ScanSettings,
    ref_chrom_list: &ChromList,
    assembly_contig_list: &ChromList,
    chrom_index: usize,
    chrom_size: u64,
    parallel_contig_data: &[ParallelContigCell],
    progress_reporter: &ProgressReporter,
) {
    let (target_region_offset, target_region_size) = (0u64, chrom_size);

    let target_region_segments = get_region_segments_with_offset(
        target_region_offset,
        target_region_size,
        scan_settings.segment_size,
    );

    rayon::scope(move |scope| {
        for (begin, end) in target_region_segments {
            let worker_thread_dataset = worker_thread_dataset.clone();

            scope.spawn(move |_| {
                let worker_id = rayon::current_thread_index().unwrap();

                let bam_reader = &mut worker_thread_dataset[worker_id].lock().unwrap().bam_reader;

                scan_chromosome_segment(
                    scan_settings,
                    bam_reader,
                    ref_chrom_list,
                    assembly_contig_list,
                    chrom_index,
                    begin as i64,
                    end as i64,
                    parallel_contig_data,
                    progress_reporter,
                );
            });
        }
    });
}

fn post_process_contig_cell(
    contig_index: usize,
    pcc: ParallelContigCell,
    ref_chrom_list: &ChromList,
    assembly_contig_list: &ChromList,
    target_region: Option<&GenomeSegment>,
    missing_supplementary_match_count: &mut i32,
    contigs_with_missing_supp_match_count: &mut HashSet<String>,
) -> ContigMappingInfo {
    let mut unlocked_pcc = pcc.lock().unwrap();
    let mut contig_mapping_info = unlocked_pcc.contig_mapping_info.take().unwrap_or_default();

    let mut missing_match_failure = false;

    for contig_segment in contig_mapping_info
        .ordered_contig_segment_info
        .iter_mut()
        .filter(|x| !x.seq_order_segment.from_primary_bam_record)
    {
        let split_read_key = {
            let (leading_cigar_clip, trailing_cigar_clip) = {
                let ignore_hard_clip = false;
                let cigar_vec = &contig_segment.seq_order_segment.cigar.0;
                let (read_start, read_end, read_size) =
                    get_read_clip_positions(cigar_vec, ignore_hard_clip);
                (read_start as u32, (read_size - read_end) as u32)
            };
            SplitReadKey {
                chrom_index: contig_segment.seq_order_segment.chrom_index,
                pos: contig_segment.seq_order_segment.pos,
                is_fwd_strand: contig_segment.seq_order_segment.is_fwd_strand,
                leading_cigar_clip,
                trailing_cigar_clip,
            }
        };

        match unlocked_pcc.supp_cigars.get(&split_read_key) {
            Some(x) => {
                contig_segment.seq_order_segment.cigar = x.0.clone();
                contig_segment.contig_to_ref_map = x.1.clone();
            }
            None => {
                // Getting here means that we never found the supplementary read bam record corresponding to the segment described in the SA tag
                //
                // We expect this to happen in targeted runs, but it shouldn't happen in a WGS run.
                //
                let allow_missing_supplementary_matches_in_wgs = false;
                if target_region.is_none() {
                    let contig_name = &assembly_contig_list.data[contig_index].label;
                    if allow_missing_supplementary_matches_in_wgs {
                        *missing_supplementary_match_count += 1;
                        contigs_with_missing_supp_match_count.insert(contig_name.to_string());
                    } else {
                        let chrom_name = &ref_chrom_list.data[split_read_key.chrom_index].label;
                        eprintln!(
                            "Can't find supplementary alignment record corresponding to segment reported in SA tag for contig '{contig_name}'"
                        );
                        eprintln!(
                            "Missing segment maps to tid: {} chrom: {chrom_name} pos: {} fwd_strand?: {}",
                            split_read_key.chrom_index,
                            split_read_key.pos,
                            split_read_key.is_fwd_strand
                        );
                        missing_match_failure = true;
                    }
                }
            }
        }
    }

    if missing_match_failure {
        eprintln!("All split alignments for contig:");
        for (contig_segment_index, seq_order_segment) in contig_mapping_info
            .ordered_contig_segment_info
            .iter()
            .map(|x| &x.seq_order_segment)
            .enumerate()
        {
            let segment_range = IntRange::from_pair(
                seq_order_segment.seq_order_read_start as i64,
                seq_order_segment.seq_order_read_end as i64,
            );
            eprintln!(
                "contig split index: {contig_segment_index} range: {segment_range:?} chrom_index: {} pos: {} mapq: {} primary?: {} is_fwd: {} cigar: {}",
                seq_order_segment.chrom_index,
                seq_order_segment.pos,
                seq_order_segment.mapq as i64,
                seq_order_segment.from_primary_bam_record,
                seq_order_segment.is_fwd_strand,
                seq_order_segment.cigar,
            );
        }
        panic!(
            "Can't find supplementary alignment record corresponding to segment reported in SA tag"
        );
    }

    contig_mapping_info
}

/// Convert distributed contig data in each contig_cell into complete split read mapping data structure which can be
/// used downstream in read liftover procedure
///
fn post_process_contig_data(
    parallel_contig_data: Vec<ParallelContigCell>,
    ref_chrom_list: &ChromList,
    assembly_contig_list: &ChromList,
    target_region: Option<&GenomeSegment>,
) -> Vec<ContigMappingInfo> {
    let mut missing_supplementary_match_count = 0;
    let mut contigs_with_missing_supp_match_count: HashSet<String> = Default::default();

    let processed_contig_data = parallel_contig_data
        .into_iter()
        .enumerate()
        .map(|(i, x)| {
            post_process_contig_cell(
                i,
                x,
                ref_chrom_list,
                assembly_contig_list,
                target_region,
                &mut missing_supplementary_match_count,
                &mut contigs_with_missing_supp_match_count,
            )
        })
        .collect();

    if missing_supplementary_match_count > 0 {
        log::warn!(
            "Couldn't match {} supplementary alignments from {} different contigs back to their primary record. These contig alignment segments will be lost.",
            missing_supplementary_match_count,
            contigs_with_missing_supp_match_count.len(),
        );
    }

    processed_contig_data
}

/// After all contig_to_ref updates have been made for each contig, create an inverted ref_to_contig map
///
/// The inverted map is only used for phasing
///
fn add_ref_to_contig_map(result: &mut AllContigMappingInfo) {
    for contig_mapping_info in result
        .iter_mut()
        .filter(|x| !x.ordered_contig_segment_info.is_empty())
    {
        for contig_segment in contig_mapping_info.ordered_contig_segment_info.iter_mut() {
            contig_segment.ref_to_contig_map =
                RefToReadTreeMap::from_read_to_ref_tree_map(&contig_segment.contig_to_ref_map);
        }
    }
}

struct ScanSettings {
    /// This defines how large of a chromosome segment should be processed by a single thread
    segment_size: u64,

    /// If assembly contigs are partially phased then we need to store them all for downstream phasing process
    are_input_contigs_partially_phased: bool,
}

/// Scan the alignment file containing contig to ref mappings, and extract all information needed for downstream read liftover
///
/// Returns a 2-tuple of:
/// 1. A vector with one entry for each contig, providing all contig mapping information for that contig
/// 2. XXX
///
pub fn scan_contig_bam(
    settings: &cli::Settings,
    thread_count: usize,
    ref_chrom_list: &ChromList,
    assembly_contig_list: &ChromList,
    target_region: Option<&GenomeSegment>,
) -> (AllContigMappingInfo, Vec<ChromPhaseSetMap>) {
    let are_input_contigs_partially_phased =
        settings.input_assembly_mode == InputAssemblyMode::partially_phased;

    let scan_settings = &ScanSettings {
        segment_size: 20_000_000,
        are_input_contigs_partially_phased,
    };

    assert!(thread_count > 0);

    let bam_filename = &settings.assembly_to_ref_bam;

    info!("Processing contig-to-ref alignment file '{bam_filename}'");

    // Setup shared worker thread data structures:
    let worker_thread_dataset = get_bam_reader_worker_thread_data(thread_count, bam_filename);

    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .unwrap();

    let chrom_count = ref_chrom_list.data.len();

    let genome_kb = ref_chrom_list.data.iter().map(|x| x.length).sum::<u64>() / 1000;
    let progress_reporter = ProgressReporter::new(
        genome_kb,
        "Processed assembly contig alignments on",
        "ref genome kb",
        false,
    );

    let progress_reporter = &progress_reporter;

    // Create array of contig results to be shared across all threads, the array is fixed size and each element has its
    // own lock to minimize contention. The locks are still needed for the case where different split read entries for the
    // same contig are being processed by different threads.
    let contig_count = assembly_contig_list.data.len();
    let parallel_contig_data =
        std::iter::repeat_with(|| Arc::new(Mutex::new(ContigCell::default())))
            .take(contig_count)
            .collect::<Vec<_>>();
    let parallel_contig_data_ref = &parallel_contig_data;

    worker_pool.scope(move |scope| {
        for chrom_index in 0..chrom_count {
            let chrom_info = &ref_chrom_list.data[chrom_index];
            let chrom_size = chrom_info.length;

            let worker_thread_dataset = worker_thread_dataset.clone();
            scope.spawn(move |_| {
                scan_chromosome_segments(
                    worker_thread_dataset,
                    scan_settings,
                    ref_chrom_list,
                    assembly_contig_list,
                    chrom_index,
                    chrom_size,
                    parallel_contig_data_ref,
                    progress_reporter,
                );
            });
        }
    });

    // Convert consolidated contig data in each contig_cell into final complete split read mapping data structure
    let mut processed_contig_data = post_process_contig_data(
        parallel_contig_data,
        ref_chrom_list,
        assembly_contig_list,
        target_region,
    );

    // The strategy for target region is to ignore it for the geneome scan, but then filter out any read segments that
    // are not targetted here. Although this is little less effecient, it greatly simplifies the target regions code
    // path, and givent that this is intended for debug anyway that should be good enough.
    //
    filter_non_targeted_segments(target_region, &mut processed_contig_data);

    clip_repeated_contig_matches(&mut processed_contig_data);

    join_colinear_contig_segments(&mut processed_contig_data);

    if settings.input_assembly_mode == InputAssemblyMode::partially_phased {
        add_ref_to_contig_map(&mut processed_contig_data);
    }

    add_contig_segment_haplotype_labels(chrom_count, &mut processed_contig_data);

    let annotated_phase_set_ranges =
        add_contig_segment_phase_set_labels(ref_chrom_list, &mut processed_contig_data);

    (processed_contig_data, annotated_phase_set_ranges)
}
