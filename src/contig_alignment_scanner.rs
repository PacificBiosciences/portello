use std::collections::{BTreeMap, HashSet};
use std::sync::{Arc, Mutex};

use camino::Utf8Path;
use log::info;
use rust_htslib::bam::record::CigarString;
use rust_htslib::bam::{self, Read};
use rust_vc_utils::cigar::get_read_clip_positions;
use rust_vc_utils::{
    ChromList, GenomeSegment, IntRange, ProgressReporter, ReadToRefTreeMap,
    SeqOrderSplitReadSegment, get_read_segment_to_ref_pos_tree_map,
    get_region_segments_with_offset, get_seq_order_read_split_segments, rev_comp_in_place,
};
use unwrap::unwrap;

use crate::worker_thread_data::{BamReaderWorkerThreadDataSet, get_bam_reader_worker_thread_data};

pub struct ContigMappingSegmentInfo {
    pub seq_order_segment: SeqOrderSplitReadSegment,
    pub contig_to_ref_map: ReadToRefTreeMap,
}

/// All contig info that will ultimately be shared, including all split alignments, and the reverse contig sequence for
/// any contigs with reverse mappings to the reference
///
#[derive(Default)]
pub struct ContigMappingInfo {
    pub ordered_contig_segment_info: Vec<ContigMappingSegmentInfo>,
    pub rev_contig_seq: Option<Vec<u8>>,
}

#[derive(Eq, Ord, PartialEq, PartialOrd)]
struct SplitReadZipCode {
    pub chrom_index: usize,
    pub pos: i64,
    pub is_fwd_strand: bool,
    pub leading_cigar_clip: u32,
    pub trailing_cigar_clip: u32,
}

/// All data pertaining to a contig accumulated during the bam scan
///
#[derive(Default)]
struct ContigCell {
    pub contig_mapping_info: Option<ContigMappingInfo>,

    /// Accurate cigar strings taken directly from the split read alignments
    pub supp_cigars: BTreeMap<SplitReadZipCode, (CigarString, ReadToRefTreeMap)>,
}

type ParallelContigCell = Arc<Mutex<ContigCell>>;

/// Public format distributed for read re-mapping
///
/// The outer vector is ordered by assembly contig index
///
pub type AllContigMappingInfo = Vec<ContigMappingInfo>;

#[allow(clippy::too_many_arguments)]
fn scan_chromosome_segment(
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

        // Requiring the contig to start in the worker region prevents double entries (when running in WGS mode).
        //
        let contig_starts_in_worker_region = worker_region.intersect_pos(record.pos());
        if !contig_starts_in_worker_region {
            continue;
        }

        let contig_id =
            assembly_contig_list.label_to_index[std::str::from_utf8(record.qname()).unwrap()];

        if !record.is_supplementary() {
            let ordered_contig_segments =
                get_seq_order_read_split_segments(ref_chrom_list, &record);

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
                    }
                })
                .collect::<Vec<_>>();

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

            let contig_mapping_info = ContigMappingInfo {
                ordered_contig_segment_info,
                rev_contig_seq,
            };

            parallel_contig_data[contig_id]
                .lock()
                .unwrap()
                .contig_mapping_info = Some(contig_mapping_info);
        } else {
            let zip = {
                let (leading_cigar_clip, trailing_cigar_clip) = {
                    let ignore_hard_clip = false;
                    let cigar_vec = &record.cigar().0;
                    let (read_start, read_end, read_size) =
                        get_read_clip_positions(cigar_vec, ignore_hard_clip);
                    (read_start as u32, (read_size - read_end) as u32)
                };
                SplitReadZipCode {
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
                .insert(zip, (record.cigar().take(), read_to_ref_tree_map));

            if ret.is_some() {
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
                let (a, _) = ret.unwrap();
                eprintln!("old insert cigar string: {a}");
                panic!("Can't uniquely identify split read alignment info in contig '{qname}'");
            }
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
    target_region: Option<&GenomeSegment>,
    progress_reporter: &ProgressReporter,
) {
    let (target_region_offset, target_region_size) = match target_region {
        Some(x) => (x.range.start as u64, x.range.size() as u64),
        None => (0u64, chrom_size),
    };

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

struct ScanSettings {
    /// This defines how large of a chromosome segment should be processed by a single thread
    segment_size: u64,
}

pub fn scan_contig_bam(
    bam_filename: &Utf8Path,
    thread_count: usize,
    ref_chrom_list: &ChromList,
    assembly_contig_list: &ChromList,
    target_region: Option<&GenomeSegment>,
) -> AllContigMappingInfo {
    let scan_settings = &ScanSettings {
        segment_size: 20_000_000,
    };

    assert!(thread_count > 0);

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

    // Create array of contig results to be shared across all threads, the array is fixed size and each element has it's own lock:
    let contig_count = assembly_contig_list.data.len();
    let parallel_contig_data =
        std::iter::repeat_with(|| Arc::new(Mutex::new(ContigCell::default())))
            .take(contig_count)
            .collect::<Vec<_>>();
    let parallel_contig_data_ref = &parallel_contig_data;

    // We make 2 passes on the contig file to work around the minimap2 SA format, wherein the split read cigar strings minimap writes out
    // are appoximations, which we must fill in form supplementary reads.
    //

    // Pass 1 - gather data from primary reads
    worker_pool.scope(move |scope| {
        for chrom_index in 0..chrom_count {
            let chrom_info = &ref_chrom_list.data[chrom_index];
            let chrom_size = chrom_info.length;

            if let Some(target_region) = target_region {
                if chrom_index != target_region.chrom_index {
                    continue;
                }
            }

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
                    target_region,
                    progress_reporter,
                );
            });
        }
    });

    let mut missing_supplementary_match_count = 0;
    let mut contigs_with_missing_supp_match_count: HashSet<String> = Default::default();

    // Convert consolidated contig data in each contig_cell into final complete split read mapping data structure
    let result = parallel_contig_data
        .into_iter().enumerate()
        .map(|(contig_index, x)| {
            let mut ulx = x.lock().unwrap();
            let mut contig_mapping_info = ulx
                .contig_mapping_info
                .take()
                .unwrap_or_default();

            let mut missing_match_failure = false;

            for contig_segment in contig_mapping_info.ordered_contig_segment_info
                .iter_mut()
                .filter(|x| !x.seq_order_segment.from_primary_bam_record)
            {
                let zip = {
                    let (leading_cigar_clip, trailing_cigar_clip) = {
                        let ignore_hard_clip = false;
                        let cigar_vec = &contig_segment.seq_order_segment.cigar.0;
                        let (read_start, read_end, read_size) = get_read_clip_positions(cigar_vec, ignore_hard_clip);
                        (read_start as u32, (read_size-read_end) as u32)
                    };
                    SplitReadZipCode {
                        chrom_index: contig_segment.seq_order_segment.chrom_index,
                        pos: contig_segment.seq_order_segment.pos,
                        is_fwd_strand: contig_segment.seq_order_segment.is_fwd_strand,
                        leading_cigar_clip,
                        trailing_cigar_clip,
                    }
                };

                match ulx.supp_cigars.get(&zip) {
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
                                missing_supplementary_match_count += 1;
                                contigs_with_missing_supp_match_count.insert(contig_name.to_string());
                            } else {
                                let chrom_name = &ref_chrom_list.data[zip.chrom_index].label;
                                eprintln!("Can't find supplementary alignment record corresponding to segment reported in SA tag for contig '{contig_name}'");
                                eprintln!("Missing segment maps to tid: {} chrom: {chrom_name} pos: {} fwd_strand?: {}", zip.chrom_index, zip.pos, zip.is_fwd_strand);
                                missing_match_failure = true;
                            }
                        }
                    }
                }
            }

            if missing_match_failure {
                eprintln!("All split alignments for contig:");
                for (contig_segment_index, seq_order_segment) in contig_mapping_info.ordered_contig_segment_info.iter().map(|x| &x.seq_order_segment).enumerate() {
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
                panic!("Can't find supplementary alignment record corresponding to segment reported in SA tag");
            }
            contig_mapping_info
        })
        .collect();

    if missing_supplementary_match_count > 0 {
        log::warn!(
            "Couldn't match {} supplementary alignments from {} different contigs back to their primary record. These contig alignment segments will be lost.",
            missing_supplementary_match_count,
            contigs_with_missing_supp_match_count.len(),
        );
    }

    result
}
