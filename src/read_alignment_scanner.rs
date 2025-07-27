use std::sync::{Arc, Mutex};

use camino::Utf8Path;
use log::info;
use rust_htslib::bam::record::{Aux, CigarString, Record};
use rust_htslib::bam::{self, CompressionLevel, Read};
use rust_htslib::htslib;
use rust_vc_utils::aux::remove_aux_if_found;
use rust_vc_utils::cigar::get_cigar_ref_offset;
use rust_vc_utils::{
    ChromList, IntRange, ProgressReporter, SeqOrderSplitReadSegment, bam_reg2bin,
    get_alignment_end, get_region_segments, get_seq_order_read_split_segments, rev_comp_in_place,
};
use unwrap::unwrap;

use crate::cli;
use crate::contig_alignment_scanner::{AllContigMappingInfo, ContigMappingSegmentInfo};
use crate::globals::{PROGRAM_NAME, PROGRAM_VERSION};
use crate::left_shift_alignment::left_shift_alignment;
use crate::liftover_read_alignment::liftover_read_alignment;
use crate::simplify_alignment_indels::simplify_alignment_indels;
use crate::worker_thread_data::{BamReaderWorkerThreadDataSet, get_bam_reader_worker_thread_data};

type SharedBamWriter = Arc<Mutex<bam::Writer>>;

const NM_AUX_TAG: &[u8] = b"NM";
const SA_AUX_TAG: &[u8] = b"SA";
const PS_AUX_TAG: &[u8] = b"PS";
const ZM_AUX_TAG: &[u8] = b"ZM";

/// Create bam header for the portello alignment output files
///
/// This is a simple header containing just the contig info, and the portello commandline as a "PG" entry
///
fn get_alignment_file_header(chrom_list: &ChromList) -> bam::header::Header {
    let mut new_header = bam::header::Header::new();

    let mut hd_record = bam::header::HeaderRecord::new(b"HD");
    hd_record.push_tag(b"VN", "1.6");
    hd_record.push_tag(b"SO", "unsorted");
    new_header.push_record(&hd_record);

    for chrom_info in chrom_list.data.iter() {
        let mut sq_record = bam::header::HeaderRecord::new(b"SQ");
        sq_record.push_tag(b"SN", &chrom_info.label);
        sq_record.push_tag(b"LN", chrom_info.length);
        new_header.push_record(&sq_record);
    }

    let cmdline = std::env::args().collect::<Vec<_>>().join(" ");
    let mut pg_record = bam::header::HeaderRecord::new(b"PG");
    pg_record.push_tag(b"PN", PROGRAM_NAME);
    pg_record.push_tag(b"ID", format!("{PROGRAM_NAME}-{PROGRAM_VERSION}"));
    pg_record.push_tag(b"VN", PROGRAM_VERSION);
    pg_record.push_tag(b"CL", &cmdline);

    new_header.push_record(&pg_record);
    new_header
}

fn get_shared_bam_writer(
    chrom_list: &ChromList,
    filename: &Utf8Path,
    thread_count: usize,
) -> SharedBamWriter {
    let output_bam_header = get_alignment_file_header(chrom_list);
    let writer = if filename.as_str() == "-" {
        let mut x = bam::Writer::from_stdout(&output_bam_header, bam::Format::Bam).unwrap();
        x.set_compression_level(CompressionLevel::Uncompressed)
            .unwrap();
        x
    } else {
        let mut x = bam::Writer::from_path(filename, &output_bam_header, bam::Format::Bam).unwrap();
        x.set_threads(thread_count).unwrap();
        x
    };
    Arc::new(Mutex::new(writer))
}

fn get_contig_split_segments_from_read_mapping(
    read_segment: &SeqOrderSplitReadSegment,
    contig_segments: &[ContigMappingSegmentInfo],
) -> Vec<usize> {
    let mut ret = Vec::new();

    // Start and end positions of the read segment's alignment, in contig coordinates:
    let read_segment_contig_map_range = IntRange::from_pair(
        read_segment.pos,
        read_segment.pos + get_cigar_ref_offset(&read_segment.cigar),
    );
    for (contig_segment_index, contig_segment) in contig_segments.iter().enumerate() {
        let seq_order_segment = &contig_segment.seq_order_segment;
        let segment_range = IntRange::from_pair(
            seq_order_segment.seq_order_read_start as i64,
            seq_order_segment.seq_order_read_end as i64,
        );
        if segment_range.intersect_range(&read_segment_contig_map_range) {
            ret.push(contig_segment_index);
        }
    }

    ret
}

fn clone_record(record: &bam::Record) -> bam::Record {
    let mut remapped_record = record.clone();

    // remove aux tags that will be invalidated
    remove_aux_if_found(&mut remapped_record, NM_AUX_TAG);

    // remove all aux tags that will be added/modified:
    remove_aux_if_found(&mut remapped_record, SA_AUX_TAG);
    remove_aux_if_found(&mut remapped_record, PS_AUX_TAG);
    remove_aux_if_found(&mut remapped_record, ZM_AUX_TAG);

    remapped_record
}

/// Reverse a bam record, except the cigar string:
///
/// 1. Flip the reverse bit
/// 2. Revcomp seq
/// 3. Reverse qual
///
fn reverse_alignment_seq_and_qual(record: &mut bam::Record) {
    record.set_flags(record.flags() ^ htslib::BAM_FREVERSE as u16);
    let mut rev_seq = record.seq().as_bytes();
    rev_comp_in_place(&mut rev_seq);
    let rev_qual = record.qual().iter().rev().cloned().collect::<Vec<_>>();
    let qname = record.qname().to_vec();
    let cigar = record.cigar();
    record.set(&qname, Some(&cigar), &rev_seq, &rev_qual);
}

#[allow(clippy::too_many_arguments)]
fn get_liftover_alignment_for_read_and_contig_segment(
    reference: &[Vec<u8>],
    contig_list: &ChromList,
    record: &Record,
    read_segment: &SeqOrderSplitReadSegment,
    contig_segment_index: usize,
    contig_to_ref_mapping_segment_info: &ContigMappingSegmentInfo,
    rev_contig_seq: Option<&[u8]>,
    debug: bool,
) -> Option<Record> {
    // Get contig to ref map for this contig segment:
    let contig_to_ref_map = &contig_to_ref_mapping_segment_info.contig_to_ref_map;

    let contig_is_fwd_strand = contig_to_ref_mapping_segment_info
        .seq_order_segment
        .is_fwd_strand;

    let need_flipped_read_alignment = {
        let read_segment_changes_strand_from_primary =
            record.is_reverse() == read_segment.is_fwd_strand;
        (!contig_is_fwd_strand) ^ read_segment_changes_strand_from_primary
    };

    let (read_to_contig_pos_on_ref_strand, read_to_contig_cigar_on_ref_strand) =
        if contig_is_fwd_strand {
            (read_segment.pos, read_segment.cigar.to_vec())
        } else {
            let contig_index = read_segment.chrom_index;
            let contig_length = contig_list.data[contig_index].length as i64;
            let read_segment_end = read_segment.pos + get_cigar_ref_offset(&read_segment.cigar);
            let rev_pos = contig_length - read_segment_end;
            let rev_cigar = read_segment.cigar.iter().rev().cloned().collect::<Vec<_>>();

            // In addition to reversing the cigar elements, we need correct indel left-shift as well
            let mut read_seq = record.seq().as_bytes();
            if need_flipped_read_alignment {
                rev_comp_in_place(&mut read_seq);
            };
            let rev_contig_seq = rev_contig_seq.unwrap();
            left_shift_alignment(rev_pos, &rev_cigar, rev_contig_seq, &read_seq)
        };

    // recompute ref_pos and cigar string for this read
    let liftover_alignment = liftover_read_alignment(
        contig_to_ref_map,
        read_to_contig_pos_on_ref_strand,
        &read_to_contig_cigar_on_ref_strand,
    );

    if let Some((ref2_pos_orig, ref2_cigar_orig)) = liftover_alignment {
        if debug {
            eprintln!("Liftover result:");
            eprintln!("  Input read_segment pos: {}", read_segment.pos);
            eprintln!("Input read_segment cigar: {}", read_segment.cigar);
            eprintln!("  Remapped pos: {ref2_pos_orig}");
            eprintln!("Remapped cigar: {}", CigarString(ref2_cigar_orig.clone()));
            let contig_index = read_segment.chrom_index;
            let contig_name = &contig_list.data[contig_index].label;
            eprintln!(
                "target read_segment mapped to contig: {contig_name} at 1-index pos: {}",
                read_segment.pos + 1
            );
            eprintln!("contig mapping to ref:");
            for (k, v) in contig_to_ref_map.get_map().iter() {
                eprintln!("({k},{v:?})");
            }
        }

        // TEMP verify the length of the remapped cigar string:
        use rust_vc_utils::cigar::get_cigar_read_offset;
        let cigar_read_len = get_cigar_read_offset(&ref2_cigar_orig, false);
        if record.seq_len() != cigar_read_len {
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();

            eprintln!("Failed to remap qname: {qname} record: {record:?}");
            eprintln!(
                "Seq len: {} new cigar len: {}",
                record.seq_len(),
                cigar_read_len
            );
            eprintln!(".  Input read_segment cigar: {}", read_segment.cigar);
            eprintln!("Remapped cigar: {}", CigarString(ref2_cigar_orig.clone()));
            let contig_index = read_segment.chrom_index;
            let contig_name = &contig_list.data[contig_index].label;
            eprintln!(
                "target read_segment mapped to contig: {contig_name} at 1-index pos: {}",
                read_segment.pos + 1
            );
            eprintln!("contig mapping to ref:");
            for (k, v) in contig_to_ref_map.get_map().iter() {
                eprintln!("({k},{v:?})");
            }
            panic!("Aborting...");
        }

        // Add additional refinement step canonicalize indel clusters in the cigar string alignments:
        let chrom_index = contig_to_ref_mapping_segment_info
            .seq_order_segment
            .chrom_index;

        let (ref2_pos, ref2_cigar) = {
            let chrom_ref = &reference[chrom_index];
            let mut read_seq = record.seq().as_bytes();
            if need_flipped_read_alignment {
                rev_comp_in_place(&mut read_seq);
            };
            simplify_alignment_indels(ref2_pos_orig, &ref2_cigar_orig, chrom_ref, &read_seq)
        };

        let mut remapped_record = clone_record(record);

        remapped_record.set_tid(chrom_index as i32);

        // Adopt MAPQ value from contig:
        let contig_mapq = contig_to_ref_mapping_segment_info.seq_order_segment.mapq;
        let original_read_mapq = remapped_record.mapq();
        remapped_record.set_mapq(contig_mapq);

        // Set PS tag based on the contig id and split segment number
        let contig_index = read_segment.chrom_index;
        let contig_name = &contig_list.data[contig_index].label;
        let ps_tag = format!(
            "{}_split{}{}",
            contig_name,
            contig_segment_index,
            if contig_is_fwd_strand { "+" } else { "-" }
        );
        remapped_record
            .push_aux(PS_AUX_TAG, Aux::String(&ps_tag))
            .unwrap();

        remapped_record
            .push_aux(ZM_AUX_TAG, Aux::U8(original_read_mapq))
            .unwrap();

        remapped_record.set_pos(ref2_pos);
        remapped_record.set_cigar(Some(&CigarString(ref2_cigar)));

        if need_flipped_read_alignment {
            reverse_alignment_seq_and_qual(&mut remapped_record)
        };

        let ref2_end = get_alignment_end(&remapped_record) as usize;
        remapped_record.set_bin(bam_reg2bin(remapped_record.pos() as usize, ref2_end));

        // Set all reads to supplementary until we determine the primary alignment
        remapped_record.set_supplementary();

        Some(remapped_record)
    } else {
        None
    }
}

/// Produce one split read segment string formatted for the SA aux tag
///
fn get_sa_tag_segment(chrom_list: &ChromList, record: &Record) -> String {
    let chrom = &chrom_list.data[record.tid() as usize].label;
    let schar = if record.is_reverse() { '-' } else { '+' };
    format!(
        "{chrom},{},{schar},{},{},0;",
        record.pos() + 1,
        record.cigar(),
        record.mapq(),
    )
}

/// Resolve full set of split alignments lifted-over from one primary read-to-assembly alignment:
///
/// - If no reads lifted over to the reference, then generate an unmapped read in the output bam.
/// - Otherwise:
///   - Determine which split read will be primary and updates tags
///   - Generate SA tags for all reads
///
fn finish_remapped_alignment_set(
    ref_chrom_list: &ChromList,
    orig_primary_record: &bam::Record,
    mut remapped_records: Vec<bam::Record>,
    is_target_region: bool,
) -> Vec<bam::Record> {
    let remapped_record_count = remapped_records.len();
    if remapped_record_count == 0 {
        if is_target_region {
            vec![]
        } else {
            // Copy primary read-to-contig alignment into an unmapped alignment:
            let mut unmapped_record = clone_record(orig_primary_record);
            unmapped_record.set_unmapped();
            unmapped_record.unset_supplementary();
            unmapped_record.set_cigar(None);
            unmapped_record.set_mapq(255);
            unmapped_record.set_tid(-1);
            unmapped_record.set_pos(-1);

            if unmapped_record.is_reverse() {
                reverse_alignment_seq_and_qual(&mut unmapped_record);
            };

            vec![unmapped_record]
        }
    } else {
        // Determine which alignment is primary and flag it
        let mut primary_record_index = 0;
        for remapped_record_index in 1..remapped_record_count {
            if remapped_records[primary_record_index].mapq()
                < remapped_records[remapped_record_index].mapq()
            {
                primary_record_index = remapped_record_index;
            }
        }
        remapped_records[primary_record_index].unset_supplementary();

        // Update SA tags on entire read set:
        for remapped_record_index in 0..remapped_record_count {
            let mut aux_str = String::new();
            for remapped_record_index2 in
                (0..remapped_record_count).filter(|&x| x != remapped_record_index)
            {
                aux_str +=
                    get_sa_tag_segment(ref_chrom_list, &remapped_records[remapped_record_index2])
                        .as_str();
            }
            if !aux_str.is_empty() {
                remapped_records[remapped_record_index]
                    .push_aux(SA_AUX_TAG, Aux::String(&aux_str))
                    .unwrap();
            }
        }
        remapped_records
    }
}

#[allow(clippy::too_many_arguments)]
fn scan_chromosome_segment(
    bam_reader: &mut bam::IndexedReader,
    reference: &[Vec<u8>],
    ref_chrom_list: &ChromList,
    contig_list: &ChromList,
    contig_index: usize,
    begin: i64,
    end: i64,
    all_contig_mapping_info: &AllContigMappingInfo,
    is_target_region: bool,
    remapped_bam_writer: &mut SharedBamWriter,
    progress_reporter: &ProgressReporter,
) {
    bam_reader
        .fetch(bam::FetchDefinition::Region(
            contig_index as i32,
            begin,
            end,
        ))
        .unwrap();

    let worker_region = IntRange::from_pair(begin, end);

    let mut record = bam::Record::new();
    while let Some(r) = bam_reader.read(&mut record) {
        unwrap!(r, "Failed to parse alignment record");

        assert!(!record.is_unmapped());

        // Determine if this is the primary or supplemental alignment, we want to parse all the breakpoints from a read
        // once, so working exclusively with the primary alignment should allow this.
        //
        // Requiring the read to start in the worker region prevents duplicated read output.
        //
        let read_starts_in_worker_region = worker_region.intersect_pos(record.pos());
        if record.is_supplementary() || (!read_starts_in_worker_region) {
            continue;
        }

        let mut debug = false;
        if debug {
            let target_qname = "m84011_220902_175841_s1/197591380/ccs";
            let qname = std::str::from_utf8(record.qname()).unwrap().to_string();
            debug = qname == target_qname;

            if debug {
                eprintln!("Target qname found: {qname}");
            }
        }

        let mut remapped_records = Vec::new();

        let ordered_splits = get_seq_order_read_split_segments(contig_list, &record);

        if debug {
            eprintln!("target qname read-to-contig alignment segments:");
            for (index, read_segment) in ordered_splits.iter().enumerate() {
                eprintln!("Index: {index} seg: {read_segment:?}");
            }
        }

        for (read_segment_index, read_segment) in ordered_splits.iter().enumerate() {
            if debug {
                eprintln!("Starting to process read_segment_index {read_segment_index}");
            }

            let contig_mapping_info = &all_contig_mapping_info[read_segment.chrom_index];
            let contig_segments = &contig_mapping_info.ordered_contig_segment_info;

            // Find the range of contig segments this read segment is mapped to:
            let contig_segment_indexes =
                get_contig_split_segments_from_read_mapping(read_segment, contig_segments);

            if debug {
                eprintln!(
                    "read segment associated with {} contig_segments.",
                    contig_segment_indexes.len()
                );
                for x in contig_segment_indexes.iter() {
                    let contig_to_ref_mapping_segment_info = &contig_segments[*x];
                    eprintln!(
                        "contig_segment_index {} details: {:?}",
                        *x, contig_to_ref_mapping_segment_info.seq_order_segment
                    );
                }
            }

            for contig_segment_index in contig_segment_indexes {
                let contig_to_ref_mapping_segment_info = &contig_segments[contig_segment_index];
                let remapped_record = get_liftover_alignment_for_read_and_contig_segment(
                    reference,
                    contig_list,
                    &record,
                    read_segment,
                    contig_segment_index,
                    contig_to_ref_mapping_segment_info,
                    contig_mapping_info.rev_contig_seq.as_deref(),
                    debug,
                );
                if let Some(x) = remapped_record {
                    remapped_records.push(x);
                }
            }
        }

        let remapped_records = finish_remapped_alignment_set(
            ref_chrom_list,
            &record,
            remapped_records,
            is_target_region,
        );

        // Output full primary+supplementary set of remapped alignments
        {
            let mut unlocked_writer = remapped_bam_writer.lock().unwrap();
            for remapped_record in remapped_records {
                unlocked_writer.write(&remapped_record).unwrap();
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
    reference: &[Vec<u8>],
    ref_chrom_list: &ChromList,
    contig_list: &ChromList,
    contig_index: usize,
    contig_size: u64,
    all_contig_mapping_info: &AllContigMappingInfo,
    is_target_region: bool,
    remapped_bam_writer: &mut SharedBamWriter,
    progress_reporter: &ProgressReporter,
) {
    let chrom_segments = get_region_segments(contig_size, scan_settings.segment_size);

    rayon::scope(move |scope| {
        for (begin, end) in chrom_segments {
            let worker_thread_dataset = worker_thread_dataset.clone();
            let mut remapped_bam_writer = remapped_bam_writer.clone();

            scope.spawn(move |_| {
                let worker_id = rayon::current_thread_index().unwrap();
                let bam_reader = &mut worker_thread_dataset[worker_id].lock().unwrap().bam_reader;

                scan_chromosome_segment(
                    bam_reader,
                    reference,
                    ref_chrom_list,
                    contig_list,
                    contig_index,
                    begin as i64,
                    end as i64,
                    all_contig_mapping_info,
                    is_target_region,
                    &mut remapped_bam_writer,
                    progress_reporter,
                );
            });
        }
    });
}

fn scan_unmapped_reads(
    bam_reader: &mut bam::IndexedReader,
    unmapped_bam_writer: &mut SharedBamWriter,
) {
    // bam de-facto convention is to put the unmapped reads at either the front or the back of the mapped&sorted
    // content, so ideally we would check both regions. Since htslib how has an unmapped fetch lets see if that works
    // instead.

    bam_reader.fetch(bam::FetchDefinition::Unmapped).unwrap();

    //bam_reader.seek(0).unwrap();

    let mut record = bam::Record::new();
    while let Some(r) = bam_reader.read(&mut record) {
        unwrap!(r, "Failed to parse alignment record");

        if !record.is_unmapped() {
            continue;
        }

        unmapped_bam_writer.lock().unwrap().write(&record).unwrap();
    }
}

struct ScanSettings {
    /// This defines how large of a chromosome segment should be processed by a single thread
    segment_size: u64,
}

pub fn scan_and_remap_reads(
    settings: &cli::Settings,
    thread_count: usize,
    reference: &[Vec<u8>],
    ref_chrom_list: &ChromList,
    all_contig_mapping_info: &AllContigMappingInfo,
    is_target_region: bool,
) {
    let scan_settings = &ScanSettings {
        segment_size: 20_000_000,
    };

    assert!(thread_count > 0);

    let bam_filename = &settings.read_to_assembly_bam;

    info!("Processing read-to-contig alignment file '{bam_filename}'");

    // Immediately take reference-value for use in the worker_pool below
    let assembly_contig_list = &ChromList::from_bam_filename(bam_filename.as_str());

    // Open up bam output structures
    //
    let writer_thread_count = std::cmp::max(1, thread_count / 2);

    let remapped_bam_writer = get_shared_bam_writer(
        ref_chrom_list,
        &settings.remapped_read_output,
        writer_thread_count,
    );

    let mut unmapped_bam_writer = get_shared_bam_writer(
        ref_chrom_list,
        &settings.unassembled_read_output,
        writer_thread_count,
    );

    // Setup shared worker thread data structures:
    let worker_thread_dataset = get_bam_reader_worker_thread_data(thread_count, bam_filename);

    let worker_pool = rayon::ThreadPoolBuilder::new()
        .num_threads(thread_count)
        .build()
        .unwrap();

    let contig_count = assembly_contig_list.data.len();

    let genome_kb = assembly_contig_list
        .data
        .iter()
        .map(|x| x.length)
        .sum::<u64>()
        / 1000;
    let progress_reporter = &ProgressReporter::new(
        genome_kb,
        "Remapped read alignments from",
        "assembly contig kb",
        false,
    );

    worker_pool.scope(move |scope| {
        {
            let worker_thread_dataset = worker_thread_dataset.clone();

            scope.spawn(move |_| {
                let worker_id = rayon::current_thread_index().unwrap();
                let bam_reader = &mut worker_thread_dataset[worker_id].lock().unwrap().bam_reader;

                scan_unmapped_reads(bam_reader, &mut unmapped_bam_writer);
            });
        }

        for contig_index in 0..contig_count {
            let chrom_info = &assembly_contig_list.data[contig_index];
            let contig_size = chrom_info.length;

            let worker_thread_dataset = worker_thread_dataset.clone();
            let mut remapped_bam_writer = remapped_bam_writer.clone();
            scope.spawn(move |_| {
                scan_chromosome_segments(
                    worker_thread_dataset,
                    scan_settings,
                    reference,
                    ref_chrom_list,
                    assembly_contig_list,
                    contig_index,
                    contig_size,
                    all_contig_mapping_info,
                    is_target_region,
                    &mut remapped_bam_writer,
                    progress_reporter,
                );
            });
        }
    });
}
