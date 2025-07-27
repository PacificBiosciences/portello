mod cli;
mod contig_alignment_scanner;
mod globals;
mod left_shift_alignment;
mod liftover_read_alignment;
mod logger;
mod read_alignment_scanner;
mod simplify_alignment_indels;
mod worker_thread_data;

use std::{error, process};

use hhmmss::Hhmmss;
use log::info;
use rust_vc_utils::{ChromList, GenomeSegment, get_genome_ref_from_fasta};

use crate::cli::{DerivedSettings, Settings, validate_and_fix_settings, validate_settings_data};
use crate::contig_alignment_scanner::scan_contig_bam;
use crate::globals::{PROGRAM_NAME, PROGRAM_VERSION};
use crate::logger::setup_logger;
use crate::read_alignment_scanner::scan_and_remap_reads;

fn run(
    settings: &Settings,
    derived_settings: &DerivedSettings,
) -> Result<(), Box<dyn error::Error>> {
    info!("Starting {PROGRAM_NAME} {PROGRAM_VERSION}");
    info!(
        "cmdline: {}",
        std::env::args().collect::<Vec<_>>().join(" ")
    );
    info!("Running on {} threads", derived_settings.thread_count);

    let start = std::time::Instant::now();

    let ref_chrom_list = ChromList::from_bam_filename(settings.assembly_to_ref_bam.as_str());
    let assembly_contig_list = ChromList::from_bam_filename(settings.read_to_assembly_bam.as_str());

    let target_region = settings
        .target_region
        .as_ref()
        .map(|x| GenomeSegment::from_region_str(&ref_chrom_list, x));

    // Load reference sequence but move it from chrom name hash to chrom index lookup array:
    let reference = {
        let mut r = get_genome_ref_from_fasta(settings.ref_filename.as_str())
            .chroms
            .into_iter()
            .map(|(label, val)| (*ref_chrom_list.label_to_index.get(&label).unwrap(), val))
            .collect::<Vec<_>>();
        r.sort_by_key(|(x, _)| *x);
        r.into_iter().map(|(_, x)| x).collect::<Vec<_>>()
    };

    let all_contig_mapping_info = scan_contig_bam(
        &settings.assembly_to_ref_bam,
        derived_settings.thread_count,
        &ref_chrom_list,
        &assembly_contig_list,
        target_region.as_ref(),
    );

    scan_and_remap_reads(
        settings,
        derived_settings.thread_count,
        &reference,
        &ref_chrom_list,
        &all_contig_mapping_info,
        target_region.is_some(),
    );

    info!(
        "{PROGRAM_NAME} completed. Total Runtime: {}",
        start.elapsed().hhmmssxxx()
    );
    Ok(())
}

fn main() {
    let settings = cli::parse_settings();
    let derived_settings = validate_and_fix_settings(&settings);
    validate_settings_data(&settings);

    setup_logger(false).unwrap();

    if let Err(err) = run(&settings, &derived_settings) {
        eprintln!("{err}");
        process::exit(2);
    }
}
