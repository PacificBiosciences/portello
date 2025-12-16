mod cli;
mod contig_alignment_scanner;
mod globals;
mod liftover_read_alignment;
mod logger;
mod read_alignment_scanner;
mod simplify_alignment_indels;
mod worker_thread_data;

use std::process;

use hhmmss::Hhmmss;
use log::info;
use rust_vc_utils::{ChromList, GenomeSegment, get_genome_ref_from_fasta};

use crate::cli::{DerivedSettings, Settings, validate_and_fix_settings, validate_settings_data};
use crate::contig_alignment_scanner::scan_contig_bam;
use crate::globals::{PROGRAM_NAME, PROGRAM_VERSION};
use crate::logger::setup_logger;
use crate::read_alignment_scanner::scan_and_remap_reads;

/// Load reference seqeunce from fasta, stored as an array in ref_chrom_list order
///
fn get_chrom_array(ref_filename: &str, ref_chrom_list: &ChromList) -> Vec<Vec<u8>> {
    let mut genome_ref = get_genome_ref_from_fasta(ref_filename);

    // For each chrom listed in the asssembly-to-ref bam index, check that it is found in genome_ref,
    // and transfer the chrom sequence from genome_ref to a chrom array
    //
    let mut r = Vec::new();
    let mut reference_consistency_error = false;
    for chrom_info in ref_chrom_list.data.iter() {
        let chrom_label = &chrom_info.label;
        let chrom_len = chrom_info.length;
        match genome_ref.chroms.remove(chrom_label) {
            Some(x) => {
                if x.len() != chrom_len as usize {
                    log::error!(
                        "Chromosome \"{chrom_label}\" specified with inconsistent length: {chrom_len} in the assembly-to-ref alignment file, and {} in the reference fasta",
                        x.len()
                    );
                    reference_consistency_error = true;
                } else {
                    r.push(x);
                }
            }
            None => {
                log::error!(
                    "Chromosome \"{chrom_label}\" specified in the assembly-to-ref alignment file, but not in the reference fasta"
                );
                reference_consistency_error = true;
            }
        };
    }

    if reference_consistency_error {
        log::error!("Exiting due to one or more reference consistency issues");
        std::process::exit(exitcode::DATAERR);
    }

    r
}

fn run(
    settings: &Settings,
    derived_settings: &DerivedSettings,
) -> Result<(), Box<dyn std::error::Error>> {
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

    let reference = get_chrom_array(settings.ref_filename.as_str(), &ref_chrom_list);

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
