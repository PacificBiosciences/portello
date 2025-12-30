use camino::{Utf8Path, Utf8PathBuf};
use chrono::Datelike;
use clap::{Parser, ValueEnum};
use simple_error::{SimpleResult, bail};

use crate::globals::PROGRAM_VERSION;

/// Input assembly mode options
///
/// Used to select the strategy to annotate phase information on reads
///
/// Note that enum values and their doc comments appear in the cmdline help.
///
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, strum::Display, ValueEnum)]
#[allow(non_camel_case_types)]
pub enum InputAssemblyMode {
    /// Set this option if input assembly contigs are fully phased (by trio-binning, Hi-C, etc..) such that the assembly haplotypes
    /// should be trusted when making read phasing haplotags.
    ///
    fully_phased,

    /// Set this option if input assembly assembly contigs are psuedo-phased/dual-assembly output
    /// that we assume contain many haplotype swich errors. In this case read support will be used
    /// to segment smaller phase sets for read phasing haplotags.
    ///
    partially_phased,
}

#[derive(Parser)]
#[command(
    author,
    version = PROGRAM_VERSION,
    about,
    after_help = format!("Copyright (C) 2004-{}     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.", chrono::Utc::now().year()),
    help_template = "\
{before-help}{name} {version}
{author-with-newline}{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}"
)]
#[clap(rename_all = "kebab_case")]
pub struct Settings {
    /// Assembly contig to reference genome alignment file in BAM/CRAM format
    ///
    /// Alignment file must be sorted and indexed. Minimap2 or pbmm2 alignments can be used.
    ///
    #[arg(long = "assembly-to-ref", value_name = "FILE")]
    pub assembly_to_ref_bam: Utf8PathBuf,

    /// Read to assembly alignment file in BAM/CRAM format
    ///
    /// Alignment file must be sorted and indexed. Only pbmm2 alignments are supported for these alignments.
    ///
    #[arg(long = "read-to-assembly", value_name = "FILE")]
    pub read_to_assembly_bam: Utf8PathBuf,

    /// Filename to use for remapped read output, our '-' for stdout
    ///
    /// This file is written in unsorted BAM format. If written to stdout, this output will be uncompressed to optimize
    /// piping into samtools sort or similar method.
    ///
    /// Unmapped reads in this file have a good alignment to an assembly contig region which is not mappable to the
    /// reference, they should not be remapped to the reference sequence.
    ///
    #[arg(long, value_name = "FILE")]
    pub remapped_read_output: Utf8PathBuf,

    /// Filename to use for unmapped reads which are not (well) mapped to any assembly contig
    ///
    /// This file is written in BAM format
    ///
    #[arg(long, value_name = "FILE")]
    pub unassembled_read_output: Utf8PathBuf,

    /// Genome reference in FASTA format
    #[arg(long = "ref", value_name = "FILE")]
    pub ref_filename: Utf8PathBuf,

    /// The input assembly mode changes portello's strategy for phasing annotation in the reampped read
    /// output.
    ///
    #[arg(long, value_enum, value_name = "MODE")]
    pub input_assembly_mode: InputAssemblyMode,

    /// Specify a target region for conversion
    ///
    /// This option is for debugging only, to allow fast test conversion of smaller regions. This will not give the same
    /// result as a non-targeted run with region selection on the sorted remapped read bam.
    ///
    #[arg(long)]
    pub target_region: Option<String>,

    /// Specify an output VCF prefix for phased heterozygous small variant calls from the input assembly contigs
    ///
    /// This option is only valid for the "partially-phased" input-assembly-mode, and is only intended for debugging and
    /// benchmarking the accuracy of the read haplotagging. Per this goal, only phase-eligible heterozygous variants are
    /// included in the output. It is intended to be used with a benchmarking tool such as whatshap compare.
    ///
    /// The output files produced by this option will be {PREFIX}.vcf.gz and {PREFIX}.vcf.gz.tbi
    ///
    #[arg(long, value_name = "PREFIX")]
    pub phased_het_vcf_prefix: Option<Utf8PathBuf>,

    /// Number of threads to use. Defaults to all logical cpus detected.
    #[arg(long = "threads", value_name = "THREAD_COUNT")]
    thread_count_option: Option<usize>,
}

/// Values immediately computed from the user settings, but not part of direct user inputs
///
pub struct DerivedSettings {
    /// Global thread count for pb-CpG-tools to use
    pub thread_count: usize,
}

/// Validate settings and use these to produce derived settings
///
fn validate_and_fix_settings_impl(settings: &Settings) -> SimpleResult<DerivedSettings> {
    fn check_required_filename(filename: &Utf8Path, label: &str) -> SimpleResult<()> {
        if filename.as_str().is_empty() {
            bail!("Must specify {label} file");
        }
        if !filename.exists() {
            bail!("Can't find specified {label} file: '{filename}'");
        }
        Ok(())
    }

    check_required_filename(&settings.assembly_to_ref_bam, "contig-to-ref bam")?;

    check_required_filename(&settings.read_to_assembly_bam, "read-to-contig bam")?;

    check_required_filename(&settings.ref_filename, "reference fasta")?;

    fn validate_output_dir_exists(filename: &Utf8Path, label: &str) -> SimpleResult<()> {
        if filename.as_str().is_empty() {
            bail!("Must specify {label} file");
        }
        let parent = filename.parent().unwrap();
        if !parent.as_str().is_empty() && !parent.exists() {
            bail!("Can't find existing directory for {label} file: '{filename}'");
        }
        Ok(())
    }

    if settings.remapped_read_output.as_str() != "-" {
        validate_output_dir_exists(&settings.remapped_read_output, "remapped read output")?;
    }

    validate_output_dir_exists(&settings.unassembled_read_output, "unassembled read output")?;

    if let Some(phased_het_vcf_prefix) = settings.phased_het_vcf_prefix.as_ref() {
        if settings.input_assembly_mode != InputAssemblyMode::partially_phased {
            bail!(
                "phased-het-vcf-prefix can only be specified for the partially-phased input-assembly-mode"
            );
        } else {
            validate_output_dir_exists(phased_het_vcf_prefix, "phased het output vcf")?;
        }
    }

    let thread_count = match settings.thread_count_option {
        Some(count) => {
            if count == 0 {
                bail!("--threads argument must be greater than 0");
            }
            count
        }
        None => num_cpus::get(),
    };

    Ok(DerivedSettings { thread_count })
}

pub fn validate_and_fix_settings(settings: &Settings) -> DerivedSettings {
    match validate_and_fix_settings_impl(settings) {
        Ok(x) => x,
        Err(msg) => {
            eprintln!("Invalid command-line setting: {msg}");
            std::process::exit(exitcode::USAGE);
        }
    }
}

fn assert_mapped_and_indexed_bam(filename: &Utf8Path) {
    use rust_htslib::bam::{self, Read};
    use rust_vc_utils::{ChromList, assert_bam_eof};

    // Note that IndexedReader is used here just to check that the index file is present
    //
    let bam_reader = match bam::IndexedReader::from_path(filename) {
        Ok(x) => x,
        Err(error) => {
            panic!("Failed to open input alignment file: {error}");
        }
    };
    assert_bam_eof(&bam_reader);

    let chrom_list = ChromList::from_bam_header(bam_reader.header());

    // Check for unmapped input
    if chrom_list.data.is_empty() {
        panic!("Input alignment file is not mapped: '{filename}'");
    }
}

/// Extended input data/settings validation that's too complex/slow to put in the cmdline parser
///
pub fn validate_settings_data(settings: &Settings) {
    assert_mapped_and_indexed_bam(&settings.assembly_to_ref_bam);
    assert_mapped_and_indexed_bam(&settings.read_to_assembly_bam);
}

pub fn parse_settings() -> Settings {
    Settings::parse()
}
