use camino::Utf8Path;
use log::info;
use rust_htslib::bcf;
use rust_vc_utils::ChromList;

use crate::vcf_utils;

#[allow(clippy::vec_init_then_push)]
fn get_vcf_header(ref_filename: &Utf8Path, chrom_list: &ChromList) -> bcf::Header {
    let sample_names = vec!["sample"];
    let mut header =
        vcf_utils::get_basic_vcf_header(ref_filename.as_ref(), chrom_list, &sample_names);

    let mut records = Vec::<&[u8]>::new();

    // Note that the `PASS` and `.` FILTER records below are not typically included in the header,
    // but there's an oddity in programatically creating vcfs in htslib that forces these to be present,
    // and specifically ordered before INFO/FORMAT.
    //
    records.push(br#"##FILTER=<ID=PASS,Description="All filters passed">"#);
    records.push(br#"##FILTER=<ID=.,Description="Unknown filtration status">"#);

    records.push(br#"##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"#);
    records.push(br#"##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">"#);

    for x in records.into_iter() {
        header.push_record(x);
    }

    header
}
/// Reformat het-site details to be more pertinant to VCF output:
pub struct VcfHetSite {
    pub chrom_index: usize,
    pub pos: i64,

    pub ref_seq: Vec<u8>,

    pub alt_seqs: Vec<Vec<u8>>,

    /// Indicate if the first alt corresponds to the HP1 or HP2 haplotype
    pub is_alt0_hp0: bool,

    /// This item is filled in later:
    pub phase_set: Option<i64>,
}

pub fn write_phased_het_vcf(
    ref_filename: &Utf8Path,
    chrom_list: &ChromList,
    thread_count: usize,
    phased_het_vcf: &str,
    mut vcf_het_sites: Vec<VcfHetSite>,
) {
    info!("Writing phased het VCF output to {}", phased_het_vcf);

    // sort vcf het sites
    vcf_het_sites.sort_by(|a, b| a.chrom_index.cmp(&b.chrom_index).then(a.pos.cmp(&b.pos)));

    {
        let header = get_vcf_header(ref_filename, chrom_list);
        let mut vcf =
            bcf::Writer::from_path(phased_het_vcf, &header, false, bcf::Format::Vcf).unwrap();

        for vcf_het_site in vcf_het_sites {
            let mut record = vcf.empty_record();
            record.set_rid(Some(vcf_het_site.chrom_index as u32));
            record.set_pos(vcf_het_site.pos);

            {
                let mut alleles = Vec::new();

                alleles.push(vcf_het_site.ref_seq.as_ref());

                for alt_seq in vcf_het_site.alt_seqs.iter() {
                    alleles.push(alt_seq.as_ref());
                }
                record.set_alleles(&alleles).unwrap();
            }
            // Add sample info:
            {
                use bcf::record::GenotypeAllele::*;

                let is_phased = vcf_het_site.phase_set.is_some();
                let gt = if is_phased {
                    if vcf_het_site.alt_seqs.len() == 1 {
                        if vcf_het_site.is_alt0_hp0 {
                            vec![i32::from(Phased(1)), i32::from(Phased(0))]
                        } else {
                            vec![i32::from(Phased(0)), i32::from(Phased(1))]
                        }
                    } else {
                        //
                        if vcf_het_site.is_alt0_hp0 {
                            vec![i32::from(Phased(1)), i32::from(Phased(2))]
                        } else {
                            vec![i32::from(Phased(2)), i32::from(Phased(1))]
                        }
                    }
                } else {
                    //
                    if vcf_het_site.alt_seqs.len() == 1 {
                        vec![i32::from(Unphased(0)), i32::from(Unphased(1))]
                    } else {
                        vec![i32::from(Unphased(1)), i32::from(Unphased(2))]
                    }
                };

                record.push_format_integer(b"GT", &gt).unwrap();

                let ps = if let Some(phase_set) = vcf_het_site.phase_set {
                    (phase_set + 1) as i32
                } else {
                    bcf::record::Numeric::missing()
                };
                record.push_format_integer(b"PS", &[ps]).unwrap();
            }

            vcf.write(&record).unwrap();
        }
    }

    let build_tbi = true;
    vcf_utils::build_bcf_index(phased_het_vcf, thread_count, build_tbi).unwrap();
}
