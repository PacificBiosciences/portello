# Portello user guide

## Overview

Portello is a method to liftover reads from sample assembly to reference
- Reads are lifted **from** their mapping to the denovo-assembly contigs for the sample (it is assumed that the contigs
  were assembled from this same set of reads).
- Reads are lifted **to** a standard reference genome, like GRCh38.

The expected workflow to use it is:
- Run diploid assembly on sample reads
- Consolidate all assembly contigs (from both haplotypes) to fasta
- Align sample reads directly to consolidated assembly contigs using pbmm2
- Align assembly contigs to standard reference (e.g. GRCh38) using minimap2
- Run portello on the two mapping outputs to lift sample read mappings to the standard reference.

## Usage Example:

Portello requires a reference genome and two alignment file inputs:
1. Assembly contig to reference alignments
2. Read to assembly contig alignments

See the inputs section below for more details on recommendations for these alignment files.

Given these two input alignment files represented as `$asm_to_ref_bam` and `$read_to_asm_bam`, a typical liftover
command is demonstrated below. Note that all direct BAM output from portello is unsorted, so piping the liftover reads
directly into a tool like `samtools sort` is a best practice to efficiently obtain an indexed liftover bam file.

```
threads=16
portello \
  --threads $threads \
  --ref $ref_fasta \
  --assembly-to-ref $asm_to_ref_bam \
  --read-to-assembly $read_to_asm_bam \
  --remapped-read-output - \
  --unassembled-read-output unasm.bam |\
samtools sort -@$threads - --write-index -o remapped.sort.bam
```

When viewing the output `remapped.sort.bam` in IGV, it is recommended to group on the `PS` tag, to see a pseudo-phased
view of the read alignments by grouping the reads to their respective assembly contigs. See the [Tags](#tags) section
below for details of the `PS` tag string.

## Inputs

Best practice and expectations for input alignment files are provided below. Note that secondary reads in either input
file will be ignored. It is assumed that both input alignment files are left-shifting alignment indels, portello will
take steps to preserve this left-shifting in the remapped output.

### Read-to-assembly alignments

For the read-to-assembly alignment input, all sequencing reads should be mapped to a reference genome comprised of all
high-quality assembly contigs for the given sample. Portello is designed and tested for diploid assemblies, so for this
case the contigs from both haplotypes should be concatenated together as the reference sequence.

The read-to-assembly alignment must be from pbmm2 so that portello can process the split read alignments correctly (in
particular the 'compressed' CIGAR strings in minimap2's `SA` tag output cannot be used).

If the reads are being re-mapped to the assembly contigs from an already-mapped bam, you must either use pbmm2 v1.17.0+
or else separately remove old `SA` aux tags from the bam file before providing this as input to pbmm2 to prevent stale
`SA` aux tags from interrupting the liftover process.

### Assembly-to-ref alignments

The best practice assembly-to-ref alignment uses a recent version of minimap2 (v2.26+) with the parameters from the
following example command:

```
minimap2 \
  -t $threads \
  -L \
  --secondary=no \
  -a \
  -x asm5 \
  -R "@RG\tID:HG002_hifiasm\tSM:HG002" \
  $ref \
  $reads |\
samtools sort -@$threads - --write-index -O BAM -o HG002.asm.GRCh38.bam
```

This strategy uses the `asm5` minimap2 preset to map sample-specific assembly contigs to the reference, which should be
a good choice for human sample analysis.

Note the common additional minimap2 options `--eqx`, `--cs` and `-Y` add more information to the alignments that
portello doesn't use. These won't help or hurt the process.

## Outputs

All BAM outputs from portello are unsorted. The two bam output files are for "unassembled" and "remapped" reads. Every
input read should be represented by exactly one primary read in one of the two outputs files. Both of these output files
may contain unmapped reads, but they're separated into two separate files because they may have different
interpretations relevant to downstream processing steps.

### Unassembled read output

This output separates out all of the unmapped reads from the input read-to-contig alignments. These are interpreted as
'unassembled' since they don't have a mapping to the assembly contigs, and the default portello workflow assumption is
that the input reads are the same ones used to generate the assembly.

Followup usage of this file may depend on workflow details. At higher depth, these reads are likely to correspond to
regions that are more difficult to assemble or have QC issues. At lower depth, there may be more reads in this set that
simply couldn't be assembled due to coverage, so conventional re-mapping of these reads could be useful to supplement
the liftover reads.

### Remapped read output

This output contains all reads that mapped to the assembly contigs, and have been lifted over to the target reference
genome. Unmapped reads in this file should correspond to assembly contig regions which did not map to the reference, so
may represent population-specific sequence or other content missing from the reference.


## Other Notes

### MAPQ

The MAPQ value of liftover read alignment segments is taken form the contig split alignment segment that that read
segment was mapped to. The read's original MAPQ to the contig alignments is stored in the record under `ZM`. This will
often be zero, which is most often a reflection of ambiguous alignment between the hap1 and hap2 contigs, so it seems
like this should be used for phasing quality at some point.

### Tags

Temporary diagnostic tags added to the remapped read output, these are likely to change in future:

`PS` - This tag can be used for read grouping in IGV to create a pseudo-phased view of the alignments. It is
`{contig_name}_split{split_alignment_no}{strand}`, where `contig_name` is the contig the read aligned to, `split
alignment no` refers to which of that contigs split alignments the read aligned to, and `strand` is "+" or "-" for the
contig alignments orientation to the reference.

`ZM` - original MAPQ of the read alignments
