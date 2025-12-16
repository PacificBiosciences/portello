# Portello user guide

## Overview

Portello is a method to liftover reads from sample assembly to reference
- Reads are lifted **from** their mapping to the denovo-assembly contigs for the sample (it is assumed that the contigs
  were assembled from this same set of reads).
- Reads are lifted **to** a standard reference genome, like GRCh38.

The expected workflow to use it is:
- Run diploid assembly on sample reads
- Consolidate haplotype-resolved assembly contigs from both haplotypes to a single fasta
- Align assembly contigs to standard reference (e.g. GRCh38) using minimap2
- Align sample reads directly to consolidated assembly contigs using pbmm2
- Run portello on the above two alignment outputs to lift sample read mappings from the assembly to the reference.

## Usage Example:

Portello requires a reference genome and two alignment file inputs:
1. Assembly contig to reference alignments
2. Read to assembly contig alignments

See the inputs section below for more details on recommendations for these alignment files.

Given these two input alignment files represented as `$asm_to_ref_bam` and `$read_to_asm_bam`, a typical liftover
command is demonstrated below. Note that all direct BAM output from portello is unsorted, so piping the liftover reads
into a tool like `samtools sort` is a best practice to efficiently obtain an indexed liftover bam file, as illustrated
in the example below:

```
threads=16
portello \
  --threads $threads \
  --ref $ref_fasta \
  --assembly-to-ref $asm_to_ref_bam \
  --read-to-assembly $read_to_asm_bam \
  --remapped-read-output - \
  --unassembled-read-output unassembled.bam |\
samtools sort -@$threads - --write-index -o remapped.sort.bam
```

When viewing the output `remapped.sort.bam` in IGV, it is recommended to group on the `PS` tag, to see a pseudo-phased
view of the read alignments by grouping the reads to their respective assembly contigs. See the [Tags](#tags) section
below for details of the `PS` tag string.

## Inputs

Requirements and best practice for generation of the two input alignment files used by portello are provided below.

### Diploid assembly

Haplotype specific contigs from a diploid assembler are required as a first step to generate portello inputs. These
contigs can be either partially-phased/dual-assembly haplotypes or fully-phased haplotypes from a process such as
trio-binning. Portello has been tested with such inputs from both [hifiasm](https://github.com/chhylp123/hifiasm) and
[verkko](https://github.com/marbl/verkko).

For either assembler, the contigs for each haplotype should be extracted to fasta format, and the two sets of
haplotype specific contigs should be concatenated into a single fasta for use in the alignment steps below.

### Read-to-assembly alignments

For the read-to-assembly alignment input, all sequencing reads should be mapped to the concatenated diploid assembly
haplotypes (see [diploid assembly](#diploid-assembly) section above).

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
  --eqx \
  -x asm5 \
  -R "@RG\tID:HG002_hifiasm\tSM:HG002" \
  $ref \
  $contigs |\
samtools sort -@$threads - --write-index -O BAM -o HG002.asm.GRCh38.bam
```

Here, `$contigs` should be the concatenated diploid assembly haplotypes for the sample (see [diploid
assembly](#diploid-assembly) section above).

The example uses the `asm5` minimap2 preset to map the assembly contigs to the reference. This is the recommended
setting for human sample analysis. Note that the `--eqx` option to provide CIGAR strings with `=` and `X` match values
is required.

The common additional minimap2 options `--cs` and `-Y` add more information to the alignments that portello doesn't use.
These won't help or hurt the process.

### Additional input alignment notes

For both input alignments, secondary reads will be ignored. It is assumed that both input files have alignments with
left-shifted indels. Portello will take steps to preserve this left-shifting in the remapped output.

## Outputs

All BAM outputs from portello are unsorted. The two bam output files are for "unassembled" and "remapped" reads. Every
input read should be represented by exactly one primary read in one of the two outputs files. Both of these output files
may contain unmapped reads, but they're separated into two separate files because the unmapped reads in each file have
different interpretations relevant to downstream processing steps, as discussed below.

### Unassembled read output

This output contains all unmapped reads from the input read-to-contig alignment file. These reads are interpreted as
'unassembled' since they don't have a mapping to any assembly contig, and the default portello workflow assumption is
that the input reads are the same ones used to generate the assembly.

Followup usage of this file may depend on workflow details. At higher depth, these reads are likely to correspond to
regions that are more difficult to assemble or have QC issues. At lower depth, there may be more reads in this set that
simply couldn't be assembled due to coverage, so conventional re-mapping of these reads could be useful to supplement
the liftover reads.

### Remapped read output

This output contains all reads that mapped to the assembly contigs, and have been lifted over to the target reference
genome. Unmapped reads in this file should correspond to assembly contig regions which did not map to the reference, so
may represent population-specific and large insertion sequences.

## Secondary analysis on portello remapped output

Portello remapped read output has been tested on DeepVariant and some other secondary analysis tools designed for
conventional pbmm2-mapped BAM inputs. Results generally seem very good. DeepVariant, run with current best-practice
models for conventionally mapped BAM inputs, typically shows slightly higher accuracy on portello BAMs. Although the
primary benefits of assembly-based liftover are for structural variants and larger-scale events rather than small
variant calling, this DeepVariant performance trend demonstrates that portello alignments are properly aligned in
other contexts as well.

Note that although we recommend viewing portello alignments in IGV with reads grouped by the `PS` tag, we recommend that
you **Interpret the PS tag with care**. The PS tag in the portello output shows which contig each read was mapped to,
and more specifically which split read segment of the contig it was mapped to. When visualizing portello BAMs in IGV it
is extremely useful to show a pseudo-phased view of the reads by grouping them on this tag. Note, however that each
assembly algorithm has a different approach to contig contiguity through LOH regions, so a read segment with one PS tag
is not equivalent to a single haplotype phase block, since it could continue through large LOH regions with possible
phase switching.

## Other Notes

### Processing of contig alignments

Assembly contig alignments are processed through two major steps prior to being used for read liftover. This may lead to some
important differences if viewing both the contig-to-reference alignments and portello liftover read alignments together.

1. Repeated match trimming

In the repeated match trimming step, we identify regions of the assembly contig which have been repeated aligned to
different parts of the genome by minimap2. For each such region, we test which instance of the repeated region alignment
has the highest gap compressed identity and trim this segment out of all other alignments. This ensures that each base
of the assembly contig is mapped to no more than one base in the reference genome, and therefor, that each read going
through the liftover process preserves this same property, so long as the read-to-assembly input alignment does as well.
When following the portello best-practice protocols for input alignments, this will be the case, since pbmm2 is the
recommended mapper for all read-to-assembly input.

2. Joining colinear segments

After repeated match trimming, portello identifies adjacent split alignment segments of the assembly contigs which are
colinear, and joins these together. Joining such segments into one continuous alignment makes it easier to view the
alignments in IGV and recognize that the phasing relationship across the split alignment segment gap. Such gaps are
typically created by minimap2 at regions meeting its Z-drop criteria indicating punctuated sample divergence from the
reference. Note that due to this joining step, the 'split read index' shown in the portello `PS` tag output refers to
the index of the joined split alignment segments.


### MAPQ

The MAPQ value of liftover read alignment segments is taken form the contig split alignment segment that that read
segment was mapped to. The read's original MAPQ to the contig alignments is stored in the record under `ZM`. This will
often be zero, which is most often a reflection of ambiguous alignment between the 2 parental contigs.

### Tags

The following diagnostic tags added to the remapped read output, these are likely to change in future:

`PS` - This tag can be used for read grouping in IGV to create a pseudo-phased view of the alignments. It is
`{contig_name}_split{split_alignment_no}{strand}`, where `contig_name` is the contig the read aligned to, `split
alignment no` refers to which of that contig's (colinear-joined) split alignment segments the read aligned to, and
`strand` is "+" or "-" for the contig alignments orientation to the reference.

`ZM` - original MAPQ of the read alignments

### Determinism

All bam outputs from portello are unsorted, and with a non-deterministic output order. However, if fully sorted, the bam
output contents should be identical over repeated runs with matching inputs.
