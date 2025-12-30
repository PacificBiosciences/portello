# Portello user guide

## Table of Contents

* [Overview](#overview)
* [Getting started](#getting-started)
* [Inputs](#inputs)
* [Outputs](#outputs)
* [Other Notes](#other-notes)

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

Note that the choice of mappers and mapping parameters can have a substantial impact on downstream small variant calling
quality. It is recommended to carefully review the corresponding [Inputs](#inputs) section before starting.

## Getting Started

### Installation

Portello binaries are available for 64-bit Linux platforms. These can be installed either directly from the GitHub
release tarball, or via conda as described below.

#### Install from GitHub

To install portello from github, download the latest release tarball compiled for 64-bit Linux on the [github release
channel](https://github.com/PacificBiosciences/portello/releases/latest), then unpack the tar file. Using v0.7.0 as an
example, the tar file can be downloaded and unpacked as follows:

    wget https://github.com/PacificBiosciences/portello/releases/download/v0.7.0/portello-v0.7.0-x86_64-unknown-linux-gnu.tar.gz
    tar -xzf portello-v0.7.0-x86_64-unknown-linux-gnu.tar.gz

The portello binary is found in the `bin/` directory of the unpacked file distribution. This can be run with the help
option to test the binary and review latest usage details:

    portello-v0.7.0-x86_64-unknown-linux-gnu/bin/portello --help

#### Install from conda

For [conda](https://github.com/conda/conda) users, installing portello on conda may be a more convenient option. Portello
is available for conda on Linux from the `bioconda` channel. A new conda environment with the latest portello release can
be created as follows:

    conda create -n portello -c bioconda portello

### Usage Example

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
  --input-assembly-mode partially-phased \
  --remapped-read-output - \
  --unassembled-read-output unassembled.bam |\
samtools sort -@$threads - --write-index -o remapped.sort.bam
```

The output `remapped.sort.bam` will include standard haplotagging annotations using the `PS` and `HP` tags, which should
help to visualize typical diploid regions. For more complicated regions, it is best to group reads in IGV by the `ZS`
tag, which annotates the corresponding contig alignment segment for each read. See the [Tags](#tags) section below for
details of the `ZS` tag string.

## Inputs

Requirements and best practice for generation of the two input alignment files used by portello are provided below.

### Diploid assembly

Haplotype specific contigs from a diploid assembler are required as a first step to generate portello inputs. These
contigs can be either partially-phased/dual-assembly haplotypes or fully-phased haplotypes from a process such as
trio-binning. Portello has been tested with such inputs from both [hifiasm](https://github.com/chhylp123/hifiasm) and
[verkko](https://github.com/marbl/verkko).

For either assembler, the contigs for each haplotype should be extracted to fasta format, and the two sets of
haplotype specific contigs should be concatenated into a single fasta for use in the alignment steps below.

Note that the type of input contigs need to be conveyed to portello as well so that its phasing/read-haplotagging
routines can be appropriately set. Use `--input-assembly-mode partially-phased` for partially-phased/dual assembly
output and `--input-assembly-mode fully-phased` otherwise. See the [Phasing](#phasing) section below for more details.

Note that while portello has been primarily tested with assemlby contigs, it can work with unitig input as well. Phasing
and small variant calling quality appear to be lower from unitig input compared to the partially-phased contig output
from the same assembly.

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
input read should be represented by exactly one primary read in one of the two output files. Both of these output files
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
may represent population-specific and large non-templated insertion sequences.

## Secondary analysis on portello remapped output

### Small variant calling

Portello remapped read output has been tested on DeepVariant and some other secondary analysis tools designed for
conventional pbmm2-mapped BAM inputs. Results generally seem very good. DeepVariant, run with current best-practice
models for conventionally mapped BAM inputs, typically shows slightly higher accuracy on portello BAMs. Although the
primary benefits of assembly-based liftover are for structural variants and larger-scale events rather than small
variant calling, this DeepVariant performance trend demonstrates that portello alignments are properly aligned in
other contexts as well.

### Phasing / IGV visualization

For viewing in IGV, portello provides haplotagged BAM output automatically, using the read to contig mappings to
implement a variation on read-backed phasing, with this logic adjusted appropriately for either partially-phased or
fully-phased contig input. This means that typical phased read viewing patterns relying on `PS` and `HP` tag annotations
on the reads will work in typical diploid regions.

For viewing the reads in IGV in more complex regions, the reads can be grouped by the `ZS` tag, which annotates each read
according to the assembly contig alignment segment it was aligned to. This creates a phasing-like view of the reads, but the viewer
should note the potential for phasing switch errors when these annotation are derived from partially-phased contigs.

See the [Phasing](#phasing) and [Tags](#tags) sections below for more details.

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

### Phasing

Portello will automatically add haplotagging annotations (`HP` and `PS` tags) to the remapped read output based on the
read-to-contig alignments. These annotations are currently restricted to simple diploid regions where exactly two contig
alignment segments span the reference. The phasing strategy will change based on whether the input assembly contigs are
fully-phased or partially-phased (as conveyed to portello by the command-line `--input-assembly-mode` setting).

For fully-phased contigs, the phasing annotation will create phase-set blocks spanning each region where two contig
alignment segments span the reference, assuming the contigs themselves are free of switch errors.

For partially-phased contigs, the phasing annotation method will start with the larger phase-set annotations from the
fully-phased annotation mode, but then further subdivide phase-set blocks according to whether any reads span
heterozygous variant sites. The method is similar to mapping-based read backed phasing except that the read support
is evaluated based on alignments to the contigs (the phasing evaluation occurs before read-to-reference re-mapping).
The partially-phased haplotagging mode will slightly increase portello'w runtime and memory usage.

#### Phased heterozygous VCF output

When run with the `partiall-phased` input assembly mode, an optional VCF file can be produced containing all
phase-eligible heterozygous small variants with their corresponding phasing annotations. To produce this VCF, provide a
file prefix using the `--phased-het-vcf-prefix` argument. The VCF output produced by this procedure is only intended for
debugging and benchmarking the BAM haplotag output -- the heterozygous variants provided in this output file (and used
for the phasing procedure) are produced directly from differences in the assembly contig sequencers, so should not be as
accurate as the output of a dedicated small variant calling method such as DeepVariant.

### Tags

The following standard haplotagging values are added to the remapped read output:

`HP` - The `HP` tag with value 1 or 2 is added in regions with coverage by two contig alignment segments.

`PS` - The `PS` tag provides the phase set id in diploid regions. The PS value corresponds to the start position of the
first heterozygous variant in the phase set block.

The following custom diagnostic tags added to the remapped read output, these could change in future updates:

`ZS` - This tag can be used for read grouping in IGV to create a pseudo-phased view of the alignments. It is
`{contig_name}_split{split_alignment_no}{strand}`, where `contig_name` is the contig the read aligned to, `split
alignment no` refers to which of that contig's (colinear-joined) split alignment segments the read aligned to, and
`strand` is "+" or "-" for the contig alignments orientation to the reference.

`ZM` - original MAPQ of the read alignments

### Determinism

All bam outputs from portello are unsorted, and with a non-deterministic output order. However, if fully sorted, the bam
output contents should be identical over repeated runs with matching inputs.
