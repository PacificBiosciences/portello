# Portello methods

## Method Overview

For both contig and read alignment inputs, all secondary mappings are ignored.

### Contig alignment processing

Contig alignments are scanned to store the full set of split read alignments for each contig. By the bam spec, it should
be possible to do this from the primary alignment record alone, but minimap2 provides approximate CIGAR strings in the
SA tag for split reads, so the full bam must be traversed to store the full CIGAR from all supplementary alignment
records and fill these into the contig split alignment structures after the bam scan is complete.

After scanning, contig alignments are reviewed for repeated matches to the reference, meaning segments of the contig
which are mapped to two different genomic locations across two split read alignments. These repeats are removed from
the contig mapping information by selecting the segment match location with high gap-compressed identity and clipping
out any other mappings of that segment in other split reads.

### Read alignment scan

Read alignments are scanned for primary alignments only. Unlike the contig alignments, portello requires that these
reads must have complete `SA` tag CIGAR strings, because it would be problematic to store a data structure for each read
and fill everything in at the end.

## Read alignment liftover

In the primary read liftover loop, each primary read mapping is analyzed to get its full split alignment structure, then
lifted over according to the the relevant contig-to-ref mapping patterns.

For reads that map to contigs which themselves have a reverse mapping to the reference, the conventional alignment liftover
would result in right-shifted indels (assuming the input read mapping process left-shifts indels). For these reads, prior to
liftover, all indels are right-shifted, so that they will result in a lift-over alignment with left-shifted indels.

## Selection of primary alignment

The read split alignment segment selected as primary is based on the following criteria:
  1. The read segment mapping to the contig with the highest contig-to-ref MAPQ score
  2. The read segment starting first in read sequencing order

