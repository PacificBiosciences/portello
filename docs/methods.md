# Portello methods

## Method Overview

Contig alignments are scanned to store the full set of split read alignments for each contig. By the bam spec, it should
be possible to do this from the primary alignment record alone, but minimap2 provides approximate CIGAR strings in the
SA tag for split reads, so the full bam must be traversed to store the full CIGAR from all supplementary alignment
records and fill these into the contig split alignment structures.

The read alignments are scanned for primary alignments only. Unlike the contig alignments, portello requires that these
reads must have complete `SA` tag CIGAR strings, because it would be problematic to store a data structure for each read
and fill everything in at the end.

In the primary read liftover loop, each primary read mapping is analyzed to get its full split alignment structure, then
lifted over according to the the relevant contig-to-ref mapping patterns.

## Selection of primary alignment

The read split alignment segment selected as primary is based on the following criteria:
  1. The read segment mapping to the contig with the highest contig-to-ref MAPQ score
  2. The read segment starting first in read sequencing order

