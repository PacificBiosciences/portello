# Change Log

## v0.6.1 - 2025-12-16

- Fix error message to describe chromosome consistency issues between assembly-to-ref alignments and the reference fasta [PacificBiosciences/portello#1]

## v0.6.0 - 2025-09-15

- CR-532 Join colinear contig alignments
  - These are splits in the contig alignment triggered by minimap2 Z-drop criteria.
- Fix repeated match trimming for assembly contigs
- Fix target region option for debugging

## v0.5.0 - 2025-08-07

- CR-499 Add repeated match trimming to clip repeated assembly contig mapping regions down to a single copy
  - As a result of this update, each read base should now map to no more than one reference site in the remapped bam
    output (if following best-practice protocols for the input alignments)

## v0.4.0 - 2025-07-27

- Simplify complex indels resulting from alignment liftover
  - The liftover process can create combined deletion/insertion events in the alignment. These are now locally simplified to minimize indel size and/or edit distance.
- CR-501 Improve indel left-shift consistency.
  - Reads aligning to assembly contigs which have a reverse mapping to the reference need to be updated from right-shift to left-shift when the alignment orientation flips.
- Add support for cram file inputs which include bzip2 and lzma codec blocks

## v0.3.0 - 2025-07-09

- Add initial documentation
- Complete output scheme so that all input reads are output to one of the two alignment output files

## v0.2.0 - 2025-07-07

- Update recommended input alignment methods
- Add several validity checks to input alignment files
- Fix to work with alignment files already containing 'PS'/'ZM' flags

## v0.1.0 - 2025-06-26

Initial proto release

