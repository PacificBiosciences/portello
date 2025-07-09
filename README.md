# Portello

Portello is a method to liftover reads from sample assembly to reference

- Reads are lifted **from** their mapping to the denovo-assembly contigs for the sample (it is assumed that the contigs were assembled from this same set of reads).

- Reads are lifted **to** a standard reference genome, like GRCh38.

Example IGV view showing portello liftover mappings vs standard pbmm2 mappings from the same set of reads:

<p align="center">
  <img src="img/portello_igv_example.png" width="80%" />
</p>

The example image show portello alignments (top) vs. pbmm2 (bottom) for the same set of HG002 HiFi reads at region `chr4:40,294,833-40,295,706`, the portello
reads are grouped by the `PS` tag to cluster them by their associated assembly contig, to create a phased-like view of the liftover output.

## Getting started

See the [user guide](docs/user_guide.md) for details on getting started with portello.

## DISCLAIMER

THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
