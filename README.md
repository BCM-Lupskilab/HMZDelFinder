# HMZDelFinder
CNV calling algorithm for detection of homozygous and hemizygous deletions from whole exome sequencing data



## Prerequisites
Followin R libraries are required to run HMZDelFinder:
 * gdata
 * data.table
 * DNAcopy
 * GenomicRanges
 * parallel
 * Hmisc

## Running HMZDelFinder
1. Download HMZDelFinder.R and example_run.R to your working directory
2. Collect BAMs and VCFs for all samples you want to include in the analysis
3. Generate RPKM files from BAMs (e.g. using featureCount method from Rsubread R package)
4. Collect RPKM files in the single directory
5. Obtain your design bed file 
6. In the example_run.R, replace paths (PATHS section) according to location of your data
7. execute example_run.R
