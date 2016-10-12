# HMZDelFinder
CNV calling algorithm for detection of rare, homozygous and hemizygous deletions from whole exome sequencing data



## Prerequisites
* R in version >= 3.0.1 
Following R libraries are required to run HMZDelFinder:
 * RCurl (version >= 1.95.4.7)
 * gdata (version >= 2.17.0)
 * data.table (version >= 1.9.6)
 * DNAcopy (version >= 1.36.0)
 * GenomicRanges (version >= 1.14.4)
 * parallel (version >= 3.0.1)
 * Hmisc (version >= 3.16.0)
 * matrixStats (version >= 0.50.2)
 * Rsubread (version >= 1.20.3)

To install missing packages, run the code from the appropriate sections ('install missing packages from ...') at  example/example_run.R

## Running HMZDelFinder

* Working example that runs HMZDelFinder on 50 samples from 1000 genomes is available at example/example_run.R


