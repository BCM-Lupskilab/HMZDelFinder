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

## Format of input files

### BED file

Tab delimited file without header and four columns: 
* Chromosome
* Start
* Stop
* Gene symbol

### RPKM files

Tab delimited file with a header and two columns:
* count <i>// the number of reads overlapping with each capture target</i>
* RPKM  <i>// the RPKM value for each capture target</i>

IMPORTANT: The number of rows and the order of capture targets have to correspond to the number of rows and the order defined in the BED file. 

To generate RPKM files from BAM files, see comments at example/example_run.R.


### VCF files

VCF files are required for AOH analysis and further filtering of identfied deletion calls. 




## Format of output files

### Object returned by runHMZDelFinder(...), contains the following items:
* filteredCalls   <i>// the list of calls after AOH and deletion size filtering</i>
* allCalls        <i>// the list of calls before AOH and deletion size filtering </i>
* bedOrdered      <i>// the data.table containing ordered coordinates of coverage targets</i>
* rpkmDtOrdered   <i>// the data.table containing RPKM data for all samples; rowNames corresponds to sample identifiers and columns to coverage targets.</i>


### Format of filteredCalls/allCalls

Both objects are data.frames with the following columns:

* Chr         <i>// Chromosome </i>
* Start       <i>// Start position of deletion call</i>
* Stop        <i>// End position of deletion call</i>
* Genes       <i>// Comma separated list of genes encompassed by deletion</i>
* Start_idx   <i>// Index of first target</i>
* Mark_num    <i>// Number of targets that indicate deletion</i>
* Exon_num    <i>// Number of exons afftected (total number of targets encompassed by deletion)</i>
* FID         <i>// Sample identifier</i>
* Length      <i>// Length of deletion</i>
* BAB         <i>// Internal sample name (used only for BHCMG samples)</i>
* project     <i>// Project name (used only for BHCMG samples)</i>
* PoorSample  <i>// TRUE if the number of calls in the sample > 98 quantile</i>
* posKey      <i>// (chr+start+stop)</i>
* key         <i>// (sampleId+chr+start+stop)</i>
* inAOH_1000  <i>//  TRUE if deletion overlap with any AOH region greater than 1000bp</i>
* ZScore      <i>// z-score</i>
* OverlapCnt  <i>// number of overlapping calls in other samples</i>
* PerSampleNr <i>// number of calls in this sample</i>


