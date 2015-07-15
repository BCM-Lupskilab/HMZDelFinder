
########################################
# THRESHOLDS
#
# See description of HMZDelFinder function for details
########################################
is_cmg <- FALSE 		# only for CMG project - otherwhise use FALSE
lowRPKMthreshold <- 0.5 # RPKM threshold 
maxFrequency <- 0.005	# max frequncy of HMZ deletion
minAOHsize <- 1000		# min AOH size
minAOHsig <- 0.45		# min AOH signal threshold
mc.cores=8 				# number of cores

########################################
#  PATHS
#
# Replace paths in this section
########################################

mainDir <- "/mnt/bigData/cmg/"											# main working directory
bedFile <- paste(mainDir, 
		"SRC/CNV/conifer/conifer_v0.2.2/vcrome2.1_hg19_gn.bed",sep="")	# path to the design bed file
aohRDataOut <- paste(mainDir, "AOH/extAOH_small.RData", sep="")			# temporary output file to store AOH data
outputDir <- paste(mainDir, "HOMDEL_OUT/" , sep="")						# output direcory
snps <- dir (paste(mainDir, "INPUT/",sep=""),
		"SNPs_Annotated.vcf.bz2$", recursive=T, 
		include.dirs=FALSE)[1:200]
snpFids <- sapply(strsplit(snps,"[/\\.]"), function(x){x[2]})			# List of sample names in the same 
																		# order as the list of paths to SNP data
snpPaths<- paste(paste(mainDir, "INPUT/",sep=""),
		snpFids,"/", 
		sapply(strsplit(snps,"/"), 
		function(x){x[2]}), sep="")										# List of the paths to SNP data
rpkmFiles <- dir(paste(mainDir, "RPKM/",sep=""), ".txt$")
rpkmFids <- gsub(".rpkm.txt", "", rpkmFiles) 
fidsSubset <- which(rpkmFids %in% snpFids) 
rpkmFids <- rpkmFids[fidsSubset] 										# List of sample names in the same 
																		# order as the list of paths to RPKM data
rpkmPaths <- paste(paste(mainDir, "RPKM/",sep=""), 
					rpkmFiles[fidsSubset], sep="") 						#  List of the paths to RPKM data


########################################
# MAIN
#
# Example run of HMZDelFinder
########################################
	# running HMZ del calling algorithm
	source("HMZDelFinder.R")
	results <- runHMZDelFinder (snpPaths, snpFids,
							rpkmPaths, rpkmFids,
							mc.cores, aohRDataOut,
							bedFile, lowRPKMthreshold,
							minAOHsize, minAOHsig, is_cmg)
					
	## saving the list of calls
	write.csv(results$filterdCalls,paste(outputDir, "/hom_dels_small.csv", sep=""))
	## plots
	dir.create(paste(outputDir, "plots_example/", sep=""))
	lapply(1:nrow(results$filteredCalls), function(i){
				print(i); 
				plotDeletion(results$filteredCalls, i, 
							results$bedOrdered, 
							results$rpkmDtOrdered , 
							paste(outputDir, "plots_example/", sep=""))})

