## functions for reading vcf...
getKVdf <- function(xx, split.char=",",removequotes=F){
	yy <- strsplit(xx, split.char)[[1]]
	kv <- do.call(rbind, lapply ( (strsplit(yy, "=")), function(ll){cbind(ll[1], paste(ll[2:length(ll)], collapse="="))}))
	colnames(kv) <- c("key", "value")
	if (removequotes){
		kv <- gsub ( "\"", "", kv)
	}
	kv <- data.frame(kv, stringsAsFactors=F)
	colnames(kv) <- c("key", "value")
	kv
}
read.vcf.header <- function(file){
	line <- scan(file=pipe(paste ("bzcat", file)), n=1, skip=0, sep="\n", what="character", quiet=TRUE)
	skip <- 1
	header <- c()
	while(substr(line, 0 ,1) == "#" ){
		header <- c(header, gsub("#", "", line) )
		line <- scan(file=file, n=1, skip=skip, sep="\n", what="character", quiet=TRUE)
		skip <- skip + 1
	}
	lines <- header
	getSubHeader <- function(headerTitle, lines){
		starts <- regexpr(headerTitle, lines) + nchar(headerTitle) 
		stops <- regexpr( ">", lines) -1

		hd <- lapply(1:length(starts), function(i){
					start <- starts[i]
					stop <- stops[i]
					if (start > nchar(headerTitle) )
					{
						xx <- substr(lines[i], start, stop)
						kv <- getKVdf(xx, removequotes=T)
					}
					else {NULL}
				})
		hdf <- Filter(function(x) !is.null(x), hd)
		names(hdf) <- lapply(hdf, function(x)x$value[which(x$key=="ID")])
		hdf
	}
	info <- getSubHeader( "INFO=<", lines)
	format <-  getSubHeader( "FORMAT=<", lines)
	list (header=strsplit(lines[length(lines)], "\t")[[1]], INFO=info,FORMAT =format, nlines=length(lines) + 1)
}
read.vcf.quick.noinfo <- function(file){
	header <- read.vcf.header(file)
	#system (paste("bzcat ",file," | tail -n +",header$nlines," | cut -f1-7,9,10 > " , file, ".noheader",  sep="" ))
	#data <- fread(paste(file , ".noheader", sep=""), header=FALSE, stringsAsFactors=F, sep="\t")
	data <- fread(paste0(" bzcat ",file," | tail -n +",header$nlines," | cut -f1-7,9,10"), header=FALSE, stringsAsFactors=F, sep="\t")
	file.remove(paste (file, ".noheader", sep=""))
	setnames(data, header$header[-9])
	data
}

## end of functions for reading vcf

## download file from dropbox
dl_from_dropbox <- function(x, key, out) {	
	print (paste0("Downloading from dropbox:" , x))
	bin <- getBinaryURL(paste0("https://dl.dropboxusercontent.com/s/", key, "/", x),ssl.verifypeer = FALSE)
	con <- file(x, open = "wb")
	writeBin(bin, con)
	close(con)
}


##------------------------------------------------------------------------------
##' Wrapper around mclapply to track progress
##' 
##' Based on http://stackoverflow.com/questions/10984556
##' 
##' @param X         a vector (atomic or list) or an expressions vector. Other
##'                  objects (including classed objects) will be coerced by
##'                  ‘as.list’
##' @param FUN       the function to be applied to
##' @param ...       optional arguments to ‘FUN’
##' @param mc.preschedule see mclapply
##' @param mc.set.seed see mclapply
##' @param mc.silent see mclapply
##' @param mc.cores see mclapply
##' @param mc.cleanup see mclapply
##' @param mc.allow.recursive see mclapply
##' @param mc.progress track progress?
##' @param mc.style    style of progress bar (see txtProgressBar)
##'
##' @examples
##' x <- mclapply2(1:1000, function(i, y) Sys.sleep(0.01))
##' x <- mclapply2(1:10, function(i){Sys.sleep(1)}, mc.cores=4)
##------------------------------------------------------------------------------
mclapply2 <- function(X, FUN, ..., 
		mc.preschedule = TRUE, mc.set.seed = TRUE,
		mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
		mc.cleanup = TRUE, mc.allow.recursive = TRUE,
		mc.progress=TRUE, mc.style=3) 
{
	if (!is.vector(X) || is.object(X)) X <- as.list(X)
	
	if (mc.progress) {
		f <- fifo(tempfile(), open="w+b", blocking=T)
		p <- parallel:::mcfork()
		#print(p$pid)
		pb <- txtProgressBar(0, length(X), style=mc.style)
		setTxtProgressBar(pb, 0) 
		progress <- 0
		if (inherits(p, "masterProcess")) {
			while (progress < length(X)) {
				readBin(f, "double")
				progress <- progress + 1
				setTxtProgressBar(pb, progress) 
			}
			cat("\n")
			#print("ending ...")
			parallel:::mcexit()
		}
	}
	tryCatch({
				result <- mclapply(X, function(...) {
							res <- FUN(...)
							if (mc.progress) writeBin(1, f)
							res
						}, 
						mc.preschedule = mc.preschedule, mc.set.seed = mc.set.seed,
						mc.silent = mc.silent, mc.cores = mc.cores,
						mc.cleanup = mc.cleanup, mc.allow.recursive = mc.allow.recursive
				)
				
			}, finally = {
				if (mc.progress) close(f)
				mccollect(p)
			})
	result
}


##------------------------------------------------------------------------------
##' Reads VCF files and generates AOH data using CBS
##' 
##' 
##' @param fileNames	a vector of paths to VCF files
##' @param fids			a vector of  sample names of the same length as fileName
##' @param mc.cores		number of cores (see mclapply)
##------------------------------------------------------------------------------
prepareAOHData <- function(fileNames, fids , mc.cores)
{
	
	print ("Reading VCFs and generating AOH data using CBS ...")
	allAOH <- mclapply2(1:length(fileNames), function(i){
				res <- NULL
				tryCatch({
							fileSNP<- fileNames[i]
							dataSNP <- read.vcf.quick.noinfo(fileSNP)
							# parsing genotype
							gt<- do.call(rbind, strsplit(dataSNP[[9]], ":"))
							vR <- gt[,2] # number of alt allele reads
							tR <- gt[,4] # total number of reads
							dataSNP[,vR:=as.numeric(vR)]
							dataSNP[,tR:=as.numeric(tR)]
							dataSNP[,af:=(as.numeric(vR)/ as.numeric(tR))]
							# Circular Binary Segmentation
							dataSNP <- dataSNP[which(dataSNP$FILTER == "PASS"),]
							dataSNP[,baf:=abs(af - 0.5)]				
							CNA.obj = CNA(dataSNP$baf, dataSNP$CHROM, dataSNP$POS, data.type = "binary")
							segment.obj = segment(CNA.obj,   verbose = 0)
							out <- segment.obj$output
							out$Length <- out$loc.end - out$loc.start
							out$Name <-fids[i]
							res <- data.table(out)
						}, error=function(e){},finally={})
				res
			},mc.cores=mc.cores)
	
	aoh <- rbindlist(allAOH)
	print ("Extending AOHs by length of uncertain regions ...")
	#browser()
	# extends AOH by length of uncertain regions; note that new segments will overlap
	extAOH <- rbindlist(mclapply2(1:length(unique(aoh$Name)), function(i){ nm <- unique(aoh$Name)[i]
						tmpAOH <- aoh[aoh$Name == nm,]
						tmpAOH <- tmpAOH[order(tmpAOH$chrom, tmpAOH$loc.start),]
						newAOH <- rbindlist(lapply(c(1:22,"X","Y"), function(chr){
											chrAOH <- tmpAOH[tmpAOH$chrom == chr, ]
											chrAOH$loc.start.ext <- 0
											if (nrow(chrAOH) > 1) chrAOH$loc.start.ext[2:(nrow(chrAOH))] <- chrAOH$loc.end[1:(nrow(chrAOH)-1)]
											chrAOH$loc.end.ext <- chrAOH$loc.end
											if (nrow(chrAOH) > 1) chrAOH$loc.end.ext[1:(nrow(chrAOH)-1)] <- chrAOH$loc.start[2:nrow(chrAOH)]
											chrAOH
										}))
						newAOH
					}, mc.cores=mc.cores))
	extAOH$Length <- extAOH$loc.end.ext - extAOH$loc.start.ext
	extAOH
}


##------------------------------------------------------------------------------
##' Reads RPKM files and creates rpkmDt data.table with
##' samples in rows and exons in columns
##' 
##' 
##' @param fileNames	a vector of paths to RPKM files
##' @param fids			a vector of sample names of the same length as fileNames
##' @param mc.cores		number of cores (see mclapply)
##------------------------------------------------------------------------------
prepareRPKMData <- function(fileNames, fids, mc.cores)
{
	print ("Reading RPKM files ...")
	rpkmList <- mclapply2(1:length(fileNames), function(i){
				file <- fileNames[i]; fid <- fids[i]
				if (file.info(file)$size ==0) {return(NULL)}
				t <- fread(file)
				t$File <- fid
				t
			}, mc.cores=mc.cores)
	
	names(rpkmList) <- fids
	print ("Removing empty elements ...")
	rpkmList2 <- Filter(function(x){!is.null(x)}, rpkmList)
	print ("Creating matrix ...")
	rpkmDf <- do.call(rbind,mclapply2(rpkmList, function(x){x$RPKM}, mc.cores=mc.cores))
	rownames(rpkmDf) <- names(rpkmList2)
	rm(rpkmList);rm(rpkmList2);gc();gc();
	sn <- rownames(rpkmDf)
	print ("Creating data.table (may take a while)...")
	rpkmDt <- data.table(rpkmDf)
	rm(rpkmDf)
	rownames(rpkmDt) <- sn
	rpkmDt
}

##------------------------------------------------------------------------------
##' Reorder BED and RPKM matrices
##' 
##' 
##' @param bedFile    tab separeated bed file containing capture design probes with the following 4 columns:
##'                    - chromosme
##'                    - start
##'                    - stop
##'                    - gene name
##' @param rpkmDt	 object returned by prepareRPKMData  
##------------------------------------------------------------------------------
reorderBedAndRpkmDt <- function(bedFile, rpkmDt)
{
	bed <- read.csv(bedFile, sep="\t", stringsAsFactors=F, header=F)
	ord <- order(bed$V1, bed$V2)
	bedOrdered <- bed[ord,]
	rpkmDtOrdered <- rpkmDt[,ord, with=F]
	toRemIdx <- which(duplicated(paste(bedOrdered$V1,"_", bedOrdered$V2, "_",bedOrdered$V3)))
	if (length(toRemIdx)>0){
		bedOrdered <- bedOrdered[-(toRemIdx),]
		rpkmDtOrdered<- rpkmDtOrdered[,-(toRemIdx),with=F]
	}
	rownames(rpkmDtOrdered) <- rownames(rpkmDt)
	list(bedOrdered=bedOrdered, rpkmDtOrdered=rpkmDtOrdered)
}


##------------------------------------------------------------------------------
##' Selects exons with potential HMZ deletions
##' 
##' 
##' @param rpkmDtOrdered	rpkmDtOrdered object returned by reorderBedAndRpkmDt
##' @param bedOrdered		bedOrdered object returned by reorderBedAndRpkmDt
##' @param mc.cores			number of cores to be used
##' @param lowRPKMthreshold	RPKM threshold
##' @param maxFrequency 	max frequency of HMZ deletion
##------------------------------------------------------------------------------
processRPKM <- function(rpkmDtOrdered, bedOrdered, mc.cores, lowRPKMthreshold ,exonsToExclude,  maxFrequency=0.005)
{
	print("Computing initial statistics for each exon...")
	a <- unlist(mclapply2(rpkmDtOrdered,function(x){sum(x<lowRPKMthreshold)}, mc.cores=mc.cores))
	gc();gc();
	maxHoms <- ceiling(nrow(rpkmDtOrdered)*maxFrequency)
	selectedExons <- setdiff(which( a >0 & a <= maxHoms ) , exonsToExclude)
	
	print("Selecting a subset of exons with potentiall deletions (may take a while)...")
	tmp <- rpkmDtOrdered[, selectedExons,with=F]
	bedTmp <- bedOrdered[selectedExons,]
	print("Calculating number of potential deletions for each sample (may take a while)...")
	nrOfPotDeletions <- apply(tmp,1,  function(x){length( which(x<lowRPKMthreshold & bedTmp$V1 != "X"))})
	print("Detetermining 5% of samples with the highest number of deletions")	
	toRemSamples <- which(nrOfPotDeletions> quantile(nrOfPotDeletions, 0.98))
	print("Removing these samples (may take a while)...")
	rpkmTmp <- rpkmDtOrdered[-toRemSamples,]
	print("Recomputing per-exon statistics without low quality samples ...")
	a1 <- unlist(mclapply2(rpkmTmp,function(x){sum(x<lowRPKMthreshold)}, mc.cores=mc.cores))
	rm(rpkmTmp)
	gc();gc();
	selectedExonsFinal <- setdiff((which( a1 >0 & a1 <= maxHoms )), exonsToExclude)
	list(selectedExonsFinal=selectedExonsFinal, a=a, a1=a1, toRemSamples=toRemSamples)
}

##------------------------------------------------------------------------------
##' Retrieves potential HMZ deletions
##' 
##' 
##' @param rpkmDtOrdered		rpkmDtOrdered object returned by reorderBedAndRpkmDt
##' @param bedOrdered			bedOrdered object returned by reorderBedAndRpkmDt
##' @param selectedExonsFinal	object returned by processRPKM
##' @param mc.cores				number of cores (see mclapply2)
##' @param lowRPKMthreshold		RPKM threshold
##------------------------------------------------------------------------------
getCandidateExonCalls <- function(rpkmDtOrdered, bedOrdered, selectedExonsFinal,  mc.cores, lowRPKMthreshold)
{
	selectedExonsDT<- rpkmDtOrdered[,selectedExonsFinal,with=F]
	print("Calling deletions")
	calls1 <- mclapply2((1:nrow(rpkmDtOrdered)), function(i){
				candidateDels <- data.table(bedOrdered[selectedExonsFinal,][which(selectedExonsDT[i,]<lowRPKMthreshold),])
				gc()
				if (nrow(candidateDels)==0){return(NULL)}
				candidateDels
			}, mc.cores=mc.cores)
	
	calls2 <- lapply((1:nrow(rpkmDtOrdered)), function(i){
				res <- calls1[[i]]
				if (!is.null(res))res$Sample <- rownames(rpkmDtOrdered)[i]
				res
			})
	candidates <- rbindlist(calls2)
	candidates
}


##------------------------------------------------------------------------------
##' Merges calls from consecutive exons
##' 
##' 
##' @param candidateExonCalls	object returned by getCandidateExonCalls
##' @param bedOrdered			bedOrdered object returned by reorderBedAndRpkmDt
##' @param maxGap				maximal distance between deleted exons to be merged
##------------------------------------------------------------------------------
mergeCandidates <- function(candidateExonCalls, bedOrdered, maxGap = 10)
{
	
	resTmp<- by(candidateExonCalls, candidateExonCalls$Sample, function(x){
				x$mark_num<-1
				x$exon_num<-1
				x$key <-paste(x$V1,"_",x$V2, sep="")
				bedOrderedIdx<- match(x$key, paste(bedOrdered$V1,"_", bedOrdered$V2, sep=""))
				df <- data.frame()
				j <- 1
				v  <-  bedOrdered[bedOrderedIdx[1],]
				v$mark_num<- 1
				v$exon_num <- 1
				v$start_idx<- bedOrderedIdx[1]
				diffIdx <- diff(bedOrderedIdx)
				for (i in diffIdx){
					j <- j + 1
					if (i >= maxGap){
						df <- rbind(df, v)
						v <-  bedOrdered[bedOrderedIdx[j],]
						v$mark_num <- 1
						v$exon_num <- 1
						v$start_idx<- bedOrderedIdx[j]
						
					}else{
						v$mark_num <- v$mark_num+1 
						v$V3 <-bedOrdered[bedOrderedIdx[j],"V3"] # replace stop
						gn <- bedOrdered[bedOrderedIdx[j],"V4"] # gene name
						if (gdata::trim(gn) != "" && !(gn %in% strsplit(v$V4,",")[[1]])){ v$V4 <- paste(v$V4, gn, sep=",")}
						v$exon_num <- bedOrderedIdx[j] - v$start_idx + 1
						z <- bedOrdered[bedOrderedIdx[j],]
					}
				}
				df <- rbind(df, v)
				df$Sample<-x$Sample[1]
				x <- as.data.frame(df)
				data.table(x[,c("V1", "V2", "V3", "V4","start_idx", "mark_num","exon_num","Sample")])
			})
	
	res <- rbindlist(resTmp)
	res <- res[order(res$V1, res$V2),]	
	colnames(res) <- c("Chr", "Start","Stop","Genes","Start_idx","Mark_num","Exon_num","FID")
	res$Length <- res$Stop - res$Start
	res[order(res$Length, decreasing=T), ]
	res
}


##------------------------------------------------------------------------------
##' Adds information about internal sample name and project (BHCMG specific)
##' 
##' 
##' @param candidatesMerged   object returned by mergeCandidates
##------------------------------------------------------------------------------
annotateCandidates <- function(candidatesMerged, is_cmg=FALSE)
{
	candidatesMerged$BAB <- candidatesMerged$FID;
	candidatesMerged$project <- ""
	if (is_cmg)
	{
		source("/mnt/bigData/cmg/SRC/CNV/data_handler/load_meta_func.R")
		report <- downloadReport()
		report$Fid <- sapply(strsplit(report$Path, "_"), function(x){if (length(x)==0 | !isProperFid(x[length(x)])){return (NA)}; x[length(x)]})
		report <- report[-which(is.na(report$Fid)),]
		candidatesMerged$BAB <- report$Internal.Processing.Sample.ID[match(candidatesMerged$FID, report$Fid)]
		candidatesMerged$project <- report$Metaproject[match(candidatesMerged$FID, report$Fid)]
		candidatesMerged$BAB[is.na(candidatesMerged$BAB)] <- candidatesMerged$FID[is.na(candidatesMerged$BAB)]
	}
	candidatesMerged
}


##------------------------------------------------------------------------------
##' Finds overlap between calls and AOH regions
##' 
##' 
##' @param candidatesMergedAnnotated		object returned by annotateCandidates
##' @param extAOH							object returned by prepareAOHData
##' @param aohSize							min AOH size 
##' @param minAOHsig						AOH signal threshold
##' @param mc.cores							number of cores (see mclapply2)
##------------------------------------------------------------------------------

annotateAOH <- function(candidatesMergedAnnotated, extAOH, aohSize, minAOHsig, mc.cores){
	if (is.null(extAOH)){
		candidatesMergedAnnotated[,paste("inAOH", "_", format(aohSize, scientific=F), sep="")] <- TRUE
		return (candidatesMergedAnnotated)
	}
	tmpRes <- mclapply2((1:length(unique(candidatesMergedAnnotated$FID))), function(i){ fid <- unique(candidatesMergedAnnotated$FID)[i]
				tmpCand <- candidatesMergedAnnotated[candidatesMergedAnnotated$FID == fid,]
				tmpCand [,paste("inAOH", "_", format(aohSize, scientific=F), sep="")] <- FALSE
				starts <- tmpCand$Start
				stops <- tmpCand$Stop
				stops [which(starts > stops)] <- starts[which(starts > stops)]
				finalCandGR <- GRanges(tmpCand$Chr, IRanges(starts, stops))
				selAOH_tmp <- extAOH[extAOH$Name == fid & extAOH$Length > aohSize & extAOH$seg.mean > minAOHsig,]
				selAOH_tmp_gr <- GRanges(selAOH_tmp$chrom, IRanges(selAOH_tmp$loc.start.ext, selAOH_tmp$loc.end.ext))
				mm_tmp <- as.matrix(findOverlaps(finalCandGR, selAOH_tmp_gr))
				tmpCand[unique(mm_tmp[,1]), paste("inAOH", "_", format(aohSize, scientific=F), sep="")] <- TRUE
				tmpCand		
			},mc.cores=mc.cores)
	res <- rbindlist(tmpRes)
	res
	
}


##------------------------------------------------------------------------------
##' Adds the information about overlapping calls (OverlapCnt column) to the list of candidates
##' 
##' 
##' @param filteredCalls   object returned by annotateAOH
##------------------------------------------------------------------------------
annotateFreq <- function(filteredCalls)
{
	candidates_ir <- IRanges(filteredCalls$Start_idx , filteredCalls$Start_idx + filteredCalls$Exon_num - 1)
	mm <- (table(as.matrix(findOverlaps(candidates_ir, candidates_ir))[,1]))
	filteredCalls$OverlapCnt <- 0 
	filteredCalls$OverlapCnt[as.numeric(names(mm))] <- mm - 1
	filteredCalls 
}

##------------------------------------------------------------------------------
##' Calculates average z-RPKM for each call	
##------------------------------------------------------------------------------
calculateZscores <- function(filteredCalls,rpkmDtOrdered){
	mm <- as.matrix(rpkmDtOrdered)
	zScores <- sapply(1:nrow(filteredCalls), function(i){
				#print(i)
				fid <- filteredCalls$FID[i]
				sampleIdx <- which(rownames(rpkmDtOrdered) == fid)
				position <- filteredCalls$Start_idx[i]
				
				mean(sapply(1:filteredCalls$Exon_num[i], function(j){
									#print (paste0("j=", j))
									val <- mm[sampleIdx, position + j -1 ]
									vv <- sort(mm[-sampleIdx, position + j -1])
									x<- c(val,vv)
									(x[1] - mean(x))/sd(x)			
								}))
				
			})
	zScores
}

##------------------------------------------------------------------------------
##' Adds legend to the plot (used internally by plotDeletion)
##------------------------------------------------------------------------------
addLegend <- function(cand, trBlack, mainText="")
{
	par(xpd=NA)
	tmp <- cnvrt.coords(.7,0.9, 'tdev')$usr
	mainText <- paste(cand$BAB, " - ","chr",cand$Chr, ":" , cand$Start,"-", cand$Stop, " (",cand$Exon_num, " exons)" ,mainText, sep="")
	title(mainText, outer=TRUE)
	legend(tmp,
			c("sample with deletion","lower and upper thresholds", "other samples", "samples with no AOH overlap" ),
			col=c("red", "blue", trBlack, "dimgrey"),
			lwd=c(3,3,1,3), 
			lty=c(1,2,1,3), cex=0.5)
	
}
##------------------------------------------------------------------------------
##' Plots read depth for selected candidate of HMZ deletion
##' 
##' 
##' @param calls			calls object returned by runHMZDelFinder
##' @param i				index of the deletion to plot
##' @param bedOrdered		bedOrdered object returned by runHMZDelFinder
##' @param rpkmDtOrdered	rpkmDtOrdered object returned by runHMZDelFinder
##' @param outputDir		output directory
##------------------------------------------------------------------------------
plotDeletion <- function(calls, i, bedOrdered, rpkmDtOrdered,  outputDir, mainText=""  ){
	library(Hmisc)
	cand <- calls[i,]
	alpha=0.5
	window <-5
	idx <- max(1,(cand$Start_idx-window)):min(ncol(rpkmDtOrdered),(cand$Start_idx + cand$Exon_num + window - 1))
	replaceInf <- function(x) {
		dm <- data.matrix(x)
		dm[!is.finite(dm)] <- 0
		res <- data.table(dm)
		rownames(res) <- rownames(x)
		res
	}
	trBlack <- rgb(0,0,0,alpha=0.3) 
	calls_ir <- IRanges(calls$Start_idx , calls$Start_idx + calls$Exon_num - 1)
	cand_ir <- IRanges(cand$Start_idx , cand$Start_idx + cand$Exon_num - 1)
	other_idx <- as.matrix(findOverlaps(cand_ir, calls_ir))[,2]
	noAOH <- vector()
	if (length(other_idx)>0){
		noAOH <- calls$FID[other_idx][which(!calls$inAOH_1000[other_idx])]
	}
	ll <- replaceInf(log(rpkmDtOrdered[,idx,with=F ] + 1, 10) )
	ll2<- ll[-which(rownames(rpkmDtOrdered)%in% union(union(calls$FID[calls$PoorSample], cand$FID), noAOH)),] # removing poor quality samples
	ll3 <- ll[which(rownames(rpkmDtOrdered)%in% noAOH),] # samples not overlapping with AOH are shown as dashed line
	png(paste(outputDir,cand$BAB,"_chr",cand$Chr, "_",cand$Start,"_", cand$Stop ,"_",gsub(",","_",cand$Genes),  ".png",sep=""), width=1000,height=1000, pointsize=25)
	maxV <- max(t(ll))
	vspace <- 0.2* maxV
	par(  oma=c(2,2,3,2))
	plot(c(min(idx), max(idx) ), c(-vspace, 1.05*maxV),type="n",xlim=c(min(idx), max(idx) ),xlab="Probe (exon) number",ylab="log(RPKM + 1)")
	matplot(matrix(rep(min(idx):max(idx), nrow(ll2)), ncol=nrow(ll2)), t(ll2), type="l",lty=1, col=trBlack, add=T)# ,
	if (nrow(ll3) > 0){
		matplot(matrix(rep(min(idx):max(idx), nrow(ll3)), ncol=nrow(ll3)), t(ll3), type="l",lty=3, col="dimgrey", add=T,  lwd=3)
	}
	lines(min(idx):max(idx), ll[which(rownames(rpkmDtOrdered)== cand$FID),], col="red",lwd=3, alpha=alpha)
	matplot(t(matrix(rep(idx,2),ncol=2)), t(cbind(rep(0, length(idx)),as.numeric(ll[which(rownames(rpkmDtOrdered)== cand$FID),]))),add=T,type="l", col="red",lwd=3, lty=1)
	abline(h=log(1.5, 10), lwd=3, lty=2, col="blue")
	abline(h=log(2, 10), lwd=3, lty=2, col="blue")
	genes <- unique(bedOrdered$V4[idx])
	i <- 1
	## plotting genes and exons
	for (gene in genes){
		geneIdx <- which(gene == bedOrdered$V4[idx])
		voffset <- -vspace * 0.15 *i
		rect(idx[geneIdx]-0.2,rep(voffset + 0.03*vspace, length(idx[geneIdx])), idx[geneIdx] + 0.2, rep(voffset - 0.03*vspace, length(idx[geneIdx])), col="darkblue", border="darkblue")
		lines(c(min(idx[geneIdx]),max(idx[geneIdx])), c(voffset, voffset), col="darkblue", lwd=3)		
		text(mean(idx[geneIdx]), voffset-0.16*vspace, gene, font=3, cex=0.7)
		i <- i +1
		if (i==5){i =1}
	}
	voffset <- -vspace * 0.9
	rect( cand$Start_idx  ,voffset + 0.05*vspace, cand$Start_idx  + cand$Exon_num - 1, voffset - 0.05*vspace, col="red", border="darkblue")
	abline(v=cand$Start_idx , lty=2)
	abline(v=cand$Start_idx + cand$Exon_num - 1  , lty=2)
	addLegend(cand, trBlack, mainText)
	dev.off()
}




printBanner <- function()
{
	cat("\n",
			" _   _ __  __ _________       _ _____ _           _\n",           
			"| | | |  \\/  |__  /  _ \\  ___| |  ___(_)_ __   __| | ___ _ __\n", 
			"| |_| | |\\/| | / /| | | |/ _ \\ | |_  | | '_ \\ / _` |/ _ \\ '__|\n",
			"|  _  | |  | |/ /_| |_| |  __/ |  _| | | | | | (_| |  __/ |\n",
			"|_| |_|_|  |_/____|____/ \\___|_|_|   |_|_| |_|\\__,_|\\___|_|\n\n")  
}	


##------------------------------------------------------------------------------
##' Main function: return HMZ deletion calls 
##' 
##' 
##' @param snpPaths				list of paths to VCFs with SNP data; if NULL algorithm skips AOH analysis
##' @param snpFids				list of sample identifiers (in the same order as in snpFileNames); if NULL algorithm skips AOH analysis
##' @param rpkmPath				list of paths to RPKM data
##' @param rpkmFids				list of sample identifiers (in the same order as in rpkmPaths)
##' @param mc.cores				number of cores 
##' @param aohRDataOut			temporary file that store AOH data
##' @param bedFile				target design file
##' @param lowRPKMthreshold		RPKM threshold used in the algorithm (default=0.5)
##' @param minAOHsize			minimal size of AOH region (default=1000)
##' @param minAOHsig			threshold for calling AOH (default=0.45)
##' @param is_cmg				CMG specific annotations (default=FALSE)
##------------------------------------------------------------------------------

runHMZDelFinder <- function(snpPaths, snpFids,
		rpkmPaths, rpkmFids,
		mc.cores, aohRDataOut,
		bedFile, lowRPKMthreshold,
		minAOHsize, minAOHsig, is_cmg)
{
	
	
	## checking input parameters
	if (!is.null(snpPaths) && any(!file.exists(snpPaths))){print("[ERROR]: One or more paths to VCF file does not exist."); return (NULL)}
	if ( any(!file.exists(snpPaths))){print("[ERROR]: One or more paths to RPKM file does not exist."); return (NULL)}
	if (!file.exists(bedFile)){print("[ERROR]: BED file does not exist.");return (NULL)}
	if (length(rpkmPaths) != length(rpkmFids)){print("[ERROR]: Number of rpkmPaths differ from the number rpkmFids.");return (NULL)}
	if (!is.null(snpPaths) && length(snpPaths) != length(snpFids)){print("[ERROR]: Number of vcfPaths differ from the number vcfFids.");return (NULL)}
	
	library(gdata)
	library(data.table)
	library(GenomicRanges)
	library(parallel)
	library(matrixStats)
	
	printBanner()
	print("[step 1 out of 7] ******  AOH data ******")
	
	extAOH <- NULL
	if (is.null(snpPaths) || is.null(snpFilds)){
		print("Skipping AOH analysis ...")
		
	}
	else{
		#if (file.exists(aohRDataOut)){load(aohRDataOut)}else
		#{
			print("[step 1.5 out of 7] ****** Preparing AOH data ******")
			library(DNAcopy)
			extAOH <- prepareAOHData(snpPaths, snpFids, mc.cores)
			save(extAOH, file=aohRDataOut)
		#}
		
	}
		
	print("[step 2 out of 7] ****** Preparing RPKM data ******")
	selectedFidsIdx <- 1:length(rpkmFids)
	if (!is.null(extAOH)) {
		selectedFidsIdx <- which(rpkmFids %in% extAOH$Name)
	}
	rpkmDt <- prepareRPKMData(rpkmPaths[selectedFidsIdx], rpkmFids[selectedFidsIdx], 1)
	tmp <- reorderBedAndRpkmDt(bedFile, rpkmDt)
	bedOrdered <- tmp$bedOrdered
	rpkmDtOrdered <- tmp$rpkmDtOrdered
	
	print("[step 3 out of 7] ****** SELECTING CANDIDATE EXONS ******")
	exonMedianRpkms <- colMedians(as.matrix(rpkmDtOrdered[,,with=F]))
	exonsToExclude <- which (exonMedianRpkms < 7)
	gc()
	processRPKMResults <- processRPKM(rpkmDtOrdered, bedOrdered, mc.cores,lowRPKMthreshold,exonsToExclude,maxFrequency=maxFrequency)
	selectedExonsFinal <- processRPKMResults[["selectedExonsFinal"]]
	
	
	print("[step 4 out of 7] ****** DELETION CALLING ******")
	candidateExonCalls <- getCandidateExonCalls (rpkmDtOrdered, bedOrdered, selectedExonsFinal, mc.cores,lowRPKMthreshold)
	
	print("[step 5 out of 7] ****** MERGING CALLS FROM CONSECUTIVE EXONS ******")
	candidatesMerged <- mergeCandidates(candidateExonCalls, bedOrdered)
	
	print("[step 6 out of 7] ****** ANNOTATING CALLS ******")
	candidatesMergedAnnotated <- annotateCandidates (candidatesMerged, is_cmg) 
	poorSamplesNames <- rownames(rpkmDtOrdered)[processRPKMResults$toRemSamples]
	finalCandToRemoveIdx <- which(candidatesMergedAnnotated$FID %in% poorSamplesNames)
	candidatesMergedAnnotated$PoorSample <- FALSE
	candidatesMergedAnnotated$PoorSample[finalCandToRemoveIdx] <- TRUE
	candidatesMergedAnnotated$posKey <- paste(candidatesMergedAnnotated$Chr, ":", candidatesMergedAnnotated$Start, "_", candidatesMergedAnnotated$Stop, sep="")
	candidatesMergedAnnotated[,key:=paste(FID,"_",posKey,sep="")]
	
	print("[step 7 out of 7] ****** OVERLAPPING WITH AOH REGIONS AND FILTERING ******")
	allCalls <- annotateAOH(candidatesMergedAnnotated, extAOH, minAOHsize, minAOHsig, mc.cores)
	allCalls$ZScore <- calculateZscores(allCalls,rpkmDtOrdered)
	gc();gc();
	allCalls <- allCalls [order(allCalls$Length, decreasing=T),]
	
	## filtering out calls from low quality samples and calls that do not overlap with any AOH region
	filteredCalls <- allCalls[which(allCalls[,paste("inAOH_",format(minAOHsize, scientific=F),sep=""),with=F] & !allCalls$PoorSample & allCalls$Length>50),]
	filteredCalls <- annotateFreq(filteredCalls)
	filteredCalls[,PerSampleNr:=.N,by=BAB]
	filteredCalls$ZScore <- calculateZscores(filteredCalls,rpkmDtOrdered)
	
	results <- list (filteredCalls = filteredCalls, allCalls=allCalls, bedOrdered = bedOrdered, rpkmDtOrdered = rpkmDtOrdered)
	
	return(results)
}

