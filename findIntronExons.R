source("~/scripts/R/dna.R")

library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library('GenomicRanges')

exon<-exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,use.names=TRUE)
intron<-intronsByTranscript(TxDb.Hsapiens.UCSC.hg38.knownGene,use.names=TRUE)
onlyExon<-setdiff(stack(exon),stack(intron),ignore.strand=FALSE)
onlyIntron<-setdiff(stack(intron),stack(exon),ignore.strand=FALSE)

transcriptsByOverlapsByRegion<-function (x, ranges, maxgap = 0L, minoverlap = 1L, type = c("any", "start", "end"), ...){
	trans<-transcripts(x, columns = c("tx_id", "tx_name"))
	overlaps<-findOverlaps(ranges, trans, maxgap = maxgap, minoverlap = minoverlap, type = match.arg(type))
	out<-tapply(subjectHits(overlaps),queryHits(overlaps),function(x)trans[x])
	return(out)
}

overlappingIntronGenes<-cacheOperation('work/overlappingIntronGenes.Rdat',transcriptsByOverlapsByRegion,TxDb.Hsapiens.UCSC.hg38.knownGene,onlyIntron,EXCLUDE='TxDb.Hsapiens.UCSC.hg38.knownGene')
overlappingExonGenes<-cacheOperation('work/overlappingExonGenes.Rdat',transcriptsByOverlapsByRegion,TxDb.Hsapiens.UCSC.hg38.knownGene,onlyExon,EXCLUDE='TxDb.Hsapiens.UCSC.hg38.knownGene')

intronRegions<-do.call(rbind,cacheOperation('work/intronRegions.Rdat',mclapply,overlappingIntronGenes,function(x){if(runif(1)<.001)cat('.');Reduce(union,x)},mc.cores=10,EXCLUDE='overlappingIntronGenes'))
exonRegions<-do.call(rbind,cacheOperation('work/exonRegions.Rdat',mclapply,overlappingExonGenes,function(x){if(runif(1)<.001)cat('.');Reduce(union,x)},mc.cores=10,EXCLUDE='overlappingIntronGenes'))

if(FALSE){
	exonOut<-data.frame('chr'=exonRanges$chr,'start'=exonRanges$start-1,'end'=exonRanges$end,'name'=sprintf('exon%09d',1:nrow(exonRanges)),'dummy'=0,'strand'='+',stringsAsFactors=FALSE)
	intronOut<-data.frame('chr'=intronRanges$chr,'start'=intronRanges$start-1,'end'=intronRanges$end,'name'=intronRanges$id,'dummy'=0,'strand'='+',stringsAsFactors=FALSE)
	intronGenesOut<-data.frame('chr'=genesForIntrons$chr,'start'=genesForIntrons$start-1,'end'=genesForIntrons$end,'name'=genesForIntrons$id,'dummy'=0,'strand'='+',stringsAsFactors=FALSE)
	exonGenesOut<-data.frame('chr'=genesForExons$chr,'start'=genesForExons$start-1,'end'=genesForExons$end,'name'=genesForExons$id,'dummy'=0,'strand'='+',stringsAsFactors=FALSE)

	options(scipen=20)
	write.table(exonOut,'work/exonOnly_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
	write.table(intronOut,'work/intronOnly_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
	write.table(intronGenesOut,'work/intronGenes_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
	write.table(exonGenesOut,'work/exonGenes_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')

	for(ii in list.files('work','regions.bed')){
		#cmd<-sprintf("%s -t 12 -b work/%s %s >work/%s","/media/RAID/shescott/splicing/oldBedCount",ii,paste(bamFiles,collapse=' '),sub('bed$','count2',ii))
		outFile<-sprintf('work/%s',sub('bed$','count',ii))
		bamFileOrder<-c()
		if(!file.exists(outFile)){
			thisOuts<-c()
			for(study in unique(runInfo$study)){
				studyOutFile<-sub('\\.bed$',sprintf('_%s.count',study),ii)
				#Katze is single end so pass -s option
				thisBams<-runInfo[runInfo$study==study,'bamFile']
				cmd<-sprintf("%s -t 12 -b work/%s %s %s >%s",bedCountBin,ii,paste(thisBams,collapse=' '),ifelse(study=='Katze','-s',''),studyOutFile)
				message(cmd)
				system(cmd)
				thisOuts<-c(thisOuts,studyOutFile)
				bamFileOrder<-c(bamFileOrder,thisBams)
			}
			thisCounts<-lapply(thisOuts,read.table,stringsAsFactors=FALSE)
			if(any(sapply(thisCounts,nrow)!=nrow(thisCounts[[1]])))stop(simpleError('Missing counts in study specific files'))
			if(any(apply(do.call(cbind,lapply(thisCounts,function(x)x[,1])),1,function(x)any(x!=x[1]))))stop(simpleError('Mismatched regions in study specific files'))
			if(any(bamFileOrder!=bamFiles))stop(simpleError('Order of bam files changed. Program this more robustly'))
			outCounts<-cbind(thisCounts[[1]][,1:2],do.call(cbind,lapply(thisCounts,function(x)x[,-1:-2])))
			write.table(outCounts,outFile)
		}else{
			message(outFile,' already exists')
		}
	}


}

