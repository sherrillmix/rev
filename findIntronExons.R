
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library('GenomicRanges')
library(parallel)

exon<-exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,use.names=TRUE)
intron<-intronsByTranscript(TxDb.Hsapiens.UCSC.hg38.knownGene,use.names=TRUE)
onlyExon<-setdiff(stack(exon),stack(intron),ignore.strand=TRUE)
onlyIntron<-setdiff(stack(intron),stack(exon),ignore.strand=TRUE)

transcriptByOverlapByRegion<-function (x, ranges, maxgap = 0L, minoverlap = 1L, type = c("any", "start", "end"), ...){
	.local <- function (x, ranges, maxgap = 0L, minoverlap = 1L, type = c("any", "start", "end"), columns = c("tx_id", "tx_name")){ 
		findOverlaps(transcripts(x, columns = columns), ranges, maxgap = maxgap, minoverlap = minoverlap, type = match.arg(type))
		browser()
	}
	.local(x, ranges, maxgap, minoverlap, type, ...)
}

overlappingIntronGenes<-lapply(split(onlyIntron,1:length(onlyIntron))[1:1000],function(intron){if(runif(1)<.001)cat('.');transcriptsByOverlaps(TxDb.Hsapiens.UCSC.hg38.knownGene,intron)},mc.cores=12)
overlappingExonGenes<-mclapply(split(onlyExon,1:length(onlyExon)),function(exon){if(runif(1)<.001)cat('.');transcriptsByOverlaps(TxDb.Hsapiens.UCSC.hg38.knownGene,exon)},mc.cores=12)


	overlappingGenes<-findOverlaps(onlyIntron,geneRanges)
	overlappingExonGenes<-findOverlaps(onlyExon,geneRanges)
#library("Homo.sapiens")
