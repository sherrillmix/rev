source("~/scripts/R/dna.R")

library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library('GenomicRanges')
exon<-exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,use.names=TRUE)
intron<-intronsByTranscript(TxDb.Hsapiens.UCSC.hg38.knownGene,use.names=TRUE)
onlyExon<-setdiff(stack(exon),stack(intron),ignore.strand=TRUE)
onlyIntron<-setdiff(stack(intron),stack(exon),ignore.strand=TRUE)

transcriptsByOverlapsByRegion<-function (x, ranges, maxgap = 0L, minoverlap = 1L, type = c("any", "start", "end"), ...){
	trans<-transcripts(x, columns = c("tx_id", "tx_name"))
	overlaps<-findOverlaps(trans, ranges, maxgap = maxgap, minoverlap = minoverlap, type = match.arg(type))
	out<-tapply(subjectHits(overlaps)[1:1000],queryHits(overlaps)[1:1000],function(x)trans[x])
	return(out)
}

overlappingIntronGenes<-cacheOperation('work/overlappingIntronGenes.Rdat',transcriptsByOverlapsByRegion,TxDb.Hsapiens.UCSC.hg38.knownGene,onlyIntron,EXCLUDE='TxDb.Hsapiens.UCSC.hg38.knownGene')
overlappingExonGenes<-cacheOperation('work/overlappingExonGenes.Rdat',transcriptsByOverlapsByRegion,TxDb.Hsapiens.UCSC.hg38.knownGene,onlyExon,EXCLUDE='TxDb.Hsapiens.UCSC.hg38.knownGene')

