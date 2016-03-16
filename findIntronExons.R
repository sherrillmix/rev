source("~/scripts/R/dna.R")

library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library('GenomicRanges')
library('org.Hs.eg.db')
bedCountBin<-'~/installs/bedCount/bedCount'

exon<-exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene,use.names=TRUE)
intron<-intronsByTranscript(TxDb.Hsapiens.UCSC.hg38.knownGene,use.names=TRUE)
onlyExon<-setdiff(stack(exon),stack(intron),ignore.strand=FALSE)
onlyIntron<-setdiff(stack(intron),stack(exon),ignore.strand=FALSE)
onlyAlt<-intersect(stack(intron),stack(exon),ignore.strand=FALSE)

overlappingIntronGenes<-cacheOperation('work/overlappingIntronGenes.Rdat',transcriptsByOverlapsByRegion,TxDb.Hsapiens.UCSC.hg38.knownGene,onlyIntron,EXCLUDE='TxDb.Hsapiens.UCSC.hg38.knownGene')
overlappingExonGenes<-cacheOperation('work/overlappingExonGenes.Rdat',transcriptsByOverlapsByRegion,TxDb.Hsapiens.UCSC.hg38.knownGene,onlyExon,EXCLUDE='TxDb.Hsapiens.UCSC.hg38.knownGene')
overlappingAltGenes<-cacheOperation('work/overlappingAltGenes.Rdat',transcriptsByOverlapsByRegion,TxDb.Hsapiens.UCSC.hg38.knownGene,onlyAlt,EXCLUDE='TxDb.Hsapiens.UCSC.hg38.knownGene')

reduceFunc<-function(x){if(runif(1)<.001)cat('.');out<-Reduce(union,x);elementMetadata(out)<-NULL;return(out)}
intronRegions<-cacheOperation('work/intronRegions.Rdat',mclapply,overlappingIntronGenes,reduceFunc,mc.cores=4,EXCLUDE='overlappingIntronGenes')
if(any(sapply(intronRegions,length)>1))stop(simpleError('Multiple overlap intron region found'))
intronRegions<-do.call(c,unname(intronRegions))
exonRegions<-cacheOperation('work/exonRegions.Rdat',mclapply,overlappingExonGenes,reduceFunc,mc.cores=4,EXCLUDE='overlappingExonGenes')
if(any(sapply(exonRegions,length)>1))stop(simpleError('Multiple overlap exon region found'))
exonRegions<-do.call(c,unname(exonRegions))
altRegions<-cacheOperation('work/altRegions.Rdat',mclapply,overlappingAltGenes,reduceFunc,mc.cores=4,EXCLUDE='overlappingAltGenes')
if(any(sapply(altRegions,length)>1))stop(simpleError('Multiple overlap alt region found'))
altRegions<-do.call(c,unname(altRegions))

exonOut<-data.frame('chr'=seqnames(onlyExon),'start'=start(onlyExon)-1,'end'=end(onlyExon),'name'=sprintf('exon%09d',1:length(onlyExon)),'strand'=strand(onlyExon),stringsAsFactors=FALSE)
intronOut<-data.frame('chr'=seqnames(onlyIntron),'start'=start(onlyIntron)-1,'end'=end(onlyIntron),'name'=sprintf('intron%09d',1:length(onlyIntron)),'strand'=strand(onlyIntron),stringsAsFactors=FALSE)
altOut<-data.frame('chr'=seqnames(onlyAlt),'start'=start(onlyAlt)-1,'end'=end(onlyAlt),'name'=sprintf('intron%09d',1:length(onlyAlt)),'strand'=strand(onlyAlt),stringsAsFactors=FALSE)

exonGenesOut<-data.frame('chr'=seqnames(exonRegions),'start'=start(exonRegions)-1,'end'=end(exonRegions),'name'=sprintf('exon%09d',1:length(exonRegions)),'strand'=strand(exonRegions),stringsAsFactors=FALSE)
intronGenesOut<-data.frame('chr'=seqnames(intronRegions),'start'=start(intronRegions)-1,'end'=end(intronRegions),'name'=sprintf('intron%09d',1:length(intronRegions)),'strand'=strand(intronRegions),stringsAsFactors=FALSE)
altGenesOut<-data.frame('chr'=seqnames(altRegions),'start'=start(altRegions)-1,'end'=end(altRegions),'name'=sprintf('alt%09d',1:length(altRegions)),'strand'=strand(altRegions),stringsAsFactors=FALSE)

options(scipen=20)
write.table(exonOut,'work/exonOnly_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
write.table(intronOut,'work/intronOnly_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
write.table(altOut,'work/altOnly_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
write.table(intronGenesOut,'work/intronGenes_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
write.table(exonGenesOut,'work/exonGenes_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
write.table(altGenesOut,'work/altGenes_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')




bamDir<-'work/align'
bamFiles<-list.files(bamDir,'\\.bam$',full.name=TRUE)
names(bamFiles)<-sub('\\.bam','',basename(bamFiles))

for(ii in list.files('work','regions.bed')){
	outFile<-sprintf('work/%s',sub('bed$','count',ii))
	if(!file.exists(outFile)){
		cmd<-sprintf("%s -t 12 -b work/%s %s -s >%s",bedCountBin,ii,paste(bamFiles,collapse=' '),outFile)
		message(cmd)
		system(cmd)
	}else{
		message(outFile,' already exists')
	}
}




