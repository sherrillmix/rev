source("~/scripts/R/dna.R")

library("TxDb.Hsapiens.UCSC.hg38.knownGene")
library('GenomicRanges')
bedCountBin<-'~/installs/bedCount/bedCount'

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

reduceFunc<-function(x){if(runif(1)<.001)cat('.');out<-Reduce(union,x);elementMetadata(out)<-NULL;return(out)}
intronRegions<-cacheOperation('work/intronRegions.Rdat',mclapply,overlappingIntronGenes,reduceFunc,mc.cores=4,EXCLUDE='overlappingIntronGenes')
if(any(sapply(intronRegions,length)>1))stop(simpleError('Multiple overlap intron region found'))
intronRegions<-do.call(c,unname(intronRegions))
exonRegions<-cacheOperation('work/exonRegions.Rdat',mclapply,overlappingExonGenes,reduceFunc,mc.cores=4,EXCLUDE='overlappingExonGenes')
if(any(sapply(exonRegions,length)>1))stop(simpleError('Multiple overlap exon region found'))
exonRegions<-do.call(c,unname(exonRegions))

exonOut<-data.frame('chr'=seqnames(onlyExon),'start'=start(onlyExon)-1,'end'=end(onlyExon),'name'=sprintf('exon%09d',1:length(onlyExon)),'strand'=strand(onlyExon),stringsAsFactors=FALSE)
intronOut<-data.frame('chr'=seqnames(onlyIntron),'start'=start(onlyIntron)-1,'end'=end(onlyIntron),'name'=sprintf('intron%09d',1:length(onlyIntron)),'strand'=strand(onlyIntron),stringsAsFactors=FALSE)

exonGenesOut<-data.frame('chr'=seqnames(exonRegions),'start'=start(exonRegions)-1,'end'=end(exonRegions),'name'=sprintf('exon%09d',1:length(exonRegions)),'strand'=strand(exonRegions),stringsAsFactors=FALSE)
intronGenesOut<-data.frame('chr'=seqnames(intronRegions),'start'=start(intronRegions)-1,'end'=end(intronRegions),'name'=sprintf('intron%09d',1:length(intronRegions)),'strand'=strand(intronRegions),stringsAsFactors=FALSE)

options(scipen=20)
write.table(exonOut,'work/exonOnly_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
write.table(intronOut,'work/intronOnly_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
write.table(intronGenesOut,'work/intronGenes_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')
write.table(exonGenesOut,'work/exonGenes_regions.bed',col.names=FALSE,row.names=FALSE,quote=FALSE,sep='\t')




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

readCounts<-lapply(list.files('work','regions.count'), function(x)read.table(file.path('work',x)))
names(readCounts)<-sub('(Only)?_regions.count','',list.files('work','regions.count'))

intronExon<-mapply(function(exon,gene){names(exon)<-c('reg_reg','reg_cut',sprintf('reg_%s',names(bamFiles)));names(gene)<-c('gene_reg','gene_cut',sprintf('gene_%s',names(bamFiles)));return(cbind(exon,gene))},readCounts[c('exon','intron')],readCounts[c('exonGenes','intronGenes')],SIMPLIFY=FALSE)
names(intronExon)<-c('exon','intron')

bayesBinomDiff<-function(y1,y2,n1,n2,priors1=c(1,1),priors2=priors1,reps=100000){
	betas1<-rbeta(reps,sum(y1)+priors1[1],sum(n1-y1)+priors1[2])
	betas2<-rbeta(reps,sum(y2)+priors2[1],sum(n2-y2)+priors2[2])
	logitDiff<-log(betas1)-log(1-betas1)-log(betas2)+log(1-betas2)
	return(logitDiff)
}

#could automate this
colNames<-list("regionRev"=c('reg_rev_1','reg_rev_2'),
	"geneRev"=c('gene_rev_1','gene_rev_2'),
	"regionControl"=c('reg_control_1','reg_control_2'),
	"geneControl"=c('gene_control_1','gene_control_2')
)

diffs<-cacheOperation('work/diffs.Rdat',mclapply,intronExon,function(xx,colNames){
	apply(xx[,unlist(colNames)],1,function(y){
		quantile(bayesBinomDiff(y[colNames[["regionControl"]]],y[colNames[["regionRev"]]],y[colNames[["geneControl"]]],y[colNames[["geneRev"]]],reps=10000),c(.005,.025,.5,.975,.995))
	})
},colNames,mc.cores=2)
diffs1<-cacheOperation('work/diffs1.Rdat',mclapply,intronExon,function(xx,colNames){
	apply(xx[,unlist(colNames)],1,function(y){
		quantile(bayesBinomDiff(y[colNames[["regionControl"]][1]],y[colNames[["regionRev"]][1]],y[colNames[["geneControl"]][1]],y[colNames[["geneRev"]][1]],reps=10000),c(.005,.025,.5,.975,.995))
	})
},colNames,mc.cores=2)
diffs2<-cacheOperation('work/diffs2.Rdat',mclapply,intronExon,function(xx,colNames){
	apply(xx[,unlist(colNames)],1,function(y){
		quantile(bayesBinomDiff(y[colNames[["regionControl"]][2]],y[colNames[["regionRev"]][2]],y[colNames[["geneControl"]][2]],y[colNames[["geneRev"]][2]],reps=10000),c(.005,.025,.5,.975,.995))
	})
},colNames,mc.cores=2)

conservBound<-lapply(diffs,function(x)apply(x[c('2.5%','97.5%'),],2,conservativeBoundary))
counts<-lapply(intronExon,function(x)apply(x[unlist(colNames)],1,sum))
mapply(function(x,y)summary(x[y>20]),conservBound,counts)

