source("~/scripts/R/dna.R")
source('functions.R')
library(parallel)

bamDir<-'work/align'
bamFiles<-list.files(bamDir,'\\.bam$',full.name=TRUE)
names(bamFiles)<-sub('\\.bam','',basename(bamFiles))

readCounts<-lapply(list.files('work','regions.count'), function(x)read.table(file.path('work',x)),stringsAsFactors=FALSE)
names(readCounts)<-sub('(Only)?_regions.count','',list.files('work','regions.count'))

beds<-lapply(list.files('work','regions.bed'), function(x)read.table(file.path('work',x)),stringsAsFactors=FALSE)
names(beds)<-sub('(Only)?_regions.bed','',list.files('work','bed'))

intronExon<-mapply(function(exon,gene,exonBed,geneBed){
	names(exon)<-c('region_reg','region_cut',sprintf('region_%s',names(bamFiles)));
	names(gene)<-c('gene_reg','gene_cut',sprintf('gene_%s',names(bamFiles)));
	colnames(exonBed)<-colnames(geneBed)<-c('chr','start','end','name','strand')
	if(any(pasteRegion(exonBed$chr,exonBed$start+1,exonBed$end)!=exon$region_reg))stop(simpleError('Mismatch in bed and counts in region'))
	if(any(pasteRegion(geneBed$chr,geneBed$start+1,geneBed$end)!=gene$region_reg))stop(simpleError('Mismatch in bed and counts in genes'))
	if(any(geneBed$strand!=exonBed$strand))stop(simpleError('Strand mismatch'))
	return(cbind(exon,gene,'strand'=exonBed$strand,'name'=exonBed$name))
},readCounts[c('exon','intron','alt')],readCounts[c('exonGenes','intronGenes','altGenes')],beds[c('exon','intron','alt')],beds[c('exonGenes','intronGenes','altGenes')],SIMPLIFY=FALSE)
names(intronExon)<-c('exon','intron','alt')


#could automate this
colNames<-list("regionRev"=c('region_rev_1','region_rev_2'),
  "geneRev"=c('gene_rev_1','gene_rev_2'),
  "regionControl"=c('region_control_1','region_control_2'),
  "geneControl"=c('gene_control_1','gene_control_2')
)

diffs<-cacheOperation('work/diffs.Rdat',mclapply,lapply(intronExon,function(x)x[,unlist(colNames)]),function(xx,colNames){
  apply(xx[,unlist(colNames)],1,function(y){
	 diffs<-bayesBinomDiff(y[colNames[["regionControl"]]],y[colNames[["regionRev"]]],y[colNames[["geneControl"]]],y[colNames[["geneRev"]]],reps=10000)
    return(c(quantile(diffs,c(.005,.025,.5,.975,.995)),'p'=mean(diffs>0)))
  })
},colNames,mc.cores=8)
diffs1<-cacheOperation('work/diffs1.Rdat',mclapply,lapply(intronExon,function(x)x[,unlist(colNames)]),function(xx,colNames){
  apply(xx[,unlist(colNames)],1,function(y){
	 diffs<-bayesBinomDiff(y[colNames[["regionControl"]][1]],y[colNames[["regionRev"]][1]],y[colNames[["geneControl"]][1]],y[colNames[["geneRev"]][1]],reps=10000)
    return(c(quantile(diffs,c(.005,.025,.5,.975,.995)),'p'=mean(diffs>0)))
  })
},colNames,mc.cores=8)
diffs2<-cacheOperation('work/diffs2.Rdat',mclapply,lapply(intronExon,function(x)x[,unlist(colNames)]),function(xx,colNames){
  apply(xx[,unlist(colNames)],1,function(y){
    diffs<-bayesBinomDiff(y[colNames[["regionControl"]][2]],y[colNames[["regionRev"]][2]],y[colNames[["geneControl"]][2]],y[colNames[["geneRev"]][2]],reps=10000)
    return(c(quantile(diffs,c(.005,.025,.5,.975,.995)),'p'=mean(diffs>0)))
  })
},colNames,mc.cores=8)

enoughReads<-lapply(intronExon,function(x)apply(x[,unlist(colNames[c('regionRev','regionControl')])],1,sum)>100)


conservBound<-lapply(diffs,function(x)apply(x[c('2.5%','97.5%'),],2,conservativeBoundary))
conservBound1<-lapply(diffs1,function(x)apply(x[c('2.5%','97.5%'),],2,conservativeBoundary))
conservBound2<-lapply(diffs2,function(x)apply(x[c('2.5%','97.5%'),],2,conservativeBoundary))
counts<-lapply(intronExon,function(x)apply(x[unlist(colNames)],1,sum))
mapply(function(x,y)summary(x[y>20]),conservBound,counts)
median1<-lapply(diffs1,function(x)x['50%',])
median2<-lapply(diffs2,function(x)x['50%',])

niceNames<-c('exon'='Exclusively exonic regions','intron'='Exclusively intronic regions','alt'='Alternatively spliced regions')

for(ii in names(diffs1)){
  png(sprintf('out/reproduce_%s.png',ii),width=2000,height=2000,res=300)
    lims<-range(c(conservBound1[[ii]][enoughReads[[ii]]],conservBound2[[ii]][enoughReads[[ii]]]))
    plot(conservBound1[[ii]][enoughReads[[ii]]],conservBound2[[ii]][enoughReads[[ii]]],main=niceNames[ii],ylab='Replicate 2 rev/control expression',xlab='Replicate 1 rev/control expression',bg='#FF000022',pch=21,cex=.4,col=NA,xlim=lims,ylim=lims)
    abline(h=0,v=0,lty=2)
  dev.off()
}

for(ii in names(diffs1)){
  png(sprintf('out/median_%s.png',ii),width=2000,height=2000,res=300)
    plot(median1[[ii]][enoughReads[[ii]]],median2[[ii]][enoughReads[[ii]]],main=niceNames[ii],ylab='Replicate 2 rev/control expression',xlab='Replicate 1 rev/control expression',bg='#FF000033',pch=21,cex=.4,col=NA)
    abline(h=0,v=0,lty=2)
  dev.off()
}

for(ii in names(diffs1)){
  message(ii)
  print(table(ifelse(conservBound1[[ii]][enoughReads[[ii]]]==0,'equal',ifelse(conservBound1[[ii]][enoughReads[[ii]]]<0,'less','greater')),ifelse(conservBound2[[ii]][enoughReads[[ii]]]==0,'equal',ifelse(conservBound2[[ii]][enoughReads[[ii]]]<0,'less','greater'))))
}


pdf('out/intronExon.pdf')
  layout(matrix(c(1,2),nrow=2,byrow=FALSE),heights=c(4.25,5))
  par(mar=c(1,3.5,0.3,.3))
  xlim<-range(c(conservBound[['intron']],conservBound[['exon']]),na.rm=TRUE)
  xlim<-c(-max(abs(xlim)),max(abs(xlim)))
  breaks<-seq(xlim[1]-1e-9,xlim[2]+1e-9,length.out=100)
  xlim<-c(-3,3) #MAGIC NUMBER
  prettyX<-pretty(xlim)
  xLabels<-sapply(prettyX,function(x)as.expression(bquote(2^.(x))))
  xLabels[prettyX==0]<-1
  prettyX<-prettyX/log2(exp(1))
  xlim<-xlim/log2(exp(1))
  intronHistCounts<-hist(conservBound[['intron']][enoughReads[['intron']]],breaks=breaks,plot=FALSE)$counts
  intronHistCounts<-hist(conservBound[['intron']][enoughReads[['intron']]],las=1,xlab='',ylab='',main='',xlim=xlim,ylim=c(0,max(intronHistCounts)*1.02),breaks=breaks,xaxt='n',yaxs='i',yaxt='n',mgp=c(0,.7,0),xpd=NA)
  prettyY<-pretty(intronHistCounts$counts)
  prettyY<-prettyY[prettyY<par('usr')[4]]
  axis(2,prettyY,prettyY/1000,las=1,mgp=c(0,.5,0),tcl=-.2)
  axis(1,prettyX,rep('',length(prettyX)),mgp=c(2.7,.5,0),tcl=-.2)
  abline(v=median(conservBound[['intron']][enoughReads[['intron']]]),col='#FF000099',lwd=2)
  intronBottom<-grconvertY(par('usr')[1],to='device')
  maxIntron<-grconvertY(par('usr')[4],to='device')
  par(lheight=.9)
  mtext('Intron frequency (x1000)',2,1.8)#xpd=NA)#,at=mean(par('usr')[3:4])*.9)
  par(mar=c(2.9,3.5,0.1,.3))
  exonHistCounts<-hist(conservBound[['exon']][enoughReads[['exon']]],breaks=breaks,plot=FALSE)$counts
  exonHistCounts<-hist(conservBound[['exon']][enoughReads[['exon']]],las=1,xlab='',ylab='',main='',xlim=xlim,ylim=c(0,max(exonHistCounts)*1.02),breaks=breaks,xaxt='n',yaxs='i',yaxt='n',mgp=c(0,.7,0),xpd=NA)
  prettyY<-pretty(exonHistCounts$counts)
  prettyY<-prettyY[prettyY<par('usr')[4]]
  axis(2,prettyY,prettyY/1000,las=1,mgp=c(0,.5,0),tcl=-.2)
  axis(1,prettyX,xLabels,mgp=c(2.7,.4,0),tcl=-.2)
  abline(v=median(conservBound[['exon']][enoughReads[['exon']]]),col='#FF000099',lwd=2)
  mtext('Ratio of rev/control expression',1,1.5)
  mtext('Exon frequency (x1000)',2,1.8)#,xpd=NA,at=mean(par('usr')[3:4])*.9,adj=.5)
dev.off()

#output a csv
out<-mapply(function(coordCounts,diff,diff1,diff2,bound1,bound2,enoughRead){
	cbind(coordCounts,t(diff),'p1'=diff1['p',],'p2'=diff2['p',],'rep1bound'=bound1,'rep2bound'=bound2)[enoughRead,]
},intronExon,diffs,diffs1,diffs2,conservBound1,conservBound2,enoughReads,SIMPLIFY=FALSE)
out<-lapply(names(out),function(x){tmp<-out[[x]];tmp$type<-x;tmp})
names(out)<-names(intronExon)
out<-rbind(out[['exon']],out[['intron']])
write.csv(out,'out/summary.csv')

mostInterestingRegions<-mapply(function(readCount,bound1,bound2,enoughRead){
	readCount$bound1<-bound1
	readCount$bound2<-bound1
	readCount[enoughRead&abs(bound1)>log(2)&abs(bound2)>log(2)&grepl('chr[0-9XY][0-9]?:',readCount$gene_reg),]
},intronExon,conservBound1,conservBound2,enoughReads,SIMPLIFY=FALSE)

pdf('out/mostDifferent.pdf')
  apply(do.call(rbind,mostInterestingRegions)[,c('gene_reg','region_reg')],1,function(x){
    plotRegion(x[1],bamFiles,bam2depthBinary='~/installs/bedCount/bam2depth',normalize=FALSE,ylab='Read coverage')
    abline(v=parseRegion(x[2])[,c('start','end')],lty=2)
  })
dev.off()
