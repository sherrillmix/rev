readCounts<-lapply(list.files('work','regions.count'), function(x)read.table(file.path('work',x)))
names(readCounts)<-sub('(Only)?_regions.count','',list.files('work','regions.count'))

intronExon<-mapply(function(exon,gene){names(exon)<-c('reg_reg','reg_cut',sprintf('reg_%s',names(bamFiles)));names(gene)<-c('gene_reg','gene_cut',sprintf('gene_%s',names(bamFiles)));return(cbind(exon,gene))},readCounts[c('exon','intron','alt')],readCounts[c('exonGenes','intronGenes','altGenes')],SIMPLIFY=FALSE)
names(intronExon)<-c('exon','intron','alt')

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

enoughReads<-lapply(intronExon,function(x)apply(x[,unlist(colNames[c('regionRev','regionControl')])],1,sum)>20)


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
    plot(conservBound1[[ii]][enoughReads[[ii]]],conservBound2[[ii]][enoughReads[[ii]]],main=niceNames[ii],ylab='Replicate 2 rev/control expression',xlab='Replicate 1 rev/control expression',bg='#FF000033',pch=21,cex=.4,col=NA)
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
  xlim<-range(c(conservBound[['intron']],conservBound[['exon']]),na.rm=TRUE)
  xlim<-c(-max(abs(xlim)),max(abs(xlim)))
  breaks<-seq(xlim[1]-1e-9,xlim[2]+1e-9,length.out=100)
  xlim<-c(-7,7) #MAGIC NUMBER
  prettyX<-pretty(xlim)
  xLabels<-sapply(prettyX,function(x)as.expression(bquote(2^.(x))))
  xLabels[prettyX==0]<-1
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
  mtext('Intron frequency (x1000)',2,1.6,xpd=NA,cex=.7,at=mean(par('usr')[3:4])*.7)
  exonHistCounts<-hist(conservBound[['exon']][enoughReads[['exon']]],breaks=breaks,plot=FALSE)$counts
  exonHistCounts<-hist(conservBound[['exon']][enoughReads[['exon']]],las=1,xlab='',ylab='',main='',xlim=xlim,ylim=c(0,max(exonHistCounts)*1.02),breaks=breaks,xaxt='n',yaxs='i',yaxt='n',mgp=c(0,.7,0),xpd=NA)
  prettyY<-pretty(exonHistCounts$counts)
  prettyY<-prettyY[prettyY<par('usr')[4]]
  axis(2,prettyY,prettyY/1000,las=1,mgp=c(0,.5,0),tcl=-.2)
  axis(1,prettyX,xLabels,mgp=c(2.7,.4,0),tcl=-.2)
  abline(v=median(conservBound[['exon']][enoughReads[['exon']]]),col='#FF000099',lwd=2)
  mtext('Ratio of infected/uninfected expression',1,1.5,cex=.7)
  mtext('Exon frequency (x1000)',2,1.6,xpd=NA,cex=.7,at=mean(par('usr')[3:4])*.8,adj=.5)
dev.off()

mostInterestingRegions<-mapply(function(readCount,bound,enoughRead)readCount[enoughRead&abs(bound)>4&grepl('chr[0-9XY][0-9]?:',readCount$gene_reg),],intronExon,conservBound,enoughReads,SIMPLIFY=FALSE)

pdf('out/mostDifferent.pdf')
  apply(do.call(rbind,mostInterestingRegions)[,c('gene_reg','reg_reg')],1,function(x){
    plotRegion(x[1],bamFiles,bam2depthBinary='~/installs/bedCount/bam2depth',normalize=FALSE,ylab='Read coverage')
    abline(v=parseRegion(x[2])[,c('start','end')],lty=2)
  })
dev.off()


