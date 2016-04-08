library(GenomicFeatures)
library("TxDb.Hsapiens.UCSC.hg38.knownGene")

reg2range<-function(reg){
  chr<-strsplit(reg,':')[[1]]
  coords<-as.numeric(strsplit(sub('[-+*]$','',chr[2]),'-')[[1]])
  chr<-chr[[1]]
  GRanges(seqnames=chr,ranges=IRanges(coords[1],end=coords[2]))
}

bayesBinomDiff<-function(y1,y2,n1,n2,priors1=c(1,1),priors2=priors1,reps=100000){
  betas1<-rbeta(reps,sum(y1)+priors1[1],sum(n1-y1)+priors1[2])
  betas2<-rbeta(reps,sum(y2)+priors2[1],sum(n2-y2)+priors2[2])
  logitDiff<-log(betas1)-log(1-betas1)-log(betas2)+log(1-betas2)
  return(logitDiff)
}

transcriptsByOverlapsByRegion<-function (x, ranges, maxgap = 0L, minoverlap = 1L, type = c("any", "start", "end"), ...){
  trans<-transcripts(x, columns = c("tx_id", "tx_name"))
  overlaps<-findOverlaps(ranges, trans, maxgap = maxgap, minoverlap = minoverlap, type = match.arg(type))
  out<-tapply(subjectHits(overlaps),queryHits(overlaps),function(x)trans[x])
  return(out)
}



plotRegion<-function(reg,files,cols=rainbow(length(files)),normalize=TRUE,counts=NULL,justCounts=FALSE,bam2depthBinary='./bam2depth',logY=FALSE,coverMax=NULL,ylab=ifelse(justCounts,'Proportion of reads segments observed','Proportion of positions observed'),legendPos='top',...){
  thisGenes<-transcriptsByOverlaps(TxDb.Hsapiens.UCSC.hg38.knownGene,reg2range(reg),columns=c('tx_id','tx_name','gene_id'))
  thisGenes$symbol<-thisGenes$tx_name
  thisGenes$symbol[sapply(thisGenes$gene_id,length)>0]<-sapply(unlist(thisGenes$gene_id),function(x)tryCatch(select(org.Hs.eg.db,keys=x,columns='SYMBOL',keytype='ENTREZID')$SYMBOL,error=function(e)return(NA)))
  thisGenes[is.na(thisGenes$symbol),'symbol']<-thisGenes[is.na(thisGenes$symbol),'tx_name']
  if(length(thisGenes)==0){
    exons<-data.frame()[0,]
  }else{
    names(thisGenes)<-thisGenes$tx_name
    exons<-select(TxDb.Hsapiens.UCSC.hg38.knownGene,keys=thisGenes$tx_name,keytype='TXNAME',columns=c('TXNAME','EXONCHROM','EXONSTART','EXONEND','EXONSTRAND'))
    exons$symbol<-sapply(exons$TXNAME,function(x)thisGenes[x]$symbol)
  }
  if(is.null(names(files)))names(files)<-1:length(files)
  message(reg)
  par(mar=c(3.8,4.2,1.5,ifelse(is.null(counts)||justCounts,.2,4.2)),las=1)
  cover<-pullRegion(reg,files,bam2depthBinary)
  countCols<-2+(1:length(files))
  if(normalize)cover[,countCols]<-apply(cover[,countCols],2,function(x)if(any(x>0)) x/sum(x) else x)
  if(is.null(coverMax))coverMax<-max(cover[,countCols])
  nExon<-nrow(exons)
  nGene<-length(unique(exons$TXNAME))
  exonSpace<-.15*coverMax
  if(nGene>0){
    exonWidth<-exonSpace/nGene*.9
    if(logY){
      y<-10^seq(-1.9,-.1,length.out=nGene)
      names(y)<-unique(exons$TXNAME)
      exons$y<-y[exons$TXNAME]
      exons$top<-exons$y*10^(.5/nGene)
      exons$bottom<-exons$y*10^(-.5/nGene)
    } else {
      y<- -exonSpace/nGene*((1:nGene)-.5)
      names(y)<-unique(exons$TXNAME)
      exons$y<-y[exons$TXNAME]
      exons$top<- exons$y+exonWidth/2
      exons$bottom<- exons$y-exonWidth/2
    }
  }
  logArg<-ifelse(logY,'y','')
  cover[,countCols]<-cover[,countCols]+ifelse(logY,1,0)
  plot(cover[,2],cover[,3],type='n',ylim=c(ifelse(logY,.01,-exonSpace),coverMax),xlab='',ylab=ylab,mgp=c(3.4,1,0),yaxt='n',log=logArg,...)
  if(!justCounts){
    if(!logY)axis(2,pretty(c(0,coverMax)))
    else axis(2,10^pretty(log10(c(1,coverMax)))+1,10^pretty(log10(c(1,coverMax))))
  }
  mtext(sprintf('Chromosome %s position',sub('chr','',sub(':.*','',reg))),1,line=2.5)
  if(!justCounts)for(jj in countCols)lines(cover[,2],cover[,jj],col=cols[jj-2]) 
  legend('top',names(files),col=cols,lwd=1,inset=-.01,bty='n')
  if(nGene>0){
    for(txid in unique(exons$TXNAME)){
      thisExons<-exons[exons$TXNAME==txid,]
      widestExon<-thisExons[which.max(thisExons$EXONEND-thisExons$EXONSTART),]
      xPos<-widestExon$EXONSTART+(widestExon$EXONEND-widestExon$EXONSTART)/2
      usrPos<-par('usr')[1:2]
      if(xPos<usrPos[1]|xPos>usrPos[2])xPos<-mean(usrPos)
      text(xPos,thisExons$y[1],thisExons$symbol[1],col='grey',cex=.7)
      rect(thisExons$EXONSTART,thisExons$top,thisExons$EXONEND,thisExons$bottom)
      if(nrow(thisExons)>1)segments(thisExons$EXONEND[-nrow(thisExons)],thisExons$y[1],thisExons$EXONSTART[-1],thisExons$y[1])
    }
  }
  if(!is.null(counts)&&nrow(counts)>0){
    countProps<-apply(counts[,countCols],2,function(x)x/sum(x))
    maxCountProp<-max(countProps)
    convertProps<-countProps/maxCountProp*coverMax
    if(is.null(nrow(convertProps)))convertProps<-matrix(convertProps,nrow=1)
    for(i in 1:nrow(convertProps))segments(counts[i,'start'],convertProps[i,],counts[i,'end'],convertProps[i,],col=cols,lwd=2)
    prettyCountProps<-pretty(c(0,countProps))
    axis(ifelse(justCounts,2,4),prettyCountProps/maxCountProp*coverMax,prettyCountProps)
    if(!justCounts)mtext('Proportion of reads',4,3,las=0)
  }
}
