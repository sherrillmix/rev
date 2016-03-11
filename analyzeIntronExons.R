readCounts<-lapply(list.files('work','regions.count'), function(x)read.table(file.path('work',x)))
names(readCounts)<-sub('(Only)?_regions.count','',list.files('work','regions.count'))

intronExon<-mapply(function(exon,gene){names(exon)<-c('reg_reg','reg_cut',sprintf('reg_%s',names(bamFiles)));names(gene)<-c('gene_reg','gene_cut',sprintf('gene_%s',names(bamFiles)));return(cbind(exon,gene))},readCounts[c('exon','intron')],readCounts[c('exonGenes','intronGenes')],SIMPLIFY=FALSE)
names(intronExon)<-c('exon','intron')

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


