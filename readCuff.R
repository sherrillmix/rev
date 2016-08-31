dataDir<-'work/diff'
dataType<-'gene'
geneDiffFile<-sprintf('%s/%s_exp.diff',dataDir,dataType)
outFile<-sprintf('out/%s_exp.diff.gz',dataType)
geneDiff<-read.table(geneDiffFile,header=TRUE,stringsAsFactors=FALSE)
kgXref<-read.table('~/installs/star/kgXRef.tsv.gz',sep='\t',stringsAsFactors=FALSE,quote='',comment='')
rownames(kgXref)<-kgXref$V1
geneDiff$name<-kgXref[geneDiff$gene_id,'V5']
kgAlias<-read.table('~/installs/star/kgAlias.tsv.gz',sep='\t',stringsAsFactors=FALSE,quote='',comment='')
kgAlias<-kgAlias[kgAlias$V1!=kgAlias$V2,]
aliases<-tapply(kgAlias$V2,kgAlias$V1,paste,collapse='|')
geneDiff$name[is.na(geneDiff$name)]<-aliases[geneDiff$gene_id[is.na(geneDiff$name)]]
write.csv(geneDiff,gzfile(outFile))

