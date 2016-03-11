
reg<-'chr7:5527000-5530800'
pdf('test.pdf')
	plotRegion(reg,bamFiles,bam2depthBinary='~/installs/bedCount/bam2depth',normalize=FALSE)
dev.off()

