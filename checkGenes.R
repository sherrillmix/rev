
source('functions.R')
regs<-c('ACTB'='chr7:5527000-5530800','GAPDH'='chr12:6533900-6538400','HSP90AB1'='chr6:44246100-44253900','PPIH'='chr1:42658400-42676800','LDHA'='chr11:18394300-18408500','NONO'='chrX:71283100-71301400')
pdf('out/housekeeping.pdf')
  for(reg in regs){
    plotRegion(reg,bamFiles,bam2depthBinary='~/installs/bedCount/bam2depth',normalize=FALSE,ylab='Read coverage')
  }
dev.off()

