
source('functions.R')
regs<-c('ACTB'='chr7:5527000-5530800-','GAPDH'='chr12:6533900-6538400+','HSP90AB1'='chr6:44246100-44253900+','PPIH'='chr1:42658400-42676800+','LDHA'='chr11:18394300-18408500+','NONO'='chrX:71283100-71301400+')
pdf('out/housekeeping.pdf')
  for(regName in names(regs)){
    reg<-regs[regName]
    plotRegion(reg,bamFiles,bam2depthBinary='~/installs/bedCount/bam2depth',normalize=FALSE,ylab='Read coverage')
    title(sprintf('%s',regName))
    substring(reg,nchar(reg))<-chartr('+-','-+',substring(reg,nchar(reg)))
    plotRegion(reg,bamFiles,bam2depthBinary='~/installs/bedCount/bam2depth',normalize=FALSE,ylab='Read coverage (reverse)')
    title(sprintf('%s reverse strand',regName))
  }
dev.off()

