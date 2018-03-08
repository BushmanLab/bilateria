library(dnar)
rareN<-100
if(!dir.exists('work'))dir.create('work')

source('readInfo.R')

otu<-read.csv('../om.csv.gz',stringsAsFactors=FALSE)
rareN<-tapply(otu$count,otu$SampleID,rareEquation,rareN)
write.csv(cbind(names(rareN),'rare'=rareN),'work/rareN.csv',row.names=FALSE)
speciesAbund<-tapply(otu$count,otu$SampleID,function(xx)sort(xx[xx>0]))
save(speciesAbund,file='work/speciesAbund.Rdat')


otuFiles<-list.files('../validation','^om.*\\.csv\\.gz$',full.names=TRUE)
otu2<-do.call(rbind,lapply(otuFiles,read.csv,stringsAsFactors=FALSE))
rareN2<-tapply(otu2$count,otu2$SampleID,rareEquation,rareN)
write.csv(cbind(names(rareN2),'rare'=rareN2),'work/rareN2.csv',row.names=FALSE)


dada<-read.csv('../om.dada.csv.gz',stringsAsFactors=FALSE)
dadaRareN<-tapply(dada$count,dada$SampleID,rareEquation,500)
write.csv(cbind(names(dadaRareN),'rare'=dadaRareN),'work/dadaRareN.csv',row.names=FALSE)

dadaAbund<-tapply(dada$count,dada$SampleID,function(xx)sort(xx[xx>0]))
save(dadaAbund,file='work/dadaAbund.Rdat')


