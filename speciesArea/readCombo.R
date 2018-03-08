if(!file.exists('work/rareN.csv'))source('readData.R')
rare<-read.csv('work/rareN.csv',row.names=1)
rare2<-read.csv('work/rareN2.csv',row.names=1)
source('readInfo.R')

info$otu<-rare[rownames(info),'rare']
info2$otu<-rare2[rownames(info2),'rare']
sharedCols<-intersect(colnames(info),colnames(info2))
combo<-rbind(info[,sharedCols],info2[,sharedCols])
combo<-combo[combo$readCount>200&!is.na(combo$weight)&!is.na(combo$otu),]
combo$speciesName<-tolower(paste(combo$genus,combo$species))
combo$log.weight<-log10(combo$weight)
combo$log.otus<-log10(combo$otu)
weightSpeciesOrder<-tapply(combo$log.weight,combo$speciesName,mean)
combo$speciesId<-as.numeric(factor(combo$speciesName,levels=names(sort(weightSpeciesOrder))))
combo$classId<-as.numeric(as.factor(combo$class))
combo$studyId<-as.numeric(factor(combo$study,levels=c('bushman',unique(combo$study[combo$study!='bushman']))))


