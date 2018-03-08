info<-read.table('map.tsv',row.names=1,header=TRUE,stringsAsFactors=FALSE,sep='\t',comment='')
info$weight<-info$Weight_to_use
info$study<-'bushman'
goodIds<-rownames(info)[info$toKeep=='keep']
dadaGoodIds<-info[info$toKeep=='keep','Description']

maps<-list.files('validation','^map\\.',full.names=TRUE)
info2<-do.call(rbind,lapply(maps,read.table,stringsAsFactors=FALSE,row.names=1,header=TRUE,sep='\t',comment=''))
