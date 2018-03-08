source('readInfo.R')
if(!dir.exists('out'))dir.create('out')

dadaRareN<-read.csv('work/dadaRareN.csv',row.names=1)
dadaRareN<-structure(dadaRareN$rare,.Names=rownames(dadaRareN))
dadaRareN<-dadaRareN[goodIds]
dadaWeights<-info[names(dadaRareN),'weight']
dadaSpecies<-sub('unknown_or_none $',' sp.',apply(info[names(dadaRareN),c('genus','species','subspecies')],1,paste,collapse=' '))
dadaClass<-info[names(dadaRareN),'class']
dadaSpecies[grepl('^unknown_or_none',dadaSpecies)]<-NA
dadaCommon<-info[names(dadaRareN),'common']
dadaDat<-data.frame(
  'weight'=log10(tapply(dadaWeights,dadaSpecies,mean)),
  'otu'=log10(tapply(dadaRareN,dadaSpecies,mean)),
  'class'=tapply(dadaClass,dadaSpecies,unique),
  'common'=tapply(dadaCommon,dadaSpecies,unique),
  stringsAsFactors=FALSE
)
summary(lm(otu~weight,dadaDat))


library(ggrepel)
library(scales)
scientific_10 <- function(x) {
    parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}
classColors<-structure(
  c('#8dd3c7','#ffed6f','#bc80bd','#fccde5','#80b1d3','#fdb462','#b3de69','#fb8072','#bebada','#ccebc5'),
  .Names = c("Actinopteri", "Aves", "Chondrichthyes", "Diplopoda", "Holothuroidea", "Insecta", "Malacostraca", "Mammalia", "Polychaeta", "Reptilia")
)
pdf('out/dadaSpeciesArea.pdf',width=7,height=7,useDingbats=FALSE)
print(
  ggplot(dadaDat, aes(10^weight, 10^otu, label = sub(' ','\n',common))) +
    geom_smooth(method=lm, color='#00000033')+
    geom_text_repel(size=2.5,box.padding=.12,point.padding=.3,lineheight=.7,min.segment.length=.1,max.iter=3e4,nudge_y=.02,nudge_x=.1,color='#00000099') +
    #scale_x_continuous(trans='log10',label=scientific_10)+
    ylab('Number of sequence variants (rarefied to 500 reads)')+
    xlab('Animal weight (g)')+
    theme_classic(base_size = 16)+
    scale_x_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x,n=3),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    annotation_logticks(sides = "lb",short = unit(1,"mm"), mid = unit(1,"mm"), long = unit(2,"mm"))+
    geom_point(fill = classColors[dadaDat$class],pch=21,size=4)
)
dev.off()

