library(parallel)
library(sads)
library(vegan)
library(digest)
library(dnar)
source('functions.R')
source('readInfo.R')
if(!dir.exists('out'))dir.create('out')
if(!file.exists('work/speciesAbund.Rdat'))source('readData.R')
load('work/speciesAbund.Rdat')

fitters<-c('Broken stick'=fitbs2,'Geometric'=fitgeom2,'Log series'=fitls2,'Neutral (MZSM)'=fitmzsm2,'Power'=fitpower,'Power bend'=fitpowbend2,'Log normal'=fitlnorm,'Poisson lognormal'=fitpoilog2,'Gamma'=fitgamma2,'Weibull'=fitweibull) 
fits<-mclapply(speciesAbund[goodIds],function(xx){
  cat('.')
  out<-lapply(fitters,function(func,xx){
    tryCatch(func(xx),error=function(e)return(NULL))
  },xx)
  return(out)
},mc.cores=50,mc.preschedule=FALSE)

plotBics<-function(fits,speciesAbund,info,outFile,speciesOrder=NULL){
  bics<-do.call(rbind,lapply(fits,function(xx){
    sapply(xx,function(yy){
      if(is.null(yy))return(NA)
      else return(AIC(yy))
    })
  }))
  bicDiff<-t(apply(bics,1,function(xx)xx-min(xx,na.rm=TRUE)))
  bicDiff[is.infinite(bicDiff)]<-NA
  bicDiff<-bicDiff[do.call(order,cbind(info[rownames(bics),c('phylum','supersuperclass','superclass','class','superorder','clade','order','family','genus','species')])),]
  colnames(bicDiff)<-names(fitters)
  #need to reorder to match bicDiff reordering
  speciesAbund<-speciesAbund[rownames(bicDiff)]
  fits<-fits[rownames(bicDiff)]
  bicCondense<-do.call(rbind,by(bicDiff,info[rownames(bicDiff),'common'],function(xx)apply(xx,2,mean,na.rm=TRUE)))
  if(is.null(speciesOrder))bicCondense<-bicCondense[orderIn(rownames(bicCondense),info[rownames(bicDiff),'common']),]
  else bicCondense<-bicCondense[orderIn(rownames(bicCondense),speciesOrder),]
  rownames(bicCondense)<-sub('TigerSha','Tiger sha',sub('SSBShark','Sandbar shark',sub('SBUShark','Bull shark',sub('Vietnames ','Vietnamese ',rownames(bicCondense)))))
  bicCondense<-bicCondense[rownames(bicCondense)!='unknown_or_none',]
  #color brewer
  fitCols<-c('#a6cee3AA','#1f78b4AA','#b2df8aAA','#33a02cAA','#fb9a99AA','#e31a1cAA','#fdbf6fAA','#ff7f00AA','#cab2d6AA','#6a3d9aAA','#ffff99AA')[1:ncol(bicDiff)]
  names(fitCols)<-colnames(bicDiff)
  targets<-c('Cricket'='Cricket','Human'='Human','Tiger Shark'='Tiger Shark','Right Whale'='Right Whale')
  ids<-sapply(names(targets),function(xx)which(info[rownames(bicDiff),'common']==xx)[1])
  names(ids)<-targets
  preds<-mclapply(targets,function(target,...){
    ii<-ids[target]
    message(target)
    out<-mclapply(names(fits[[ii]]),function(jj){
      message(jj)
      #time out after .5 hour
      tryCatch(R.utils::withTimeout(lazyRadPred(fits[[ii]][[jj]]),timeout=3600*.5),TimeoutException=function(ex){warning('Time out ',jj);NULL})
    },mc.cores=10)
    names(out)<-names(fits[[ii]])
    return(out)
  },mc.cores=5)
  names(preds)<-targets
  pdf(outFile,width=8,height=10,useDingbats=FALSE)
    layout(matrix(c(rep(5,4),1:4),ncol=2),width=c(.7,.3))
    par(mar=c(3,3.5,1,1))
    targetTops<-targetBottoms<-targetLefts<-c()
    sortTargets<-targets[orderIn(targets,rownames(bicCondense),decreasing=TRUE)]
    for(target in sortTargets){
      ii<-ids[target]
      thisRad<-rad(speciesAbund[[ii]])
      plot(thisRad$rank,thisRad$abund,main=target,las=1,log='y',xlab='',ylab='OTU abundance',mgp=c(2.5,.7,0),yaxt='n')
      logAxis(2,las=1)
      title(xlab='OTU rank',mgp=c(1.6,1,0))
      isNA<-sapply(names(fitCols),function(xx){
        if(is.null(fits[[ii]][[xx]]))return(TRUE)
        if(is.null(preds[[target]][[xx]]))return(TRUE)
        thisPred<-preds[[target]][[xx]]
        thisPred[is.infinite(thisPred[,'abund']),'abund']<-sum(thisRad)*2
        lines(thisPred,col=fitCols[xx],lwd=2)
        return(FALSE)
      })
      targetLefts[target]<-grconvertX(par('usr')[1],to='ndc')
      targetTops[target]<-grconvertY(10^par('usr')[4],to='ndc')
      targetBottoms[target]<-grconvertY(10^par('usr')[3],to='ndc')
      if(target==tail(sortTargets,1))legend('topright',sub(' ',' ',names(fitCols)),col=fitCols,lwd=2,bty='n',cex=.85,ncol=1,y.intersp=.8)
    }
    par(mar=c(7,10,1,1.5))
    cols<-colorRampPalette(c('red','blue'))(200)
    breaks<-seq(min(log10(1+bicCondense),na.rm=TRUE),max(log10(1+bicCondense),na.rm=TRUE),length.out=201)
    image(1:ncol(bicCondense),1:nrow(bicCondense),t(log10(1+bicCondense)),xaxt='n',col=cols,breaks=breaks,ylab='',xlab='',yaxt='n')
    box()
    axis(2,1:nrow(bicCondense),rownames(bicCondense),las=2,cex.axis=.7)
    slantAxis(1,1:ncol(bicCondense),colnames(bicCondense),srt=-30,adj=c(.1,.5))
    for(target in targets){
      y1<-which(rownames(bicCondense)==target)+c(.5,-.5)
      x1<-rep(par('usr')[2],2)
      x2<-rep(grconvertX(targetLefts[target],'ndc','user'),2)
      y2<-grconvertY(c(targetTops[target],targetBottoms[target]),'ndc','user')
      polygon(c(x1,rev(x2)),c(y1,rev(y2)),border='#00000011',col='#00000011',xpd=NA)
    }
    ticks<-1:floor(max(log10(bicCondense),na.rm=TRUE))
    insetScale(breaks,cols,c(.015,.015,.025,.3),main='Difference from minimum AIC',at=log10(c(0,10^(ticks)+1)),labels=c(0,sapply(ticks,function(xx)as.expression(bquote(10^.(xx))))))
  dev.off()
}

commonOrder<-c("Encrusted Sand Tubed Worm", "Giant Millipede", "Woodlouse", "Pink Shrimp", "Crab", "Hermit Crab", "Cricket", "Hissing Cockroach", "Praying Mantis", "Boxelder Bug", "Bedbug", "Eastern Yellowjacket", "European Wool Carder Bee", "Two-spotted Bumble Bee", "Long-horned Bee", "Western Honey Bee", "Mealworm", "Rhinoceros Beetle", "Hornworm", "Dagger Moth", "Mosquito", "Drosophila", "Fly", "Striped Sea Cucumber", "Skate", "Tiger Shark", "Lemon Shark", "Sandbar Shark", "Bull Shark", "Sea Robin", "Fluke", "Bearded Dragon", "Eagle Owl", "Cockatiel", "Parakeet", "Horse", "Cat", "Ferret", "Dog", "Alpaca", "Pig", "Cow", "Sheep", "Goat", "Right Whale", "Humpback Whale", "Fin Whale", "Rabbit", "Guinea Pig", "Hamster", "Gerbil", "Mouse", "Rat", "Red Colobus", "Macaque", "Sooty Mangabey", "Mandrill", "Gorilla", "Bonobo", "Western Chimpanzee", "Eastern Chimpanzee", "Nigeria-Cameroon Chimpanzee", "Central Chimpanzee", "Human")
if(any(!commonOrder %in% info$common))stop('Problem ordering commons')

#dada2
load('work/dadaAbund.Rdat')
dadaFits<-mclapply(dadaAbund[goodIds],function(xx){
  cat('.')
  xx<-xx[xx>1]
  out<-lapply(fitters,function(func,xx){
    tryCatch(func(xx,trunc=1),error=function(e)return(NULL))
  },xx)
  return(out)
},mc.cores=50,mc.preschedule=FALSE)


plotBics(fits,lapply(speciesAbund,function(xx)xx[xx>1]),info,'out/speciesAbundanceFit.pdf',speciesOrder=rev(commonOrder))

plotBics(dadaFits,lapply(dadaAbund,function(xx)xx[xx>1]),info,'out/dadaSpeciesAbundanceFit.pdf',speciesOrder=rev(commonOrder))
