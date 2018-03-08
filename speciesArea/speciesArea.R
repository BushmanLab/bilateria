library(rstan)
library(vipor)
library(dnar)
library(vioplot)
#use multiple cores (careful with big models on small computers)
options(mc.cores = parallel::detectCores())
#cache compilations (doesn't actually work that well, at least with model_code=)
rstan_options(auto_write = TRUE)

source('readCombo.R')
stanCode<-'
  functions{
    int[] p_greater(vector x, real y,int n){
      int out[n];
      for(ii in 1:n){
        if(x[ii]>y)out[ii]=1;
        else out[ii]=2;
      }
      return(out);
    }
  }
  data{
    int<lower=0> n;
    real otus[n]; //log
    int<lower=1> ids[n]; //assumes id is sorted by weigh
    int<lower=0> nIds;
    vector[nIds] weights; //log
    int<lower=1> classes[nIds];
    int<lower=0> nClasses;
    int<lower=1> studies[n];
    int<lower=0> nStudies;
  }
  parameters {
    real interceptMu;
    real<lower=0> interceptSd;
    vector[nStudies-1] rawIntercepts;
    vector<lower=0>[nClasses] classSigmas;
    vector<lower=0>[nStudies] sigmas;
    vector[nIds] rawMus;
    real beta;
  }
  transformed parameters {
    vector[n] mus;
    vector[nIds] speciesMus;
    vector[nStudies] intercepts;
    speciesMus= rawMus .* classSigmas[classes];
    intercepts=interceptMu + append_row(0.0,rawIntercepts)*interceptSd;
    mus = intercepts[studies] + beta* weights[ids] + speciesMus[ids];
  }
  model {
    classSigmas ~ gamma(1,.1);
    sigmas ~ gamma(1,.1);
    rawIntercepts ~ normal(0,1);
    rawMus ~ normal(0,1);
    otus ~ normal(mus,sigmas[studies]);
  }
'
fit<-stan(
  model_code=stanCode,
  data=list(
    otus=combo$log.otus,
    n=nrow(combo),
    ids=combo$speciesId,
    nIds=max(combo$speciesId),
    weights=tapply(combo$log.weight,combo$speciesId,mean),
    classes=tapply(combo$classId,combo$speciesId,'[[',1),
    nClasses=max(combo$classId),
    studies=combo$studyId,
    nStudies=max(combo$studyId)
  ),
  control=list(max_treedepth=11,adapt_delta=.99),
  chains=32,iter=100000,thin=40
)


plotFit<-function(fit,dat){
  sims<-extract(fit)
  studyMeans<-apply(sims[['intercepts']],2,mean)
  beta<-mean(sims[['beta']])
  interceptMu<-mean(sims[['interceptMu']])
  weights<-tapply(dat$log.weight,dat$speciesId,mean)
  speciesCols<-rainbow(max(dat$speciesId),alpha=.8)
  studyCols<-rainbow.lab(max(dat$studyId),alpha=.8)
  names(studyCols)<-sapply(1:length(studyCols),function(xx)dat[dat$studyId==xx,'study'][1])
  xlim<-range(dat$log.weight)
  ylim<-range(dat$log.otus)
  dat$adjusted<-dat$log.otus-studyMeans[dat$studyId]+interceptMu
  #plotting
  plot(dat$log.weight,dat$log.otus,pch=21,bg=studyCols[dat$studyId],col=NA,lwd=2,las=1,xlab='Average weight of species',ylab='Number of rarefied OTUs',xlim=xlim,ylim=ylim,xaxt='n',yaxt='n',main='Raw',mgp=c(2.4,.8,0),cex=.7)
  prettyY<-c(1,2,5,10,20,50,100)
  axis(2,log10(prettyY),prettyY,las=1)
  yTick<-c(1:10,(2:10)*10,(2:10)*100)
  axis(2,log10(yTick),rep('',length(yTick)),tcl=-.25)
  prettyX<-10^seq(-4,8,2)
  axis(1,log10(prettyX),sapply(log10(prettyX),function(xx)as.expression(bquote(10^.(xx)))),las=1)
  #adjusted
  plot(dat$log.weight,dat$adjusted,pch=21,bg=studyCols[dat$studyId],col=NA,lwd=2,las=1,xlab='Average weight of species',ylab='Number of rarefied OTUs',xlim=xlim,ylim=ylim,main='Adjusted',xaxt='n',yaxt='n',mgp=c(2.4,.8,0),cex=.7)
  axis(2,log10(prettyY),prettyY,las=1)
  axis(2,log10(yTick),rep('',length(yTick)),tcl=-.25)
  axis(1,log10(prettyX),sapply(log10(prettyX),function(xx)as.expression(bquote(10^.(xx)))),las=1)
  fakeX<-seq(min(dat$log.weight),max(dat$log.weight),length.out=100)
  predX<-interceptMu+beta*fakeX
  multiPred<-do.call(rbind,mapply(function(a,b)a+fakeX*b,sims[['interceptMu']],sims[['beta']],SIMPLIFY=FALSE))
  cis<-apply(multiPred,2,quantile,c(.025,.975))
  lines(fakeX,predX)
  polygon(c(fakeX,rev(fakeX)),c(cis[1,],rev(cis[2,])),col='#00000022',border=NA)
  legend('bottomleft',names(studyCols),pt.bg=studyCols,pch=21,inset=c(-.8,-.46),ncol=3,xpd=NA,col=NA)
}
pdf('out/bayesSpeciesArea.pdf',width=8,height=6)
  par(mfrow=c(1,2),mar=c(9.5,4,1,.3))
  plotFit(fit,combo)
dev.off()


