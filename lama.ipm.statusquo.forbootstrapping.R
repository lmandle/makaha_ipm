############################################################################################
########################### LAMA IPM ###############################################
##status Quo restoration

rm(list=ls(all=TRUE)) ## Clear the workspace
library(car)
library(lattice)
library(MASS)
library(msm)
library(popbio)
library(gridBase)
library(lattice)
library(glmmTMB)
library(MuMIn)
library(ggplot2)

### Vital rate models *******************

##Adult survival:calculate from 2004-2017; does not vary w size

##For bootstrapping added file:
lama.ad.surv<-read.table("lama_adult.survival.forbootstrapping.txt", header=T, sep="\t")

#Calculated survival as: (nt/n0)1/t
percent<-(sum(lama.ad.surv$alive.2017)/sum(lama.ad.surv$alive.2004))
surv<-percent^(1/17)

#calculate intercept
surv.int<-log(surv/(1-surv))

##Adult growth 
lama.ag<-read.table("lama.adult.growth2.txt", header=T, sep="\t")
lama.ag$Site<-as.factor(lama.ag$Site)

mod1<-glmmTMB(log(final.size2)~log(initial.size)+(1|Site), data=lama.ag)
summary(mod1)

g1a<-fixef(mod1)

g.sig2a<-summary(mod1)$sigma^2 		##overall variance


##Seedling growth
lama.sg<-read.table("lama.seedling.growth.txt", header=T, sep="\t")
str(lama.sg)
lama.sg$year<-as.factor(lama.sg$year)

#This is status quo - using only yr 3 "1x1 plots and MG and assuming no
#negative growth - no pigs)

lama.sg3<-subset(lama.sg, year=="yr3")
mod2a.3<-glmmTMB(log(final.c)~log(initial)+(1|year), data=lama.sg3)
summary(mod2a.3)

g1s<-fixef(mod2a.3)

g.sig2s<-summary(mod2a.3)$sigma^2 		##overall variance


##Fecundity
lama.fec<-read.table("lama.fecundity2.txt", header=T, sep="\t")
lama.fec$year<-as.factor(lama.fec$year)

mod3<-glmmTMB(log(fruit)~log(dbh) +(1|tag), data=lama.fec)
summary(mod3)

f1<-fixef(mod3)


# Survival seedlings
lama.surv<-read.table("lama.seedling.survival2.txt", header=T, sep="\t")
str(lama.surv)

lama.surv$tag<-as.factor(lama.surv$tag)
lama.surv$year<-as.factor(lama.surv$year)

#We use yrs 3 and4; high survival, we assume same w and w.out weeding
lama.surv34<-subset(lama.surv,year=="yr3" | year=="yr4")
mod4b<-glmmTMB(survival~log(size) +(1|year), family="binomial", data=lama.surv34)
summary(mod4b)

s1s<-fixef(mod4b)#model without size is same AIC so we keep it in.


#######COEFFICIENTS###########################

nstate<-2
ncoef<-8
p.vec<-array(0,c(3,ncoef,nstate))#state is adult or seedling

#SURVIVAL ADULTS 
p.vec[1,1,1]<-surv.int #intercept for survival for adults, slope is zero. 
#This is based on survival from 2004/5-2017/20.

#exp(x)/(1+exp(x))

#GROWTH ADULTS
p.vec[2,1,1]<-g1a$cond[1]#intercept for growth for adults  
p.vec[2,2,1]<-g1a$cond[2]#slope for growth for adults 

p.vec[2,4,1]<-g.sig2a#g.sigma2 overall variation


##FECUNDITY
p.vec[3,1,1]<-f1$cond[1]#intercept
p.vec[3,2,1]<-f1$cond[2]#slope for fecundity for adults 


#status quo 
p.vec[3,6,1]<-(1-0.89)*0.64*0.31#prob of pred*germ *survival open, whole, all cover


##seedling size class distribution
p.vec[3,7,1]<-log(5.4) #scd MEAN size class distribution
p.vec[3,8,1]<-1.57 #standard deviation of log of size class distribution

#SURVIVAL OF SEEDLINGS
p.vec[1,1,2]<-s1s$cond[1]#intercept
p.vec[1,2,2]<-s1s$cond[2]#slope for sdlg survival

#GROWTH SEEDLINGS
p.vec[2,1,2]<- g1s$cond[1]# intercept for growth for seedlings 
p.vec[2,2,2]<- g1s$cond[2]#slope for growth for seedlings

p.vec[2,4,2]<-g.sig2s 

p.vec


## The IPM FUNCTIONS: s(x), g(y,x), f(y,x), p(y,x), K(y,x)####### 
####ADULTS
####A.- ADULT SURVIVAL function s(x)
sx.adults<-function(x, pvec){
  surv<-ifelse(x<log(5),0.996,1)#or 0.986 if 98%, this is at 99%
  xbeta.adults<-pvec[1,1,1]+pvec[1,2,1]*x; #slope is 0
  s<-(exp(xbeta.adults)/(1+exp(xbeta.adults)))*surv
  return(s);
}



#check function with graph
plot(sx.adults(seq(1,3.5,length.out=50), p.vec)~(seq(1,3.5,length.out=50)), type="l",col="black",lty=1,xlab="Adult size (log DBH)", ylab="Survival",ylim=c(0.8,1))
#


### B. ADULT GROWTH function g(y,x) 
gyx.adults<-function(y, x, pvec) {
  mux.adults<-pvec[2,1,1] + pvec[2,2,1]*x 
  sigmax2.adu<-pvec[2,4,1]#*exp(2*pvec[2,5,1]*mux.adults)
  sigmax.adu<-sqrt(sigmax2.adu)
  g<-dnorm(y, mux.adults, sigmax.adu)
  return(g)
}

##Check function with graph
gx<-function(x, pvec, year) {
  mux<-pvec[2,1,1] + pvec[2,2,1]*x 
  return(mux)
}

plot(gx(seq(0,3.2,length.out=50),p.vec)~(seq(0,3.2,length.out=50)), type="l",col="grey",lty=1,lwd=2,xlab="Adult size a t (log DBH)", ylab="Adult size at time (t+1)")
points(log(lama.ag$initial.size), log(lama.ag$final.size), col="black",xlab="DBH at t (log)", ylab="DBH at t+1 (log)")
#legend(2.6, 3.5, c("Year1", "Year2", "Year3"), bg="white", bty="n", lwd=c(3,3,3),lty=c(1,2,3))


####C.The ADULT SURVIVAL-GROWTH function P(y, x)
pyx.adults<-function(y,x, pvec) { 
  p<-sx.adults(x, pvec)*gyx.adults(y, x, pvec)
  return(p) 
}


###D.FERTILITY FUNCTION
fyx.adults<-function(y, x, pvec) {
  p.cap<-ifelse(x<0.7,0,1)#put in the function #same as log(2)
  no.cap<-pvec[3,1,1] + pvec[3,2,1]*x;
  surv<-pvec[3,6,1]
  scd<-dnorm(y, pvec[3,7,1], pvec[3,8,1])
  f<-p.cap*no.cap*surv*scd    
  
  return(f)
}

#check graph wrt fruit production
fx<-function(x, pvec) {
  p.cap<-ifelse(x<0.7,0,1)
  no.cap<-pvec[3,1,1] + pvec[3,2,1]*x;
  f<-p.cap*no.cap
  return(f)
}

plot(fx(seq(1,3.2,length.out=50),p.vec)~(seq(1,3.2,length.out=50)), type="l",col="black",lty=1,xlab="Adult size (log DBH)", ylab="Log(Seeds produced)", ylim=c(4.5,8.5))
points(log(lama.fec$dbh),log(lama.fec$fruit),type="p",lwd=1)

#Check graph - this is mean number of seedlings produced
fx<-function(x, pvec) {
  p.cap<-ifelse(x<0.7,0,1)
  no.cap<-pvec[3,1,1] + pvec[3,2,1]*x;
  surv<-pvec[3,6,1]
  f<-p.cap*no.cap*surv    #around 2 seedlings per adult
  return(f)
}

plot(fx(seq(2,3.2,length.out=50),p.vec)~(seq(2,3.2,length.out=50)), type="l",col="black",lty=1,ylab="No.seedlings produced",xlab="Adult size (log DBH)")

##Looking at probability density of new seedlings
plot(fyx.adults(seq(0.5,3,length.out=50),seq(2,3.2,length.out=50),p.vec)~(seq(0.5,3,length.out=50)), type="l",col="black")


###SEEDLINGS
###A.- SEEDLINGS SURVIVAL function s(x)
sx.seedlings<-function(x, pvec) {
  xbeta.seedling<-pvec[1,1,2]+pvec[1,2,2]*x; 
  s<-exp(xbeta.seedling)/(1+exp(xbeta.seedling))
  return(s);
}

plot(sx.seedlings(seq(0.6,5,length.out=50), p.vec)~(seq(0.6,5,length.out=50)), type="l",col="black",lty=1,xlab="Seedling size (log height)", ylab="Survival",ylim=c(0,1))
points(log(lama.surv$size),lama.surv$survival,type="p",lwd=1)


### B. SEEDLING GROWTH function g(y,x)
gyx.seedlings<-function(y, x,pvec) {
  mux.seedlings<-pvec[2,1,2] + pvec[2,2,2]*x;
  sigmax2.seedlings<-pvec[2,4,2]#*exp(2*pvec[2,5,2]*mux.seedlings)
  sigmax.seedlings<-sqrt(sigmax2.seedlings)
  g<-dnorm(y, mux.seedlings, sigmax.seedlings)
  return(g)
}



#check function with graph
gx.s<-function(x, pvec) {
  mux<-pvec[2,1,2] + pvec[2,2,2]*x 
  return(mux)
}

plot(gx.s(seq(0.6,5,length.out=50),p.vec)~(seq(0.6,5,length.out=50)), type="l",col="grey",lty=1,lwd=2,xlab="Seedling size a t (log height)", ylab="Seedling size at time (t+1)")
points(log(lama.sg$initial),log(lama.sg$final),col="black",xlab="Height at t (log)", ylab="Height at t+1 (log)")
#legend(-2, 0.75, c("Year1", "Year2", "Year3"), bg="white", bty="n", lwd=c(3,3,3),lty=c(1,2,3))


### C.SEEDLING SURVIVAL-GROWTH function P(y, x)
pyx.seedlings<-function(y,x, pvec) { 
  p<-sx.seedlings(x, pvec)*gyx.seedlings(y, x, pvec)
  return(p) 
}



#########KERNEL K + F###

#ADULTS.yx (Growth-survival+fecundity)
## The (master) KERNEL: K(y,x)= p(y,x) + f(y,x)
Kyx.adults<-function(y, x, pvec) {
  adult.k<-pyx.adults(y, x, pvec)+fyx.adults(y, x, pvec)#But we will separate these below
  return(adult.k) 
}

##SEEDLINGS.yx (growth plus survival)
Kys.seedlings<-function(y, x, pvec){
  seedling.k<-pyx.seedlings(y, x, pvec)
  return(seedling.k) 
}


##########################BIGMAT FUNCTION################################

bigmat<-function(bigM, pvec){
  
  ad.min.size<-0.9*log(1) #adults are dbh                       
  ad.max.size<-1.1*log(24)
  
  sd.min.size<-0.99*log(2) #seedlings are height
  sd.max.size<-1.01*log(150) 
  
  h.adult=(ad.max.size-ad.min.size)/(bigM+1);#bins
  y.adult = (seq(ad.min.size, ad.max.size, length=bigM) + seq(ad.min.size + h.adult, 
                                                              ad.max.size+ h.adult, length=bigM))/2; #midpoints ** OG FIXED "ad.max.size.size"
  
  h.seedling=(sd.max.size-sd.min.size)/(bigM+1);
  y.seedling=(seq(sd.min.size, sd.max.size, length=(bigM)) + seq(sd.min.size + h.seedling,
                                                                 sd.max.size + h.seedling, length=(bigM)))/2
  
  
  ## COMPLETE KERNEL REQUIRES 4 components 'bolted' together (See Zuidema et al 2010      supplementary material)
  ## (kss  kst)  where kss - seedlings staying sdlings, kts - sdlings transitioning to trees
  ## (kts  ktt)  kst - trees to seedlings (fecundity), ktt - trees staying trees
  
  ##FULL seedling matrix from which to make submatrices for
  ## kss (seedlings->seedlings) and kts (seedlings->adults)
  
  K.seedling.full=outer(y.seedling,y.seedling,pyx.seedlings, pvec);
  KD.seedling.full=h.seedling*K.seedling.full;
  
  ##ADULT MATRIX
  K.adult=outer(y.adult, y.adult, pyx.adults, pvec);
  KD.adult=h.adult*K.adult
  
  ###ADULTS TO SEEDLINGS ##### FECUNDITY MATRIX
  F.sd=outer(y.seedling, y.adult, fyx.adults, pvec);
  FD.sd=h.adult*F.sd
  
  
  ####seedlings to adults    #Code from Clay  OG FIXED a few parenthesis in L412-414
  seedling_behavior<-rbind(KD.seedling.full[1:length(y.seedling[y.seedling<=sd.max.size]),
                                            1:length(y.seedling[y.seedling<=sd.max.size])],
                           rbind(colSums(KD.seedling.full[length(y.seedling[y.seedling<=sd.max.size]):(length(y.seedling[y.seedling<=sd.max.size])+1), 1:length(y.seedling[y.seedling<=sd.max.size])]),  #extract and sum probabilities for seedlings > max.seedling.size
                                 matrix(0,(bigM-1),length(y.seedling[y.seedling<=sd.max.size]))))
  
  
  #grad<-colSums(KD.seedling.full[(length(y.seedling[y.seedling<=see.max.sz]):(length(y.seedling)+1)), 1:length(y.seedling[y.seedling<=see.max.sz])])
  
  full.matrix<-cbind(seedling_behavior,rbind(FD.sd[1:length(y.seedling[y.seedling<=sd.max.size]), 1:bigM],KD.adult))
  return(full.matrix);
  #return(grad)#to check if the colSum is working
}



#########################Eigenvalues and eigenvectors analysis, sensitivity,...
#Plot lambda as a function of bigm
# bigm<-c(500,600,700,800)
# lam<-array(0, c(length(bigm)))
# 
# for(i in 1:length(bigm)){
#       library(popbio)
#       lam[i]<-eigen.analysis(bigmat(bigm[i], pvec=p.vec))$lambda1
#     }
#   
# 
# plot(lam[]~bigm, type="l", col="blue")#lambda is decreasing, need to use larger matrix
#use 800x800, but really 400x400 is sufficient

lambda.orig<-eigen.analysis(bigmat(400, pvec=p.vec))$lambda1



#########################################################
### Bootstrap lambda (code adapted from Kuss et al. 2008)
#########################################################

n.boot=1000 #Kuss used 5000
lambda.boot=data.frame(lambda=rep(NA,n.boot))#just want to collect lambdas

for(b.samp in 1:n.boot){
  
  #Adult survival
  sample.boot=c(sample(1:nrow(lama.ad.surv),replace=T)) #generate bootstrapped sample
  percent.boot<-(sum(lama.ad.surv$alive.2017[sample.boot])/sum(lama.ad.surv$alive.2004[sample.boot])) #recalc with bootstrap sample
  surv.boot<-percent.boot^(1/17)
  if(surv.boot==1){ #add to avoid Inf problem when surv.boot is 1
    surv.boot=0.999
    } 
  surv.int.boot<-log(surv.boot/(1-surv.boot))
  
  #Adult growth
  sample.boot=c(sample(1:nrow(lama.ag),replace=T)) #generate bootstrapped sample
  lama.ag.boot<-data.frame(Site=lama.ag$Site[sample.boot], #create bootstrapped dataset
                           initial.size=lama.ag$initial.size[sample.boot],
                           final.size2=lama.ag$final.size2[sample.boot])
  mod1.boot<-update(mod1, data=lama.ag.boot) #refit model
  mod1<-glmmTMB(log(final.size2)~log(initial.size)+(1|Site), data=lama.ag)
  g1a.boot<-fixef(mod1.boot)
  g.sig2a.boot<-summary(mod1.boot)$sigma^2 
  
  #Seedling growth
  sample.boot=c(sample(1:nrow(lama.sg3),replace=T))
  lama.sg3.boot<-data.frame(initial=lama.sg3$initial[sample.boot],
                            final.c=lama.sg3$final.c[sample.boot],
                            year=lama.sg3$year[sample.boot])
  mod2a.3.boot<-update(mod2a.3, data=lama.sg3.boot)
  g1s.boot<-fixef(mod2a.3.boot)
  g.sig2s.boot<-summary(mod2a.3.boot)$sigma^2 	
  
  #Fecundity
  sample.boot=c(sample(1:nrow(lama.fec),replace=T))
  lama.fec.boot<-data.frame(fruit=lama.fec$fruit[sample.boot],
                            dbh=lama.fec$dbh[sample.boot],
                            tag=lama.fec$tag[sample.boot])
  mod3.boot<-update(mod3, data=lama.fec.boot)
  f1.boot<-fixef(mod3.boot)

  #Seedling survival
  sample.boot=c(sample(1:nrow(lama.surv34),replace=T))
  lama.surv34.boot<-data.frame(survival=lama.surv34$survival[sample.boot],
                               size=lama.surv34$size[sample.boot],
                              year=lama.surv34$year[sample.boot])
  mod4b.boot<-update(mod4b, data=lama.surv34.boot)
  s1s.boot<-fixef(mod4b.boot)
  
  #rebuild p.vec from bootstrapped params
  p.vec.boot<-array(0,c(3,ncoef,nstate))#state is adult or seedling
  
  p.vec.boot[1,1,1]<-surv.int.boot
  p.vec.boot[2,1,1]<-g1a.boot$cond[1]#intercept for growth for adults  
  p.vec.boot[2,2,1]<-g1a.boot$cond[2]#slope for growth for adults 
  p.vec.boot[2,4,1]<-g.sig2a.boot#g.sigma2 overall variation
  p.vec.boot[3,1,1]<-f1.boot$cond[1]#intercept
  p.vec.boot[3,2,1]<-f1.boot$cond[2]#slope for fecundity for adults 
  p.vec.boot[1,1,2]<-s1s.boot$cond[1]#intercept
  p.vec.boot[1,2,2]<-s1s.boot$cond[2]#slope for sdlg survival
  p.vec.boot[2,1,2]<- g1s.boot$cond[1]# intercept for growth for seedlings 
  p.vec.boot[2,2,2]<- g1s.boot$cond[2]#slope for growth for seedlings
  p.vec.boot[2,4,2]<-g.sig2s.boot 
  
  #fixed params:
  p.vec.boot[3,6,1]<-(1-0.89)*0.64*0.31#prob of pred*germ *survival open, whole, all cover
  p.vec.boot[3,7,1]<-log(5.4) #scd MEAN size class distribution
  p.vec.boot[3,8,1]<-1.57 #standard deviation of log of size class distribution
  

  lambda.boot$lambda[b.samp]<-eigen.analysis(bigmat(400, pvec=p.vec.boot))$lambda1
}

date = gsub(":","-",Sys.time()) #get date and time to append to filename  
date = gsub(" ","_",date)
write.csv(lambda.boot, file=paste0("lambda.boot", "_", date, ".csv"))

ci.normal.app=c(mean(lambda.boot$lambda)-1.96*sd(lambda.boot$lambda),mean(lambda.boot$lambda)+1.96*sd(lambda.boot$lambda))
res=c(mean(lambda.boot$lambda,na.rm=T),quantile(lambda.boot$lambda,p=c(0.025,0.5,0.975),na.rm=T),ci.normal.app)

ggplot(data=lambda.boot, aes(lambda))+  
  geom_histogram(bins=50)+
  geom_vline(xintercept=mean(lambda.boot$lambda,na.rm=T))+
  geom_vline(xintercept=median(lambda.boot$lambda, na.rm=T),color="blue")+
  geom_vline(xintercept=quantile(lambda.boot$lambda,p=0.025,na.rm=T), color="red")+
  geom_vline(xintercept=quantile(lambda.boot$lambda,p=0.975,na.rm=T), color="red")+
  geom_vline(xintercept=lambda.orig, color="green")
