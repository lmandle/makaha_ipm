############################################################################################
########################### OLOPUA IPM ###############################################

###Status quo restoration

rm(list=ls(all=TRUE)) ## Clear the workspace
library(car)
library(lattice)
library(MASS)
library(msm)
library(popbio)
library(gridBase)
library(lattice)
library(glmmTMB)
library(ggplot2)


### Vital rate models *******************

##Adult survival:calculate from 2004-2017; does not vary w size

##For bootstrapping added file:
olopua.ad.surv<-read.table("olopua_adult.survival.forbootstrapping.txt", header=T, sep="\t")

#Calculated survival as: (nt/n0)1/t
percent<-(sum(olopua.ad.surv$alive.2017)/sum(olopua.ad.surv$alive.2004))
surv<-percent^(1/16.5)

#calculate intercept
surv.int<-log(surv/(1-surv))

#Adult growth final model
olopua.ag<-read.table("olopua.adult.growth2.txt", header=T, sep="\t")

mod1<-glmmTMB(log(final)~log(initial)+(1|Site), data=olopua.ag)
summary(mod1)

g1a<-fixef(mod1)

g.sig2a<-summary(mod1)$sigma^2 		##overall variance

#Seedling growth
sdlg<-read.table("olopua.seedling.growth.txt", header=T, sep="\t")

str(sdlg)
sdlg$site<-as.factor(sdlg$site)
sdlg$year<-as.factor(sdlg$year)
sdlg$tag<-as.factor(sdlg$tag)
sdlg$yearsum<-as.factor(sdlg$yearsum)

grow<-subset(sdlg, final>0)

##Status quo, just using 2017-2020
sdlg.h<-subset(sdlg,yearsum=="2017.2")
sdlg.l<-subset(sdlg,yearsum=="2003.07")

g3.sq<-glmmTMB(log(final)~log(initial)+(1|tag), data=sdlg.h)
summary(g3.sq)
g1s<-fixef(g3.sq)
g.sig2s<-summary(g3.sq)$sigma^2 		##overall variance


#survival
#This is later years (2017.20)- assume some weeding
s4.h<-glmmTMB(survival~log(initial)+(1|tag),family="binomial", data=sdlg.h)
summary(s4.h)# intercept model only is same AIC so we keep it.

s1s<-fixef(s4.h)


##Fecundity
olopua.fec<-read.table("olopua.fecundity.txt", header=T, sep="\t")
olopua.fec$tree<-as.factor(olopua.fec$tree)
str(olopua.fec)

#Subset to those producing fruit
v.fec<-subset(olopua.fec, viable.fruit>0)
#t.fec<-subset(olopua.fec, fruit.plus40>0)
#v.fec.40<-subset(olopua.fec, viable.fruit.plus.40>0)

##This assume rat predation in canopy, count only viable fruit
mod3e<-glmmTMB(log(viable.fruit)~log(dbh)+I(log(dbh)^2)+(1|tree),data=v.fec)
summary(mod3e)
f1<-fixef(mod3e)


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
#Each year 3/11 trees produced no fruit, =73% produced

p.vec[3,1,1]<-f1$cond[1]#intercept
p.vec[3,2,1]<-f1$cond[2]#slope for fecundity for adults 
p.vec[3,3,1]<-f1$cond[3]#quadratic term

#(1-prob ground pred)*%germ*%surv to 6 months
#Status quo: 40% were predated across all cover trts since cover mean 20%
p.vec[3,6,1]<-(1-0.4)*0.64*0.5#

p.vec[3,7,1]<-log(7.6) #scd MEAN size class distribution
p.vec[3,8,1]<-log(3.2) #standard deviation size class distribution

#SURVIVAL OF SEEDLINGS
p.vec[1,1,2]<-s1s$cond[1]#intercept
p.vec[1,2,2]<-s1s$cond[2]#slope for sdlg survival

#GROWTH SEEDLINGS
p.vec[2,1,2]<-g1s$cond[1]# intercept for growth for seedlings 
p.vec[2,2,2]<-g1s$cond[2]#slope for growth for seedlings

p.vec[2,4,2]<-g.sig2s 

p.vec


## The IPM FUNCTIONS: s(x), g(y,x), f(y,x), p(y,x), K(y,x)####### 
####ADULTS
####A.- ADULT SURVIVAL function s(x)

#we use < or >10,  since so few individuals 1-5cm dbh
sx.adults<-function(x, pvec){
  surv<-ifelse(x<log(10),0.998,1)# this is 0.982 mean vs 0.98 #olopua <10
  xbeta.adults<-pvec[1,1,1]+pvec[1,2,1]*x; #slope is 0
  s<-(exp(xbeta.adults)/(1+exp(xbeta.adults)))*surv
  return(s);
}


#check function with graph
plot(sx.adults(seq(1,3.5,length.out=50), p.vec)~(seq(1,3.5,length.out=50)), type="l",col="blue",lty=1,xlab="Adult size (log DBH)", ylab="Survival",ylim=c(0,1))



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
points(log(olopua.ag$initial), log(olopua.ag$final), col="black",xlab="DBH at t (log)", ylab="DBH at t+1 (log)")
#legend(2.6, 3.5, c("Year1", "Year2", "Year3"), bg="white", bty="n", lwd=c(3,3,3),lty=c(1,2,3))


####C.The ADULT SURVIVAL-GROWTH function P(y, x)
pyx.adults<-function(y,x, pvec) { 
  p<-sx.adults(x, pvec)*gyx.adults(y, x, pvec)
  return(p) 
}


###D.FERTILITY FUNCTION
fyx.adults<-function(y, x, pvec) {
  p.cap<-ifelse(x<0.7,0,0.73)#prob of fruiting is 73%
  no.cap<-pvec[3,1,1] + pvec[3,2,1]*x+pvec[3,3,1]*I(x^2);
  surv<-pvec[3,6,1]#this is viable*pred*germ*survival
  scd<-dnorm(y, pvec[3,7,1], pvec[3,8,1])
  f<-p.cap*no.cap*surv*scd    
  
  
  return(f)
}

#check graph wrt fruit production
fx<-function(x, pvec) {
  p.cap<-ifelse(x<0.7,0,0.8)
  no.cap<-pvec[3,1,1] + pvec[3,2,1]*x+pvec[3,3,1]*I(x^2);
  f<-p.cap*no.cap
  return(f)
}


plot(fx(seq(1,3.2,length.out=50),p.vec)~(seq(1,3.2,length.out=50)), type="l",col="blue",lty=1,xlab="Adult size (log DBH)", ylab="Log(Seeds produced)")
points(log(fec$dbh),log(fec$fruit),type="p",lwd=1, col="blue")

points.default()#Check graph - this is mean number of seedlings produced
fx<-function(x, pvec) {
  p.cap<-ifelse(x<0.7,0,0.8)
  no.cap<-pvec[3,1,1] + pvec[3,2,1]*x+pvec[3,3,1]*I(x^2);
  surv<-pvec[3,6,1]
  f<-p.cap*no.cap*surv    
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
points(sx.seedlings(seq(0.6,5,length.out=50), p.vec)~(seq(0.6,5,length.out=50)), type="l",col="black",lty=1,xlab="Seedling size (log height)", ylab="Survival",ylim=c(0,1))


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

plot(gx.s(seq(0.6,5,length.out=50),p.vec)~(seq(0.6,5,length.out=50)), type="l",col="grey",lty=1,lwd=2,xlab="Adult size a t (log DBH)", ylab="Adult size at time (t+1)")
points(log(grow$initial),log(grow$final),col="black",xlab="Height at t (log)", ylab="Height at t+1 (log)")


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
  sd.max.size<-1.01*log(180) 
  
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
#library(popbio)
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

lambda.orig<-eigen.analysis(bigmat(400, pvec=p.vec))$lambda1 #you can change site and year to get the different values

#########################################################
### Bootstrap lambda (code adapted from Kuss et al. 2008)
#########################################################

n.boot=1100 #Kuss used 5000
lambda.boot=data.frame(lambda=rep(NA,n.boot))#just want to collect lambdas

for(b.samp in 1:n.boot){
  
  #Adult survival
  sample.boot=c(sample(1:nrow(olopua.ad.surv),replace=T)) #generate bootstrapped sample
  percent.boot<-(sum(olopua.ad.surv$alive.2017[sample.boot])/sum(olopua.ad.surv$alive.2004[sample.boot])) #recalc with bootstrap sample
  surv.boot<-percent.boot^(1/16.5)
  if(surv.boot==1){ #add to avoid Inf problem when surv.boot is 1
    surv.boot=0.999
  } 
  surv.int.boot<-log(surv.boot/(1-surv.boot))
  
  #Adult growth
  sample.boot=c(sample(1:nrow(olopua.ag),replace=T)) #generate bootstrapped sample
  olopua.ag.boot<-data.frame(Site=olopua.ag$Site[sample.boot], #create bootstrapped dataset
                           initial=olopua.ag$initial[sample.boot],
                           final=olopua.ag$final[sample.boot])
  mod1.boot<-update(mod1, data=olopua.ag.boot) #refit model
  g1a.boot<-fixef(mod1.boot)
  g.sig2a.boot<-summary(mod1.boot)$sigma^2 
  
  #Seedling growth
  sample.boot=c(sample(1:nrow(sdlg.h),replace=T))
  olopua.sdlg.g.boot<-data.frame(initial=sdlg.h$initial[sample.boot],
                            final=sdlg.h$final[sample.boot],
                            tag=sdlg.h$tag[sample.boot])
  #g3.sq.boot<-update(g3.sq, data=olopua.sdlg.g.boot)
  g3.sq.boot<-glm(log(final)~log(initial), data=olopua.sdlg.g.boot) #remove individual-level random effect from bootstrapped model
  g1s.boot<-coef(g3.sq.boot)
  g.sig2s.boot<-summary(g3.sq.boot)$dispersion	
  
  #Fecundity
  sample.boot=c(sample(1:nrow(v.fec),replace=T))
  v.fec.boot<-data.frame(viable.fruit=v.fec$viable.fruit[sample.boot],
                            dbh=v.fec$dbh[sample.boot],
                            tree=v.fec$tree[sample.boot])
  #mod3e.boot<-update(mod3e, data=v.fec.boot)
  mod3e.boot<-glm(log(viable.fruit)~log(dbh)+I(log(dbh)^2),data=v.fec.boot) #remove tree random eff for bootstrapping
  f1.boot<-coef(mod3e.boot)
  
  #Seedling survival
  sample.boot=c(sample(1:nrow(sdlg.h),replace=T))
  olopua.sdlg.s.boot<-data.frame(survival=sdlg.h$survival[sample.boot],
                               initial=sdlg.h$initial[sample.boot],
                               tag=sdlg.h$tag[sample.boot])
  #s4.h.boot<-update(s4.h, data=olopua.sdlg.s.boot) #fitting problems depending on resampling, b/c of tag being a factor
  s4.h.boot<-glm(survival~log(initial),family="binomial", data=olopua.sdlg.s.boot) #drop tag from random effect?
  s1s.boot<-coef(s4.h.boot)
  
  #rebuild p.vec from bootstrapped params
  p.vec.boot<-array(0,c(3,ncoef,nstate))#state is adult or seedling
  
  p.vec.boot[1,1,1]<-surv.int.boot #intercept for survival for adults, slope is zero. 
  p.vec.boot[2,1,1]<-g1a.boot$cond[1]#intercept for growth for adults  
  p.vec.boot[2,2,1]<-g1a.boot$cond[2]#slope for growth for adults 
  p.vec.boot[2,4,1]<-g.sig2a.boot#g.sigma2 overall variation
  p.vec.boot[3,1,1]<-f1.boot[1]#intercept
  p.vec.boot[3,2,1]<-f1.boot[2]#slope for fecundity for adults 
  p.vec.boot[3,3,1]<-f1.boot[3]#quadratic term
  p.vec.boot[1,1,2]<-s1s.boot[1]#intercept
  p.vec.boot[1,2,2]<-s1s.boot[2]#slope for sdlg survival
  p.vec.boot[2,1,2]<-g1s.boot[1]# intercept for growth for seedlings 
  p.vec.boot[2,2,2]<-g1s.boot[2]#slope for growth for seedlings
  p.vec.boot[2,4,2]<-g.sig2s.boot
 
  #fixed
  p.vec.boot[3,6,1]<-(1-0.4)*0.64*0.5 #(1-prob ground pred)*%germ*%surv to 6 months
  p.vec.boot[3,7,1]<-log(7.6) #scd MEAN size class distribution
  p.vec.boot[3,8,1]<-log(3.2) #standard deviation size class distribution
  
  #bootstrapped lambda
  lambda.boot$lambda[b.samp]<-tryCatch(
    expr=eigen.analysis(bigmat(400, pvec=p.vec.boot))$lambda1,
    error=function(cond){
      message(cond)
      return(NA)},
   warning=function(cond){
     message(cond)
     return(NULL)})
}

date = gsub(":","-",Sys.time()) #get date and time to append to filename  
date = gsub(" ","_",date)
write.csv(lambda.boot, file=paste0("olopua.lambda.boot", "_", date, ".csv"))

ci.normal.app=c(mean(lambda.boot$lambda)-1.96*sd(lambda.boot$lambda),mean(lambda.boot$lambda)+1.96*sd(lambda.boot$lambda))
res=c(mean(lambda.boot$lambda,na.rm=T),quantile(lambda.boot$lambda,p=c(0.025,0.5,0.975),na.rm=T),ci.normal.app)

ggplot(data=lambda.boot, aes(lambda))+  
  geom_histogram(bins=50)+
  geom_vline(xintercept=mean(lambda.boot$lambda,na.rm=T))+
  geom_vline(xintercept=median(lambda.boot$lambda, na.rm=T),color="blue")+
  geom_vline(xintercept=quantile(lambda.boot$lambda,p=0.025,na.rm=T), color="red")+
  geom_vline(xintercept=quantile(lambda.boot$lambda,p=0.975,na.rm=T), color="red")+
  geom_vline(xintercept=lambda.orig, color="green")

