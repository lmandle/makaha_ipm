############################################################################################
########################### ALAHE'E IPM ###############################################
##alahee status quo

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
##### Vital Rate Regressions *******************

####Survival adults
#No variation with size. Across sites, over 13-17 years we get a consistent 0.03 annual mortality.

##For bootstrapping added file:
alahee.ad.surv<-read.table("alahee_adult.survival.forbootstrapping.txt", header=T, sep="\t")

#Calculated survival as: (nt/n0)1/t
percent<-(sum(alahee.ad.surv$alive.2017)/sum(alahee.ad.surv$alive.2004))
surv<-percent^(1/12.5)

#calculate intercept
surv.int<-log(surv/(1-surv))
surv.int

#####Adult growth
a.grow<-read.table("alahee.growth.txt", header=T, sep="\t")
a.grow$plot<-as.factor(a.grow$plot)
a.grow$treatment<-as.factor(a.grow$treatment)
a.grow$tag<-as.factor(a.grow$tag)
a.grow$site<-as.factor(a.grow$site)

str(a.grow)

a.grow<-subset(a.grow,initial>0 &final.year>0)

mod.g1<-glmmTMB(log(final.year)~log(initial)+(1|site/plot),data=a.grow)
summary(mod.g1)

g1a<-fixef(mod.g1)

g.sig2a<-summary(mod.g1)$sigma^2 		##overall variance

###Seedling growth
alahee.sd<-read.table("alahee.sldg.txt", header=T, sep="\t")
str(alahee.sd)
alahee.sd$site<-as.factor(alahee.sd$site)
alahee.sg<-subset(alahee.sd, final.year>0)

g2s<-glmmTMB(log(final.year)~log(initial)+ (1|site), data=alahee.sg)
summary(g2s)

g1s<-fixef(g2s)

g.sig2s<-summary(g2s)$sigma^2 		##overall variance

# Survival seedlings
alahee.s<-subset(alahee.sd, survival!="NA")
head(alahee.s)
str(alahee.s)

mod1a<-glmmTMB(survival~log(initial)+(1|site), family="binomial", data=alahee.s)#tried also w Site, no diffs across sites
summary(mod1a)

s1s<-fixef(mod1a)


##clonal survival and clonal growth
alahee.cl<-read.table("alahee.clonal.txt", header=T, sep="\t")
str(alahee.cl)
alahee.cl$tag<-as.factor(alahee.cl$tag)
alahee.cl$year<-as.factor(alahee.cl$year)

#clonal growth - we keep same across trts. 
alahee.cl<-subset(alahee.cl, final>0)

g2c<-glmmTMB(log(final)~log(initial)+ (1|tag), data=alahee.cl)
summary(g2c)

g1c<-fixef(g2c)
g.sig2c<-summary(g2c)$sigma^2 

#clonal survival - this is very high, we don't vary across trts
alahee.cl.s<-subset(alahee.cl, survival!="NA")
head(alahee.cl.s)
str(alahee.cl.s)

mod5<-glmmTMB(survival~log(initial)+ (1|year), family="binomial", data=alahee.cl.s)
summary(mod5)#intercept model same fit

s1c<-fixef(mod5)


##Fecundity
alahee.fec<-read.table("alahee.fecundity2.txt", header=T, sep="\t")
str(alahee.fec)
alahee.fec$tree<-as.factor(alahee.fec$tree)
alahee.fec$year<-as.factor(alahee.fec$year)

fec<-subset(alahee.fec, fruit>0)

mod3b<-glmmTMB(log(fruit)~log(dbh)+ (1|tree),data=fec)
summary(mod3b)

f1<-fixef(mod3b)



#######COEFFICIENTS###########################
nstate<-3#adult or seedling or clone
ncoef<-8 #these are intercept, slope etc; here we don't actually use 8 but left from before
p.vec<-array(0,c(4,ncoef,nstate)) #4 is survival, growth, fecundity, clonal rep

#SURVIVAL ADULTS 
p.vec[1,1,1]<-surv.int #intercept for survival for sizes >1 (97%), slope is zero. 

#exp(3.5)/(1+exp(3.5))


#GROWTH ADULTS 
p.vec[2,1,1]<-g1a$cond[1]#intercept for growth for adults  
p.vec[2,2,1]<-g1a$cond[2]#slope for growth for adults 

p.vec[2,4,1]<-g.sig2a#g.sigma2 overall variation


##FECUNDITY
#Predisperal predation (moths)
p.vec[3,1,1]<-f1$cond[1]#intercept
p.vec[3,2,1]<-f1$cond[2]#slope for fecundity for adults 

##Scenarios
#(1-prob pre disp pred)*(1-prob post disp pred)*prob germ*surv to census
#In this case though we have predisp pred in base model (mean 62% removed; 38% remain)
#Or high range is 82%; 

#Status quo
p.vec[3,6,1]<-(1)*(1-0.89)*0.4*0.24#average


## seedling size
p.vec[3,7,1]<-log(4) #scd MEAN size class distribution
p.vec[3,8,1]<-log(2) #standard deviation size class distribution

#clonal reproduction
#exp(-3)/(1+exp(-3))=0.04 prob of producing
p.vec[4,1,1]<-1 #number of clones/tree
p.vec[4,2,1]<-0 #doesn't vary w size
p.vec[4,7,1]<-log(25.5) #scd MEAN size class distribution
p.vec[4,8,1]<-log(2.5) #standard deviation size class distribution

#SURVIVAL OF SEEDLINGS
p.vec[1,1,2]<-s1s$cond[1]#intercept
p.vec[1,2,2]<-s1s$cond[2]##slope for sdlg survival

p.vec[1,1,3]<-s1c$cond[1]#intercept clones
p.vec[1,2,3]<-s1c$cond[2]##slope for clones survival


#GROWTH SEEDLINGS
p.vec[2,1,2]<-g1s$cond[1]# intercept for growth for seedlings 
p.vec[2,2,2]<-g1s$cond[2]#slope for growth for seedlings

p.vec[2,4,2]<-g.sig2s 

#growth clones
p.vec[2,1,3]<-g1c$cond[1]# intercept for growth for clones 
p.vec[2,2,3]<-g1c$cond[2]#slope for growth for clones

p.vec[2,4,3]<-g.sig2c 

p.vec


## The IPM FUNCTIONS: s(x), g(y,x), f(y,x), p(y,x), K(y,x)####### 
####ADULTS
####A.- ADULT SURVIVAL function s(x)
sx.adults<-function(x, pvec){
  xbeta.adults<-pvec[1,1,1]+pvec[1,2,1]*x; #slope is 0
  s<-exp(xbeta.adults)/(1+exp(xbeta.adults))
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

plot(gx(seq(0,2.5,length.out=50),p.vec)~(seq(0,2.5,length.out=50)), type="l",col="grey",lty=1,lwd=2,xlab="Adult size a t (log DBH)", ylab="Adult size at time (t+1)")
points(log(a.grow$initial), log(a.grow$final.year), col="black",xlab="DBH at t (log)", ylab="DBH at t+1 (log)")


####C.The ADULT SURVIVAL-GROWTH function P(y, x)
pyx.adults<-function(y,x, pvec) { 
  p<-sx.adults(x, pvec)*gyx.adults(y, x, pvec)
  return(p) 
}


###D.FERTILITY FUNCTION
fyx.adults<-function(y, x, pvec) {
  p.cap<-ifelse(x<log(1),0,0.722)#put in the function #log(1)# 13/18 (2 yrs) trees produced fruit
  no.cap<-pvec[3,1,1]+ pvec[3,2,1]*x#+pvec[3,3,1]*I(x^2);
  surv<-pvec[3,6,1]
  scd<-dnorm(y, pvec[3,7,1], pvec[3,8,1])
  f<-p.cap*no.cap*surv*scd    
  
  return(f)
}

#check graph wrt fruit production
fx<-function(x, pvec) {
  p.cap<-ifelse(x<log(1),0,1)
  no.cap<-pvec[3,1,1] + pvec[3,2,1]*x#+pvec[3,3,1]*I(x^2);
  f<-p.cap*no.cap
  return(f)
}


plot(fx(seq(0,3.2,length.out=50),p.vec)~(seq(0,3.2,length.out=50)), type="l",col="blue",lty=1,xlab="Adult size (log DBH)", ylab="Log(Seeds produced)")
points(log(fec$dbh),log(fec$fruit),type="p",lwd=1, col="blue")

#Check graph - this is mean number of seedlings produced
fx<-function(x, pvec) {
  p.cap<-ifelse(x<log(1),0,1)
  no.cap<-pvec[3,1,1] + pvec[3,2,1]*x#+pvec[3,3,1]*I(x^2);
  surv<-pvec[3,6,1]
  f<-p.cap*no.cap*surv    #around 1 seedlings per adult
  return(f)
}

plot(fx(seq(0,3.2,length.out=50),p.vec)~(seq(0,3.2,length.out=50)), type="l",col="black",lty=1,ylab="No.seedlings produced",xlab="Adult size (log DBH)")

#this is about right based on the field

##Looking at probability density of new seedlings
plot(fyx.adults(seq(0.5,3,length.out=50),seq(0,3.2,length.out=50),p.vec)~(seq(0.5,3,length.out=50)), type="l",col="black")


# ### D. CLONAL FERTILITY
cyx.adults<-function(y, x, pvec) {
  p.cap<-ifelse(x<log(1),0,0.1)#0.04% produced root suckers that surv to census, but we estimate 10% 
  no.cap<-pvec[4,1,1]+ pvec[4,2,1]*x#+pvec[3,3,1]*I(x^2);#1 per tree
  surv<-1#surv to census, our numbers include this so it's 1
  scd<-dnorm(y, pvec[4,7,1], pvec[4,8,1])
  c<-p.cap*no.cap*surv*scd    
  
  return(c)
}

##Looking at probability density of new seedlings
plot(cyx.adults(seq(0.5,6,length.out=50),seq(0,3.2,length.out=50),p.vec)~(seq(0.5,6,length.out=50)), type="l",col="black")


###SEEDLINGS
###A.- SEEDLINGS SURVIVAL function s(x)

#This is with clones surviving at different rates than seedlings
sx.seedlings<-function(x, pvec) {
  xbeta.seedling<-(0.47*pvec[1,1,2]+0.53*pvec[1,1,3])+(0.47*pvec[1,2,2]+0.53*pvec[1,2,3])*x;
  s<-exp(xbeta.seedling)/(1+exp(xbeta.seedling))
  return(s);
}

plot(sx.seedlings(seq(0.6,5,length.out=50), p.vec)~(seq(0.6,5,length.out=50)), type="l",col="blue",lty=1,xlab="Seedling size (log height)", ylab="Survival",ylim=c(0,1))


### B. SEEDLING GROWTH function g(y,x)

#seedlings and clones together
gyx.seedlings<-function(y, x,pvec) {
  mux.seedlings<-(0.47*pvec[2,1,2]+0.53*pvec[2,1,3])+(0.47*pvec[2,2,2]+0.53*pvec[2,2,3])*x;
  sigmax2.seedlings<-(0.47*pvec[2,4,2]+0.53*pvec[2,4,3])
  sigmax.seedlings<-sqrt(sigmax2.seedlings)
  g<-dnorm(y, mux.seedlings, sigmax.seedlings)
  return(g)
}

#check graph
gx.s<-function(x, pvec) {
  mux<-(0.47*pvec[2,1,2]+0.53*pvec[2,1,3])+(0.47*pvec[2,2,2]+0.53*pvec[2,2,3])*x 
  return(mux)
}

points(gx.s(seq(0.6,5,length.out=50),p.vec)~(seq(0.6,5,length.out=50)), type="l",col="blue",lty=1,lwd=2,xlab="seedling size a t (log ht)", ylab="Seedling size at time (t+1)")


### C.SEEDLING SURVIVAL-GROWTH function P(y, x)
pyx.seedlings<-function(y,x, pvec) { 
  p<-sx.seedlings(x, pvec)*gyx.seedlings(y, x, pvec)
  return(p) 
}


#########KERNEL K + F###

#ADULTS.yx (Growth-survival+fecundity)
## The (master) KERNEL: K(y,x)= p(y,x) + f(y,x)

Kyx.adults<-function(y, x, pvec) {
  adult.k<-pyx.adults(y, x, pvec)+fyx.adults(y, x, pvec)+cyx.adults(y, x, pvec)
  return(adult.k) 
}


#put together the fertility and clonality
tyx.adults<-function(y, x, pvec) {
  fert.k<-fyx.adults(y, x, pvec)+cyx.adults(y, x, pvec)
  return(fert.k) 
}


##SEEDLINGS.yx (growth plus survival)
Kys.seedlings<-function(y, x, pvec){
  seedling.k<-pyx.seedlings(y, x, pvec)
  return(seedling.k) 
}


##########################BIGMAT FUNCTION################################

bigmat<-function(bigM, pvec){
  
  ad.min.size<-0.9*log(1) #adults are dbh                       
  ad.max.size<-1.1*log(15)
  
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
  F.sd=outer(y.seedling, y.adult, tyx.adults,pvec);
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

n.boot=50 #Kuss used 5000
lambda.boot=data.frame(lambda=rep(NA,n.boot))#collect lambdas
p.vec.boot<-array(0,c(n.boot,4,ncoef,nstate))#save p.vecs, state is adult or seedling

for(b.samp in 1:n.boot){
  
  #Adult survival
  sample.boot=c(sample(1:nrow(alahee.ad.surv),replace=T)) #generate bootstrapped sample
  percent.boot<-(sum(alahee.ad.surv$alive.2017[sample.boot])/sum(alahee.ad.surv$alive.2004[sample.boot]))
  surv.boot<-percent.boot^(1/12.5)
  if(surv.boot==1){ #add to avoid Inf problem when surv.boot is 1
    surv.boot=0.999
  } 
  surv.int.boot<-log(surv.boot/(1-surv.boot))
  surv.int.boot

  #Adult growth
  sample.boot=c(sample(1:nrow(a.grow),replace=T)) #generate bootstrapped sample
  alahee.ag.boot<-data.frame(site=a.grow$site[sample.boot], #create bootstrapped dataset
                             plot=a.grow$plot[sample.boot], 
                             initial=a.grow$initial[sample.boot],
                             final.year=a.grow$final[sample.boot])
  mod.g1.boot<-update(mod.g1, data=alahee.ag.boot) #refit model 
  g1a.boot<-fixef(mod.g1.boot)
  g.sig2a.boot<-summary(mod.g1.boot)$sigma^2 

  
  #Seedling growth
  sample.boot=c(sample(1:nrow(alahee.sg),replace=T))
  alahee.sdlg.g.boot<-data.frame(initial=alahee.sg$initial[sample.boot],
                                 final.year=alahee.sg$final.year[sample.boot],
                                 site=alahee.sg$site[sample.boot])
  g2s.boot<-update(g2s, data=alahee.sdlg.g.boot)
  g1s.boot<-fixef(g2s.boot)
  g.sig2s.boot<-summary(g2s.boot)$sigma^2 

  
  #Fecundity
  sample.boot=c(sample(1:nrow(fec),replace=T))
  fec.boot<-data.frame(fruit=fec$fruit[sample.boot],
                         dbh=fec$dbh[sample.boot],
                         tree=fec$tree[sample.boot])
  mod3b.boot<-glm(log(fruit)~log(dbh),data=fec.boot)#refit without individual-level random effects
  f1.boot<-coef(mod3b.boot)
  
  #Seedling survival
  sample.boot=c(sample(1:nrow(alahee.s),replace=T))
  alahee.sdlg.s.boot<-data.frame(survival=alahee.s$survival[sample.boot],
                                 initial=alahee.s$initial[sample.boot],
                                 site=alahee.s$site[sample.boot])
  mod1a.boot<-update(mod1a, data=alahee.sdlg.s.boot)
  s1s.boot<-fixef(mod1a.boot)

  #Clonal growth
  sample.boot=c(sample(1:nrow(alahee.cl),replace=T))
  alahee.cg.boot<-data.frame(initial=alahee.cl$initial[sample.boot],
                                 final=alahee.cl$final[sample.boot],
                                 tag=alahee.cl$tag[sample.boot])
  g2c.boot<-glm(log(final)~log(initial), data=alahee.cg.boot)
  g1c.boot<-coef(g2c.boot)
  g.sig2c.boot<-summary(g2c.boot)$dispersion
  
  #Clonal survival
  sample.boot=c(sample(1:nrow(alahee.cl.s),replace=T))
  alahee.sg.boot<-data.frame(initial=alahee.cl.s$initial[sample.boot],
                             year=alahee.cl.s$year[sample.boot],
                             survival=alahee.cl.s$survival[sample.boot])
  mod5.boot<-update(mod5, data=alahee.sg.boot) #will get false convergence warning when 100% survival
  s1c.boot<-fixef(mod5.boot)
  
  #rebuild p.vec from bootstrapped params
  p.vec.boot[b.samp,1,1,1]<-surv.int.boot #intercept for survival for sizes >1 (97%), slope is zero. 
  p.vec.boot[b.samp,2,1,1]<-g1a.boot$cond[1]#intercept for growth for adults  
  p.vec.boot[b.samp,2,2,1]<-g1a.boot$cond[2]#slope for growth for adults 
  p.vec.boot[b.samp,2,4,1]<-g.sig2a.boot#g.sigma2 overall variation
  p.vec.boot[b.samp,3,1,1]<-f1.boot[1]#intercept
  p.vec.boot[b.samp,3,2,1]<-f1.boot[2]#slope for fecundity for adults 
  p.vec.boot[b.samp,1,1,2]<-s1s.boot$cond[1]#intercept
  p.vec.boot[b.samp,1,2,2]<-s1s.boot$cond[2]##slope for sdlg survival
  p.vec.boot[b.samp,1,1,3]<-s1c.boot$cond[1]#intercept clones
  p.vec.boot[b.samp,1,2,3]<-s1c.boot$cond[2]##slope for clones survival
  p.vec.boot[b.samp,2,1,2]<-g1s.boot$cond[1]# intercept for growth for seedlings 
  p.vec.boot[b.samp,2,2,2]<-g1s.boot$cond[2]#slope for growth for seedlings
  p.vec.boot[b.samp,2,4,2]<-g.sig2s.boot 
  p.vec.boot[b.samp,2,1,3]<-g1c.boot[1]# intercept for growth for clones 
  p.vec.boot[b.samp,2,2,3]<-g1c.boot[2]#slope for growth for clones
  p.vec.boot[b.samp,2,4,3]<-g.sig2c.boot
  
  #fixed parameters
  p.vec.boot[b.samp,3,6,1]<-(1)*(1-0.89)*0.4*0.24#average
  p.vec.boot[b.samp,3,7,1]<-log(4) #scd MEAN size class distribution
  p.vec.boot[b.samp,3,8,1]<-log(2) #standard deviation size class distribution
  p.vec.boot[b.samp,4,1,1]<-1 #number of clones/tree
  p.vec.boot[b.samp,4,2,1]<-0 #doesn't vary w size
  p.vec.boot[b.samp,4,7,1]<-log(25.5) #scd MEAN size class distribution
  p.vec.boot[b.samp,4,8,1]<-log(2.5) #standard deviation size class distribution
  
 
  #bootstrapped lambda
  lambda.boot$lambda[b.samp]<-tryCatch(
    expr=eigen.analysis(bigmat(400, pvec=p.vec.boot[b.samp,,,]))$lambda1,
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
