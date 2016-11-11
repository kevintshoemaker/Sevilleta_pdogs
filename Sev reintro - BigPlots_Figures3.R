

## clear the environment
rm(list=ls())

ANA = F
KEVIN = F
KEVINOFFICE = F
ELIZABETH = F
ELIZABETHOFFICE = T

###################################
##########  SET WORKING DIRECTORY

if(KEVIN) rootDir <- "C:\\Users\\Kevin\\Dropbox\\Sev Pdog Study\\BIG plot study\\"
if(KEVINOFFICE) rootDir <- "E:\\Dropbox\\Sev Pdog Study\\BIG plot study\\"
if(ANA) rootDir <- ""
if(ELIZABETH) rootDir <- "~/Dropbox/Nevada/Sevilleta/Sev Pdog Study/BIG plot study/"
if(ELIZABETHOFFICE) rootDir <- "C:\\Users\\ehunter\\Dropbox\\Nevada\\Sevilleta\\Sev Pdog Study\\BIG plot study\\"

ScriptDir <- paste(rootDir,"Rscript",sep="")
DataDir <- paste(rootDir,"Data",sep="")
FiguresDir <- paste(rootDir,"RawFigures",sep="")
ResultsDir <- paste(rootDir,"Results",sep="")
BUGSDir <- paste(rootDir,"BUGS",sep="")

setwd(BUGSDir)

getwd()

# setwd(BUGSDir)
# load(".RData")

##################################
##########  LOAD PACKAGES

library(R2WinBUGS)
library(lubridate)
library(Hmisc)


################################
########## LOAD MODEL RESULTS AND GLOBAL PARAMETERS
setwd(BUGSDir)
load("SevModel10_1mnth_newprecip_2016-07-16.RData")
#load("SevModel_2015-11-25.RData")   # load the latest BUGS results

setwd(DataDir)
load("GlobalParams.RData")



#################################
########  HISTOGRAMS OF FITTED PARAMS
#################################

## names(Mod$sims.list)

#hist(Mod$sims.list$pbadyr)
hist(Mod$sims.list$precipEff)    # no discernible precip effect anymore
hist(Mod$sims.list$relEff)       # very strong negative release effect (what about interactions?)
#hist(Mod$sims.list$predEff)
hist(Mod$sims.list$p0)
hist(Mod$sims.list$phi0)

hist(Mod$sims.list$maleEff)
hist(Mod$sims.list$juvEff)      # no longer has any effect- it was all a release-year effect  

#hist(Mod$sims.list$juvRelEff)   # interesting- released juveniles do MUCH better than released adults
hist(Mod$sims.list$springRelEff)  # wow- spring releases do much worse than fall releases???

hist(Mod$sims.list$soilEff)    # not much soil effect- slightly negative!

hist(Mod$sims.list$hoardEff)

#hist(Mod$sims.list$cohortEff[,6])
hist(Mod$sims.list$tagRetentionRate)

###############################
############  PlOT OUT RESULTS
###############################




#################################
############  ABUNDANCE
#################################


# NOTE: to really use the H-T estimator here, do we need to account for tag loss?? To estimate abundance,
  ##  we need to be able to know the number of individuals sampled and the probability of capture. We DONT
  ##  need to account for tag loss! So, the code below should work fine...



## should include information about the number of releases. 

setwd(FiguresDir)

## Plot out the population size estimates for each site on one graph

Grids <- c("B","D","F")

plot=1
for(plot in 1:length(Grids)){
  thisGrid <- Grids[plot]
  
  graphics.off()
  setwd(FiguresDir)
  filename=sprintf("Abundance_Time_grid%s_F%s.svg",thisGrid,Sys.Date())
  svg(file=filename,height=6,width=5)
  N.mean <- numeric(nperiodsG-1)
  N.97.5 <- numeric(nperiodsG-1)
  N.02.5 <- numeric(nperiodsG-1)
  for(t in 1:(nperiodsG-1)){
    N.mean[t] <- eval(parse(text=sprintf("mean(Mod$sims.list$N_plot%s[,t])",thisGrid)))
    N.97.5[t] <- eval(parse(text=sprintf("quantile(Mod$sims.list$N_plot%s[,t],.975)",thisGrid)))
    N.02.5[t] <- eval(parse(text=sprintf("quantile(Mod$sims.list$N_plot%s[,t],.025)",thisGrid)))
  }
     # extract the number of releases for this grid
  releases <- as.vector(as.matrix(sumReleases)[thisGrid,])
  
  par(mai=c(2,1,1,1))
  plot(c(1:(nperiodsG)),releases, type="h", lend=1, lwd=10,col=gray(0.3),
       main=sprintf("Pdog abundance, sampling grid %s",thisGrid), ylab="# prairie dogs", xlab="",
       ylim=c(0,max(releases)),xaxt="n")
  errbar(x=c(2:(nperiodsG)),y=N.mean,yplus=N.97.5,yminus=N.02.5,add=T,pch=20)
  axis(1,at=c(1:(nperiodsG)),labels=realperiodsG,las=2)
  dev.off()
}

#All plots together now (I'm not sure if these sums across plots are the correct way to do this,
#but the Ntot estimates were way too high.  But the first period seems to have an N that is way too low??)
N.mean <- numeric(nperiodsG-1)
N.97.5 <- numeric(nperiodsG-1)
N.02.5 <- numeric(nperiodsG-1)
for(t in 1:(nperiodsG-1)){
	N.mean[t] <- sum(mean(Mod$sims.list$N_plotB[,t]), mean(Mod$sims.list$N_plotD[,t]), mean(Mod$sims.list$N_plotF[,t]))
	N.97.5[t] <- sum(quantile(Mod$sims.list$N_plotB[,t],.975), quantile(Mod$sims.list$N_plotD[,t],.975), quantile(Mod$sims.list$N_plotF[,t],.975))
	N.02.5[t] <- sum(quantile(Mod$sims.list$N_plotB[,t],.025), quantile(Mod$sims.list$N_plotD[,t],.025), quantile(Mod$sims.list$N_plotF[,t],.025))
  }

releases <- colSums(as.matrix(sumReleases))

setwd(FiguresDir)
filename=paste("Abundance_Time_ALL_", Sys.Date(), ".pdf", sep="")
pdf(filename, height=6, width=7)
par(mai=c(1.5,1,0.2,0.2))
plot(c(1:(nperiodsG)+0.1),releases, type="h", lend=1, lwd=10,col=gray(0.3),
       ylab="Prairie dog abundance", xlab="", ylim=c(0,max(releases)),xaxt="n")
errbar(x=c(2:(nperiodsG))-0.1,y=N.mean,yplus=N.97.5,yminus=N.02.5,add=T,pch=19)
axis(1,at=c(1:(nperiodsG)),labels=c("Summer 2010", "Summer 2011", "Spring 2012", "Summer 2012", "Spring 2013", "Summer 2013", "Spring 2014", "Summer 2014", "Summer[1] 2015", "Summer[2] 2015", "Spring 2016", "Summer 2016"),las=2)
dev.off()

##########################
########  SURVIVAL


## Plot out the survival estimates for each site on one graph

graphics.off()
setwd(FiguresDir)
filename=sprintf("Female Survival_byTime_%s.pdf",Sys.Date())
pdf(file=filename,height=5,width=5)
phi.mean <- numeric(nperiodsG-1)
phi.97.5 <- numeric(nperiodsG-1)
phi.02.5 <- numeric(nperiodsG-1)
for(t in 1:(nperiodsG-1)){
  phi.mean[t] <- mean(Mod$sims.list$femphi[,t])     # femphi  # malephi    # reladultphi
  phi.97.5[t] <- quantile(Mod$sims.list$femphi[,t],.975)
  phi.02.5[t] <- quantile(Mod$sims.list$femphi[,t],.025)
}

par(mai=c(1.5,1,0.2,0.2))
plot(c(1:(nperiodsG-1)),phi.mean,
     ylab="Survival rate", xlab="",type="p",
     ylim=c(0,1),xaxt="n", las=1)
errbar(x=c(1:(nperiodsG-1)),y=phi.mean,yplus=phi.97.5,yminus=phi.02.5,add=T,pch="")

axis(1,at=c(1:(nperiodsG-1)),labels=c("Summer 2010", "Summer 2011", "Spring 2012", "Summer 2012", "Spring 2013", "Summer 2013", "Spring 2014", "Summer 2014", "Spring 2015", "Spring 2016", "Summer 2016"),las=2)
dev.off()

###### 
#### make CSV for Sigma plot

setwd(ResultsDir)
tempname <- paste("SurvToPlot_",Sys.Date(),".csv",sep="")
write.csv(data.frame(mean=phi.mean,upper=phi.97.5,lower=phi.02.5),tempname)


names(Mod$sims.list)

############
### visualize difference between released juv, spring released adult fem, spring released adult male, etc

### mean survival of females vs males vs juveniles

# mu.phi[i,t]   <- logit.phi + precipEff*precip[t] + soilEff*soil[t] +      
#   juvEff*isJuv[i,t] + maleEff*isMale[i] + relEff*isRel[i,t] + 
#   hoardEff*isHoard[i] + juvRelEff*isJuv[i,t]*isRel[i,t] + springRelEff*springRel[i,t]


nreps <- length(Mod$sims.list$phi0) 
dims <- 4
dimnms1 <- c("Adult Female","Adult Male","Juvenile Female","Juvenile Male")    

Surv.fem <- with(Mod$sims.list,plogis(qlogis(phi0)))
Surv.mal <- with(Mod$sims.list,plogis(qlogis(phi0)+maleEff))
Surv.juv.fem <- with(Mod$sims.list,plogis(qlogis(phi0)+juvEff))
Surv.juv.mal <- with(Mod$sims.list,plogis(qlogis(phi0)+maleEff+juvEff))

setwd(FiguresDir)
filename=sprintf("Survival_categor_%s.pdf",Sys.Date()) 
pdf(filename, width=7, height=3.5)
#SexAge plot
par(mar=c(2,4,1,1))
alldata <- c(Surv.fem,Surv.mal,Surv.juv.fem,Surv.juv.mal)
alldatanames <- rep(dimnms1,each=nreps)
alldata <- data.frame(survival=alldata,group=alldatanames)
alldata$group <- factor(as.character(alldata$group),levels=dimnms1)
plot(survival~group,data=alldata)
dev.off()


#Release survival rates
nreps <- length(Mod$sims.list$phi0) 
dims <- 2
dimnms1 <- c("Spring","Summer")    

Surv.spr.adt <- Mod$sims.list$springadrelsurv_mo
Surv.sum.adt <- Mod$sims.list$falladrelsurv_mo

setwd(FiguresDir)
filename=sprintf("Survival_releases_categor_%s.pdf",Sys.Date()) 
pdf(filename, width=4, height=3.5)
#SexAge plot
par(mar=c(2,4,1,1))
alldata <- c(Surv.spr.adt, Surv.sum.adt)
alldatanames <- rep(dimnms1,each=nreps)
alldata <- data.frame(survival=alldata,group=alldatanames)
alldata$group <- factor(as.character(alldata$group),levels=dimnms1)
plot(survival~group,data=alldata, ylab="Survival (one month)")
dev.off()


#################################
############  FECUNDITY
#################################
setwd(BUGSDir)
load("SevModel3_fert_2016-07-16_1mnth_newprecip.RData")
#Histograms of fecundity and lambda
setwd(FiguresDir)
filename=sprintf("Lambda_fert_%s.pdf",Sys.Date()) 
pdf(filename, width=7, height=5)
par(mfrow=c(2,1), mar=c(4,4,1,1))
hist(fertoutput, breaks=30, xlab="Recruitment (juveniles[F]/female/year)", main="", xlim=c(0,1.5))
segments(x0=0.88, y0=0, x1=0.88, y1=1000, lty=2, lwd=2) #Recruitment needed to reach lambda=1
text(x=1.4, y=450, labels="(a)", cex=1.5)
#Estimates from literature
#points(x=1.46, y=800, pch=19) #Haynie 2003 estimate
#arrows(x0=1.46, y0=800, x1=1.46+0.6, y1=800, angle=90, length=0.1)
#arrows(x0=1.46, y0=800, x1=1.46-0.6, y1=800, angle=90, length=0.1)
#text(x=1.46, y=700, labels="(Haynie et al. 2003)")
#points(x=1.89, y=500, pch=19) #Hoogland 2001 estimate
#arrows(x0=1.89, y0=500, x1=1.89+0.62, y1=500, angle=90, length=0.1)
#arrows(x0=1.89, y0=500, x1=1.89-0.62, y1=500, angle=90, length=0.1)
#text(x=1.89, y=400, labels="(Hoogland 2001)")
#points(x=2.14, y=200, pch=19) #Fitzgerald 1974 estimate
#arrows(x0=2.14, y0=200, x1=2.14+1.25, y1=200, angle=90, length=0.1)
#arrows(x0=2.14, y0=200, x1=2.14-1.25, y1=200, angle=90, length=0.1)
#text(x=2.14, y=100, labels="(Fitzgerald & Lechleitner 1974)")
#points(x=mean(fertoutput), y=1200, pch=19) #our estimate
#arrows(x0=mean(fertoutput), y0=1200, x1=mean(fertoutput)+sd(fertoutput), y1=1200, angle=90, length=0.1)
#arrows(x0=mean(fertoutput), y0=1200, x1=mean(fertoutput)-sd(fertoutput), y1=1200, angle=90, length=0.1)
#text(x=mean(fertoutput), y=1080, labels="(this study)")
par(mar=c(4,4,1,1))
hist(lambda, breaks=50, xlab="Population growth rate (lambda)", main="", xlim=c(0.2,1.55))
segments(x0=1, y0=0, x1=1, y1=4000, lty=2, lwd=2)
text(x=1.47, y=450, labels="(b)", cex=1.5)
dev.off()

####
#Fecundity for each year
setwd(FiguresDir)
filename=sprintf("Fert_year_%s.pdf",Sys.Date()) 
pdf(filename, width=7, height=5)
plot(density(littersize[temp.12.3]), xlim=c(0,1), xlab="Recruitment (juveniles[F]/female/year)", lwd=3, col="orangered4", main="", ylab="Probability density")
lines(density(littersize[temp.13.3]), col="orangered1", lwd=3)
lines(density(littersize[temp.14.3]), col="orange", lwd=3)
lines(density(littersize[temp.15.3]), col="gold", lwd=3)
segments(x0=mean(littersize[temp.12.3]), y0=6, y1=4.6, lty=2, lwd=3, col="orangered4")
segments(x0=mean(littersize[temp.13.3]), y0=6, y1=2.95, lty=2, lwd=3, col="orangered1")
segments(x0=mean(littersize[temp.14.3]), y0=6, y1=1.65, lty=2, lwd=3, col="orange")
segments(x0=mean(littersize[temp.15.3]), y0=6, y1=1.65, lty=2, lwd=3, col="gold")
text(x=mean(littersize[temp.12.3]), y=5, labels="2012", col="orangered4", adj=-0.1)
text(x=mean(littersize[temp.13.3]), y=5, labels="2013", col="orangered1", adj=-0.1)
text(x=mean(littersize[temp.14.3]), y=5, labels="2014", col="orange", adj=-0.1)
text(x=mean(littersize[temp.15.3]), y=5, labels="2015", col="gold", adj=-0.1)
dev.off()

####
#Lambda for each year
plot(density(lambda[1:1031]), xlab="Lambda", lwd=3, col="orangered4", main="")
lines(density(lambda[1032:2062]), col="orangered1", lwd=3)
lines(density(lambda[2063:3093]), col="orange", lwd=3)
lines(density(lambda[3094:4124]), col="gold", lwd=3)
segments(x0=1, x1=1, y0=0, y1=3, lty=2, lwd=3)

###
#Combined fecundity and lambda for each year
setwd(FiguresDir)
filename=sprintf("Lambda_fert_%s.pdf",Sys.Date()) 
pdf(filename, width=7, height=5)
par(mfrow=c(2,1), mar=c(4,4,1,1), mgp=c(2.5,1,0))
plot(density(littersize[temp.12.3]), xlim=c(0,1.057), xlab="Recruitment (juveniles[F]/female/year)", lwd=3, col="orangered4", 
	main="", ylab="Probability density", ylim=c(0, 5.8))
lines(density(littersize[temp.13.3]), col="orangered1", lwd=3)
lines(density(littersize[temp.14.3]), col="orange", lwd=3)
lines(density(littersize[temp.15.3]), col="gold", lwd=3)
segments(x0=mean(littersize[temp.12.3]), y0=6, y1=4.6, lty=2, lwd=3, col="orangered4")
segments(x0=mean(littersize[temp.13.3]), y0=6, y1=2.95, lty=2, lwd=3, col="orangered1")
segments(x0=mean(littersize[temp.14.3]), y0=6, y1=1.65, lty=2, lwd=3, col="orange")
segments(x0=mean(littersize[temp.15.3]), y0=6, y1=1.65, lty=2, lwd=3, col="gold")
text(x=mean(littersize[temp.12.3]), y=5, labels="2012", col="orangered4", adj=-0.1)
text(x=mean(littersize[temp.13.3]), y=5, labels="2013", col="orangered1", adj=-0.1)
text(x=mean(littersize[temp.14.3]), y=5, labels="2014", col="orange", adj=-0.1)
text(x=mean(littersize[temp.15.3]), y=5, labels="2015", col="gold", adj=-0.1)
segments(x0=0.88, y0=0, x1=0.88, y1=1000, lty=3, lwd=2) #Recruitment needed to reach lambda=1
text(x=1.01, y=5.07, labels="(a)", cex=1.5)
plot(density(lambda[1:1031]), ylab="Probability density", xlab="Lambda", lwd=3, col="orangered4", main="", xlim=c(0,1.2), ylim=c(0,2.85))
lines(density(lambda[1032:2062]), col="orangered1", lwd=3)
lines(density(lambda[2063:3093]), col="orange", lwd=3)
lines(density(lambda[3094:4124]), col="gold", lwd=3)
segments(x0=mean(lambda[1:1031]), y0=6, y1=2.7, lty=2, lwd=3, col="orangered4")
segments(x0=mean(lambda[1032:2062]), y0=6, y1=2.4, lty=2, lwd=3, col="orangered1")
segments(x0=mean(lambda[2063:3093]), y0=6, y1=2, lty=2, lwd=3, col="orange")
segments(x0=mean(lambda[3094:4124]), y0=6, y1=2, lty=2, lwd=3, col="gold")
segments(x0=1, x1=1, y0=0, y1=3, lty=3, lwd=2)
text(x=1.15, y=2.5, labels="(b)", cex=1.5)
dev.off()

###
#Combined fecundity and lambda for each year, black and white
setwd(FiguresDir)
filename=sprintf("Lambda_fert_bw_%s.pdf",Sys.Date()) 
pdf(filename, width=7, height=5)
par(mfrow=c(2,1), mar=c(4,4,1,1), mgp=c(2.5,1,0))
plot(density(littersize[temp.12.3]), xlim=c(0,1.057), xlab="Recruitment (juveniles[F]/female/year)", lwd=3, col="black", 
	main="", ylab="Probability density", ylim=c(0, 5.8))
lines(density(littersize[temp.13.3]), col=gray(0.3), lwd=3)
lines(density(littersize[temp.14.3]), col=gray(0.6), lwd=3)
lines(density(littersize[temp.15.3]), col=gray(0.85), lwd=3)
segments(x0=mean(littersize[temp.12.3]), y0=6, y1=4.6, lty=2, lwd=3, col="black")
segments(x0=mean(littersize[temp.13.3]), y0=6, y1=2.95, lty=2, lwd=3, col=gray(0.3))
segments(x0=mean(littersize[temp.14.3]), y0=6, y1=1.6, lty=2, lwd=3, col=gray(0.6))
segments(x0=mean(littersize[temp.15.3]), y0=6, y1=1.65, lty=2, lwd=3, col=gray(0.85))
text(x=mean(littersize[temp.12.3]), y=5, labels="2012", col="black", adj=-0.1)
text(x=mean(littersize[temp.13.3]), y=5, labels="2013", col=gray(0.3), adj=-0.1)
text(x=mean(littersize[temp.14.3]), y=5, labels="2014", col=gray(0.6), adj=-0.1)
text(x=mean(littersize[temp.15.3]), y=5, labels="2015", col=gray(0.85), adj=-0.1)
segments(x0=0.88, y0=0, x1=0.88, y1=1000, lty=3, lwd=2) #Recruitment needed to reach lambda=1
text(x=1.01, y=5.07, labels="(a)", cex=1.5)
plot(density(lambda[1:1031]), ylab="Probability density", xlab="Lambda", lwd=3, col="black", main="", xlim=c(0,1.2), ylim=c(0,2.85))
lines(density(lambda[1032:2062]), col=gray(0.3), lwd=3)
lines(density(lambda[2063:3093]), col=gray(0.6), lwd=3)
lines(density(lambda[3094:4124]), col=gray(0.85), lwd=3)
segments(x0=mean(lambda[1:1031]), y0=6, y1=2.7, lty=2, lwd=3, col="black")
segments(x0=mean(lambda[1032:2062]), y0=6, y1=2.4, lty=2, lwd=3, col=gray(0.3))
segments(x0=mean(lambda[2063:3093]), y0=6, y1=2, lty=2, lwd=3, col=gray(0.6))
segments(x0=mean(lambda[3094:4124]), y0=6, y1=2, lty=2, lwd=3, col=gray(0.85))
segments(x0=1, x1=1, y0=0, y1=3, lty=3, lwd=2)
text(x=1.15, y=2.5, labels="(b)", cex=1.5)
dev.off()

#############################
#####  write parameter estimates to file

names(Mod$sims.list)

tempdf = data.frame()   

#"phi0",      # mean probability of survival
#"pbadyr",  # probability that a sampling period has low survival (due to predation?)
#"badyr",   # probability that a given season is classified as a bad year...
#"precipEff",   # effect of precip on survival
#"maleEff",     # effect of maleness on survival
#"juvEff",      # effect of being juvenile on survival
#"predEff",     # reduction in survival rate in a bad year
#"p0",            # mean nightly prob. of capture
#"gamma.prime",    # probability of staying off site
#"gamma.dprime"     # probability of leaving from on site


cphi0 <- c(mean(Mod$sims.list$phi0),quantile(Mod$sims.list$phi0,0.025),quantile(Mod$sims.list$phi0,0.975))
cpbadyr <- c(mean(Mod$sims.list$pbadyr),quantile(Mod$sims.list$pbadyr,0.025),quantile(Mod$sims.list$pbadyr,0.975))
cprecipeff <- c(mean(Mod$sims.list$precipEff),quantile(Mod$sims.list$precipEff,0.025),quantile(Mod$sims.list$precipEff,0.975))
#cmaleeff   <- c(mean(Mod$sims.list$maleEff),quantile(Mod$sims.list$maleEff,0.025),quantile(Mod$sims.list$maleEff,0.975))
#cjuveff  <- c(mean(Mod$sims.list$juvEff),quantile(Mod$sims.list$juvEff,0.025),quantile(Mod$sims.list$juvEff,0.975))
cpredeff  <- c(mean(Mod$sims.list$predEff),quantile(Mod$sims.list$predEff,0.025),quantile(Mod$sims.list$predEff,0.975))
cp0   <- c(mean(Mod$sims.list$p0),quantile(Mod$sims.list$p0,0.025),quantile(Mod$sims.list$p0,0.975))
cgammaprime  <- c(mean(Mod$sims.list$gamma.prime),quantile(Mod$sims.list$gamma.prime,0.025),quantile(Mod$sims.list$gamma.prime,0.975))
cgammadprime  <- c(mean(Mod$sims.list$gamma.dprime),quantile(Mod$sims.list$gamma.dprime,0.025),quantile(Mod$sims.list$gamma.dprime,0.975))


paramtable1 <- round(rbind(cphi0,cpbadyr,cprecipeff,cpredeff,cp0,cgammaprime,cgammadprime) ,3)  # cmaleeff,cjuveff,
colnames(paramtable1) <- c("mean","lower limit 95% CI","upper limit 95% CI")
rownames(paramtable1) <- c("mean survival","probability of bad year","effect of precipitation on survival (logit)",
                           #"effect of juvenileness on survival (logit)",
                           "reduced survival in bad year (logit)","mean capture probability","probability of staying off sampled region",
                           "probability of leaving sampled region")   # "effect of maleness on survival (logit)",

setwd(ResultsDir)
tempname <- paste("Parameter_Table1_",Sys.Date(),".csv",sep="")
write.csv(paramtable1,tempname)

setwd(BUGSDir)
tempname <- paste("LatestWorkspace_",Sys.Date(),".RData",sep="")
save.image(file = tempname)








##############  FOR SIGMA PLOT...

## Plot out the survival estimates as a function of precipitation

graphics.off()

setwd(FiguresDir)
filename=sprintf("Surv_Precip_%s.svg",Sys.Date())
svg(file=filename,height=5,width=5)
phi.mean <- numeric(length(xvars))
phi.97.5 <- numeric(length(xvars))
phi.02.5 <- numeric(length(xvars))
t=1
for(t in 1:length(xvars)){
  phi.mean[t] <- mean(precip_matrix[,t])
  phi.97.5[t] <- quantile(precip_matrix[,t],.975)
  phi.02.5[t] <- quantile(precip_matrix[,t],.025)
}

par(mai=c(2,1,1,1))
plot(c(1:length(xvars)),phi.mean,
     main="Survival vs Precip", ylab="Survival rate", xlab="",type="p",
     ylim=c(0,1),xaxt="n")
errbar(x=c(1:length(xvars)),y=phi.mean,yplus=phi.97.5,yminus=phi.02.5,add=T,pch="")

#(log(precip)-meanlogprecip)/sdlogprecip
xlabels <- round(exp((xvars*sdlogprecip)+meanlogprecip),1)
ats <- seq(1,length(xvars),5)
axis(1,at=ats,labels=xlabels[ats],las=2)
dev.off()

###### 
#### make CSV for Sigma plot

setwd(ResultsDir)
tempname <- paste("Surv_Precip_",Sys.Date(),".csv",sep="")
write.table(data.frame(precip=xlabels,mean=phi.mean,upper=phi.97.5,lower=phi.02.5),tempname,sep=",",row.names=F)



###########################################################
###################################  
#                     ANA, you can run the code up to this point
###################
###########
#########
######
#####
####




###################################
#   Run Convergence diagnostics (Gelman-Rubin diagnostic).

#chains <- list()
#chains[[1]] <- read.coda("CCmodel_coda1.txt", "CCmodel_index1.txt")  #,start, end, thin, quiet=FALSE)
#chains[[2]] <- read.coda("CCmodel_coda2.txt", "CCmodel_index1.txt")
#chains[[3]] <- read.coda("CCmodel_coda3.txt", "CCmodel_index1.txt")

#MCMC <- mcmc.list(chains[[1]][,-c(61:69)],chains[[2]][,-c(61:69)],chains[[3]][,-c(61:69)])


#GelmanRubinTest <- gelman.diag(MCMC,transform=T)
#gelman.plot(MCMC)

#  some of the monitored variables have problems with this test, even after transformation...
# this worked after I got rid of the  







############################
########  SUMMARIZE STRANGE DATA

#ndx <- which(((caphist2[,1]==1)|(caphist2[,2]==1))&(apply(caphist2[,3:12],1,sum)>0))

#newdf <- data.frame(names[ndx],caphist2[ndx,])

#names(newdf) <- c("Left Ear Tag ID", realperiodsG) 

#newdf

#getwd()
#write.csv(newdf,"Redflag_survivors.csv",row.names=F)






#names(Mod$sims.list)



#hist(Mod$sims.list$p0)


#table(pdogT$PLOT,pdogT$PTP)







######################
####### ARCHIVE CODE


#####################
### plot out the effect of precip on survival...

nsims=5000 # temporary
meanlogprecip <- mean(log(precip))             # standardize precip variable
sdlogprecip <- sd(log(precip))
precip_std <- (log(precip)-meanlogprecip)/sdlogprecip
xvars <- seq(min(as.vector(precip_std)),max(as.vector(precip_std)),0.05)

nsamples <- 1000
precip_matrix <- matrix(nrow=nsamples,ncol=length(xvars))
for(s in 1:nsamples){
  samp <- sample(c(1:nsims),1,replace=F)
  for(i in 1:length(xvars)){
    temp <- plogis(qlogis(Mod$sims.list$phi0[samp]) + Mod$sims.list$precipEff[samp]*xvars[i])
    precip_matrix[s,i] <- temp
  }
}

boxplot(precip_matrix)


### determine which year is a "bad year"

badsurv_p <- numeric(nperiodsG)
badsurv_b <- numeric(nperiodsG)

#  i=1
for(i in 1:nperiodsG){
  ndx <- which(names(table(Mod$sims.list$badyr[,i]))=="1")
  if(length(ndx)>0) temp <- as.numeric(table(Mod$sims.list$badyr[,i]))[ndx]/nsims else temp=1
  badsurv_p[i] <- temp        #probability of bad survival this period  
  badsurv_b[i] <- ifelse(temp>0.75,1,0)  # this period designated as bad survival if designated bad more than 75% of the time
}
# write this out to file...
newdf <- data.frame(period=realperiodsG,prob_highpred_period=badsurv_p,is_highpred_period=badsurv_b)
tempname <- paste("highPredation_",Sys.Date(),".csv",sep="")

setwd(ResultsDir)
write.table(newdf, file=tempname, sep=",",row.names=F)














