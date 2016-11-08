#################### 
#  R script for running pdog *fecundity* analysis for large plots at Sevilleta NWR
#
rm(list=ls())

KEVIN = F
KEVINOFFICE = F
ELIZABETH = F
ELIZABETHOFFICE = T

###################################
##########  SET WORKING DIRECTORY

if(KEVIN) rootDir <- "C:\\Users\\Kevin\\Dropbox\\Sev Pdog Study\\BIG plot study\\"
if(KEVINOFFICE) rootDir <- "E:\\Dropbox\\Sev Pdog Study\\BIG plot study\\"
if(ELIZABETH) rootDir <- "~/Dropbox/Nevada/Sev Pdog Study/BIG plot study/"
if(ELIZABETHOFFICE) rootDir <- "C:\\Users\\ehunter\\Dropbox\\Nevada\\Sevilleta\\Sev Pdog Study\\BIG plot study\\"

ScriptDir <- paste(rootDir,"Rscript",sep="")
DataDir <- paste(rootDir,"Data",sep="")
FiguresDir <- paste(rootDir,"RawFigures",sep="")
ResultsDir <- paste(rootDir,"Results",sep="")
BUGSDir <- paste(rootDir,"BUGS",sep="")

setwd(DataDir)

getwd()

  # setwd(BUGSDir)
  # load(".RData")

##################################
##########  LOAD PACKAGES

library(lubridate)
library(Hmisc)
library(vcdExtra)

###################################
##########   READ IN DATA

#Read in posterior distributions
setwd(BUGSDir)
load("SevModel10_1mnth_newprecip_2016-07-16.RData") 
names(Mod$sims.list)
hist(Mod$sims.list$phi0)

###############Everything else - taken from BigPlots_script5, unnecessary parts removed
setwd(DataDir)
pdog <- read.csv("pdog_captures_16ha_9-12-15 AD_KTS_2016.csv", header=TRUE)  ## KTS: updated with 2015 data and fixed a couple errors in the data
## EAH: updated with 2016 data through June

nCaptures <- nrow(pdog)

############
#  Process the capture data 

names(pdog)
pdog$DATE2 <- mdy(pdog$DATE)
pdog$PLOT <- as.character(pdog$PLOT)   
pdog$PTP <- paste(pdog$SEASON,pdog$YEAR,sep="")   # Primary trapping period

unique(pdog$PTP)

pdog <- subset(pdog,!PLOT=="G")    # remove plot G for now (some pdogs released but not surveyed as of 2015)

pdog <- subset(pdog,!TRAPDAY%in%c(4,5))    # remove observations on "experimental" fourth and fifth capture days

#pdog$TRAPDAY

###############################
####### DEFINE "PERIOD" FOR THIS STUDY

# Note: some plots were trapped at different times within the same season- so trap day 1 isn't necessarily always the same date...

setwd(DataDir)
filename <- "PeriodDefinition.csv"
PeriodDef <- read.csv(filename)

PeriodDef$StartDate <- mdy(PeriodDef$Start)
PeriodDef$EndDate <- mdy(PeriodDef$End)
#PeriodDef$Interval <- new_interval(PeriodDef$StartDate, PeriodDef$EndDate)

nPeriods <- nrow(PeriodDef)

  #####  DETERMINE PERIOD BETWEEN EACH PRIMARY SURVEY (periods among which population is open)

 ## convert to datetime object
pdate <- mdy(PeriodDef$Median)

 ## 
counter=1
intervals <- numeric(nPeriods-1)
for(i in 2:nPeriods){
  intervals[counter]=as.numeric(difftime(pdate[i],pdate[i-1],units="days"))/365
  counter =counter+1
}

PeriodDef$Intervals <- c(intervals,NA)


###############################
#######  READ IN ENV CONDITIONS- PRECIP AND SOIL MOISTURE
#Using the precip interval that would affect that period
#For summer periods: May-Oct
#For spring periods: Nov-Apr

setwd(DataDir)
PrecipPeriodDef <- read.csv(filename)

PrecipPeriodDef$StartDate <- mdy(PrecipPeriodDef$Start)
PrecipPeriodDef$EndDate <- mdy(PrecipPeriodDef$End)

envCond_df <- read.csv("Month_sums_50 -sev 2010-2016_ppt.csv",header=T) #Now includes 2016 through May
names(envCond_df)

  ## add a "Date" field
temp <- paste(envCond_df$mon,1,envCond_df$year,sep="/")
envCond_df$Date <- mdy(temp)

#########
#  summarize environmental conditions by period
##JUST PRECIP and SOILH2O
envCond_df <- envCond_df[,c(1,2,6,12,13)]
envCondNames <- names(envCond_df)

newline <- envCond_df[1,-which(envCondNames%in%c("Date","mon"))]
names(newline)[1] <- "Period"

envCond <- newline[-1,]

p=1
for(p in 1:nPeriods){
  ndx <- which((envCond_df$Date>=PrecipPeriodDef$StartDate[p])&(envCond_df$Date<PrecipPeriodDef$EndDate[p]))
  temp <- newline
  temp[1,1] <- as.character(PrecipPeriodDef$Period[p])
  i=2
  for(i in 2:ncol(temp)){
    colname <- names(temp)[i]
    ndx2 <- which(names(envCond_df)==colname)
    temp[1,i] <- mean(envCond_df[ndx,ndx2])
  }
  envCond <- rbind(envCond,temp)
}

precip <- envCond$ppt

       ## and soil moisture data...

soilH2O <- envCond$sH2O



###################################
   ####        DEVELOP UNIQUE ID VAR ("indiv")
   ####  KTS 11/17/15: made this more thorough- now record number of tags and type of tags for each individual


pdog$TagR <- as.character(pdog$RIGHTTAG)      
pdog$TagL <- as.character(pdog$LEFTTAG)
pdog$PIT.TAG <- as.character(pdog$PITTAG)

ndx <- which(pdog$TagR=="")                # change "" to NA
if(length(ndx)>0) pdog$TagR[ndx] <- NA 
ndx <- which(pdog$TagL=="")                # change "" to NA
if(length(ndx)>0) pdog$TagL[ndx] <- NA 
ndx <- which(pdog$PIT.TAG=="")                # change "" to NA
if(length(ndx)>0) pdog$PIT.TAG[ndx] <- NA 


### new strategy: 
       # first make a list with all single tags 
       # then loop through this vector
           # IF this ID was ever observed on the same individual as another ID later in the list, 
           # THEN remove that other ID from the list and associate it with the ID in question 
             
### NOTE: this takes a long time, totally inefficient, but whatever! 

uniqueR <- unique(pdog$TagR)
uniqueR <- uniqueR[!is.na(uniqueR)]
uniqueL <- unique(pdog$TagL)
uniqueL <- uniqueL[!is.na(uniqueL)]
uniquePIT <- unique(pdog$PIT.TAG)
uniquePIT <- uniquePIT[!is.na(uniquePIT)]

freeIDs <- c(uniqueL,uniqueR,uniquePIT)  # all unique ID tags
usedIDs <- numeric(0)

indivList <- list()
IDarray <- with(pdog,cbind(TagL,TagR,PIT.TAG) )

counter=1
t=1
for(t in 1:length(freeIDs)){      #loop through all ID tags
  thisID <- freeIDs[t]
  if(!thisID%in%usedIDs){   # if this tag hasn't been used before
    indivList[[counter]] <- thisID    # this must be a new individual 
    usedIDs <- c(usedIDs,thisID)
    temp <- IDarray[apply(IDarray,1,function(t) thisID%in%t),]   # observations involving this tag
    temp <- unique(as.vector(temp))
    temp <- temp[!is.na(temp)]
    temp <- temp[temp!=thisID]
    if(length(temp)>0){
      indivList[[counter]] <- c(indivList[[counter]],temp)  # add all other associated IDs
      usedIDs <- c(usedIDs,temp)

      while(length(temp)>0){  # go through the other IDs- make sure that these IDs are not associated with other IDs
        thisID2 <- temp[1]     #take the first element of temp
        temp2 <- IDarray[apply(IDarray,1,function(t) thisID2%in%t),]
        temp2 <- unique(as.vector(temp2))
        temp2 <- temp2[!is.na(temp2)]
        temp2 <- temp2[!temp2%in%usedIDs]
        indivList[[counter]] <- c(indivList[[counter]],temp2)    # add any new ID tags associated with this individual
        usedIDs <- c(usedIDs,temp2)
        temp<-temp[-1]
      }  
    }  
    counter <- counter + 1 
  }
}


 #####  Good- this seems to work...

#######################################
### NOW: make a new "indiv" field... 

### NOTE: again, this takes a long time, but F#@# it :) 

pdog$indiv <- 0

nan <- length(indivList)

nrow(IDarray)     # temp
nrow(pdog)


errcheck <- numeric(0)   # which individual IDs, if any, are associated with multiple individuals (error check)
noobsndx <- numeric(0)   # which observations are associated with no ID at all??

i=32
for(i in 1:nrow(pdog)){
  tempID <- IDarray[i,]
  tempID <- unique(as.vector(tempID))
  tempID <- tempID[!is.na(tempID)]
  tempID <- tempID[1]  # just pick the first of the IDs...

  temp2 <- (1:nan)[sapply(indivList,function(t) tempID%in%t)]   # which individual is this??

  if(length(temp2)>1) errcheck <- c(errcheck,i)   # if this corresponds to more than one individual, there is a problem!!
  if(length(temp2)==0) noobsndx <- c(noobsndx,i)
  pdog$indiv[i] <- ifelse(length(temp2)==1, temp2, NA)  
}

errcheck   # Good- no ambiguous identifications

noobsndx    # 9 captures with no IDs 

pdog[noobsndx,]   # checks out

##### NOW: remove any individuals that had no ID

pdog <- pdog[-noobsndx,]

pdog$indiv    # looks OK


#### NOW: generate capture histories as a list- one dataframe of observations per individual (CHList)

uniqueNames <- unique(pdog$indiv)
nan <- length(uniqueNames)
i=1122
CHList <- list()   
for(i in 1:nan){        # remember to sort by date
  eval(parse(text=sprintf("CHList$indiv%s <- subset(pdog,indiv==uniqueNames[i])",uniqueNames[i])))
  eval(parse(text=sprintf("CHList$indiv%s <- CHList$indiv%s[order(CHList$indiv%s$DATE2),]",uniqueNames[i],uniqueNames[i],uniqueNames[i])))
}

  
####################################
##########  SUMMARY DATA (global: suffix "G")

names(pdog)

nobsG <- nrow(pdog)       # 3272 observations

nindG <- length(unique(pdog$indiv))     # 2505 individuals 

realindG <- unique(pdog$indiv)  

nyearsG <- length(unique(pdog$YEAR))     # 6 years

realyearsG <- sort(unique(pdog$YEAR))

nperiodsG <- length(unique(pdog$PTP))    # 10 primary trapping periods

realperiodsG <- unique(pdog$PTP)

nplotsG <- length(unique(pdog$PLOT))     # 3 plots

realplotsG <- as.character(unique(pdog$PLOT)) 

ndatesG <- length(unique(pdog$DATE))     # 73 survey occasions / release events

realsurveydatesG <- format(sort(unique(pdog$DATE2[which(!is.na(pdog$TRAPDAY))])),format="%m/%d/%Y")

nsurveysG <- length(realsurveydatesG)   # 38 trapping surveys

realreleasedatesG <- format(sort(unique(pdog$DATE2[which(is.na(pdog$TRAPDAY))])),format="%m/%d/%Y")

nreleasesG <- length(realreleasedatesG)   # 31 release dates

realdatesG <- format(sort(unique(pdog$DATE2)),format="%m/%d/%Y")

realperiodyearsG <- as.numeric(na.omit(as.numeric(unlist(strsplit(unlist(realperiodsG), "[^0-9]+")))))

   # look at num caps each year by plot

table(pdog$PLOT,pdog$YEAR)   

pdog$PTP2 <- factor(pdog$PTP,levels=realperiodsG)   # make the periods are ordered properly

captures <- subset(pdog,!RECAPTYPE=="RELEASE")  

table(captures$PLOT,captures$YEAR)

releases <- subset(pdog,RECAPTYPE=="RELEASE")

sumReleases <- table(releases$PLOT,releases$PTP2)

################################
############## SUMMARIZE DATA FOR WINBUGS

   ## loop through data and summarize the data for each dimension of interest

ns <- nPeriods
nss <- array(0,dim=c(nPeriods))   # number of surveys (suboccasions)
nan <- nindG   # number of animals 
names <- realindG  #  list()

y=4
for(y in 1:nPeriods){
  ndx <- which((pdog$PTP==realperiodsG[y])&(!is.na(pdog$TRAPDAY)))     
  nss[y] <- length(unique(pdog$TRAPDAY[ndx]))
}
nss[1] <- 1    # allow the first season to have a dummy trapping event (these are all releases... )

maxnan <- max(nan)
maxnss <- max(nss)      # max of 3 suboccasions... 

##########################
#######       MAKE CAPTURE HISTORY DATA
######## INITIALIZE DATA ARRAYS

caphist <- array(0,dim=c(maxnan,nperiodsG,maxnss))     # Master 3-D capture history array
caphist_2D <- array(0,dim=c(maxnan,nperiodsG))
NMarks <- array(NA,dim=c(maxnan,nperiodsG))     # choose the min number of marks in each PTP
LeftEarTag <- array(NA,dim=c(maxnan,nperiodsG))
RightEarTag <- array(NA,dim=c(maxnan,nperiodsG))

Released <- array(0,dim=c(maxnan,nperiodsG))    # if release happens after trapping, code as "1" in final suboccasion. If released before trapping, put "1" in initial suboccasion.
SpringRelease <- array(0,dim=c(maxnan,nperiodsG)) # NOTE: actually, releases ALWAYS happened after the final trapping event (this is good!)
 #Juvenile <- array(0,dim=c(maxnan,nperiodsG))
PITTag <- array(0,dim=c(maxnan,nperiodsG))     # if it has a pit tag, can't lose tags
AddLeftTag <- array(0,dim=c(maxnan,nperiodsG))
AddRightTag <- array(0,dim=c(maxnan,nperiodsG)) 

firsts <- numeric(maxnan)
first2 <- numeric(maxnan)

addTagNdx <- character(0)   # observations where tag was added...
  

ind=285#2510#2327#285#1122#1#
for(ind in 1:nan){
  df <- CHList[[ind]] 
              # make sure to sort by date!
  df <- df[order(df$DATE2),]
  tags <- with(df,cbind(TagL,TagR,PIT.TAG))
                                                         ### ensure that once it has a pit tag, it always has that pit tag... 
  pitpresent <- !is.na(tags[,3])   # is there a PIT tag?
  firstpit <- ifelse(any(pitpresent),which(pitpresent)[1],NA)
  if(!is.na(firstpit)) tags[min(nrow(df),(firstpit+1)):nrow(df),3] <- tags[firstpit,3]    # make sure that pit tag is not lost... since pit tags cant be lost in this model

  df$nMarks <-  apply(tags,1,function(t) length(t[!is.na(t)])) 
  
  obs <- 1
  for(obs in 1:nrow(df)){   # loop through by observation.. 
    thisObs <- df[obs,]
    thisTags <- tags[obs,]
    thisPeriod <- which(PeriodDef$Period==thisObs$PTP)
    thisSub <- ifelse(!is.na(thisObs$TRAPDAY),thisObs$TRAPDAY,max(nss[thisPeriod]))    # if it's just released, put it at the last suboccasion

    if(obs==1){    # if this is the first capture
      firsts[ind] <- thisPeriod            # record period of first capture
      first2[ind] <- thisSub
    }

    if((!is.na(firstpit))&(obs==firstpit)){  # if has pit and this is the first observation with pit tag
      PITTag[ind,thisPeriod:nPeriods] <- 1      # then we assume this individual has a pit tag forever, as long as it is alive. can't lose this mark
    }
      
    caphist[ind,thisPeriod,thisSub] <- 1                              # captured this suboccasion?
    caphist_2D[ind,thisPeriod] <- 1                                   # capture this period? (standard 2D capture history)

    if(thisObs$RECAPTYPE=="RELEASE"){
      Released[ind,thisPeriod] <- 1    # released this period?   
	  if(month(thisObs$DATE2)%in%c(4,5,6)) SpringRelease[ind,thisPeriod] <- 1                             # is it a spring release?
	}

    NMarks[ind,thisPeriod] <- df$nMarks[obs]  # number of tags observed on this individual
    LeftEarTag[ind,thisPeriod] <- ifelse(!is.na(tags[,"TagL"][obs]),1,0) 
    RightEarTag[ind,thisPeriod] <- ifelse(!is.na(tags[,"TagR"][obs]),1,0)   
    tagsadded <- (!is.na(thisTags))&((!thisTags==tags[max(1,obs-1),])|(is.na(tags[max(1,obs-1),])))
    AddLeftTag[ind,thisPeriod] <- ifelse(tagsadded[1],1,0)
    AddRightTag[ind,thisPeriod] <- ifelse(tagsadded[2],1,0)
    if(any(tagsadded)) addTagNdx <- c(addTagNdx,sprintf("ind%s_obs%s",ind,obs))
  }  
}

   # make sure once indiv has PIT, it is not lost... 
tail(PITTag,500)


##############################
##        STORE SEX AND AGE FOR EACH INDIVIDUAL

isMale <- numeric(nan)

for(i in 1:nan){
  ndx <- which(pdog$indiv==names[i])
  sexvec <- pdog$SEX[ndx]
  sexvec <- sexvec[which(!is.na(sexvec))]
  sex <- as.character(sexvec[length(sexvec)])    # final determined sex is the true sex...
  isMale[i] <- ifelse(length(sex>0),ifelse(sex=="M",1,0),0)
   #  pdogT[ndx,]
}

 
isJuv  <- array(1,dim=c(nan,nperiodsG))

i=1   #363
p=firsts[i]
for(i in 1:nan){
  firstyear <- realperiodyearsG[firsts[i]]
  for(p in firsts[i]:nperiodsG){
    if(realperiodyearsG[p]==firstyear){   # if within the first year of capture
      ndx <- which((pdog$indiv==names[i]) &
                 (pdog$YEAR==firstyear) )
      if(length(ndx)>0){
        tempage <- as.character(pdog$AGE[ndx])[1]
        isJuv[i,p] <- ifelse(tempage=="J",1,0)
      }
    } else{
      isJuv[i,p] <- 0
    } 
  }
}

isJuv <- apply(isJuv,c(1,2),function(t) ifelse(is.na(t),0,t))

tail(isJuv,500)


########################################
#########    standardize covariates

meanlogprecip <- mean(log(precip))             # standardize precip variable
sdlogprecip <- sd(log(precip))
precip_std <- (log(precip)-meanlogprecip)/sdlogprecip

meanlogsoil <- mean(log(soilH2O))             # standardize soil variable
sdlogsoil <- sd(log(soilH2O))
soil_std <- (log(soilH2O)-meanlogsoil)/sdlogsoil


######################
#Input data for simulation (releases):
#Collapse "releases" into table divided up by age, sex (just females), year, and season
ReleasedAge <- table(releases$AGE, releases$YEAR, releases$SEX, releases$SEASON, exclude=c("", "M"), dnn=c("Age", "Year", "Sex", "Season"))
ReleasedAge <- as.data.frame(collapse.table(ReleasedAge, Sex=c("F", "F", "F", "F", "F"), Age=c("A", "J", "J")))
ReleasedAge$Year <- as.numeric(as.character(ReleasedAge$Year))
ReleasedAge <- rbind(ReleasedAge, c(as.character("A"), as.numeric(2014), as.character("F"), as.character("SPRING"), as.numeric(0)), 
	c(as.character("J"), as.numeric(2014), as.character("F"), as.character("SPRING"), as.numeric(0)), 
	c(as.character("A"), as.numeric(2014), as.character("F"), as.character("SUMMER"), as.numeric(0)), 
	c(as.character("J"), as.numeric(2014), as.character("F"), as.character("SUMMER"), as.numeric(0)),
	c(as.character("A"), as.numeric(2016), as.character("F"), as.character("SPRING"), as.numeric(0)), 
	c(as.character("J"), as.numeric(2016), as.character("F"), as.character("SPRING"), as.numeric(0)))
ReleasedAge$Year <- as.numeric(ReleasedAge$Year); ReleasedAge$Freq <- as.numeric(ReleasedAge$Freq)
ReleasedAge <- ReleasedAge[order(ReleasedAge$Year),]
names(ReleasedAge) <- c("Age", "Year", "Sex", "Season", "N")

#Need "ReleasedSpring" and "ReleasedSummer" dataframes, separating out the females released in spring in 2012 and 2013
ReleasedSpring <- matrix(c(seq(2010,2016, length=nyearsG), rep(0,times=nyearsG)), nrow=nyearsG, ncol=2)
ReleasedSpring <- as.data.frame(ReleasedSpring); names(ReleasedSpring) <- c("Year", "N")
ReleasedSpring$N[ReleasedSpring$Year==2012] <- ReleasedAge$N[ReleasedAge$Year==2012 & ReleasedAge$Season=="SPRING" & ReleasedAge$Age=="A"]
ReleasedSpring$N[ReleasedSpring$Year==2013] <- ReleasedAge$N[ReleasedAge$Year==2013 & ReleasedAge$Season=="SPRING" & ReleasedAge$Age=="A"]
ReleasedSpring <- as.vector(ReleasedSpring$N)

#Remove spring releases from ReleasedSummer, condense to single season, break into A/J
temp <- ReleasedAge
temp$N[temp$Year==2012 & temp$Season=="SPRING" & temp$Age=="A"] <- 0
temp$N[temp$Year==2013 & temp$Season=="SPRING" & temp$Age=="A"] <- 0
temp <- aggregate(temp$N, by=list(temp$Year, temp$Age), FUN=sum)

ReleasedSummer <- array(0, dim=c(nyearsG, 2)) #2 ages
ReleasedSummer[,1] <- temp[8:14, 3] #Juveniles
ReleasedSummer[,2] <- temp[1:7, 3]  #Adults

##########################
#Matching projections to proportion "nativo"
#Getting "nativos" from pdog dataframe (N in RECAPTYPE)
pdognativos <- pdog[pdog$RECAPTYPE=="N",] #104 individuals, but there are 193 individuals with no release date in "Released" dataframe?
#"Released" has tag losses, so 104 is the more accurate count.
pdognativost <- table(pdognativos$SEX, pdognativos$YEAR, dnn=c("Sex", "Year"))
pdognativost <- collapse.table(pdognativost, Sex=c("F", "F", "F", "F", "F", "F", "M"))

#Summary statistic: slope of proportion of nativos in the population over time
#Need total number of observed pdogs for each year
temp <- table(pdog$indiv, pdog$YEAR, pdog$SEX, dnn=c("indiv", "YEAR", "SEX"))
temp <- collapse.table(temp, SEX=c("F", "F", "F", "F", "F", "F", "M"))
temp2 <- ifelse(temp>0, 1, 0)
totalobserved <- t(colSums(temp2))
totalobserved <- totalobserved[,3:7]

obsfracnativos <- as.data.frame(pdognativost/totalobserved)
obsfracnativos$Year <- as.numeric(as.character(obsfracnativos$Year))
plot(obsfracnativos$Freq ~ obsfracnativos$Year, col=obsfracnativos$Sex, type="p", pch=19, xaxp=c(2012,2015,3), xlab="Year", ylab="Fraction Nativos")
legend("bottomright", c("Female", "Male"), col=c("black", "red"), pch=19)
obsfracnativosF <- obsfracnativos[obsfracnativos$Sex=="F",]  #Just use females for pattern matching
obsfracnativosF <- obsfracnativosF[1:4,] #Just 2012-2015
obsfracnativosF <- rbind(obsfracnativosF, c("F", as.numeric(2011), as.numeric(0))) #Add in 2011, which had 0 nativos
obsfracnativosF$Year <- as.numeric(obsfracnativosF$Year); obsfracnativosF$Freq <- as.numeric(obsfracnativosF$Freq)
obsfracnativosF <- obsfracnativosF[order(obsfracnativosF$Year),]
#obsslope <- summary(lm(obsfracnativosF$Freq ~ obsfracnativosF$Year))$coefficients[2,1]
#This slope indicates number of NEW nativos per year (either juvenile or adult), once they are tagged, no longer a "nativo"

#######################
###Functions:
nMCMC <- length(Mod$sims.list$phi0)

basesurvival <- function(draw) {
	eval(Mod$sims.list$phi0)[draw]
	}

captureprob <- function(draw) {
	eval(Mod$sims.list$p0)[draw]
	}
	
logit <- function(p) {
	log(p/(1-p))
	}
	
precipB <- function(draw) {
	eval(Mod$sims.list$precipEff)[draw]
	}
	
soilB <- function(draw) {
	eval(Mod$sims.list$soilEff)[draw]
	}
	
juvB <- function(draw) {
	eval(Mod$sims.list$juvEff)[draw]
	}
	
relB <- function(draw) {
	eval(Mod$sims.list$relEff)[draw]
	}
	
juvRelB <- function(draw) {
	eval(Mod$sims.list$juvRelEff)[draw]
	}
	
springRelB <- function(draw) {
	eval(Mod$sims.list$springRelEff)[draw]
	}

juvsurvival <- function(draw, periods) {
	temp <- logit(basesurvival(draw)) + precipB(draw)*precip_std[periods] + soilB(draw)*soil_std[periods] + juvB(draw)
	return(1/(1 + exp(-1*temp)))
	}  
	
adultsurvival <- function(draw, periods) {
	temp <- logit(basesurvival(draw)) + precipB(draw)*precip_std[periods] + soilB(draw)*soil_std[periods]
	return(1/(1 + exp(-1*temp)))
	}

survivalprob <- function(draw, periods) {
	temp <- c(juvsurvival(draw, periods), adultsurvival(draw, periods))
	return(temp)
	}

Reljuvsurvival <- function(draw, periods, Season) {
	temp <- logit(basesurvival(draw)) + juvB(draw) + relB(draw) + 
		 springRelB(draw)*Season   #+ juvRelB(draw) + precipB(draw)*precip_std[periods] + soilB(draw)*soil_std[periods] 
	return(1/(1 + exp(-1*temp)))
	} 
	
Reladultsurvival <- function(draw, periods, Season) {
	temp <- logit(basesurvival(draw)) + relB(draw) + 
		springRelB(draw)*Season   #+ precipB(draw)*precip_std[periods] + soilB(draw)*soil_std[periods] 
	return(1/(1 + exp(-1*temp)))
	} 

Relsurvivalprob <- function(draw, periods, Season) {
	temp <- c(Reljuvsurvival(draw, periods, Season), Reladultsurvival(draw, periods, Season))
	return(temp)
	}
	#Season=1 is spring

##############################
#SIMULATION SCENARIOS

nscenarios <- 500
nsims <- 10000
nages <- 2
nperiods <- c(1, 1, 2, 2, 2, 2) #Number of periods in each year

#Make littersize vector to draw from
littersize <- seq(0.01, 2.5, length=nscenarios) #Only a single litter per year, includes underground pup survival rates, prob of pregnancy

ff_samp <-sample(seq(0.4, 0.6, length=100), nMCMC, replace=T) #Need to check this range
fracfem <- function(draw) ff_samp[draw]
	
fertility <- function(draw) {
	temp <- littersize[s] * fracfem(draw)
	return(temp)
	}

###############################
#START SIMULATION

N <- array(0, dim=c(nscenarios, nsims, nyearsG))
NSpring <- array(0, dim=c(nscenarios, nsims, nyearsG))
RelSpring <- array(0, dim=c(nscenarios, nsims, nyearsG))
RelSpringTemp <- array(0, dim=c(nscenarios, nsims, nyearsG))
RelSummer <- array(0, dim=c(nscenarios, nsims, nages, nyearsG))
RelSummerTemp <- array(0, dim=c(nscenarios, nsims, nages, nyearsG))
Nativo <- array(0, dim=c(nscenarios, nsims, nyearsG))
NativoSpring <- array(0, dim=c(nscenarios, nsims, nyearsG))
fracnativo <- array(0, dim=c(nscenarios,nsims, nyearsG))
TotN <- array(0, dim=c(nscenarios, nsims, nyearsG))
CapturedNativos <- array(0, dim=c(nscenarios, nsims, nyearsG))
releaseEff_interval <- 1/12 #Release effect only lasts for 1 months
MCMCsample <- array(0, dim=c(nscenarios, nsims))

for (s in 1:nscenarios){
for (i in 1:nsims){
	randnum <- sample(c(1:nMCMC),1)
	MCMCsample[s,i] <- randnum
	periods <- 1
for (y in 1:(nyearsG-1)){  
	#SPRING
		#Releases:
			rel_interval <- min(releaseEff_interval, intervals[periods]) #Release effect is 2 months or shorter, actual interval
			Relsurv <- Relsurvivalprob(draw=randnum, periods=periods, Season=1)^rel_interval
			RelSpringTemp[s,i,y] <- rbinom(1, ReleasedSpring[y], Relsurv[2])
			#Released individuals then survive the rest of the interval at the "normal" rate
			surv_interval <- intervals[periods] - rel_interval
			surv <- survivalprob(draw=randnum, periods=periods)^surv_interval 
			RelSpring[s,i,y] <- rbinom(1, RelSpringTemp[s,i,y], surv[2])
		#Survival (of non-released population):
			if (y > 1) {  #This only occurs after the first year
			surv <- survivalprob(draw=randnum, periods=periods)^intervals[periods]  
			NSpring[s,i,y] <- sum(rbinom(1, N[s,i,y-1], surv[2]),
						RelSpring[s,i,y],  #Add in released, survived adults
						round(rbinom(1, Nativo[s,i,y-1], surv[2])*(captureprob(randnum)))) #Captured nativos are marked
			NativoSpring[s,i,y] <- round(rbinom(1, Nativo[s,i,y-1], surv[2])*(1-captureprob(randnum))) 
			}
		if (nperiods[y] == 2){periods <- periods + 1}
			

	#REPRODUCTION
		#Adult N (and nativos) produce offspring as a function of fecundity parameters 
			Noffspring <- round(sum(NSpring[s,i,y], NativoSpring[s,i,y]) * fertility(randnum))

	#SUMMER
		#Releases:
		for (a in 1:nages){
			rel_interval <- min(releaseEff_interval, intervals[periods]) #Release effect is 2 months or shorter, actual interval
			Relsurv <- Relsurvivalprob(draw=randnum, periods=periods, Season=0)^intervals[periods]
			RelSummerTemp[s,i,a,y] <- rbinom(1, ReleasedSummer[y,a], Relsurv[a])
			#Released individuals then survive the rest of the interval at the "normal" rate
			surv_interval <- intervals[periods] - rel_interval
			surv <- survivalprob(draw=randnum, periods=periods)^surv_interval 
			RelSummer[s,i,a,y] <- rbinom(1, RelSummerTemp[s,i,a,y], surv[a])
		}
		#Survival:
			surv <- survivalprob(draw=randnum, periods=periods)^intervals[periods] 
			N[s,i,y] <- sum(rbinom(1, NSpring[s,i,y], surv[2]), #Survival of adult to adult
						sum(RelSummer[s,i,,y]), #Add in the survived, released summer females (summed juv and adult)
						round(rbinom(1, Noffspring, surv[1])*(captureprob(randnum))), #Captured nativos are marked
						round(rbinom(1, NativoSpring[s,i,y], surv[2])*(captureprob(randnum)))) 
		#Female nativos survive, and those that are not captured stay "nativo"
			Nativo[s,i,y] <- sum(round(rbinom(1, Noffspring, surv[1])*(1-captureprob(randnum))), #Survival of juv to adult 
						round(rbinom(1, NativoSpring[s,i,y], surv[2])*(1-captureprob(randnum)))) #Survival of adult to adult 
		#Captured nativos for pattern matching
			CapturedNativos[s,i,y] <- sum(round(rbinom(1, Nativo[s,i,y], surv[2])*(captureprob(randnum))),
						round(rbinom(1, Noffspring, surv[1])*(captureprob(randnum))))
		periods <- periods + 1

	#Get the fraction nativo for each year (need to sum captured and uncaptured nativos) out of total population size
		fracnativo[s,i,y] <- sum(Nativo[s,i,y], CapturedNativos[s,i,y]) / sum(N[s,i,y], Nativo[s,i,y])
		TotN[s,i,y] <- sum(N[s,i,y], Nativo[s,i,y])
}
}
}

#Find the fracnativo's that are within 5% of obsfracnativosF for each year 2012-2015
temp.12 <- fracnativo[,,3] <= (0.05*obsfracnativosF$Freq[2] + obsfracnativosF$Freq[2]) & fracnativo[,,3] >= (obsfracnativosF$Freq[2] - 0.05*obsfracnativosF$Freq[2]) 
temp.13 <- fracnativo[,,4] <= (0.05*obsfracnativosF$Freq[3] + obsfracnativosF$Freq[3]) & fracnativo[,,4] >= (obsfracnativosF$Freq[3] - 0.05*obsfracnativosF$Freq[3]) 
temp.14 <- fracnativo[,,5] <= (0.05*obsfracnativosF$Freq[4] + obsfracnativosF$Freq[4]) & fracnativo[,,5] >= (obsfracnativosF$Freq[4] - 0.05*obsfracnativosF$Freq[4]) 
temp.15 <- fracnativo[,,6] <= (0.05*obsfracnativosF$Freq[5] + obsfracnativosF$Freq[5]) & fracnativo[,,6] >= (obsfracnativosF$Freq[5] - 0.05*obsfracnativosF$Freq[5]) 
temp.rows <- row(fracnativo[,,3])
temp.12.2 <- temp.rows[which(temp.12==TRUE)]
temp.13.2 <- temp.rows[which(temp.13==TRUE)]
temp.14.2 <- temp.rows[which(temp.14==TRUE)]
temp.15.2 <- temp.rows[which(temp.15==TRUE)]
#Sample from each year (so that no year is overrepresented).  
#Doing this because practically impossible to get a single fertility that matches all three years.
samps <- min(c(length(temp.12.2), length(temp.13.2), length(temp.14.2), length(temp.15.2)))
samp.12 <- sample(1:length(temp.12.2), size=samps)
temp.12.3 <- temp.12.2[samp.12]
samp.13 <- sample(1:length(temp.13.2), size=samps)
temp.13.3 <- temp.13.2[samp.13]
samp.14 <- sample(1:length(temp.14.2), size=samps)
temp.14.3 <- temp.14.2[samp.14]
samp.15 <- sample(1:length(temp.15.2), size=samps)
temp.15.3 <- temp.15.2[samp.15]
fertoutput <- littersize[c(temp.12.3, temp.13.3, temp.14.3, temp.15.3)]
hist(fertoutput, breaks=20)
quantile(fertoutput, probs=c(0.025, 0.975))

plot(density(littersize[temp.12.3]), xlim=c(0,1))
lines(density(littersize[temp.13.3]), col="red")
lines(density(littersize[temp.14.3]), col="gray")
lines(density(littersize[temp.15.3]), col="blue")

#Get the MCMC run from which each fecundity value was taken - to account for linkage b/w fertility and survival
#Get the columns
temp.cols <- col(fracnativo[,,3])
cols.12.2 <- temp.cols[which(temp.12==TRUE)]
cols.12.3 <- cols.12.2[samp.12]
cols.13.2 <- temp.cols[which(temp.13==TRUE)]
cols.13.3 <- cols.13.2[samp.13]
cols.14.2 <- temp.cols[which(temp.14==TRUE)]
cols.14.3 <- cols.14.2[samp.14]
cols.15.2 <- temp.cols[which(temp.15==TRUE)]
cols.15.3 <- cols.15.2[samp.15]

MCMCfert <- NULL
for (i in 1:samps){
	x.12 <- MCMCsample[temp.12.3[i], cols.12.3[i]]
	x.13 <- MCMCsample[temp.13.3[i], cols.13.3[i]]
	x.14 <- MCMCsample[temp.14.3[i], cols.14.3[i]]
	x.15 <- MCMCsample[temp.15.3[i], cols.15.3[i]]
	MCMCfert <- c(MCMCfert, x.12, x.13, x.14, x.15)
}


###Calculate lambda using popbio package
library(popbio)
nsims <- length(MCMCfert)
lambda <- 0
for (i in 1:nsims){
	randnum <- MCMCfert[i]
	p <- sample(c(1:periods), 1)
	m = matrix(c(0, fertoutput[i], juvsurvival(randnum, p), adultsurvival(randnum, p)), nrow=2, ncol=2, byrow=TRUE)
	lambda[i] <- lambda(m)
}
hist(lambda)
quantile(lambda, probs=c(0.025, 0.975))

setwd(BUGSDir)
save.image("SevModel3_fert_2016-07-16_1mnth_newprecip.RData")

#Caveat: this is kind of a "prospective" lambda in that survival and fertility are not linked,
#So these are lambdas that we might expect to see in the future under similar conditions

#######
#Looking at effect of precip on fecundity
setwd(BUGSDir)
#load("SevModel3_fert_2016-05-11.RData")
#Compare each year to the Nov-April precip starting in the previous year
pN11.A12 <- sum(envCond_df[envCond_df$Date > "2011-11-01" & envCond_df$Date < "2012-05-01",]$ppt)
pN12.A13 <- sum(envCond_df[envCond_df$Date > "2012-11-01" & envCond_df$Date < "2013-05-01",]$ppt)
pN13.A14 <- sum(envCond_df[envCond_df$Date > "2013-11-01" & envCond_df$Date < "2014-05-01",]$ppt)
pN14.A15 <- sum(envCond_df[envCond_df$Date > "2014-11-01" & envCond_df$Date < "2015-05-01",]$ppt)


