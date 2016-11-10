#################### 
#  R script for running pdog survival analysis for large plots at Sevilleta NWR
#


#TODO:   

# plot potential weight loss/gain over time for individuals recaptured multiple times. Correlate with environmental conditions?

# NOTE: trapping days 4 and 5 in 2015 (plots b and f) were experimental, and should probably not be used

rm(list=ls())

ANA = F
KEVIN = F
KEVINOFFICE = F
ELIZABETHOFFICE = T

###################################
##########  SET WORKING DIRECTORY

if(KEVIN) rootDir <- "C:\\Users\\Kevin\\Dropbox\\Sev Pdog Study\\BIG plot study\\"
if(KEVINOFFICE) rootDir <- "E:\\Dropbox\\Sev Pdog Study\\BIG plot study\\"
if(ANA) rootDir <- ""
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

library(R2WinBUGS)
library(lubridate)
library(Hmisc)

##################################
######### SET MCMC Params

ni = 5000
nc <- 1
nb <- 1000
nt <- 1

nsims <- (ni-nb)/nt


###################################
##########   READ IN DATA

############
#  "pdog" : Prairie dog capture data. 

setwd(DataDir)
pdog <- read.csv("pdog_captures_16ha_9-12-15 AD_KTS_2016.csv", header=TRUE)  ## KTS: updated with 2015 data and fixed a couple errors in the data
## EAH: updated with 2016 data through June

 #  head(pdog)
 #  nrow(pdog)   # 3423 observations

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

pdog$TRAPDAY

###############################
####### DEFINE "PERIOD" FOR THIS STUDY

# Note: some plots were trapped at different times within the same season- so trap day 1 isn't necessarily always the same date...

setwd(DataDir)
filename <- "PeriodDefinition.csv"
sink(filename)
cat("Period,Start,End,StartTrap,EndTrap,TrapDays,StartRelease,EndRelease,Median
SUMMER2010,6/1/2010,5/28/2011,NA,NA,0,6/18/2010,8/19/2010,7/1/2010
SUMMER2011,6/1/2011,2/28/2012,6/13/2011,6/15/2011,3,6/28/2011,9/22/2011,6/14/2011
SPRING2012,3/1/2012,5/28/2012,4/2/2012,4/2/2012,1,4/9/2012,4/9/2012,4/2/2012
SUMMER2012,6/1/2012,2/28/2013,6/13/2012,6/15/2012,3,6/18/2012,9/14/2012,6/14/2012
SPRING2013,3/1/2013,5/28/2013,3/25/2013,4/5/2013,3,4/9/2013,4/9/2013,4/1/2013
SUMMER2013,6/1/2013,2/28/2014,6/10/2013,6/22/2013,3,NA,NA,6/15/2013
SPRING2014,3/1/2014,5/28/2014,3/17/2014,3/27/2014,3,NA,NA,3/20/2014
SUMMER2014,6/1/2014,5/28/2015,6/9/2014,6/18/2014,3,NA,NA,6/13/2014
SPRING2015,6/1/2015,6/30/2015,6/1/2015,6/5/2015,3,6/23/2015,6/25/2015,6/3/2015
SUMMER2015,7/1/2015,2/28/2016,7/21/2015,7/23/2015,3,8/13/2015,8/13/2015,7/22/2015
SPRING2016,3/1/2016,4/28/2016,4/11/2016,4/13/2016,3,NA,NA,4/12/2016
SUMMER2016,5/1/2016,6/7/2016,6/5/2016,6/7/2016,3,NA,NA,6/6/2016")
sink()

setwd(DataDir)
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


#######PRECIPITATION INTERVALS
#For each period, define the precip interval that would affect that period
#For summer periods: May-Oct
#For spring periods: Nov-Apr
setwd(DataDir)
filename <- "PeriodDefinition_Precip.csv"
sink(filename)
cat("Period,Start,End
SUMMER2010,5/1/2010,4/30/2011
SUMMER2011,5/1/2011,10/31/2011
SPRING2012,11/1/2011,4/30/2012
SUMMER2012,5/1/2012,10/31/2012
SPRING2013,11/1/2012,4/30/2013
SUMMER2013,5/1/2013,10/31/2013
SPRING2014,11/1/2013,4/30/2014
SUMMER2014,5/1/2014,10/31/2014
SPRING2015,11/1/2014,4/30/2015
SUMMER2015,5/1/2015,10/31/2015
SPRING2016,11/1/2015,4/30/2016
SUMMER2016,5/1/2016,6/7/2016")
sink()

setwd(DataDir)
PrecipPeriodDef <- read.csv(filename)

PrecipPeriodDef$StartDate <- mdy(PrecipPeriodDef$Start)
PrecipPeriodDef$EndDate <- mdy(PrecipPeriodDef$End)

###############################
#######  READ IN ENV CONDITIONS- PRECIP AND SOIL MOISTURE

setwd(DataDir)
      # KTS: now we read in monthly data and summarize in R
 #precip_df <- read.csv("Precip data 2010-2014 - new - summarized for actual periods.csv",header=T)   # actually includes 2015  
 #precip_df$Period <- as.character(precip_df$Period)
 #names(precip_df)


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


setwd(DataDir)
save(nobsG,nindG,realperiodyearsG,realdatesG,nreleasesG,realreleasedatesG,nsurveysG,realsurveydatesG,ndatesG,realplotsG,nplotsG,
     realperiodsG,nperiodsG,realyearsG,nyearsG,realindG,sumReleases,precip,soilH2O,
     nsims,file="GlobalParams.RData")

#####################################
##############  Verify survey effort across sites

for(t in 1:nyearsG){
  ndx <- which(pdog$YEAR==realyearsG[t])
  print(table(pdog$PLOT[ndx]))
}

    ### there was survey effort conducted at plots B and C every year. Then at F as soon as it came online


###############################
############   SUMMARIZE RELEASES  

 ## NOTE: TRAPDAY signifies real trapping occasions/suboccasions. When TRAPDAY is NA, is was a release occasion


names(pdog)
table(pdog$PTP,pdog$TRAPDAY)   # no trapping in 2010- just releases


table(pdog$RELEASE)    # 2315 pdogs released. 957 pdogs captured

        ## KTS: note: do we have recruitment information in this data set? That is, individuals recruited into
               # the population are those that were not released! 


################################
############## SUMMARIZE DATA FOR WINBUGS

names(pdog)
head(pdog)

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
#######       MAKE CAPTURE HISTORY DATA FOR WINBUGS

      # note: if possible, capture histories should include the releases/initial captures... #KTS: Why didn't they before????


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


###############################
############ MAKE SITE VARIABLE  (history of plots each individual was captured on) (check: did any move from plot to plot?)


plothist <- array(NA,dim=c(maxnan,nperiodsG)) 
i=1

factorplot <- as.factor(pdog$PLOT) 
for(i in 1:nan){
  p=1
  for(p in 1:nperiodsG){
      ndx <- which((pdog$PTP==realperiodsG[p])&
                   (pdog$indiv==realindG[i]))
      if(length(ndx)>0) plothist[i,p] <- as.numeric(factorplot[ndx[1]]) 
  }
}

head(plothist,500)
tail(plothist,500)

   ## interesting- several pdogs actually do switch plots...
   ## Actually- not true- several ear tags are doubles! notably, 321-324 were re-used in 2015

plothist2 <- apply(plothist,c(1,2),function(t) ifelse(is.na(t),0,t))
temp <- apply(plothist2,1,function(t) length(unique(t)))

length(which(temp==3))  # total of 5 pdogs moved from one plot to another.
realindG[which(temp==3)]   # the following individuals "moved": 563  646  730 2109 2214

CHList$indiv563
CHList$indiv646
CHList$indiv1507    # one individual moved twice!!

# NOTE: KTS: fixed the obvious errors here with duplicate tags... 
# NOTE: two individuals released in Plot "f" moved to different plots, where they were subsequently recaptured. There is little evidence of any other such movements


############### MAKE DUMMY MATRIX FOR EACH PLOT

plothist3 <- plothist2*(1-Released)     # only those individuals that were captured

isPlotB <- apply(plothist3,c(1,2),function(t) ifelse(t==1,1,0) )
isPlotD <- apply(plothist3,c(1,2),function(t) ifelse(t==2,1,0) )
isPlotF <- apply(plothist3,c(1,2),function(t) ifelse(t==3,1,0) )

head(isPlotB,100)

############################
##    MAKE COHORT VARIABLE

       ### loop through years. all individuals released on a given year get assigned a unique cohort... 
       ### need an "is.rel" variable also. Should we model a release effect?

cohort <- numeric(nan)
for(i in 1:nan){
  for(p in 1:length(realperiodsG)){
    ndx <- which((pdog$PTP==realperiodsG[p])&
                 (pdog$RELEASE=="Y")&
                 (pdog$indiv==realindG[i]))
    if(length(ndx)>0) cohort[i] <- p  
  }
}    
ndx <- which(cohort==0)           # NOTE: The last cohort represents native born individuals
cohort[ndx] <- max(cohort)+1
ncohorts <- max(cohort)

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


##################
# determine if the pdog was from taos (hoarded individuals)

hoarded <- numeric(nan)

i=1
for(i in 1:nan){
  thisch <- CHList[[i]]
  thisrel <- which(thisch$RELEASE=="Y")
  if(length(thisrel)>0){
    temp <- grep("aos",thisch$SITE)
    if(length(temp)>0) hoarded[i] <- 1
  }  
}

hoardndx <- which(hoarded==1)

CHList[[hoardndx[1]]]

##############################
#######   WINBUGS CODE


######################
#     SET BUGS OR JAGS DIRECTORY
#####################

if(KEVIN){
  bugs.directory="C:\\Users\\Kevin\\Documents\\Employment\\ESF\\Bog Turtle\\DATA\\software\\BUGS\\WinBUGS14"
    #if(HRAPC) bugs.directory= "C:\\Users\\Kevin\\Desktop\\Kevin\\WinBUGS14"
}

if(KEVINOFFICE){
  bugs.directory="C:\\WinBUGS\\winbugs14\\WinBUGS14"
  #if(HRAPC) bugs.directory= "C:\\Users\\Kevin\\Desktop\\Kevin\\WinBUGS14"
}

if(ELIZABETHOFFICE){
  bugs.directory= "C:\\Users\\eahunter\\Documents\\WinBUGS14\\"
}

working.directory= BUGSDir   
setwd(working.directory)

################################
#      # WRITE BUGS OR JAGS MODEL TO FILE
################################

filename <- "SevPDModel_20160617_oneMonthRelEff"  

model.file <- paste(filename,".bug",sep="")


sink(model.file)

cat("
    
# Cormack Jolly Seber model, robust design 
#  K Shoemaker Nov 2015  

    ##########################################
    #  BEGIN MODEL
    model {
    
    ###############################################
    #  SET PRIORS
    
    #####################
    # survival terms 
    
    phi0 ~ dunif(0.01,0.99)           # mean survival rate
    logit.phi <- log(phi0/(1-phi0))   # transform to logit scale... 
    precipEff ~ dunif(-2,2)           # beta term for effect of (normalized) precip on survival  
    soilEff ~ dunif(-2,2)             # beta term for effect of (normalized) soil moisture on survival
    maleEff ~ dunif(-2,2)             # beta term for sex effect on survival   
    juvEff ~ dunif(-2,2)              # beta term for juvenile/age effect on survival
    relEff ~ dunif(-28,0.5) #-4            # beta term for release effect on survival# NOTE: effect lasts for only one month!
    #juvRelEff ~ dunif(-4,4)           # beta term for juvenile/release interaction 
    springRelEff ~ dunif(-8,8)        # beta term for spring release effect        
    hoardEff ~ dunif(-2,2)            # beta term for effect of hoarded individuals (taos 2012) on survival
    
 
    ##########################
    # capture probability terms    

    p0 ~ dunif(0.01,0.75)           # mean probability of capture per secondary trapping occasion    
    logit.p <- log(p0/(1-p0))       # transform to logit scale

     sd.ind   ~ dunif(0.1,1)         #  random effect of individual on capture probability (due to overlap of home range)  # KTS: added this back
     prec.ind <- 1/pow(sd.ind,2) 

    for(i in 1:nan){
      indeff[i] ~ dnorm(0,prec.ind)   # logit-normal random intercept (KTS: added this back in when I removed the temp migration term)
    }

    ##################
    # tag loss terms (KTS: added 11/16/15) 
             # note: we assume that PIT tags are never lost
             #  todo: account for loss of both tags in the same time interval??????

    tagRetentionRate ~ dunif(0.4,0.99)         # annual rate of ear tag retention (per tag) 

    for(i in 1:nan){
      tagged[i,first[i]] <- 1                        # no chance of losing any tags at release 
      for(t in (first[i]+1):ns){
        pHasLeftEarTag[i,t] <- (1-leftTagAdded[i,t])*(hasLeftEarTag[i,t-1] * pow(tagRetentionRate,interval[t-1])) + leftTagAdded[i,t]     # probability of retaining a tag in the interval since last survey (to lose a tag, it had to have at least one!).  NOTE: assumes only one tag can be lost each period
        hasLeftEarTag[i,t] ~ dbern(pHasLeftEarTag[i,t])              # latent state/DATA NODE- Does it have a left ear tag? This is sometimes observed.     
        pHasRightEarTag[i,t] <- (1-rightTagAdded[i,t])*(hasRightEarTag[i,t-1] * pow(tagRetentionRate,interval[t-1])) + rightTagAdded[i,t]     # probability of retaining a tag in the interval since last survey (to lose a tag, it had to have at least one!).  NOTE: assumes only one tag can be lost each period
        hasRightEarTag[i,t] ~ dbern(pHasRightEarTag[i,t])              # latent state/DATA NODE- does it have a right ear tag? This is sometimes observed.
        nTags[i,t] <- hasLeftEarTag[i,t] + hasRightEarTag[i,t] + PITTag[i,t]            
        tagged[i,t] <- step(nTags[i,t]-1)            # latent state: is the individual lost from the population??
      }
    }   
   
    #####################################################
    #  SET OVERALL PROCESS MODEL
    
    
    #  SURVIVAL MODEL (no release effect)
    
    for(i in 1:nan){
      for(t in 1:(ns-1)){
        mu.phi[i,t]   <- logit.phi + precipEff*precip[t] + soilEff*soil[t] +        
                               juvEff*isJuv[i,t] + maleEff*isMale[i] +    ### + relEff*isRel[i,t] +    ## KTS: removed release effects here
                               hoardEff*isHoard[i] #### + juvRelEff*isJuv[i,t]*isRel[i,t] + springRelEff*springRel[i,t]  # function for survival rate      
        phi[i,t]      <- 1/(1+exp(-1*mu.phi[i,t]))                     
      }
    }
    

    #  SURVIVAL MODEL (with release effect)  ## NOTE: this effect is annualized, to make sure everything is on a common currency.
                                                       ## HOWEVER, the effect only ever lasts for one month. 
    
    for(i in 1:nan){
      for(t in 1:(ns-1)){
        mu.phiRel[i,t]   <- logit.phi +		## precipEff*precip[t] + soilEff*soil[t] +   #No precip effect for release     
        juvEff*isJuv[i,t] + maleEff*isMale[i] + relEff*isRel[i,t] +    # with release-year effects 
        hoardEff*isHoard[i]  + springRelEff*springRel[i,t]       #Removed + juvRelEff*isJuv[i,t]*isRel[i,t]  - confounded w/ spring
        phiRel[i,t]      <- 1/(1+exp(-1*mu.phiRel[i,t]))                     
      }
    }
    
    ########## Latent variable: living or dead 
    
    for(i in 1:nan){
      alive[i,first[i]] ~ dbern(1)   # at first capture, the animal is alive
    
      for(t in (first[i]+1):ns){
        transition_norel[i,t] <- pow(phi[i,t-1],interval[t-1])   
        transition_rel[i,t] <- pow(phiRel[i,t-1],(1/12)) * pow(phi[i,t-1],(interval[t-1]-(1/12)))  # release effect lasts for only one month...
        transition[i,t] <- isRel[i,t-1]*transition_rel[i,t] + (1-isRel[i,t-1])*transition_norel[i,t]
        mualive[i,t] <- alive[i,(t-1)] * transition[i,t]  # probability of each ind being alive, all subsequent periods
        alive[i,t] ~ dbern(mualive[i,t])        # latent variable- is it still alive?
      }
    }
    
    ####################################
    #  OBSERVATION MODEL   (actual data likelihood)

                        # use the year/season of first capture to help refine p estimates. Assume no tags fall off within a single 3-day trapping session
    for(i in 1:nan){
      for(j in (first2[i]+1):nss[first[i]]){
        mu.p[i,first[i],j] <- logit.p + indeff[i]   # KTS: added individual effect back
        p[i,first[i],j] <- (1/(1+exp(-1*mu.p[i,first[i],j])))   #* surveyed[first[i],j]
        muy[i,first[i],j] <- alive[i,first[i]]*p[i,first[i],j] 
        y[i,first[i],j] ~ dbern(muy[i,first[i],j])    #  DATA NODE   likelihood of observed data... 
      }
    }
    
    for(i in 1:nan) { 
      for(t in (first[i]+1):ns) { 
        for(j in 1:nss[t]){                
          mu.p[i,t,j]  <- logit.p  + indeff[i]            #  logit probability of capture for this individual   [[add term for number of grids deployed]]  
          p[i,t,j]     <- (1/(1+exp(-1*mu.p[i,t,j]))) #* surveyed[t,j]                                   #  convert back to probability scale
          muy[i,t,j]   <- alive[i,t]*tagged[i,t]*p[i,t,j] # *onsite[i,t]                    # if it's alive and tagged, then it's seen with prob. p.  
          y[i,t,j] ~ dbern(muy[i,t,j])                 #  DATA NODE   likelihood of observed data...  
        }
      }       
    }		
    
    #################################################
    #  ADDITIONAL CALCULATIONS
    
    #################
    ####  CALCULATE ABUNDANCE USING HORVITZ-THOMPSON ESTIMATOR....
    # note: only applies to sampled areas (grids)
    # account for # grids deployed
    # note: should be per hectare.
     
    for(i in 1:nan){
      for(t in 2:ns){
        mu.p2[i,t,1] <- logit.p + indeff[i]
        p2[i,t,1] <- 1/(1+exp(-1*mu.p2[i,t,1]))   # back to prob. scale
        pncap[i,t,1] <- 1-p2[i,t,1]        # pncap refers to the probability of not capturing for entire 3-day interval
        for(j in 2:nss[t]){
          mu.p2[i,t,j] <- logit.p + indeff[i]
          p2[i,t,j] <- 1/(1+exp(-1*mu.p2[i,t,j]))   # back to prob. scale
          pncap[i,t,j] <- pncap[i,t,j-1]*(1-p2[i,t,j])
        }
        pcap[i,t] <- 1-pncap[i,t,nss[t]]     # pcap refers to the prob of being captured at least once for entire 3-day period
        invpcap[i,t] <- 1/pcap[i,t]    # for H-T estimator
      } 
    }


    for(t in 2:ns){
      Ntot[t] <- inprod(invpcap[1:nan,t],ch2D[1:nan,t])     #nan2[t] / pcap[t]   # estimate of total abundance within sampled region... 
      #N[t] <- N[t]/nGrids[t]    # convert to abundance per ha (assuming 4 grids were always deployed...)
      N_plotB[t] <- inprod(invpcap[1:nan,t],isB[1:nan,t])
      N_plotD[t] <- inprod(invpcap[1:nan,t],isD[1:nan,t])
      N_plotF[t] <- inprod(invpcap[1:nan,t],isF[1:nan,t])
      
    }
    
    #################
    ##### CALCULATE SURVIVAL RATE OVER TIME...
    
    # store the survival rate for each year
    for(t in 1:(ns-1)){
      femphi2[t] <- logit.phi + precipEff*precip[t] + soilEff*soil[t]   # + predEff*badyr[t]
      femphi[t] <- 1/(1+exp(-1*femphi2[t]))
      malephi2[t] <- logit.phi + precipEff*precip[t] + soilEff*soil[t] +  maleEff 
      malephi[t] <- 1/(1+exp(-1*malephi2[t]))
      juvphi2[t] <- logit.phi + precipEff*precip[t] + soilEff*soil[t] + juvEff
      juvphi[t] <- 1/(1+exp(-1*juvphi2[t]))
      #reladultphi2[t] <- logit.phi + precipEff*precip[t] + soilEff*soil[t] + relEff
      #reladultphi[t] <- 1/(1+exp(-1*reladultphi2[t]))
    } 

    ### determine average 1-month survival of released individuals...

    falljuvrelsurv_mo <-  pow(1/(1+exp(-1*(logit.phi + relEff + juvEff))),(1/12))  #Removed + juvRelEff 
    falladrelsurv_mo <-  pow(1/(1+exp(-1*(logit.phi + relEff))),(1/12))
    springjuvrelsurv_mo <-  pow(1/(1+exp(-1*(logit.phi + relEff + juvEff + springRelEff))),(1/12))  #Removed + juvRelEff 
    springadrelsurv_mo <-  pow(1/(1+exp(-1*(logit.phi + relEff + springRelEff))),(1/12))        
    
    ##########################################
    }  # END MODEL
    ##########################################


    


    ",fill = TRUE)

sink()                   ## Writes BUGS code to the appropriate file



######################################
########  summarize total individuals captured in trap surveys (per grid...)

ndx <- which((pdog$RELEASE=="Y")) 
pdogTRAPPED <- pdog[-ndx,] 

     # summarize the number captured in each year

temp <- table(pdogTRAPPED$PTP)
namesP <- names(temp)

ndx<-match(realperiodsG,namesP)

nanTRAPPED <- as.numeric(temp)[ndx]
nanTRAPPED[which(is.na(nanTRAPPED))] <- 0

table(pdogTRAPPED$PLOT,pdogTRAPPED$PTP)[,ndx]

nGrids <- c(0,2,2,2,3,3,3,3,3,3)

nanTRAPPED_perGrid <- nanTRAPPED/nGrids

nanTRAPPED_perGrid[1] <- 0
 
########################################
########  START BUILDING DATA FOR WINBUGS


     ########################################
     #########    standardize covariates

meanlogprecip <- mean(log(precip))             # standardize precip variable
sdlogprecip <- sd(log(precip))
precip_std <- (log(precip)-meanlogprecip)/sdlogprecip

meanlogsoil <- mean(log(soilH2O))             # standardize soil variable
sdlogsoil <- sd(log(soilH2O))
soil_std <- (log(soilH2O)-meanlogsoil)/sdlogsoil

Data <- list(nan = nan, 
             #nan2 = nanTRAPPED_perGrid,
             ns = nperiodsG,
             #nGrids = nGrids,
             #ncohorts = ncohorts,
             #cohort = cohort,
             soil = as.vector(soil_std),
             nss = as.vector(nss),
             first = as.vector(firsts),
             first2 = as.vector(first2),
             interval = as.vector(intervals),  
             precip=as.vector(precip_std),
             y = caphist,
             ch2D = caphist_2D,
             isB = isPlotB,
             isD = isPlotD,
             isF = isPlotF,
             isMale = isMale,
             isJuv = isJuv,
             isRel = Released,
             springRel = SpringRelease,
             PITTag = PITTag,
             leftTagAdded = AddLeftTag,
             rightTagAdded = AddRightTag,
             hasLeftEarTag = LeftEarTag,
             hasRightEarTag = RightEarTag,
             isHoard = hoarded      
)

Z <- array(1,dim=c(nan,nperiodsG))

Inits <- function() list(       #  initial values for all parameters for BUGS MCMC routine
                phi0 = runif(1,.4,.6),
                precipEff = runif(1,0.5,1),
                #predEff = runif(1,-2,-1),
                soilEff = runif(1,-2,0),
                maleEff = runif(1,-1,0),
                juvEff = runif(1,-0.5,0.5),
                relEff = runif(1,-3.5,-2.5),
                #juvRelEff = runif(1,-0.5,0.5),
                springRelEff = runif(1,-0.5,0.5),
                #sd.ind = runif(1,0.2,0.3),
		            p0 = runif(1,.2,.3),
                tagRetentionRate = runif(1,0.6,0.7),
                #pbadyr = runif(1,0.4,0.6),
                #gamma.prime = runif(1,0.4,0.6),
                #gamma.dprime = runif(1,0.4,0.6),
                #indeff = runif(nan,-0.1,0.1),
                #badyr3 = rep(1,times=nperiodsG),
                #onsite = Z,
                alive = Z               
		)  

Par <- c(          #   parameters to save the posterior samples for
         #"N",     # H-T estimate of abundance        
         "phi0",      # mean probability of survival
         #"phi",      # actual estimate of survival for each year
	   #"badyr",
         #"pbadyr",  # probability that a sampling period has low survival (due to predation?)
         "precipEff",   # effect of precip on survival
         "soilEff",
         "hoardEff",
         "maleEff",
         #"predEff",     # reduction in survival rate in a bad year
         #"cohortEff",
         "juvEff",
         "relEff",
         #"juvRelEff",
         "springRelEff",
         "tagRetentionRate",
         #"sd.ind",
         "p0",            # mean nightly prob. of capture
         #"gamma.prime",    # probability of staying off site
         #"gamma.dprime",     # probability of leaving??
         "femphi",
         "malephi",
         "juvphi",
         #"reladultphi",
		     "Ntot", 
         "N_plotB",
         "N_plotD",
         "N_plotF",
	   "falljuvrelsurv_mo", 
    		"falladrelsurv_mo",
    		"springjuvrelsurv_mo", 
   		"springadrelsurv_mo"
)   


setwd(BUGSDir)

Mod <- bugs(data=Data, inits=Inits, parameters.to.save=Par,  #  run BUGS model
                   model.file=model.file, n.chains=2, n.iter=100000,                            # nc / ni / nb / nt
                   n.burnin=50000, n.thin=5,over.relax = TRUE,bugs.directory=bugs.directory,debug=FALSE)
                                         

setwd(BUGSDir)
save(Mod,file=sprintf("SevModel10_1mnth_newprecip_%s.RData",Sys.Date()))

nsims = 2500
									   
cat(paste("SCRIPT FINISHED RUNNING SUCCESSFULLY, ",Sys.time(),", on ",Sys.Date(),"\n",sep=""))



######################
#Summary stats for table
setwd(BUGSDir)
load("SevModel10_1mnth_newprecip_2016-07-16.RData")

#Phi
mean(Mod$sims.list$phi)
quantile(Mod$sims.list$phi ,c(0.025,0.975))

#Release effect
mean(Mod$sims.list$relEff)
quantile(Mod$sims.list$relEff ,c(0.025,0.975))

#Spring Release effect
mean(Mod$sims.list$springRelEff)
quantile(Mod$sims.list$springRelEff ,c(0.025,0.975))

#Juvenile effect
mean(Mod$sims.list$juvEff)
quantile(Mod$sims.list$juvEff ,c(0.025,0.975))

#Male effect
mean(Mod$sims.list$maleEff)
quantile(Mod$sims.list$maleEff ,c(0.025,0.975))

#Precipitation effect
mean(Mod$sims.list$precipEff)
quantile(Mod$sims.list$precipEff ,c(0.025,0.975))

#Soil moisture effect
mean(Mod$sims.list$soilEff)
quantile(Mod$sims.list$soilEff ,c(0.025,0.975))

#p
mean(Mod$sims.list$p0)
quantile(Mod$sims.list$p0 ,c(0.025,0.975))

#Tag retention rate
mean(Mod$sims.list$tagRetentionRate)
quantile(Mod$sims.list$tagRetentionRate ,c(0.025,0.975))

#######Derived
#Residents
#Adult female:
mean(Mod$sims.list$femphi)
quantile(Mod$sims.list$femphi ,c(0.025,0.975))

#Adult male:
mean(Mod$sims.list$malephi)
quantile(Mod$sims.list$malephi ,c(0.025,0.975))

#Juvenile:
mean(Mod$sims.list$juvphi)
quantile(Mod$sims.list$juvphi ,c(0.025,0.975))

#Releases:
#Spring:
mean(Mod$sims.list$springadrelsurv_mo)
quantile(Mod$sims.list$springadrelsurv_mo, c(0.025,0.975))

#Summer:
mean(Mod$sims.list$falladrelsurv_mo)
quantile(Mod$sims.list$falladrelsurv_mo, c(0.025,0.975))


##############################################
############### END OF SCRIPT ################
