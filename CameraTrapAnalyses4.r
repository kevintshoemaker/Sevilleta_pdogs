#####################
### Process data for Sevilleta camera trap analysis
#####################

rm(list=ls())
##################

#####################
# LOAD PACKAGES
#####################

library(unmarked)
library(lubridate)
library(Hmisc)
library(ggplot2)

#####################
# PREPARE WORKSPACE
#####################

KEVIN_OFFICE = TRUE # FALSE # FALSE
KEVIN_LAPTOP = FALSE # TRUE # TRUE

if(KEVIN_OFFICE) BaseDir <- "E:\\Dropbox\\Camera trap study\\" 
if(KEVIN_LAPTOP) BaseDir <- "C:\\Users\\Kevin\\Dropbox\\Camera trap study\\"
ExamplesDir <- sprintf("%s%s",BaseDir,"Examples")
DataDir <- sprintf("%s%s",BaseDir,"Data")
FiguresDir <- sprintf("%s%s",BaseDir,"figures")

setwd(DataDir)
# new data!

seasons <- c("spring","summer","fall")

#####################
# READ CAMERA DEPLOYMENT DATA
#####################

setwd(DataDir)
tempdf <- read.csv("CameraDeploymentDates.csv",header = T,na.strings=c("missing","NA"),stringsAsFactors = F)
head(tempdf)

##################
### PROCESS THE CAMERA DEPLOYMENT DATA
##################

tempdf$fall <- gsub("Sept","Sep",tempdf$fall)

nsesaons <- length(seasons)

CameraDates <- list()

CameraDates[["CameraNames"]] <- tempdf$camera

sitenames <- CameraDates[["CameraNames"]]
nsites <- length(sitenames)

s=3
for(s in 1:length(seasons)){
  CameraDates[[seasons[s]]] <- list()
  CameraDates[[seasons[s]]][["START"]] <- dmy(unlist(lapply(strsplit(tempdf[,(s+1)],split="[.]"),function(t) t[1] )) ) #gsub("2014","2015",tempdf$First.date)       # ANA- is this okay?
  CameraDates[[seasons[s]]][["END"]] <- dmy(unlist(lapply(strsplit(tempdf[,(s+1)],split="[.]"),function(t) t[2] )) )
  CameraDates[[seasons[s]]][["TOTAL_DAYS"]] <- as.numeric(difftime(CameraDates[[seasons[s]]][["END"]],CameraDates[[seasons[s]]][["START"]],units="days"))
  names(CameraDates[[seasons[s]]][["START"]]) <- CameraDates[["CameraNames"]]
  names(CameraDates[[seasons[s]]][["END"]]) <- CameraDates[["CameraNames"]]
  names(CameraDates[[seasons[s]]][["TOTAL_DAYS"]]) <- CameraDates[["CameraNames"]]
}


###################
# SPECIES LIST
###################

setwd(DataDir)
speciesdf <- read.csv("SpeciesListNEW.csv",header=T,stringsAsFactors = F)    # "SpeciesListNEW-altList.csv"
head(speciesdf)

speciesnames <- as.character(speciesdf$Species)

speciescodes <- c("CYGU",as.character(speciesdf$CODE))
speciesnames <- c("Gunnison's Prairie Dog",speciesnames)

#speciesnames[which(speciesnames=="Unknown raptor")=="Raptor spp."]

cbind(speciescodes,speciesnames)


#####################
# READ CAMERA TRAP DATA
#####################

tempdf <- read.csv("Camera data for occupancy AD 7-26-16.csv",             
                         header=T,na.strings=c(".","NA"),stringsAsFactors = F)
head(tempdf)


####################
# PROCESS DATES
####################

tempdf$Sample.period <- factor(tempdf$Sample.period,levels=c("Spring","Summer","Fall"))
levels(tempdf$Sample.period) <- seasons    # change to match the deployment data
tempdf$SEASON <- tempdf$Sample.period   # season

tempdf <- tempdf[order(tempdf$Sample.period),]    # order by season...
head(tempdf)

tempdf$DATE_TAKEN <- gsub("2014","2015",tempdf$Date.taken)
tempdf$DATE_TAKEN <- mdy(tempdf$DATE_TAKEN)
   

#################
# REMOVE DUPLICATES (not sure how those got in there!!)
#################

newdf <- tempdf[-(1:nrow(tempdf)),]  
species <- "INSP"
for(species in speciescodes){
  this <- subset(tempdf,CODE==species)
  unique.pics <- unique(this$Pic.No) 
  ndx <- match(unique.pics,this$Pic.No)
  newdf <- rbind(newdf,this[ndx,])
}
#newdf
tempdf <- newdf

#################
# CORRECT THE CRAZY ISSUE WITH PHOTO DATES (as good as possible...)
#################


tempdf$DATE_TAKEN_CORRECTED <- tempdf$DATE_TAKEN
errors <- NULL
camera <- "A1"
season <- "summer"
for(camera in CameraDates$CameraNames){
  for(season in seasons){
    this <- subset(tempdf,(SEASON==season)&(Camera==camera))
    real.first.date <- CameraDates[[season]]$START[camera]
    real.end.date <- CameraDates[[season]]$END[camera]
    too.early <- which(this$DATE_TAKEN<real.first.date)
    too.late <- which(this$DATE_TAKEN>real.end.date)
    if( all(c(length(too.early)>0,length(too.late)>0)) ){
      errors <- c(errors,paste("ERROR",camera,season))
      #break
    }
    if( length(too.early)>0 ){
      difs <- as.numeric(difftime(this$DATE_TAKEN[too.early],real.first.date,units="days"))
      maxdif <- difs[which.max(abs(difs))]
      correction <- -1*maxdif
      ndx <- which(tempdf$Pic.No%in%this$Pic.No)#match(this$Pic.No,tempdf$Pic.No)
      tempdf$DATE_TAKEN_CORRECTED[ndx]<-tempdf$DATE_TAKEN[ndx]+days(correction)
    }
    if( length(too.late)>0 ){
      difs <- as.numeric(difftime(this$DATE_TAKEN[too.late],real.end.date,units="days"))
      maxdif <- difs[which.max(abs(difs))]
      correction <- -1*maxdif
      ndx <- which(tempdf$Pic.No%in%this$Pic.No)
      tempdf$DATE_TAKEN_CORRECTED[ndx]<-tempdf$DATE_TAKEN[ndx]+days(correction)
    }
  }
}
errors

#tempdf$DATE_TAKEN <- tempdf$DATE_TAKEN_CORRECTED 

tempdf$TIME_TAKEN <- hms(as.character(tempdf$Time.Taken))
tempdf$TIME_TAKEN[which(is.na(tempdf$TIME_TAKEN))] <- hms("9:00:00")    # ANA- is this okay?

tempdf$DATETIME_TAKEN <- tempdf$DATE_TAKEN_CORRECTED+tempdf$TIME_TAKEN

tempdf$MONTH <- month(tempdf$DATETIME_TAKEN)
tempdf$YEAR <- year(tempdf$DATETIME_TAKEN)
tempdf$HOUR <- hour(tempdf$DATETIME_TAKEN)
tempdf$WEEK <- week(tempdf$DATETIME_TAKEN)

table(tempdf$MONTH)
table(tempdf$WEEK)
table(tempdf$YEAR)

tempdf$SITE <- tempdf$Camera
tempdf$TREATMENT <- tempdf$Trt
tempdf$GRID <- tempdf$Plot


########
## DEBUG

tempdf$LASTWEEK <- NA
tempdf$FIRSTWEEK <- NA
tempdf$LASTDATE <- tempdf$DATE_TAKEN
tempdf$FIRSTDATE <- tempdf$DATE_TAKEN
i=1
for(i in 1:nrow(tempdf)){
  tempdf$LASTWEEK[i] <- week(CameraDates[[as.character(tempdf$SEASON[i])]]$END[tempdf$Camera[i]]) 
  tempdf$FIRSTWEEK[i] <- week(CameraDates[[as.character(tempdf$SEASON[i])]]$START[tempdf$Camera[i]])
  tempdf$LASTDATE[i] <- CameraDates[[as.character(tempdf$SEASON[i])]]$END[tempdf$Camera[i]] 
  tempdf$FIRSTDATE[i] <- CameraDates[[as.character(tempdf$SEASON[i])]]$START[tempdf$Camera[i]]
}

ndx <- which(tempdf$WEEK>tempdf$LASTWEEK)
debug <- tempdf[ndx,]
head(debug)
setwd(DataDir)
write.csv(debug,"DEBUG_for_Ana1.csv")
table(debug$Camera,debug$SEASON)

ndx <- which(tempdf$WEEK<tempdf$FIRSTWEEK)
debug2 <- tempdf[ndx,]
head(debug2)
setwd(DataDir)
write.csv(debug2,"DEBUG_for_Ana2.csv")
table(debug2$Camera,debug2$SEASON)


df <- subset(tempdf,subset=((Camera=="F2")&(SEASON=="spring")))
ndx <- grep("3898",df$Pic.No)

df[ndx,]


ndx <- grep("3898",tempdf$Pic.No)
tempdf[ndx,]


ndx <- grep("3898",this$Pic.No)   # note: there are some duplicate photos- for example, this one!
this[ndx,]


camera <- "A1"
season <- "summer"
this <- subset(tempdf,(SEASON==season)&(Camera==camera))
setwd(DataDir)
write.csv(this,"DEBUG_for_Ana3.csv")


ndx <- which(duplicated(tempdf$Pic.No))   #find duplicate pictures
setwd(DataDir)
write.csv(tempdf[ndx,],"DEBUG_for_Ana_duplicates.csv")


length(ndx)
nrow(tempdf)






#####################
# GLOBAL VARS
#####################

masterdf <- tempdf

names(masterdf)

  #occnames <- with(masterdf,  as.character(unique(PERIOD))  )
occnames <- with(masterdf,  as.character(unique(WEEK))  )        # KTS: changed to weeks
occnames <- occnames[order(as.numeric(occnames))]
nocc <- length(occnames)
  #sitenames <- with(masterdf, as.character(unique(SITE)) )
  #nsites <- length(sitenames)

masterdf$PERIOD <- masterdf$WEEK





#####################
# PROCESS DATA
#####################


#########
# site specific covariates


covariates <- data.frame(TREATMENT=rep(NA,times=nsites),GRID=NA)
rownames(covariates) = sitenames
  
r=sitenames[1]
for(r in sitenames){
  temp1 <- subset(masterdf,SITE==r)
  covariates$TREATMENT[rownames(covariates)==r] <- names(which.max(table(as.character(temp1$TREATMENT))))
  covariates$GRID[rownames(covariates)==r] <- names(which.max(table(as.character(temp1$GRID))))
}



##############
#  Fill in the template occupancy matrix with NAs for non-deployments
##############

tempoccmat <- data.frame(matrix(NA,ncol=nocc,nrow=nsites))      # make occupancy matrix
names(tempoccmat) <- occnames
rownames(tempoccmat) <- sitenames

r = "C5"
for(r in sitenames){
  #temp1 <- subset(masterdf,SITE==r)

  
  s = seasons[1]
  allweeks <- c()
  for(s in seasons){
    #temp2 <- subset(temp1,SEASON==i)
    #temp3 <- subset(masterdf,SEASON==i)
    min <- week(CameraDates[[s]][["START"]][r])
    max <- week(CameraDates[[s]][["END"]][r])
    if(!is.na(min)){
      allweeks <- unique(c(allweeks,c(min:max)))
    }
  }
  tempoccmat[r,as.character(allweeks)] <- 0
}

tempoccmat

subset(masterdf,WEEK%in%c(9:12))  # NOTE: there are lots of obs for which the date is before the start date!


      # remove unsurveyed times
tempoccmat <- tempoccmat[,-c(1)]
occnames <- occnames[-c(1)]
nocc <- nocc-1

########
# site-visit covariates

  # determine number of days of "effort"

covariates2 <- data.frame(CTDAYS=rep(NA,times=nsites*nocc))
# r=sitenames[1]
# counter=1
# for(r in sitenames){
#   c=as.numeric(occnames[1])
#   for(c in occnames){
#     temp1 <- subset(masterdf,((SITE==r)&(PERIOD==c)))
#     if(nrow(temp1)>0){
#       covariates2$CTDAYS[counter] <- as.numeric(names(which.max(table(as.character(temp1$CTDAYS)))))
#     } else{
#       temp1 <- subset(masterdf,PERIOD==c)
#       covariates2$CTDAYS[counter] <- as.numeric(names(which.max(table(as.character(temp1$CTDAYS)))))
#     }
#     counter=counter+1
#   }
# }
covariates2[] <- 1

head(covariates2)
#########
# determine occupancy

speciescode <- "XESP"
speciescode <- "CAAU"

umf.list <- list()    #set up list for unmarked package
for(speciescode in speciescodes){
  occmat <- tempoccmat
  temp1 <- subset(masterdf,CODE==speciescode)
  r=sitenames[5]
  for(r in sitenames){
    temp2 <- subset(temp1,SITE==r)
    c=occnames[1]
    for(c in occnames){
      temp3 <- subset(temp2,PERIOD==c)
      if(nrow(temp3)>0) occmat[r,c] <- 1
    }
  }

  # tempoccmat
  umf.list[[speciescode]] <- unmarkedFrameOccu(y = occmat, siteCovs = covariates, obsCovs = covariates2)
  #obsCovs(umf.list[[speciescode]])=scale(obsCovs(umf.list[[speciescode]]))
  # plot(umf)
  # summary(umf)
}

summary(umf.list[[2]])

### artificially add a burrowing owl to a control plot (aid in convergence)
umf.list[['ATCU']]@y[1,3] <- 1
summary(umf.list[['ATCU']])

# ### artificially add a white-winged dove to a pdog plot (aid in convergence)
# umf.list[['ZEAS']]@y[6,3] <- 1
# summary(umf.list[['ZEAS']])

setwd(DataDir)
filename = "OccupData4.RData"
save(umf.list,file=filename)


####################
# Total counts

speciescode <- "CYGU"

count.data <- data.frame(SPECIES=rep(speciesnames,each=2),
                         CODE=factor(rep(speciescodes,each=2),levels=levels(as.factor(masterdf$CODE))),
                         TREATMENT=rep(c("C","P"),times=length(speciesnames)),
                         COUNT=NA)
i=1
for(i in 1:nrow(count.data)){
  temp1 <- subset(masterdf,((CODE==count.data$CODE[i])&(TREATMENT==count.data$TREATMENT[i])))
  count.data$COUNT[i] <- nrow(temp1)
}


#############################
# OCCUPANCY ANALYSIS
#############################

# Note: in future try to use offset to account for number of camera deployment days.

setwd(DataDir)
load("OccupData4.RData")

fm.list <- list()

speciescode <- "CYGU"
for(speciescode in speciescodes){
   #eval(parse(text=sprintf("umf.list$%s <-
  occmod <- occuRN(~1 ~TREATMENT, umf.list[[speciescode]])
  fm.list[[speciescode]]<- occmod # CTDAYS starts=c(,)
}

setwd(DataDir)
filename = "OccupAnal4.RData"
save(umf.list,fm.list,file=filename)


#####
# Pdog model with site effect
#####

gridnames <- levels(tempdf$GRID)
pdog.grids <- c("B","D","F","G")
keep <- which(covariates$GRID%in%pdog.grids)
umf.pdog <- unmarkedFrameOccu(y=umf.list[['CYGU']]@y[keep,],siteCovs=umf.list[['CYGU']]@siteCovs[keep,],
               obsCovs=data.frame(CTDAYS=as.vector(t(matrix(as.vector(umf.list[['CYGU']]@obsCovs[,1]),nrow=nsites,byrow=T)[keep,]))))
fm.pdog <- occuRN(~1 ~GRID, umf.pdog)


#############################
# RESULTS: TREATMENT EFFECT (on-colony vs off-colony)
#############################

treatmenteff <- data.frame(mean=rep(NA,times=length(speciescodes)),lower=NA,upper=NA)
rownames(treatmenteff) <- speciesnames

s=1
for(s in 1:length(speciescodes)){
  treatmenteff$mean[s] <- as.numeric(fm.list[[speciescodes[s]]]@estimates@estimates$state@estimates['TREATMENTP'])
  temp <- confint(fm.list[[speciescodes[s]]],type="state")
  treatmenteff$lower[s] <- temp['lam(TREATMENTP)','0.025']
  treatmenteff$upper[s] <- temp['lam(TREATMENTP)','0.975']
}

treatmenteff <- treatmenteff[order(treatmenteff$mean),]    # sort in order of increasing effect

ignore <- which(treatmenteff$upper>50)     # remove crazy estimeates (non-converged)
if(length(ignore)>0) treatmenteff <- treatmenteff[-ignore,]    

### remove gunnison's prairie dog

treatmenteff <- treatmenteff[-which(rownames(treatmenteff)=="Gunnison's Prairie Dog"),]

setwd(FiguresDir)
filename = "CommunityEffect5.svg"
svg(width=4,height=4,filename=filename)
dotchart(treatmenteff$mean,labels=row.names(treatmenteff),cex=.7,
  	main="Effect of being on prairie dog colony\nvs off colony", 
   xlab="Effect size ('beta' term)",xlim=c(-4,6))

segments(treatmenteff$lower, 1:nrow(treatmenteff), treatmenteff$upper, 1:nrow(treatmenteff))
abline(v=0,col=gray(0.3),lwd=2)
dev.off()



###########################
# RESULTS: GUNNISON'S PDOG ABUNDANCE (2015)   [note: should be better once we are able to get count data...]
###########################


fm.pdog

newdata1 <- data.frame(GRID=factor(c("B","D","F"),levels=pdog.grids))
pred.ct <- predict(fm.pdog, type="state", newdata=newdata1)   # predicted abundance at all pdog grids 


   ### From mark recap analysis....    [for now, just eyeballing from figure]

pred.cmr <- data.frame(Predicted=c(5,45,33),SE=c(1.5,3,4),lower=c(3,40,31),upper=c(11,51,40))

newdf <- data.frame(means=c(pred.ct[,"Predicted"],pred.cmr[,"Predicted"]))
newdf$lowers <- c(pred.ct[,"lower"],pred.cmr[,"lower"])
newdf$uppers <- c(pred.ct[,"upper"],pred.cmr[,"upper"])
newdf$origin <- rep(c("Camera","Live Capture"),each=3)
newdf$grid <- rep(c("B","D","F"),times=2)


setwd(FiguresDir)
filename="AbundComparison2.svg"
svg(filename=filename,width=5,height=4)
ggplot(newdf, aes(x=as.factor(grid), y=means, fill=origin)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  scale_fill_manual(values=c(gray(0.2), gray(0.7))) +
  geom_errorbar(aes(ymin=lowers, ymax=uppers), width=.2,position=position_dodge(.9)) +
  ggtitle("Abundance Estimates") +
  labs(x="Trapping Grid",y="Pdog Abundance")
dev.off()
 


###########################
# RESULTS: TOTAL COUNTS BY TREATMENT (data summary)
###########################

totcount <- with(count.data, tapply(COUNT,SPECIES,sum) )
sorted <- sort(totcount,decreasing=T)

ndx <- which(count.data$TREATMENT=="P")
new <- count.data[ndx,]
newlevels <- as.character(new[order(new$COUNT,decreasing=T),]$SPECIES)
count.data$SPECIES <- factor(count.data$SPECIES,levels=newlevels)

count.data$TREATMENT <- factor(count.data$TREATMENT,levels=c("P","C"))
levels(count.data$TREATMENT) <- c("On-colony","Off-colony")

setwd(FiguresDir)
filename="RawCounts_byTreatment_with_pdog2.svg"
svg(filename=filename,width=6,height=4)

ggplot(count.data, aes(x=as.factor(SPECIES), y=COUNT, fill=TREATMENT)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  scale_fill_manual(values=c(gray(0.2), gray(0.7))) +
  #scale_y_log10() +
  ggtitle("Raw Counts (# photo captures)") +
  labs(x="Species",y="Raw Counts")

dev.off()


setwd(FiguresDir)
filename="RawCounts_byTreatment_without_pdog2.svg"
svg(filename=filename,width=6,height=4)

ggplot(count.data[-c(1,2),], aes(x=as.factor(SPECIES), y=COUNT, fill=TREATMENT)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  scale_fill_manual(values=c(gray(0.2), gray(0.7))) +
  #scale_y_log10() +
  ggtitle("Raw Counts (# photo captures)") +
  labs(x="Species",y="Raw Counts")

dev.off()




##########################
###################

# END SCRIPT















