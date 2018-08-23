########################################################################
########################################################################
# Last edited by B.Goller on 26 June 2016
#To do:
########################################################################
########################################################################
#Clear the workspace
rm(list = ls())
########################################################################
########################################################################
if(! "circular" %in% .packages(all=TRUE))
  install.packages("circular")
if(! "nlme" %in% .packages(all=TRUE))
  install.packages("nlme")
if(! "multcomp" %in% .packages(all=TRUE))
  install.packages("multcomp")
if(! "lme4" %in% .packages(all=TRUE))
  install.packages("lme4")

library(circular)
library(nlme)
library(multcomp)
library(lme4)
########################################################################
########################################################################
#DEFINE FUNCTIONS
########################################################################
########################################################################
quick.mean <- function(x){return(mean(x,na.rm = TRUE))}
########################################################################
########################################################################
quick.stdev <- function(x){return(sd(x,na.rm = TRUE))}
########################################################################
########################################################################
quick.sum <- function(x){return(sum(x,na.rm = TRUE))}
########################################################################
########################################################################
response.to.density <- function(dir, firing) {
  #Use the firing to create repeats of the direction value
  output <- NULL
  
  num.dirs <- unique(dir)
  rep.min <- min(as.numeric(lapply(num.dirs, FUN = function(x) length(which(dir == x)))))
  
  for(i in 1:(length(num.dirs)*rep.min)) {
    output <- c(output,rep(dir[i],round(10*firing[i])))
  }
  return(output)
}
########################################################################
########################################################################
extract.num.from.string <- function(s) {
  s.split <- strsplit(s, "cell")
  s.id <- as.numeric(unlist(strsplit(s.split[[1]][1], "[^[:digit:]]")))
  s.id <- s.id[!is.na(s.id)][1:3]
  
  s.cell <- as.numeric(unlist(strsplit(s.split[[1]][2], "[^[:digit:]]")))
  s.cell <- s.cell[!is.na(s.cell)][1]
  return(c(s.id, s.cell))
}
########################################################################
########################################################################
dir.by.mean <- function(dir, firing){
  v <- rep(0,8)
  h <- rep(0,8)
  
  unique.dir <- unique(dir)
  
  for(i in 1:length(unique.dir)){
    #Vertical component of firing vector
    v[i] <- sum(firing[which(dir == unique.dir[i])]*sin(-unique.dir[i]*pi/180))
    #Horizontal component of firing vector
    h[i] <- sum(firing[which(dir == unique.dir[i])]*cos(-unique.dir[i]*pi/180))
  }
  #calculate preferred direction as arctangent of the sums of vertical and horizontal components
  pd <- atan(sum(v)/sum(h))*(180/pi)
  
  ifelse(sum(v)>0, ifelse(sum(h)>0, pd.adj<-360-pd, pd.adj <- 180+abs(pd)),
         ifelse(sum(h)>0,pd.adj <- abs(pd),pd.adj <- 180-abs(pd)))
  
  return(c(sum(v),sum(h),pd, pd.adj))
}

########################################################################
########################################################################
dir.by.mean2 <- function(dir, firing){
  v <- rep(0,8)
  h <- rep(0,8)
  
  unique.dir <- unique(dir)
  rep.min <- min(as.numeric(lapply(unique.dir, FUN = function(x) length(which(dir == x)))))
  print(rep.min)
  
  for(i in 1:length(unique.dir)){
    #Vertical component of firing vector
    v[i] <- sum(firing[which(dir == unique.dir[i])[1:rep.min]]*sin(-unique.dir[i]*pi/180))
    #Horizontal component of firing vector
    h[i] <- sum(firing[which(dir == unique.dir[i])[1:rep.min]]*cos(-unique.dir[i]*pi/180))
    #rep.test[i] <- length(firing[which(dir == unique.dir[i])])
  }
  #calculate preferred direction as arctangent of the sums of vertical and horizontal components
  pd <- atan(sum(v)/sum(h))*(180/pi)
  
  ifelse(sum(v)>0, ifelse(sum(h)>0, pd.adj<-360-pd, pd.adj <- 180+abs(pd)),
         ifelse(sum(h)>0,pd.adj <- abs(pd),pd.adj <- 180-abs(pd)))
  
  return(c(sum(v),sum(h),pd, pd.adj))
}
########################################################################
########################################################################
get.preferred.direction <- function(file.list){
  percent <- seq(0.5,0.95,by = 0.05)
  cell.ids <- matrix(NA,length(file.list),4)
  prefdir.sum <- matrix(NA,length(file.list),4)
  tuning.sums <- matrix(NA,length(file.list),length(percent))
  rayleigh.test <- matrix(NA,length(file.list),2)
  
  for(i in 1:length(file.list))
  {
    print(i)
    onecell.data <- read.csv(file.list[i], header = FALSE)
    onecell.data <- cbind(apply(onecell.data[,3:22],1,quick.mean),onecell.data)
    cell.ids[i,] <- extract.num.from.string(file.list[i])
    
    onecell.meanbase <- quick.mean(onecell.data[is.na(onecell.data[,2]),1])
    
    #    if(max(abs(onecell.data[!is.na(onecell.data[,2]),1] - onecell.meanbase)) > (onecell.meanbase * 0.2)) {
    #prefdir.sum[i,] <- dir.by.mean2(onecell.data[!is.na(onecell.data[,2]),2],(onecell.data[!is.na(onecell.data[,2]),1]-onecell.meanbase))
    
    #Collect overall means by direction for tuning test
    onecell.nonbase <- onecell.data[!is.na(onecell.data[,2]),1:2]
    onecell.nonbase <- cbind(onecell.nonbase,(onecell.nonbase[,1]-onecell.meanbase))
    colnames(onecell.nonbase) <- c("mean.resp", "dir", "delta.resp")
    onecell.means <- aggregate(delta.resp ~ dir, data = as.data.frame(onecell.nonbase[,2:3]), quick.mean)
    onecell.error <- aggregate(delta.resp ~ dir, data = as.data.frame(onecell.nonbase[,2:3]), quick.stdev)
    onecell.maxresp <- max(onecell.means$delta.resp)
    
    prefdir.sum[i,] <- dir.by.mean2(onecell.means[,1], onecell.means[,2])
    
    #Test for uniformity
    #onecell.relfiring <- onecell.data[!is.na(onecell.data[,2]),1]-onecell.meanbase + abs(min(onecell.data[!is.na(onecell.data[,2]),1]-onecell.meanbase))
    onecell.relfiring <- onecell.means[,2]+abs(min(onecell.means[,2]))
    for.circ.density <- response.to.density(onecell.means[,1], onecell.relfiring)
    ray.test <- rayleigh.test(circular(for.circ.density,type = "angles",units = "degrees",rotation = "clock"))
    rayleigh.test[i,] <- c(round(ray.test$statistic,4), round(ray.test$p.value,4))
    
    #Test tuning
    if(rayleigh.test[i,2] < 0.05){
      tuning.hold <- matrix(NA,length(percent),9)
      for(j in 1:length(percent)){
        tuning.hold[j,] <- c(percent[j],simple.tuning(onecell.means[,1],onecell.means[,2],onecell.maxresp*percent[j]))
      }
      tuning.sums[i,] <- apply(tuning.hold[,2:9], 1, quick.sum)
    }
    #    }
  }
  tuning.sums[which(tuning.sums[] == 0)] <- NA  
  #col1 is the sum method, col2 is the mle method
  return(cbind(cell.ids,prefdir.sum[,4],rayleigh.test,tuning.sums))
}
########################################################################
########################################################################
quick.speed <- function(x){
  degrees <- acos(((sqrt((((x)*144*0.02768)/2)^2+30^2))^2+(sqrt((((x)*144*0.02768)/2)^2+30^2))^2-((x)*144*0.02768)^2)/(2*(sqrt((((x)*144*0.02768)/2)^2+30^2))*(sqrt((((x)*144*0.02768)/2)^2+30^2))))*180/pi
  return(round(degrees,2))
}
########################################################################
########################################################################
simple.tuning <- function(grps, response, threshold){
  result <- rep(0, length(grps))
  for(i in 1:length(grps)){
    if(response[i]>threshold){result[i] <- 1}
  }
  return(result)
}
########################################################################
########################################################################
extract.num.from.string <- function(s) {
  s.split <- strsplit(s, "cell")
  s.id <- as.numeric(unlist(strsplit(s.split[[1]][1], "[^[:digit:]]")))
  s.id <- s.id[!is.na(s.id)][1:3]
  
  s.cell <- as.numeric(unlist(strsplit(s.split[[1]][2], "[^[:digit:]]")))
  s.cell <- s.cell[!is.na(s.cell)][1]
  return(c(s.id, s.cell))
}

########################################################################
########################################################################
trapezoid.areas <- function(speed, firing)
{
  trapezoids <- matrix(0,length(speed)-1,1)
  if(length(speed) == length(firing)) {
    firing = firing - firing[1]
    for(i in 2:length(speed)) {
      width <- speed[i]-speed[i-1]  
      if(firing[i]>firing[i-1]) {
        tri <- 0.5*width*(firing[i]-firing[i-1])
        rect <- width*firing[i-1]
        area <- tri+rect
      }
      else {
        tri <- 0.5*width*(firing[i-1]-firing[i])
        rect <- width*firing[i]
        area <- tri+rect
      }
      trapezoids[i-1] <- area
    }
    return(trapezoids)
  }
  else {print("***Speeds and firing do NOT match!***")}
}
########################################################################
########################################################################
get.prefdir.for.cell <- function(cell.info){
  return(all.cell.analysis$pd.sum[which(all.cell.analysis$bird == cell.info[1] & all.cell.analysis$trk == cell.info[2] & all.cell.analysis$site == cell.info[3] & all.cell.analysis$cell == cell.info[4])])
}
########################################################################
########################################################################
angle.diff <- function(theta1, theta2){
  delta.theta <- abs(theta1-theta2)
  if(delta.theta > 180){delta.theta <- 360-delta.theta}
  return(abs(delta.theta))
}
########################################################################
########################################################################
get.preferred.speed <- function(file.list){
  percent <- seq(0.5,0.95,by = 0.05)
  cell.ids <- matrix(NA,length(file.list),4)
  pref.area <- matrix(NA,length(file.list),4)
  norm.pref.area <- matrix(NA,length(file.list),4)
  pref.speed <- rep(NA,length(file.list))
  areas <- matrix(NA,12,2)
  norm.areas <- matrix(NA,12,2)
  tuning.sums <- matrix(NA,length(file.list),length(percent))
  tuning.80 <- matrix(NA,length(file.list),12)
  
  for(i in 1:length(file.list)){
    cell.for.check <- NULL
    print(i)
    cell.for.check <- extract.num.from.string(file.list[i])
    if(check.match(cell.for.check,all.cell.analysis[,1:4])){
      onecell.data <- read.csv(file.list[i], header = FALSE)
      onecell.data <- cbind(apply(onecell.data[,3:22],1,quick.mean),onecell.data)
      cell.ids[i,] <- cell.for.check
      
      onecell.meanbase <- quick.mean(onecell.data[is.na(onecell.data[,2]),1])
      
      onecell.nonbase <- onecell.data[!is.na(onecell.data[,2]),1:3]
      onecell.max <- max(abs(onecell.nonbase[,1]-onecell.meanbase))
      onecell.nonbase <- cbind(onecell.nonbase,(onecell.nonbase[,1]-onecell.meanbase)/onecell.max)
      colnames(onecell.nonbase) <- c("mean.resp", "dir", "speed", "norm.to.max")
      onecell.nonbase[,"dir"] <- onecell.nonbase[,"dir"]%%360
      
      #Now get the area under the response using trapezoid.areas() then select a preferred direction
      onecell.means <- aggregate(norm.to.max ~ dir + speed, data = as.data.frame(onecell.nonbase[,2:4]), quick.mean)
      onecell.error <- aggregate(norm.to.max ~ dir + speed, data = as.data.frame(onecell.nonbase[,2:4]), quick.stdev)
      onecell.means$speed <- quick.speed(onecell.means$speed)
      onecell.error$speed <- quick.speed(onecell.error$speed)
      onecell.dirs <- unique(onecell.means$dir)
      
      #Get prefdir from the reference list
      onecell.prefdir <- get.prefdir.for.cell(cell.ids[i,])
      
      #Check if one of the tested directions is within 45 degrees of the prefdir
      diff.from.prefdir <- c(angle.diff(onecell.dirs[1], onecell.prefdir), angle.diff(onecell.dirs[2], onecell.prefdir))
      closest.testdir <- numeric(0)
      closest.testdir <- onecell.dirs[diff.from.prefdir < 45]
      
      if(length(closest.testdir) == 1){
        print("Cell was tested appropriately")
        #Rearrange the onecell.dirs variable (prefdir first, then antipref)
        ifelse(onecell.dirs[1] == closest.testdir, onecell.dirs, onecell.dirs <- c(onecell.dirs[2],onecell.dirs[1]))
        pref.speed[i] <- onecell.means$speed[onecell.means$norm.to.max == max(onecell.means$norm.to.max[onecell.means$dir == closest.testdir])]
        
        for(j in 1:length(onecell.dirs)){
          areas[,j] <- trapezoid.areas(c(0,onecell.means$speed[onecell.means$dir == onecell.dirs[j]]),
                                       c(0,onecell.means$norm.to.max[onecell.means$dir == onecell.dirs[j]]))
        }
        
        #Store the areas
        pref.area[i,] <- c(onecell.dirs[1]%%360,sum(areas[,1]),onecell.dirs[2]%%360,sum(areas[,2]))
        #Normalize the area by step size
        steps <- diff(c(0,unique(onecell.means$speed)))
        norm.areas <- areas*max(steps)/steps
        norm.pref.area[i,] <- c(onecell.dirs[1]%%360,sum(norm.areas[,1]),onecell.dirs[2]%%360,sum(norm.areas[,2]))
        
        #Test tuning
        tuning.hold <- matrix(NA,length(percent),13)
        onecell.max.mean <- max(onecell.means$norm.to.max[onecell.means$dir%%360 == pref.area[i,1]])
        for(j in 1:length(percent)){
          tuning.hold[j,] <- c(percent[j],simple.tuning(onecell.means$speed[onecell.means$dir%%360 == pref.area[i,1]],onecell.means$norm.to.max[onecell.means$dir%%360 == pref.area[i,1]],onecell.max.mean*percent[j]))
        }
        tuning.sums[i,] <- apply(tuning.hold[,2:13], 1, quick.sum)
        #tuning.80[i,] <- tuning.hold[which(tuning.hold[,1] == 0.8),2:13]
        #instead of getting the tuning, just get the mean firing values for each prefdir speed
        tuning.80[i,] <- onecell.means$norm.to.max[which(onecell.means$dir == onecell.dirs[1])]
        if(max(tuning.80[i,])==max(abs(tuning.80[i,]))){
          tuning.80[i,] <- tuning.80[i,]/max(tuning.80[i,])
        }
        else{
          tuning.80[i,] <- rep(NA,12)
        }
      }
      else{
        pref.area[i,] <- rep(NA,4)
        norm.pref.area[i,] <- rep(NA,4)
        pref.speed[i] <- NA
        tuning.80[i,] <- rep(NA,12)
      }
    }
  }
  tuning.sums[which(tuning.sums[] == 0)] <- NA 
  return(cbind(cell.ids,pref.area,pref.speed,norm.pref.area,tuning.sums,tuning.80))
}
########################################################################
########################################################################
check.match <- function(test.cell, ref.cell.info){
  #Given the information of a new cell, check whether it is on the reference list
  #test.cell and ref.cell.info must have the same information: bird, trk, site, cell
  match <- FALSE
  if(length(which(test.cell[1] == ref.cell.info[,1] & test.cell[2] == ref.cell.info[,2] & test.cell[3] == ref.cell.info[,3] & test.cell[4] == ref.cell.info[,4])) == 1) {match <- TRUE}
  return(match)
}
########################################################################
########################################################################
dir <- INSERT PATH TO "2015-CALAN-PSTH-direction/"
setwd(dir)
calan.file.list <- list.files(pattern = "CALAN.+\\.csv$", ignore.case = TRUE)
calan.prefdir.tuned <- get.preferred.direction(calan.file.list)
colnames(calan.prefdir.tuned) <- c("bird", "trk", "site", "cell", "pd.sum", "rayleigh.stat", "rayleigh.pval",
                                   "50tune", "55tune", "60tune", "65tune","70tune", "75tune", "80tune", "85tune","90tune", "95tune")
########################################################################
########################################################################
#Set the working directory (must contain data file)
dir <- INSERT PATH TO "2015-ZF-PSTH-direction/"
setwd(dir)
zf.file.list <- list.files(pattern = "ZF.+\\.csv$", ignore.case = TRUE)
zf.prefdir.tuned <- get.preferred.direction(zf.file.list)
colnames(zf.prefdir.tuned) <- c("bird", "trk", "site", "cell", "pd.sum", "rayleigh.stat", "rayleigh.pval",
                                "50tune", "55tune", "60tune", "65tune","70tune", "75tune", "80tune", "85tune","90tune", "95tune")
########################################################################
########################################################################
#END PREFDIR CODE

all.cell.analysis <- as.data.frame(rbind(calan.prefdir.tuned[which(calan.prefdir.tuned[,"rayleigh.pval"] < 0.05 & !is.na(calan.prefdir.tuned[,"pd.sum"])),1:5],
                                         zf.prefdir.tuned[which(zf.prefdir.tuned[,"rayleigh.pval"] < 0.05 & !is.na(zf.prefdir.tuned[,"pd.sum"])),1:5]))
all.cell.analysis["species"] <- c(rep(2,length(which(calan.prefdir.tuned[,"rayleigh.pval"] < 0.05 & !is.na(calan.prefdir.tuned[,"pd.sum"])))),
                                  rep(1,length(which(zf.prefdir.tuned[,"rayleigh.pval"] < 0.05 & !is.na(zf.prefdir.tuned[,"pd.sum"])))))
all.cell.analysis[c("test.dir1","dir1.area","test.dir2","dir2.area","pref.speed")] <- NA
all.cell.analysis[c("sp80.1", "sp80.2", "sp80.3", "sp80.4", "sp80.5", "sp80.6", "sp80.7", "sp80.8", "sp80.9", "sp80.10", "sp80.11", "sp80.12")] <- NA
all.cell.analysis[c("tune50","tune55","tune60","tune65","tune70","tune75","tune80","tune85","tune90","tune95")] <- NA
########################################################################
########################################################################
########################################################################
########################################################################
#Cell match code
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################
dir <- INSERT PATH TO "2015-CALAN-PSTH-speed/"
setwd(dir)

calan.file.list <- list.files(pattern = "CALAN.+\\.csv$", ignore.case = TRUE)
calan.prefspeeds.tuned <- get.preferred.speed(calan.file.list)
colnames(calan.prefspeeds.tuned) <- c("bird", "trk", "site", "cell",
                                      "pref.dir","pref.area", "antipref.dir", "antipref.area", "pref.speed",
                                      "norm.pref.dir","norm.pref.area", "norm.antipref.dir", "norm.antipref.area",
                                      "tune50","tune55","tune60","tune65","tune70","tune75","tune80","tune85","tune90","tune95",
                                      "80.1", "80.2", "80.3", "80.4", "80.5", "80.6", "80.7", "80.8", "80.9", "80.10", "80.11", "80.12"
                                      )
########################################################################
########################################################################
dir <- INSERT PATH TO "2015-ZF-PSTH-speed/"
setwd(dir)

zf.file.list <- list.files(pattern = "ZF.+\\.csv$", ignore.case = TRUE)
zf.prefspeeds.tuned <- get.preferred.speed(zf.file.list)
colnames(zf.prefspeeds.tuned) <- c("bird", "trk", "site", "cell",
                                   "pref.dir","pref.area", "antipref.dir", "antipref.area", "pref.speed",
                                   "norm.pref.dir","norm.pref.area", "norm.antipref.dir", "norm.antipref.area",
                                   "tune50","tune55","tune60","tune65","tune70","tune75","tune80","tune85","tune90","tune95",
                                   "80.1", "80.2", "80.3", "80.4", "80.5", "80.6", "80.7", "80.8", "80.9", "80.10", "80.11", "80.12"
                                   )
########################################################################
########################################################################
for(i in 1:dim(all.cell.analysis)[1]){
  row <- which(calan.prefspeeds.tuned[,"bird"] == all.cell.analysis$bird[i] & calan.prefspeeds.tuned[,"trk"] == all.cell.analysis$trk[i] & calan.prefspeeds.tuned[,"site"] == all.cell.analysis$site[i] & calan.prefspeeds.tuned[,"cell"] == all.cell.analysis$cell[i])
  if(length(row) == 1){
    all.cell.analysis[i,7:33] <- calan.prefspeeds.tuned[row,c(5:9,24:35,14:23)]
  }
  if(length(row) == 0){
    row <- which(zf.prefspeeds.tuned[,"bird"] == all.cell.analysis$bird[i] & zf.prefspeeds.tuned[,"trk"] == all.cell.analysis$trk[i] & zf.prefspeeds.tuned[,"site"] == all.cell.analysis$site[i] & zf.prefspeeds.tuned[,"cell"] == all.cell.analysis$cell[i])
    if(length(row) == 1){
      all.cell.analysis[i,7:33] <- zf.prefspeeds.tuned[row,c(5:9,24:35,14:23)]
    }
  }
}
########################################################################
########################################################################
all.cell.analysis <- within(all.cell.analysis, {bird <- factor(bird)
                                                trk <- factor(trk)
                                                site <- factor(site)
                                                cell <- factor(cell)
                                                species <- factor(species)
                                              })
########################################################################
########################################################################
speeds <- c(0.24,0.48,0.99,1.98,4.03,7.98,16.03,24.00,32.03,47.99,64.02,79.86)
log.speeds <- log10(speeds)
calan.speedbins <- apply(all.cell.analysis[which(all.cell.analysis$species == 2),12:23],2,quick.sum)/length(which(all.cell.analysis$species == 2 & is.na(all.cell.analysis$dir1.area) == FALSE))
zf.speedbins <- apply(all.cell.analysis[which(all.cell.analysis$species == 1),12:23],2,quick.sum)/length(which(all.cell.analysis$species == 1 & is.na(all.cell.analysis$dir1.area) == FALSE))

barplot(rbind(calan.speedbins,zf.speedbins), beside = TRUE, col = c(2,1),
        xlab = "speed bins (log speed)", ylab = "proportion of cell pop.", names.arg = round(log.speeds,3))
########################################################################
########################################################################
#Data for speed fitting
cal.prop.tested <- all.cell.analysis[which(as.numeric(all.cell.analysis$species) == 2 & !is.na(all.cell.analysis$test.dir1)),c((1:5),7,(11:23))]
cal.speed <- all.cell.analysis[which(as.numeric(all.cell.analysis$species) == 2 & !is.na(all.cell.analysis$sp80.1)),c((1:5),7,(11:23))]
zf.prop.tested <- all.cell.analysis[which(as.numeric(all.cell.analysis$species) == 1 & !is.na(all.cell.analysis$test.dir1)),c((1:5),7,(11:23))]
zf.speed <- all.cell.analysis[which(as.numeric(all.cell.analysis$species) == 1 & !is.na(all.cell.analysis$sp80.1)),c((1:5),7,(11:23))]

#create a column with all the cells
cal.rows <- dim(cal.speed)[1]*length(speeds)
cal.df <- data.frame(species = rep(2,cal.rows), bird = rep(0,cal.rows), trk = rep(0,cal.rows), site = rep(0,cal.rows), cell = rep(0,cal.rows),
                     speed = rep(0,cal.rows), fire = rep(0,cal.rows))
for(i in 1:dim(cal.speed)[1]){
  for(j in 1:length(speeds)){
    cal.df[((i-1)*12+j),2:7] <- c(cal.speed[i,1:4],speeds[j],cal.speed[i,7+j])
  }
}
cal.df <- within(cal.df, {species <- factor(species)
                          bird <- factor(bird)
                          trk <- factor(trk)
                          site <- factor(site)
                          cell <- factor(cell)
                          })





zf.rows <- dim(zf.speed)[1]*length(speeds)
zf.df <- data.frame(species = rep(1,zf.rows), bird = rep(0,zf.rows), trk = rep(0,zf.rows), site = rep(0,zf.rows), cell = rep(0,zf.rows),
                     speed = rep(0,zf.rows), fire = rep(0,zf.rows))
for(i in 1:dim(zf.speed)[1]){
  for(j in 1:length(speeds)){
      zf.df[((i-1)*12+j),2:7] <- c(zf.speed[i,1:4],speeds[j],zf.speed[i,7+j])
  }
}
zf.df <- within(zf.df, {species <- factor(species)
                        bird <- factor(bird)
                        trk <- factor(trk)
                        site <- factor(site)
                        cell <- factor(cell)
                        })

all.df <- rbind(zf.df,cal.df)
##################################################################################################
#Try logistic regression
#Need to have the binary version of the analysis for this to work of course...
#bin.glm <- glmer(fire~speed+species+(1|bird)+(1|trk)+(1|site)+(1|cell), data = all.df, family = binomial)
#summary(bin.glm)
#######################################################################################################################
#Define sigmoidal functions
firing  <- function(asym, slope,center,int, speed.bin)     {int + asym / (1 + exp(-slope * (speed.bin-center)))}
##################################################################################################
#Split the data sets for hb and zf
#Fit very simple models that just vary the slope
#NEED TUNING VERSION 2 - NORMALIZED FIRING...see the speed analysis
#CALAN analysis
cal.df.means <- aggregate(fire~speed, data = cal.df, quick.mean)
cal.asym <- max(cal.df.means$fire)
cal.int <- cal.df.means$fire[1]
calan.sig.nlme <- nlme(fire ~ firing(a,s,c,0,speed),
                       data = cal.df,
                       fixed = a+s+c ~ 1,
                       random = s ~ 1|bird/trk/site/cell,
                       start = c(a =cal.asym,s = 0.05, c = 30))
plot(fire~speed, data = cal.df.means)
points(speeds,firing(mean(coef(calan.sig.nlme)$a), mean(coef(calan.sig.nlme)$s), mean(coef(calan.sig.nlme)$c), 0,speeds), type = 'l')

#fit 2
calan.sig.nlme2 <- nlme(fire ~ firing(a,s,c,0,speed),
                        data = cal.df,
                        fixed = a+s+c ~ 1,
                        random = c ~ 1|bird/trk/site/cell,
                        start = c(a = 0.75,s = 0.01, c = 30))
points(speeds,firing(mean(coef(calan.sig.nlme2)$a), mean(coef(calan.sig.nlme2)$s), mean(coef(calan.sig.nlme2)$c), 0,speeds), type = 'l', col = 2)

#ZF analysis
zf.df.means <- aggregate(fire~speed, data = zf.df, quick.mean)
zf.asym <- max(zf.df.means$fire)
zf.sig.nlme <- nlme(fire ~ firing(a,s,c,0,speed),
                    data = zf.df,
                    fixed = a+s+c ~ 1,
                    random = s ~ 1|bird/trk/site/cell,
                    start = c(a = zf.asym,s = 0.01, c = 60))
plot(fire~speed, data = zf.df.means)
points(speeds,firing(mean(coef(zf.sig.nlme)$a), mean(coef(zf.sig.nlme)$s), mean(coef(zf.sig.nlme)$c), 0,speeds), type = 'l')

#fit 2
zf.sig.nlme2 <- nlme(fire ~ firing(a,s,c,0,speed),
                     data = zf.df,
                     fixed = a+s+c ~ 1,
                     random = c ~ 1|bird/trk/site/cell,
                     start = c(a = 0.8,s = 0.1, c=20))
points(speeds,firing(mean(coef(zf.sig.nlme2)$a), mean(coef(zf.sig.nlme2)$s), mean(coef(zf.sig.nlme2)$c), 0,speeds), type = 'l', col = 2)

zf.grab <- aggregate(fire~cell+site+trk+bird+species, data = zf.df, quick.mean)
hb.grab <- aggregate(fire~cell+site+trk+bird+species, data = cal.df, quick.mean)

zf.coef.summary <- cbind(zf.grab,coef(zf.sig.nlme)$s,coef(zf.sig.nlme2)$c)
colnames(zf.coef.summary)[7:8] <- c("slope", "center")
hb.coef.summary <- cbind(hb.grab,coef(calan.sig.nlme)$s,coef(calan.sig.nlme2)$c)
colnames(hb.coef.summary)[7:8] <- c("slope", "center")

all.coefs <- rbind(zf.coef.summary,hb.coef.summary)
center.lme <- lme(center~species,random=~1|bird/trk/site/cell,data=all.coefs)
anova(center.lme)
slope.lme <- lme(slope~species,random=~1|bird/trk/site/cell,data=all.coefs)
anova(slope.lme)
fire.lme <- lme(fire~species,random=~1|bird/trk/site/cell,data=all.coefs)
anova(fire.lme)

boxplot(center~species, data = all.coefs, main = "centerpoint", ylab = "speed (deg/s)", xlab = "species")
boxplot(slope~species, data = all.coefs, main = "slope", ylab = "spikes/deg", xlab = "species")
