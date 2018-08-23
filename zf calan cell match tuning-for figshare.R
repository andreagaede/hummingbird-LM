########################################################################
########################################################################
# Last edited by A.Gaede on 7 July 2016
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

library(circular)
library(nlme)
library(multcomp)
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
        tuning.80[i,] <- tuning.hold[which(tuning.hold[,1] == 0.8),2:13] # changed to 0.7 for 70%, 0.9 for 90%
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
dir <- INSERT PATH TO "2015-CALAN-PSTH-direction"
setwd(dir)
calan.file.list <- list.files(pattern = "CALAN.+\\.csv$", ignore.case = TRUE)
calan.prefdir.tuned <- get.preferred.direction(calan.file.list)
colnames(calan.prefdir.tuned) <- c("bird", "trk", "site", "cell", "pd.sum", "rayleigh.stat", "rayleigh.pval",
                                   "50tune", "55tune", "60tune", "65tune","70tune", "75tune", "80tune", "85tune","90tune", "95tune")
########################################################################
########################################################################
#Set the working directory (must contain data file)
dir <- INSERT PATH TO "2015-ZF-PSTH-direction"
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
dir <- INSERT PATH TO "2015-CALAN-PSTH-speed"
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
dir <- INSERT PATH TO "2015-ZF-PSTH-speed"
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
for(i in 1:nrow(all.cell.analysis)){
  row <- which(calan.prefspeeds.tuned[,"bird"] == all.cell.analysis$bird[i] & calan.prefspeeds.tuned[,"trk"] == all.cell.analysis$trk[i] & calan.prefspeeds.tuned[,"site"] == all.cell.analysis$site[i] & calan.prefspeeds.tuned[,"cell"] == all.cell.analysis$cell[i])
  print(sprintf("i = %d, matches = %d", i,length(row)))
  if(length(row) == 1){
    all.cell.analysis[i,7:33] <- calan.prefspeeds.tuned[row,c(5:9,24:35,14:23)]
  }
  if(length(row) == 0){
    row <- which(zf.prefspeeds.tuned[,"bird"] == all.cell.analysis$bird[i] & zf.prefspeeds.tuned[,"trk"] == all.cell.analysis$trk[i] & zf.prefspeeds.tuned[,"site"] == all.cell.analysis$site[i] & zf.prefspeeds.tuned[,"cell"] == all.cell.analysis$cell[i])
    print(sprintf("2: i = %d, matches = %d", i,length(row)))
    if(length(row) == 1){
      all.cell.analysis[i,7:33] <- zf.prefspeeds.tuned[row,c(5:9,24:35,14:23)]
    }
  }
}
########################################################################
########################################################################
#Generate quadrant plot

all.cell.rows <- which(is.na(all.cell.analysis$test.dir1) == FALSE) 
# excludes cells that were dropped between direction and speed test (i.e. cells that fall in top left quadrant)
excl <- which((all.cell.analysis$dir1.area < 0) & (all.cell.analysis$dir2.area > 0))
all.cell.rows <- all.cell.rows[!(all.cell.rows %in% excl)]

plot(all.cell.analysis$dir1.area[all.cell.rows],
     all.cell.analysis$dir2.area[all.cell.rows], 
     col = all.cell.analysis$species[all.cell.rows],
     ylim = c(-80,60), xlim = c(-40,80),
     xlab = "preferred direction response area",
     ylab = "anti-preferred direction response area",
     main = "cell response areas by species"
     )
lines(c(0,0), c(-80,60), lty = 2)
lines(c(-40,80), c(0,0), lty = 2)

#######################################################
# ggplot quadrant
#######################################################
# library(ggplot2)
# 
# all.cell.rows <- which(is.na(all.cell.analysis$test.dir1) == FALSE) 
# # excludes cells that were dropped between direction and speed test (i.e. cells that fall in top left quadrant)
# excl <- which((all.cell.analysis$dir1.area < 0) & (all.cell.analysis$dir2.area > 0))
# all.cell.rows <- all.cell.rows[!(all.cell.rows %in% excl)]
# all.cell.df <- data.frame(all.cell.analysis$species[all.cell.rows], all.cell.analysis$dir1.area[all.cell.rows], all.cell.analysis$dir2.area[all.cell.rows])
# 
# #all.cell.df <- data.frame(all.cell.analysis$species, all.cell.analysis$dir1.area, all.cell.analysis$dir2.area)
# colnames(all.cell.df) <- c("species", "area1", "area2")
# 
# ggplot(data = all.cell.df, aes(y=all.cell.df$area2, x=all.cell.df$area1, shape = all.cell.df$species)) + 
#   geom_point(size = 3) +
#   scale_shape_manual(values = c(1,2)) # Use open circle and triangle


########################################################################
########################################################################
all.cell.analysis <- within(all.cell.analysis, {bird <- factor(bird)
                                                trk <- factor(trk)
                                                site <- factor(site)
                                                cell <- factor(cell)
                                                species <- factor(species)
                                              })

dir1.area.lme <- lme(dir1.area ~ species, random = ~1|bird/trk/site/cell, data = all.cell.analysis, na.action = na.exclude)
anova(dir1.area.lme)

dir2.area.lme <- lme(dir2.area ~ species, random = ~1|bird/trk/site/cell, data = all.cell.analysis, na.action = na.exclude)
anova(dir2.area.lme)
########################################################################
########################################################################
#Generate circular plots for each of the species
calan.circ <- circular(all.cell.analysis$pd.sum[which(all.cell.analysis$species == 2)],type = "angles",units = "degrees",rotation = "clock")
zf.circ <- circular(all.cell.analysis$pd.sum[which(all.cell.analysis$species == 1)],type = "angles",units = "degrees",rotation = "clock")

plot(calan.circ, shrink = 1, stack = TRUE, main = sprintf("CALAN preferred direction (n = %i)",length(which(all.cell.analysis$species == 2))))
rose.diag(calan.circ, bins = 24, add = TRUE)
calan.raytest <- rayleigh.test(calan.circ)

plot(zf.circ, shrink = 1, stack = TRUE, main = sprintf("ZF preferred direction (n = %i)",length(which(all.cell.analysis$species == 1))))
rose.diag(zf.circ, bins = 24, add = TRUE)
zf.raytest <- rayleigh.test(zf.circ)
zf.mletest <- mle.vonmises(zf.circ)
########################################################################
########################################################################
# 80% max firing speed tuning 
speeds <- c(0.24,0.48,0.99,1.98,4.03,7.98,16.03,24.00,32.03,47.99,64.02,79.86)
log.speeds <- log10(speeds)
calan.speedbins <- apply(all.cell.analysis[which(all.cell.analysis$species == 2),12:23],2,quick.sum)/length(which(all.cell.analysis$species == 2 & is.na(all.cell.analysis$dir1.area) == FALSE))
zf.speedbins <- apply(all.cell.analysis[which(all.cell.analysis$species == 1),12:23],2,quick.sum)/length(which(all.cell.analysis$species == 1 & is.na(all.cell.analysis$dir1.area) == FALSE))

barplot(rbind(calan.speedbins,zf.speedbins), beside = TRUE, col = c(2,1),
        xlab = "speed bins (log speed)", ylab = "proportion of cell pop.", main = "pref speeds at 80% max firing", names.arg = round(log.speeds,3))
########################################################################
########################################################################
bin.count <- matrix(NA,2,length(speeds))
for(j in 1:length(speeds)){
  bin.count[1:2,j] <- c(length(which(all.cell.analysis$species == 2 & all.cell.analysis$pref.speed == speeds[j]))/length(which(all.cell.analysis$species == 2 & is.na(all.cell.analysis$dir1.area) == FALSE)),
                        length(which(all.cell.analysis$species == 1 & all.cell.analysis$pref.speed == speeds[j]))/length(which(all.cell.analysis$species == 1 & is.na(all.cell.analysis$dir1.area) == FALSE)))
}

barplot(bin.count, beside = TRUE, col = c(2,1),
        xlab = "speed bins (log speed)", ylab = "proportion of cell pop.", main = "pref speeds", names.arg = round(log.speeds,3))
########################################################################
########################################################################
#speed tuning in zf vs. hb
#tuning in 24:33
#cell ids in species, bird, trk, site, cell, pd.sum
plot(all.cell.analysis$tune80[which(all.cell.analysis$species == 1)] ~ all.cell.analysis$pd.sum[which(all.cell.analysis$species == 1)])
points(all.cell.analysis$tune80[which(all.cell.analysis$species == 2)] ~ all.cell.analysis$pd.sum[which(all.cell.analysis$species == 2)], col = 2)

tuning.levels <- seq(0.5,0.95, by = 0.05)
all.tuning.for.lme <- data.frame(species = numeric(0), bird = numeric(0), trk = numeric(0), site = numeric(0), 
                                 cell = numeric(0), pref.dir = numeric(0), level = numeric(0), tuning = numeric(0), level.count = numeric(0))
for(b in 1:dim(all.cell.analysis)[1]){
  all.tuning.for.lme <- rbind(all.tuning.for.lme,
                              cbind(all.cell.analysis[rep(b,10),c(6,1:5)], level = tuning.levels, tuning = as.numeric(all.cell.analysis[b,24:33]), level.count = c(1:10)))
}
all.tuning.for.lme <- within(all.tuning.for.lme, {species <- factor(species)
                                                  bird <- factor(bird)
                                                  trk <- factor(trk)
                                                  site <- factor(site)
                                                  cell <- factor(cell)
                                                  })
speed.tuning.lme <- lme(tuning ~ level*species, random = ~1|bird/trk/site/cell, data = all.tuning.for.lme, na.action = na.exclude)
anova(speed.tuning.lme)
summary(speed.tuning.lme)

speed.zf.tuning.lme <- lme(tuning ~ level, random = ~1|bird/trk/site/cell, data = all.tuning.for.lme[which(all.tuning.for.lme$species == 1),], na.action = na.exclude)
anova(speed.zf.tuning.lme)
summary(speed.zf.tuning.lme)
zf.coef <- coef(speed.zf.tuning.lme)

speed.calan.tuning.lme <- lme(tuning ~ level, random = ~1|bird/trk/site/cell, data = all.tuning.for.lme[which(all.tuning.for.lme$species == 2),], na.action = na.exclude)
anova(speed.calan.tuning.lme)
summary(speed.calan.tuning.lme)
calan.coef <- coef(speed.calan.tuning.lme)

plot(seq(0,1,by=0.05), mean(calan.coef[,2])*seq(0,1,by=0.05)+mean(calan.coef[,1]), ylim = c(0,14), col = 2)
points(seq(0,1,by=0.05), mean(zf.coef[,2])*seq(0,1,by=0.05)+mean(zf.coef[,1]), pch = 2)

##################
# ggplot boxplot #
##################
library(ggplot2)

all.tuning.df <- data.frame(all.tuning.for.lme$species, all.tuning.for.lme$level, all.tuning.for.lme$tuning)
colnames(all.tuning.df) <- c("species", "level", "tuning")
all.tuning.df$level<- as.factor(all.tuning.df$level)

all.tuning.df$species.c <- c("calan")
all.tuning.df$species.c[which(all.tuning.df$species==1)]<-c("zf")

p <- ggplot(data = all.tuning.df, aes(y=all.tuning.df$tuning, x=all.tuning.df$level, fill = all.tuning.df$species.c)) + scale_fill_hue(l=40)
p + geom_boxplot(outlier.size = 0) + geom_point(aes(colour=all.tuning.df$species.c, shape=all.tuning.df$species, size = 2), position = position_jitterdodge(jitter.height = 0.3)) +  
  stat_summary(fun.y=mean, geom = "point", shape = 23, size = 4)

#######################################################################
########### speed tuning histograms for supplementary figure ##########
#######################################################################

# make new data frame for histograms in ggplot
library(ggplot2)

sp.tuning.df <- data.frame(all.cell.analysis$species,all.cell.analysis$tune70,all.cell.analysis$tune80,all.cell.analysis$tune90)
colnames(sp.tuning.df) <- c("species","tune70","tune80","tune90")

sp.tuning.df<- na.omit(sp.tuning.df)


## histograms for speed tuning width at 70% 80% and 90%
sp.tuning.df$tune70<- as.factor(sp.tuning.df$tune70)
sp.tuning.df$tune80<- as.factor(sp.tuning.df$tune80)
sp.tuning.df$tune90<- as.factor(sp.tuning.df$tune90)

# build proportion of population plot
c.counts <- matrix(NA,10,3)

for(r in 1:10){
  c.counts[r,1]<- length(which(sp.tuning.df$tune70[which(sp.tuning.df$species==2)]==r))
  c.counts[r,2]<- length(which(sp.tuning.df$tune80[which(sp.tuning.df$species==2)]==r))
  c.counts[r,3]<- length(which(sp.tuning.df$tune90[which(sp.tuning.df$species==2)]==r))
}  

z.counts <- matrix(NA,10,3)
for(t in 1:10){
  z.counts[t,1]<- length(which(sp.tuning.df$tune70[which(sp.tuning.df$species==1)]==t))
  z.counts[t,2]<- length(which(sp.tuning.df$tune80[which(sp.tuning.df$species==1)]==t))
  z.counts[t,3]<- length(which(sp.tuning.df$tune90[which(sp.tuning.df$species==1)]==t))
}

# get proportion of population with each tuning width for each species at all 3 thresholds
z.counts.pro<-z.counts/107 #divide by number of zf cells
c.counts.pro<-c.counts/56 #divide by number of calan cells

z.counts.pro<-data.frame(z.counts.pro)
z.counts.pro<-cbind("zf", 1:10, z.counts.pro)
colnames(z.counts.pro)<-c("species","width","tune70","tune80","tune90")

c.counts.pro<-data.frame(c.counts.pro)
c.counts.pro<-cbind("calan", 1:10, c.counts.pro)
colnames(c.counts.pro)<-c("species","width","tune70","tune80","tune90")

# combine both data frames, so zf and calan data can be combined in histograms
all.counts.pro<-data.frame(rbind(c.counts.pro, z.counts.pro))

# make bar plots - side by side
ggplot(data=all.counts.pro, aes(x=width, y=tune70, fill=species)) + 
  ylim(0,0.6) +
  xlab("speed tuning width at 70% threshold") + ylab("proportion of population") + 
  geom_bar(stat="identity", position=position_dodge(), colour="black")

ggplot(data=all.counts.pro, aes(x=width, y=tune80, fill=species)) + 
  ylim(0,0.6) +
  xlab("speed tuning width at 80% threshold") + ylab("proportion of population") + 
  geom_bar(stat="identity", position=position_dodge(), colour="black")

ggplot(data=all.counts.pro, aes(x=width, y=tune90, fill=species)) + 
  ylim(0,0.6) +
  xlab("speed tuning width at 90% threshold") + ylab("proportion of population") + 
  geom_bar(stat="identity", position=position_dodge(), colour="black")


########################################################################
########################################################################
#get the counts for slow (<6deg/s) and fast (>6deg/s) cells
cell.speed.classes <- matrix(NA,2,5)

for(c in 1:2){
  cell.speed.classes[c,] <- c(c,length(which(all.cell.analysis$species == c & all.cell.analysis$pref.speed < 6)),
                              length(which(all.cell.analysis$species == c & all.cell.analysis$pref.speed > 6)),
                              length(which(all.cell.analysis$species == c & is.na(all.cell.analysis$pref.speed))),
                              length(which(all.cell.analysis$species == c)))
}
########################################################################
########################################################################

## quadrant plot ggplot

all.cell.rows <- which(is.na(all.cell.analysis$test.dir1) == FALSE) 
# excludes cells that were dropped between direction and speed test (i.e. cells that fall in top left quadrant)
excl <- which((all.cell.analysis$dir1.area < 0) & (all.cell.analysis$dir2.area > 0))
all.cell.rows <- all.cell.rows[!(all.cell.rows %in% excl)]
all.cell.df <- data.frame(all.cell.analysis$species[all.cell.rows], all.cell.analysis$dir1.area[all.cell.rows], all.cell.analysis$dir2.area[all.cell.rows])

#all.cell.df <- data.frame(all.cell.analysis$species, all.cell.analysis$dir1.area, all.cell.analysis$dir2.area)
colnames(all.cell.df) <- c("species", "area1", "area2")

ggplot(data = all.cell.df, aes(y=all.cell.df$area2, x=all.cell.df$area1, shape = all.cell.df$species)) + 
  geom_point(size = 3) +
  scale_shape_manual(values = c(1,2)) # Use open circle and triangle


########################################################################
########################################################################
# Last edited by B.Goller on 02 November 2015
########################################################################
########################################################################
#Clear the workspace
#rm(list = ls())
########################################################################
########################################################################
if(! "circular" %in% .packages(all=TRUE))
  install.packages("circular")
if(! "nlme" %in% .packages(all=TRUE))
  install.packages("nlme")
if(! "multcomp" %in% .packages(all=TRUE))
  install.packages("multcomp")

library(circular)
library(nlme)
library(multcomp)
########################################################################
########################################################################
########################################################################
########################################################################
#dir <- "C:/Users/Benny Goller/Documents/R programming/neurovis analysis/cellmatch data/2015-CALAN-PSTH-speed/"
#dir <- "C:/Users/Flightlab/Documents/Benny/R programming/neurovis analysis/wylie data/"
#dir <- "~/Documents/DATA/zf hb pg ephys"
#dir <- "~/Documents/DATA/Scripts for LM 2016 paper/scripts"
dir <- "~/Documents/DATA/Data and analysis LM 2016 paper/scripts"

setwd(dir)

pigeon.data <- read.csv("LM_dir_and_vel_pigeon_BG_edit.csv", header = TRUE)
########################################################################
########################################################################
########################################################################
########################################################################
unique.pigeon.speeds <- sort(unique(pigeon.data$velocity.from.peak[which(!is.na(pigeon.data$velocity.from.peak))]))
log.pigeon.speeds <- log10(unique.pigeon.speeds)
bin.count.pigeon <- rep(NA,length(unique.pigeon.speeds))
for(k in 1:length(unique.pigeon.speeds)){
  bin.count.pigeon[k] <- length(which(pigeon.data$velocity.from.peak == unique.pigeon.speeds[k]))/length(which(!is.na(pigeon.data$velocity.from.peak)))
}

#barplot(bin.count.pigeon, names.arg = round(log.pigeon.speeds,3))
########################################################################
########################################################################
all.speeds <- sort(unique(round(c(speeds,unique.pigeon.speeds),1)))
log.all.speeds <- log10(all.speeds)
all.bin.count <- matrix(NA,3,length(all.speeds))
for(l in 1:length(all.speeds)){
  all.bin.count[1:3,l] <- c(length(which(all.cell.analysis$species == 2 & round(all.cell.analysis$pref.speed,1) == all.speeds[l]))/length(which(all.cell.analysis$species == 2 & is.na(all.cell.analysis$dir1.area) == FALSE)),
                            length(which(all.cell.analysis$species == 1 & round(all.cell.analysis$pref.speed,1) == all.speeds[l]))/length(which(all.cell.analysis$species == 1 & is.na(all.cell.analysis$dir1.area) == FALSE)),
                            length(which(pigeon.data$velocity.from.peak == all.speeds[l]))/length(which(!is.na(pigeon.data$velocity.from.peak))))
  #length(which(round(pigeon.data$velocity.from.peak,1) == all.speeds[l]))/length(which(!is.na(pigeon.data$velocity.from.peak))))
}

barplot(all.bin.count, beside = TRUE, col = c(2,1,4),
        xlab = "speed bins (log speed)", ylab = "proportion of cell pop.", main = "pref speeds", names.arg = round(log.all.speeds,3))


# pref speed vs pref direction for hummingbirds and zf
plot(all.cell.analysis$pref.speed[which(all.cell.analysis$species == 1)] ~ all.cell.analysis$pd.sum[which(all.cell.analysis$species == 1)])
points(all.cell.analysis$pref.speed[which(all.cell.analysis$species == 2)] ~ all.cell.analysis$pd.sum[which(all.cell.analysis$species == 2)], col = 2)

