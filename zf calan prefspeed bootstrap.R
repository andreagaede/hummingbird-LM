set.seed(1)
########################################################################
########################################################################
# Last edited by B.Goller on 26 June 2016
########################################################################
########################################################################
#Clear the workspace
rm(list = ls())
########################################################################
########################################################################
if(! "circular" %in% .packages(all=TRUE))
  install.packages("circular")

library(circular)
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
simple.tuning <- function(grps, response, threshold){
  result <- rep(0, length(grps))
  for(i in 1:length(grps)){
    if(response[i]>threshold){result[i] <- 1}
  }
  return(result)
}
########################################################################
########################################################################
########################################################################
dir.by.mean <- function(dir, firing){
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
  rayleigh.test <- matrix(NA,length(file.list),2)
  
  for(i in 1:length(file.list))
  {
    print(i)
    onecell.data <- read.csv(file.list[i], header = FALSE)
    onecell.data <- cbind(apply(onecell.data[,3:22],1,quick.mean),onecell.data)
    cell.ids[i,] <- extract.num.from.string(file.list[i])
    
    onecell.meanbase <- quick.mean(onecell.data[is.na(onecell.data[,2]),1])
    
    #Collect overall means by direction for tuning test
    onecell.nonbase <- onecell.data[!is.na(onecell.data[,2]),1:2]
    onecell.nonbase <- cbind(onecell.nonbase,(onecell.nonbase[,1]-onecell.meanbase))
    colnames(onecell.nonbase) <- c("mean.resp", "dir", "delta.resp")
    onecell.means <- aggregate(delta.resp ~ dir, data = as.data.frame(onecell.nonbase[,2:3]), quick.mean)
    onecell.error <- aggregate(delta.resp ~ dir, data = as.data.frame(onecell.nonbase[,2:3]), quick.stdev)
    onecell.maxresp <- max(onecell.means$delta.resp)
    
    #Test for uniformity
    onecell.relfiring <- onecell.means[,2]+abs(min(onecell.means[,2]))
    for.circ.density <- response.to.density(onecell.means[,1], onecell.relfiring)
    ray.test <- rayleigh.test(circular(for.circ.density,type = "angles",units = "degrees",rotation = "clock"))
    rayleigh.test[i,] <- c(round(ray.test$statistic,4), round(ray.test$p.value,4))
    
    if(ray.test$p.value<0.05){
      prefdir.sum[i,] <- dir.by.mean(onecell.means[,1], onecell.means[,2])
    }
  }
  return(cbind(cell.ids,prefdir.sum[,4],rayleigh.test))
}
########################################################################
########################################################################
quick.speed <- function(x){
  degrees <- acos(((sqrt((((x)*144*0.02768)/2)^2+30^2))^2+(sqrt((((x)*144*0.02768)/2)^2+30^2))^2-((x)*144*0.02768)^2)/(2*(sqrt((((x)*144*0.02768)/2)^2+30^2))*(sqrt((((x)*144*0.02768)/2)^2+30^2))))*180/pi
  return(round(degrees,2))
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
simulate.cell.response <- function(file.name, reps){
  raw.data <- read.csv(file.name, header = FALSE)
  mean.data <- cbind(raw.data[,1:2],apply(raw.data[,3:22],1,quick.mean))
  dirs <- unique(mean.data[!is.na(mean.data[,1]),1])
  spds <- unique(mean.data[!is.na(mean.data[,2]),2])
  
  mean.base <- quick.mean(mean.data[is.na(mean.data[,2]),3])
  cell.max <- max(abs(mean.data[,3]-mean.base))
  mean.data <- cbind(mean.data,(mean.data[,3]-mean.base)/cell.max)
  
  grabber <- matrix(NA, length(dirs)*length(spds),2+reps)
  grabber[,1] <- dirs
  for(h in 1:length(spds)){grabber[c(((h-1)*length(dirs)+1):(h*length(dirs))),2] <- rep(spds[h],length(dirs))}
  grabber <- grabber[order(grabber[,2],grabber[,1]),]
  
  for(i in 1:reps){
    for(j in 1:dim(grabber)[1]){
      subset <- mean.data[which(mean.data[,1]==grabber[j,1] & mean.data[,2]==grabber[j,2]),4]
      grabber[j,i+2] <- mean(subset[sample(1:length(subset),length(subset),replace = TRUE)])
    }
  }
  return(grabber)
}
########################################################################
########################################################################
get.prefdir.for.cell <- function(cell.info,ref.info){
  match <- which(ref.info$bird == cell.info[1] & ref.info$trk == cell.info[2] & ref.info$site == cell.info[3] & ref.info$cell == cell.info[4])
  if(length(match) == 1){
    return(ref.info$pd.sum[match])
  }
  else{return(NA)}
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
check.match <- function(test.cell, ref.cell.info){
  #Given the information of a new cell, check whether it is on the reference list
  #test.cell and ref.cell.info must have the same information: bird, trk, site, cell
  match <- FALSE
  if(length(which(test.cell[1] == ref.cell.info[,1] & test.cell[2] == ref.cell.info[,2] & test.cell[3] == ref.cell.info[,3] & test.cell[4] == ref.cell.info[,4])) == 1) {match <- TRUE}
  return(match)
}
########################################################################
########################################################################
width.around.peak <- function(value,p){
  width <- 0
  test = TRUE
  
  #Get peak
  max <- max(value)
  start <- min(which(value == max))
  
  step = start
  while(test == TRUE & step <= length(value)){
    if(value[step] >= max*p){
      width <- width + 1
      step <- step + 1
    }
    else{test = FALSE}
  }
  step = start-1
  test = TRUE
  while(test == TRUE & step >= 1){
    if(value[step] >= max*p){
      width <- width + 1
      step <- step - 1
    }
    else{test = FALSE}
  }
  return(width)
}
########################################################################
########################################################################
dir <- INSERT PATH TO "2015-CALAN-PSTH-direction/"
setwd(dir)
calan.file.list <- list.files(pattern = "CALAN.+\\.csv$", ignore.case = TRUE)
calan.prefdir <- get.preferred.direction(calan.file.list)
colnames(calan.prefdir) <- c("bird", "trk", "site", "cell", "pd.sum", "rayleigh.stat", "rayleigh.pval")
########################################################################
########################################################################
#Set the working directory (must contain data file)
dir <- INSERT PATH TO "2015-ZF-PSTH-direction/"
setwd(dir)
zf.file.list <- list.files(pattern = "ZF.+\\.csv$", ignore.case = TRUE)
zf.prefdir <- get.preferred.direction(zf.file.list)
colnames(zf.prefdir) <- c("bird", "trk", "site", "cell", "pd.sum", "rayleigh.stat", "rayleigh.pval")
########################################################################
########################################################################
#END PREFDIR CODE - NEED PREFDIR TO DETERMINE WHICH CELLS WERE APPROPRIATELY TESTED FOR SPEED TUNING
########################################################################
########################################################################
#Generate a single holding data frame for all the information
all.cell.analysis <- as.data.frame(rbind(calan.prefdir[which(calan.prefdir[,"rayleigh.pval"] < 0.05 & !is.na(calan.prefdir[,"pd.sum"])),1:5],
                                         zf.prefdir[which(zf.prefdir[,"rayleigh.pval"] < 0.05 & !is.na(zf.prefdir[,"pd.sum"])),1:5]))
all.cell.analysis["species"] <- c(rep(2,length(which(calan.prefdir[,"rayleigh.pval"] < 0.05 & !is.na(calan.prefdir[,"pd.sum"])))),
                                  rep(1,length(which(zf.prefdir[,"rayleigh.pval"] < 0.05 & !is.na(zf.prefdir[,"pd.sum"])))))
all.cell.analysis[c("test.dir1","dir1.area","test.dir2","dir2.area","pref.speed")] <- NA
########################################################################
########################################################################
########################################################################
########################################################################
dir <- INSERT PATH TO "2015-CALAN-PSTH-speed/"
setwd(dir)

calan.file.list <- list.files(pattern = "CALAN.+\\.csv$", ignore.case = TRUE)
repeats <- 1000
calan.all.cells <- matrix(NA,length(calan.file.list),1+repeats)
calan.sim.areas <- matrix(NA,repeats,2)
calan.sim.prefspeed <- matrix(NA,length(calan.file.list),2+repeats)
calan.tuning.width <- matrix(NA,length(calan.file.list),2+repeats)
calan.pref.area <- matrix(NA,length(calan.file.list),2+repeats)
calan.antipref.area <- matrix(NA,length(calan.file.list),2+repeats)

for(g in 1:length(calan.file.list)){
  cell.for.check <- NULL
  cell.for.check <- extract.num.from.string(calan.file.list[g])
  if(check.match(cell.for.check,all.cell.analysis[,1:4])){
    calan.sim.cell <- simulate.cell.response(calan.file.list[g],repeats)
    calan.sim.cell[,1] <- calan.sim.cell[,1]%%360
    calan.sim.cell[,2] <- quick.speed(calan.sim.cell[,2])
    
    #Check to see if the speed test was performed in an appropriate set of directions
    test.dirs <- unique(calan.sim.cell[,1])
    cell.prefdir <- get.prefdir.for.cell(cell.for.check,all.cell.analysis)
    diff.from.prefdir <- c(angle.diff(test.dirs[1], cell.prefdir), angle.diff(test.dirs[2], cell.prefdir))
    closest.testdir <- numeric(0)
    closest.testdir <- test.dirs[diff.from.prefdir < 45]
    
    if(length(closest.testdir) == 1){
      #Rearrange the test.dirs variable (prefdir first, then antipref)
      ifelse(test.dirs[1] == closest.testdir, test.dirs, test.dirs <- c(test.dirs[2],test.dirs[1]))
      
      calan.sim.prefspeed[g,1:2] <- c(cell.prefdir,closest.testdir)
      calan.tuning.width[g,1:2] <- c(cell.prefdir,closest.testdir)
      calan.pref.area[g,1:2] <- c(cell.prefdir,test.dirs[1])
      calan.antipref.area[g,1:2] <- c(cell.prefdir,test.dirs[2])
      for(h in 1:repeats){
        max.resp <- max(calan.sim.cell[calan.sim.cell[,1] == closest.testdir,h+2])
        calan.sim.prefspeed[g,h+2] <- min(calan.sim.cell[which(calan.sim.cell[,h+2] == max.resp),2])
        calan.tuning.width[g,h+2] <- width.around.peak(calan.sim.cell[calan.sim.cell[,1] == closest.testdir,h+2],0.8)
        calan.pref.area[g,h+2] <- sum(trapezoid.areas(c(0,calan.sim.cell[which(calan.sim.cell[,1] == test.dirs[1]),2]),
                                                      c(0,calan.sim.cell[which(calan.sim.cell[,1] == test.dirs[1]),h+2])))
        calan.antipref.area[g,h+2] <- sum(trapezoid.areas(c(0,calan.sim.cell[which(calan.sim.cell[,1] == test.dirs[2]),2]),
                                                          c(0,calan.sim.cell[which(calan.sim.cell[,1] == test.dirs[2]),h+2])))
      }
    }
  }
}
  
########################################################################
########################################################################
dir <- INSERT PATH TO "2015-ZF-PSTH-speed/"
setwd(dir)

zf.file.list <- list.files(pattern = "ZF.+\\.csv$", ignore.case = TRUE)
repeats <- 1000
zf.all.cells <- matrix(NA,length(zf.file.list),1+repeats)
zf.sim.areas <- matrix(NA,repeats,2)
zf.sim.prefspeed <- matrix(NA,length(zf.file.list),2+repeats)
zf.tuning.width <- matrix(NA,length(zf.file.list),2+repeats)
zf.pref.area <- matrix(NA,length(zf.file.list),2+repeats)
zf.antipref.area <- matrix(NA,length(zf.file.list),2+repeats)

for(g in 1:length(zf.file.list)){
  cell.for.check <- NULL
  cell.for.check <- extract.num.from.string(zf.file.list[g])
  if(check.match(cell.for.check,all.cell.analysis[,1:4])){
    zf.sim.cell <- simulate.cell.response(zf.file.list[g],repeats)
    zf.sim.cell[,1] <- zf.sim.cell[,1]%%360
    zf.sim.cell[,2] <- quick.speed(zf.sim.cell[,2])
    
    #Check to see if the speed test was performed in an appropriate set of directions
    test.dirs <- unique(zf.sim.cell[,1])
    cell.prefdir <- get.prefdir.for.cell(cell.for.check,all.cell.analysis)
    diff.from.prefdir <- c(angle.diff(test.dirs[1], cell.prefdir), angle.diff(test.dirs[2], cell.prefdir))
    closest.testdir <- numeric(0)
    closest.testdir <- test.dirs[diff.from.prefdir < 45]
      
    if(length(closest.testdir) == 1){
      #Rearrange the test.dirs variable (prefdir first, then antipref)
      ifelse(test.dirs[1] == closest.testdir, test.dirs, test.dirs <- c(test.dirs[2],test.dirs[1]))
      
      zf.sim.prefspeed[g,1:2] <- c(cell.prefdir,closest.testdir)
      zf.tuning.width[g,1:2] <- c(cell.prefdir,closest.testdir)
      zf.pref.area[g,1:2] <- c(cell.prefdir,test.dirs[1])
      zf.antipref.area[g,1:2] <- c(cell.prefdir,test.dirs[2])
      for(h in 1:repeats){
        max.resp <- max(zf.sim.cell[zf.sim.cell[,1] == closest.testdir,h+2])
        zf.sim.prefspeed[g,h+2] <- min(zf.sim.cell[which(zf.sim.cell[,h+2] == max.resp),2])
        zf.tuning.width[g,h+2] <- width.around.peak(zf.sim.cell[zf.sim.cell[,1] == closest.testdir,h+2],0.8)
        zf.pref.area[g,h+2] <- sum(trapezoid.areas(c(0,zf.sim.cell[which(zf.sim.cell[,1] == test.dirs[1]),2]),
                                                   c(0,zf.sim.cell[which(zf.sim.cell[,1] == test.dirs[1]),h+2])))
        zf.antipref.area[g,h+2] <- sum(trapezoid.areas(c(0,zf.sim.cell[which(zf.sim.cell[,1] == test.dirs[2]),2]),
                                                       c(0,zf.sim.cell[which(zf.sim.cell[,1] == test.dirs[2]),h+2])))
      }
    }
  }
}
 

########################################################################
########################################################################
speeds <- c(0.24,0.48,0.99,1.98,4.03,7.98,16.03,24.00,32.03,47.99,64.02,79.86)
log.speeds <- log10(speeds)
calan.bin.count <- matrix(NA,repeats,length(speeds))
calan.bin.tuning <- matrix(NA,repeats,dim(calan.tuning.width)[1])
zf.bin.count <- matrix(NA,repeats,length(speeds))
zf.bin.tuning <- matrix(NA,repeats,dim(zf.tuning.width)[1])
for(i in 1:repeats){
  for(j in 1:length(speeds)){
    calan.bin.count[i,j] <- length(which(calan.sim.prefspeed[,i+2] == speeds[j]))/length(which(is.na(calan.sim.prefspeed[,i+2]) == FALSE))
    zf.bin.count[i,j]    <- length(which(zf.sim.prefspeed[,i+2] == speeds[j]))/length(which(is.na(zf.sim.prefspeed[,i+2]) == FALSE))
  }
}
speed.means <- rbind(apply(calan.bin.count,2,quick.mean), apply(zf.bin.count,2,quick.mean))
speed.error <- rbind(apply(calan.bin.count,2,quick.stdev), apply(zf.bin.count,2,quick.stdev))
bp1 <- barplot(speed.means, beside = TRUE, col = c(2,1), ylim = c(0,0.3), names.arg = round(log.speeds,3),
               xlab = "speed bins (log speed)", ylab = "proportion of cell pop.", main = "pref speeds")
arrows(bp1, speed.means, bp1, speed.means + speed.error, length=0.05, angle=90, code=2, col = 1)
########################################################################
########################################################################
calan.max.width <- max(calan.tuning.width[which(is.na(calan.tuning.width[,1])==FALSE),c(3:(repeats+2))])
zf.max.width <- max(zf.tuning.width[which(is.na(zf.tuning.width[,1])==FALSE),c(3:(repeats+2))])

calan.cells <- which(is.na(calan.tuning.width[,1])==FALSE)
zf.cells <- which(is.na(zf.tuning.width[,1])==FALSE)

calan.tuning <- matrix(NA,repeats,max(calan.max.width,zf.max.width))
zf.tuning <- matrix(NA,repeats,max(calan.max.width,zf.max.width))

for(i in 1:repeats){
  for(j in 1:max(calan.max.width,zf.max.width)){
    calan.tuning[i,j] <- length(which(calan.tuning.width[calan.cells,i+2] == j))/length(calan.cells)
    zf.tuning[i,j] <- length(which(zf.tuning.width[zf.cells,i+2] == j))/length(zf.cells)
  }
}

tuning.means <- rbind(apply(calan.tuning,2,quick.mean), apply(zf.tuning,2,quick.mean))
tuning.error <- rbind(apply(calan.tuning,2,quick.stdev), apply(zf.tuning,2,quick.stdev))
bp2 <- barplot(tuning.means, beside = TRUE, col = c(2,1), ylim = c(0,1), names.arg = c(1:max(calan.max.width,zf.max.width)),
               xlab = "width bins", ylab = "proportion of cell pop.", main = "width of speed tuning peak")
arrows(bp2, tuning.means, bp2, tuning.means + tuning.error, length=0.05, angle=90, code=2, col = 1)
########################################################################
########################################################################
# ggplot speed tuning width
library(ggplot2)
#library(gridExtra)

#calan.matrix <- matrix(c(tuning.means[1,],tuning.error[1,]), nrow=10, ncol=2)
#zf.matrix <- matrix(c(tuning.means[2,],tuning.error[2,]), nrow=10, ncol=2)

#zf.df<-data.frame(zf.matrix,"zf",1:10)
#colnames(zf.df)<-c("proportion","error", "species","width")
#calan.df<-data.frame(calan.matrix,"calan",1:10)
#colnames(calan.df)<-c("proportion","error", "species","width")

#gg.means<-data.frame(rbind(calan.df,zf.df))

#ggplot(gg.means, aes(x=width,y=proportion,fill=species), color=species) +  
#  stat_summary(fun.y=mean,position=position_dodge(),geom="bar")

## bootstrapped speed tuning width graphs
#ggplot(data=gg.means, aes(x=width, y=proportion, fill=species)) + 
#  ylim(0,1) +
#  geom_bar(stat="identity", position=position_dodge(), colour="black") +
#  geom_errorbar(aes(ymin=proportion, ymax=proportion+error), width=.2, position=position_dodge(width=0.9), color="black")

#ggplot(gg.means, aes(x=width, y=proportion, fill=species)) +
#  ylim(0,1) +
#  geom_histogram(alpha=0.5, stat = "identity", position = "identity")+
#  geom_errorbar(aes(ymin=proportion, ymax=proportion+error), width=.2, position= "identity", color="dark gray")


########################################################################
########################################################################
calan.areas <- cbind(apply(calan.pref.area[calan.cells,3:repeats+2],1,quick.mean),
                     apply(calan.pref.area[calan.cells,3:repeats+2],1,quick.stdev),
                     apply(calan.antipref.area[calan.cells,3:repeats+2],1,quick.mean),
                     apply(calan.antipref.area[calan.cells,3:repeats+2],1,quick.stdev))

zf.areas <- cbind(apply(zf.pref.area[zf.cells,3:repeats+2],1,quick.mean),
                  apply(zf.pref.area[zf.cells,3:repeats+2],1,quick.stdev),
                  apply(zf.antipref.area[zf.cells,3:repeats+2],1,quick.mean),
                  apply(zf.antipref.area[zf.cells,3:repeats+2],1,quick.stdev))

plot(calan.areas[,1], calan.areas[,3], col =2, pch = 24,
     ylim = c(-80,60), xlim = c(-40,80),
     xlab = "preferred direction response area",
     ylab = "anti-preferred direction response area",
     main = "bootstrapped cell response areas by species")
arrows(calan.areas[,1]-calan.areas[,2], calan.areas[,3], calan.areas[,1]+calan.areas[,2], calan.areas[,3], length=0.05, angle=90, code=3, col = 2)
arrows(calan.areas[,1], calan.areas[,3]-calan.areas[,4], calan.areas[,1], calan.areas[,3]+calan.areas[,4], length=0.05, angle=90, code=3, col = 2)
points(zf.areas[,1], zf.areas[,3], col =1, pch = 16)
arrows(zf.areas[,1]-zf.areas[,2], zf.areas[,3], zf.areas[,1]+zf.areas[,2], zf.areas[,3], length=0.05, angle=90, code=3, col = 1)
arrows(zf.areas[,1], zf.areas[,3]-zf.areas[,4], zf.areas[,1], zf.areas[,3]+zf.areas[,4], length=0.05, angle=90, code=3, col = 1)
lines(c(0,0), c(-80,60), lty = 2)
lines(c(-40,80), c(0,0), lty = 2)

########################################################################
########################################################################
calan.prefspeed.byprefdir <- cbind(calan.sim.prefspeed[,1:2],apply(calan.sim.prefspeed[,3:repeats+2],1,quick.mean))
plot(calan.prefspeed.byprefdir[,3] ~ calan.prefspeed.byprefdir[,2])

zf.prefspeed.byprefdir <- cbind(zf.sim.prefspeed[,1:2],apply(zf.sim.prefspeed[,3:repeats+2],1,quick.mean))
plot(zf.prefspeed.byprefdir[,3] ~ zf.prefspeed.byprefdir[,2])
