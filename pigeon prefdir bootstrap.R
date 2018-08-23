########################################################################
########################################################################
# Last edited by B.Goller on 06 June 2016
# To do: 
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
angle.diff <- function(theta1, theta2){
  delta.theta <- abs(theta1-theta2)
  if(delta.theta > 180){delta.theta <- 360-delta.theta}
  return(delta.theta)
}
########################################################################
########################################################################
diff.grid.mean <- function(x){
  grid <- matrix(NA,length(x),length(x))
  for(i in 1:length(x)){
    for(j in 1:length(x)){
      grid[i,j] <- angle.diff(x[i],x[j])
    }
  }
  return(mean(grid[upper.tri(grid)]))
}
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
dir.by.mean <- function(dir, firing){
  v <- rep(0,8)
  h <- rep(0,8)
  
  unique.dir <- unique(dir)
  rep.min <- min(as.numeric(lapply(unique.dir, FUN = function(x) length(which(dir == x)))))
  #print(rep.min)
  
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
simulate.cell.response <- function(file.name, reps){
  raw.data <- read.csv(file.name, header = FALSE)
  mean.data <- cbind(raw.data[,1:2],apply(raw.data[,3:22],1,quick.mean))
  dirs <- unique(mean.data[!is.na(mean.data[,1]),1])
  spds <- unique(mean.data[!is.na(mean.data[,2]),2])
  print(length(dirs))
  mean.base <- quick.mean(mean.data[is.na(mean.data[,2]),3])
  mean.data <- cbind(mean.data,mean.data[,3]-mean.base)
  
  grabber <- matrix(NA, length(dirs)*length(spds),3+reps)
  grabber[,1] <- dirs
  grabber[,2] <- spds
  
  #check the sampling function arguments in other bootstrapping scripts...
  for(i in 1:reps){
    for(j in 1:nrow(grabber)){
      subset <- mean.data[which(mean.data[,1]==grabber[j,1] & mean.data[,2]==grabber[j,2]),4]
      grabber[j,3] <- length(subset)
      grabber[j,i+3] <- mean(subset[sample(1:length(subset),length(subset),replace = TRUE)])
    }
  }
  return(grabber)
}
########################################################################
########################################################################
dir <- INSERT PATH TO "PIGEON-PSTH-direction/"
setwd(dir)

set.seed(100)

pigeon.file.list <- list.files(pattern = "\\.csv$", ignore.case = TRUE)
repeats <- 1000
pigeon.sim.prefdir <- matrix(NA,length(pigeon.file.list),repeats)
pigeon.sim.prefdir2 <- matrix(NA,length(pigeon.file.list),repeats)
pigeon.nosim.prefdir <- matrix(NA,length(pigeon.file.list),4)

for(g in 1:length(pigeon.file.list)){
  print(g)
  #Get the overall value for the cell
  raw.data <- read.csv(pigeon.file.list[g], header = FALSE)
  mean.data <- cbind(raw.data[,1:2],apply(raw.data[,3:22],1,quick.mean))
  mean.base <- quick.mean(mean.data[is.na(mean.data[,2]),3])
  mean.data <- cbind(mean.data,mean.data[,3]-mean.base)
  colnames(mean.data) <- c("dir", "pattern", "mean", "relmean")
  cellfire <- aggregate(relmean~dir,mean.data,quick.mean)
  cellfire$relfire <- cellfire$relmean + abs(min(cellfire$relmean))
  onecell.circ <- response.to.density(cellfire$dir+270, cellfire$relfire)
  ray.test <- rayleigh.test(circular(onecell.circ,type = "angles",units = "degrees",rotation = "clock"))
  pigeon.nosim.prefdir[g,1] <- ray.test$p.value
  if(ray.test$p.value<0.05){
    pigeon.nosim.prefdir[g,3:4] <- as.numeric(mle.vonmises(circular(onecell.circ,type = "angles",units = "degrees",rotation = "clock"))[c("mu","kappa")])
    pigeon.nosim.prefdir[g,2] <- dir.by.mean(cellfire$dir+270,cellfire$relmean)[4]
  }
  #####################
  #Bootstrap
  pigeon.sim.cell <- simulate.cell.response(pigeon.file.list[g],repeats)
  #The plus-270 degrees rotates D. Wylie's coordinate system (up = 0, forward = 90)
  #to match ours (forward = 0, up = 270)
  for(h in 1:repeats){
    cell.relfire <- pigeon.sim.cell[,3+h] + abs(min(pigeon.sim.cell[,3+h]))
    cell.fire.circ <- response.to.density(pigeon.sim.cell[,1]+270, cell.relfire)
    ray.test <- rayleigh.test(circular(cell.fire.circ,type = "angles",units = "degrees",rotation = "clock"))
    if(as.numeric(ray.test["p.value"])<0.05){
      pigeon.sim.prefdir2[g,h] <- as.numeric(mle.vonmises(circular(cell.fire.circ,type = "angles",units = "degrees",rotation = "clock"))["mu"])
      prefdir.by.sum <- dir.by.mean((pigeon.sim.cell[,1]+270),pigeon.sim.cell[,2+h])
      ifelse(prefdir.by.sum[4]>180, 
             pigeon.sim.prefdir[g,h] <- prefdir.by.sum[4]-360,
             pigeon.sim.prefdir[g,h] <- prefdir.by.sum[4])  
    }
  }
}
########################################################################
########################################################################
pigeon.circ <- circular(pigeon.nosim.prefdir[,2],type = "angles",units = "degrees",rotation = "clock")

pigeon.ray.test <- rayleigh.test(pigeon.circ)
nosim.pop <- mle.vonmises(circular(pigeon.circ,type = "angles",units = "degrees",rotation = "clock"))

plot(pigeon.circ, shrink = 1, stack = TRUE, main = sprintf("PIGEON nosim pref dir (n = %i)",length(pigeon.circ[!is.na(pigeon.circ)])))
rose.diag(pigeon.circ, bins = 24, add = TRUE)
########################################################################
########################################################################
pigeon.sim.pops <- matrix(NA,repeats,4)

for(k in 1:repeats){
  pigeon.sim.pops[k,] <- c(k,as.numeric(rayleigh.test(circular(pigeon.sim.prefdir[,k],type = "angles",units = "degrees",rotation = "clock"))["p.value"]),
                          as.numeric(mle.vonmises(circular(pigeon.sim.prefdir[,k],type = "angles",units = "degrees",rotation = "clock"))[c("mu","kappa")]))
}

pigeon.pop.prefdirs.boot <- circular(pigeon.sim.pops[which(pigeon.sim.pops[,2] < 0.05),3],type = "angles",units = "degrees",rotation = "clock")
pigeon.pop.mle <- mle.vonmises(pigeon.pop.prefdirs.boot)
plot(pigeon.pop.prefdirs.boot, shrink = 1.5, bins = 24, stack = TRUE, main = sprintf("pigeon Bootstrapped populations: preferred direction (%i/%i were non-uniform)",length(which(pigeon.sim.pops[,2] < 0.05)),repeats))
rose.diag(pigeon.pop.prefdirs.boot, bins = 24, add = TRUE)

sorted.prefdir <- sort(pigeon.sim.pops[,3])
sorted.prefdir[c(51,949)]
