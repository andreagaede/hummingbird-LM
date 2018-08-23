set.seed(1)
# set.seed(64)
########################################################################
########################################################################
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
quick.circ.mean <- function(x){return(mean(circular(x,type = "angles",units = "degrees",rotation = "clock"),na.rm = TRUE))}
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
########################################################################
simulate.cell.response <- function(file.name, reps){
  raw.data <- read.csv(file.name, header = FALSE)
  mean.data <- cbind(raw.data[,1:2],apply(raw.data[,3:22],1,quick.mean))
  dirs <- unique(mean.data[!is.na(mean.data[,1]),1])
  spds <- unique(mean.data[!is.na(mean.data[,2]),2])
  
  mean.base <- quick.mean(mean.data[is.na(mean.data[,2]),3])
  mean.data <- cbind(mean.data,mean.data[,3]-mean.base)
  
  grabber <- matrix(NA, length(dirs)*length(spds),2+reps)
  grabber[,1] <- dirs
  grabber[,2] <- spds
  
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
bootstrap.sites <- function(uni.sites, reps){
  samples <- matrix(NA,dim(uni.sites)[1],reps)
  balls.in <- dim(uni.sites)[1]
  if(reps<(balls.in^balls.in)){
    for(i in 1:reps){
      samples[,i] <- sample(1:balls.in,balls.in,replace = TRUE)
    }
    return(samples)
  }
  else{return(samples)}
}
########################################################################
########################################################################
site.grabber <- function(sampler, ref, data){
  resp.vector <- NULL
  for(i in 1:length(sampler)){
    resp.vector <- c(resp.vector, data[(data[,1] == ref[sampler[i],1] & data[,2] == ref[sampler[i],2] & data[,3] == ref[sampler[i],3]),5])
  }
  return(resp.vector)
}
########################################################################
########################################################################
run.prefdir.test <- function(angle.data) {
  circ.stat.out <- c(NA,NA,NA)
  #Angles must be in degrees and the rotation must be clockwise
  ray.test <- rayleigh.test(circular(angle.data,type = "angles",units = "degrees",rotation = "clock"))
  circ.stat.out[1] <- as.numeric(ray.test$p.value)
  if(circ.stat.out[1]<0.05){
    pref.dir <- mle.vonmises(circular(angle.data,type = "angles",units = "degrees",rotation = "clock"))
    circ.stat.out[2:3] <- c(as.numeric(pref.dir$mu), as.numeric(pref.dir$kappa))
  }
  return(circ.stat.out)
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
dir <- INSERT PATH TO "2015-CALAN-PSTH-direction/"
setwd(dir)

calan.file.list <- list.files(pattern = "CALAN.+\\.csv$", ignore.case = TRUE)
repeats <- 1000
calan.sim.prefdir <- matrix(NA,length(calan.file.list),repeats+4)
calan.nosim.prefdir <- matrix(NA,length(calan.file.list),4)

for(g in 1:length(calan.file.list)){
  #Get the overall value for the cell
  raw.data <- read.csv(calan.file.list[g], header = FALSE)
  mean.data <- cbind(raw.data[,1:2],apply(raw.data[,3:22],1,quick.mean))
  mean.base <- quick.mean(mean.data[is.na(mean.data[,2]),3])
  mean.data <- cbind(mean.data,mean.data[,3]-mean.base)
  colnames(mean.data) <- c("dir", "pattern", "mean", "relmean")
  cellfire <- aggregate(relmean~dir,mean.data,quick.mean)
  cellfire$relfire <- cellfire$relmean + abs(min(cellfire$relmean))
  onecell.circ <- response.to.density(cellfire$dir, cellfire$relfire)
  ray.test <- rayleigh.test(circular(onecell.circ,type = "angles",units = "degrees",rotation = "clock"))
  calan.nosim.prefdir[g,1] <- ray.test$p.value
  if(ray.test$p.value<0.05){
    calan.nosim.prefdir[g,3:4] <- as.numeric(mle.vonmises(circular(onecell.circ,type = "angles",units = "degrees",rotation = "clock"))[c("mu","kappa")])
    calan.nosim.prefdir[g,2] <- dir.by.mean(cellfire$dir,cellfire$relmean)[4]
  }
  #####################
  #Bootstrap
  calan.sim.cell <- simulate.cell.response(calan.file.list[g],repeats)
  calan.sim.prefdir[g,1:4] <- as.numeric(unlist(extract.num.from.string(calan.file.list[g])))
  for(h in 1:repeats){
    cell.relfire <- calan.sim.cell[,2+h] + abs(min(calan.sim.cell,2+h))
    cell.fire.circ <- response.to.density(calan.sim.cell[,1], cell.relfire)
    ray.test <- rayleigh.test(circular(cell.fire.circ,type = "angles",units = "degrees",rotation = "clock"))
    if(as.numeric(ray.test["p.value"])<0.05){
      prefdir.by.sum <- dir.by.mean(calan.sim.cell[,1],calan.sim.cell[,2+h])
      ifelse(prefdir.by.sum[4]>180, 
             calan.sim.prefdir[g,h+4] <- prefdir.by.sum[4]-360,
             calan.sim.prefdir[g,h+4] <- prefdir.by.sum[4])  
    }
  }
}

# Bootstrap the sites as well as the cells...look for the impact on population direction preference
calan.mean.sim <- cbind(calan.sim.prefdir[,1:4], apply(calan.sim.prefdir[,5:(repeats+4)],1,quick.circ.mean), t(apply(calan.sim.prefdir[,5:(repeats+4)],1,run.prefdir.test)))
calan.uni.sites <- unique(calan.mean.sim[,1:3])
calan.boot.sites <- bootstrap.sites(calan.uni.sites,repeats)
calan.site.bootstrap.means <- apply(calan.boot.sites,2,site.grabber,ref=calan.uni.sites, data=calan.mean.sim)
calan.site.bootstrap.dirs <- lapply(calan.site.bootstrap.means,run.prefdir.test)
calan.site.results <- t(matrix(as.numeric(unlist(calan.site.bootstrap.dirs)),3,repeats))

calan.site.circ <- circular(calan.site.results[,2],type = "angles",units = "degrees",rotation = "clock")
calan.site.ray.test <- rayleigh.test(calan.site.circ)
plot(calan.site.circ, shrink = 1, stack = TRUE, main = sprintf("CALAN site bootstrap pops: preferred direction (%i/%i were non-uniform)",length(which(calan.site.results[,1] < 0.05)),repeats))
rose.diag(calan.site.circ, bins = 24, add = TRUE)

calan.cell.info <- apply(calan.sim.prefdir[,5:(repeats+4)],1,quick.mean)
########################################################################
########################################################################
dir <- INSERT PATH TO "2015-ZF-PSTH-direction/"
setwd(dir)

zf.file.list <- list.files(pattern = "ZF.+\\.csv$", ignore.case = TRUE)
repeats <- 1000
zf.sim.prefdir <- matrix(NA,length(zf.file.list),repeats+4)
zf.nosim.prefdir <- matrix(NA,length(zf.file.list),4)

for(g in 1:length(zf.file.list)){
  #Get the overall value for the cell
  raw.data <- read.csv(zf.file.list[g], header = FALSE)
  mean.data <- cbind(raw.data[,1:2],apply(raw.data[,3:22],1,quick.mean))
  mean.base <- quick.mean(mean.data[is.na(mean.data[,2]),3])
  mean.data <- cbind(mean.data,mean.data[,3]-mean.base)
  colnames(mean.data) <- c("dir", "pattern", "mean", "relmean")
  cellfire <- aggregate(relmean~dir,mean.data,quick.mean)
  cellfire$relfire <- cellfire$relmean + abs(min(cellfire$relmean))
  onecell.circ <- response.to.density(cellfire$dir, cellfire$relfire)
  ray.test <- rayleigh.test(circular(onecell.circ,type = "angles",units = "degrees",rotation = "clock"))
  zf.nosim.prefdir[g,1] <- ray.test$p.value
  if(ray.test$p.value<0.05){
    zf.nosim.prefdir[g,3:4] <- as.numeric(mle.vonmises(circular(onecell.circ,type = "angles",units = "degrees",rotation = "clock"))[c("mu","kappa")])
    zf.nosim.prefdir[g,2] <- dir.by.mean(cellfire$dir,cellfire$relmean)[4]
  }
  #####################
  #Bootstrap
  zf.sim.cell <- simulate.cell.response(zf.file.list[g],repeats)
  zf.sim.prefdir[g,1:4] <- as.numeric(unlist(extract.num.from.string(zf.file.list[g])))
  for(h in 1:repeats){
    cell.relfire <- zf.sim.cell[,2+h] + abs(min(zf.sim.cell,2+h))
    cell.fire.circ <- response.to.density(zf.sim.cell[,1], cell.relfire)
    ray.test <- rayleigh.test(circular(cell.fire.circ,type = "angles",units = "degrees",rotation = "clock"))
    if(as.numeric(ray.test["p.value"])<0.05){
      prefdir.by.sum <- dir.by.mean(zf.sim.cell[,1],zf.sim.cell[,2+h])
      ifelse(prefdir.by.sum[4]>180, 
             zf.sim.prefdir[g,h+4] <- prefdir.by.sum[4]-360,
             zf.sim.prefdir[g,h+4] <- prefdir.by.sum[4])  
    }
  }
}

# Bootstrap the sites as well as the cells...look for the impact on population direction preference
zf.mean.sim <- cbind(zf.sim.prefdir[,1:4], apply(zf.sim.prefdir[,5:(repeats+4)],1,quick.circ.mean), t(apply(zf.sim.prefdir[,5:(repeats+4)],1,run.prefdir.test)))
zf.uni.sites <- unique(zf.mean.sim[,1:3])
zf.boot.sites <- bootstrap.sites(zf.uni.sites,repeats)
zf.site.bootstrap.means <- apply(zf.boot.sites,2,site.grabber,ref=zf.uni.sites, data=zf.mean.sim)
zf.site.bootstrap.dirs <- lapply(zf.site.bootstrap.means,run.prefdir.test)
zf.site.results <- t(matrix(as.numeric(unlist(zf.site.bootstrap.dirs)),3,repeats))

zf.site.circ <- circular(zf.site.results[,2],type = "angles",units = "degrees",rotation = "clock")
zf.site.ray.test <- rayleigh.test(zf.site.circ)
plot(zf.site.circ, shrink = 1, stack = TRUE, main = sprintf("ZF site bootstrap pops: preferred direction (%i/%i were non-uniform)",length(which(zf.site.results[,1] < 0.05)),repeats))
rose.diag(zf.site.circ, bins = 24, add = TRUE)

zf.cell.info <- apply(zf.sim.prefdir[,5:repeats+4],1,quick.mean)
########################################################################
########################################################################
calan.circ <- circular(calan.cell.info,type = "angles",units = "degrees",rotation = "clock")
zf.circ <- circular(zf.cell.info,type = "angles",units = "degrees",rotation = "clock")

calan.ray.test <- rayleigh.test(calan.circ)
zf.ray.test <- rayleigh.test(zf.circ)

plot(calan.circ, shrink = 1, stack = TRUE, main = "CALAN means")
rose.diag(calan.circ, bins = 24, add = TRUE)

plot(zf.circ, shrink = 1, stack = TRUE, main = "ZF means")
rose.diag(zf.circ, bins = 24, add = TRUE)
########################################################################
########################################################################
calan.sim.pops <- matrix(NA,repeats,4)
zf.sim.pops <- matrix(NA,repeats,4)


for(k in 1:repeats){
  calan.sim.pops[k,] <- c(k,as.numeric(rayleigh.test(circular(calan.sim.prefdir[,k+4],type = "angles",units = "degrees",rotation = "clock"))["p.value"]),
                          as.numeric(mle.vonmises(circular(calan.sim.prefdir[,k+4],type = "angles",units = "degrees",rotation = "clock"))[c("mu","kappa")]))
  zf.sim.pops[k,] <- c(k,as.numeric(rayleigh.test(circular(zf.sim.prefdir[,k+4],type = "angles",units = "degrees",rotation = "clock"))["p.value"]),
                       as.numeric(mle.vonmises(circular(zf.sim.prefdir[,k+4],type = "angles",units = "degrees",rotation = "clock"))[c("mu","kappa")]))
}

zf.pop.prefdirs.boot <- circular(zf.sim.pops[which(zf.sim.pops[,2] < 0.05),3],type = "angles",units = "degrees",rotation = "clock")
zf.pop.mle <- mle.vonmises(zf.pop.prefdirs.boot)
plot(zf.pop.prefdirs.boot, shrink = 2, stack = TRUE, main = sprintf("ZF Bootstrapped populations: preferred direction (%i/%i were non-uniform)",length(which(zf.sim.pops[,2] < 0.05)),repeats))
rose.diag(zf.pop.prefdirs.boot, bins = 24, add = TRUE)

calan.pop.prefdirs.boot <- circular(calan.sim.pops[which(calan.sim.pops[,2] < 0.05),3],type = "angles",units = "degrees", rotation = "clock")
calan.pop.mle <- mle.vonmises(calan.pop.prefdirs.boot)
plot(calan.pop.prefdirs.boot, shrink = 1.5, bins = 24, stack = TRUE, main = sprintf("CALAN Bootstrapped populations: preferred direction (%i/%i were non-uniform)",length(which(calan.sim.pops[,2] < 0.05)),repeats))
rose.diag(calan.pop.prefdirs.boot, bins = 24, add = TRUE)

#95% confidence intervals
zf.sorted.prefdir <- sort(zf.sim.pops[,3])
zf.sorted.prefdir[c(51,949)]


##### plot bootstrapped significant + non-significant
calan.pop.prefdirs.boot <- circular(calan.sim.pops[which(calan.sim.pops[,2] < 0.05),3],type = "angles",units = "degrees") #,rotation = "clock")
calan.pop.prefdirs.boot.non <- circular(calan.sim.pops[which(calan.sim.pops[,2] > 0.05),3],type = "angles",units = "degrees")

#plot(calan.pop.prefdirs.boot, shrink = 1.5, bins = 24, stack = TRUE, main = sprintf("CALAN Bootstrapped populations: preferred direction (%i/%i were non-uniform)",length(which(calan.sim.pops[,2] < 0.05)),repeats))


rayleigh.test(calan.pop.prefdirs.boot)
rayleigh.test(zf.pop.prefdirs.boot)

########################################################################
########################################################################
calan.circ <- circular(calan.nosim.prefdir[,2],type = "angles",units = "degrees",rotation = "clock")

calan.ray.test <- rayleigh.test(calan.circ)

plot(calan.circ, shrink = 1, stack = TRUE, main = sprintf("CALAN nosim pref dir (n = %i)",length(calan.circ[!is.na(calan.circ)])))
rose.diag(calan.circ, bins = 24, add = TRUE)
########################################################################
########################################################################
zf.circ <- circular(zf.nosim.prefdir[,2],type = "angles",units = "degrees",rotation = "clock")

zf.ray.test <- rayleigh.test(zf.circ)
zf.nosim.pop <- mle.vonmises(circular(zf.circ,type = "angles",units = "degrees",rotation = "clock"))

plot(zf.circ, shrink = 1, stack = TRUE, main = sprintf("ZF nosim pref dir (n = %i)",length(zf.circ[!is.na(zf.circ)])))
rose.diag(zf.circ, bins = 24, add = TRUE)
