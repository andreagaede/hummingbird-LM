## Small multiples speed tuning curves for Figure S3
## Edited for figshare on Nov 30, 2016

## may need to install tidyverse package
## install.packages("tidyverse")
library("tidyverse")

# Choose one of the directories, have to manually tell program which file to use
# dir <- INSERT PATH TO "/2015-ZF-PSTH-speed"
# dir <- INSERT PATH TO "/2015-CALAN-PSTH-speed"
setwd(dir)


## CALAN file numbers: 1, 4, 13, 21, 28, 45, 57, 64, 80, 72, 83, 86
## ZF file numbers: 115, 20, 66, 77, 80, 95, 105, 120, 32, 56, 58, 78

## CALAN file names that could be used: 
## "CALAN194-trk1-3007-sp-vary-cell1-2-0C.csv" 
## "CALAN194-trk1-3062-sp-vary-cell2-2-01.csv" 
## "CALAN198-trk3-2990-sp-vary-cell2A-2-0C.csv" 
## "CALAN199-trk1-2468-sp-vary-cell1-2-01.csv" 
## "CALAN199-trk1-2711-sp-vary-cell3-2-0F.csv" 
## "CALAN200-trk2-2974-sp-vary-cell2-2-01.csv" 
## "CALAN201-trk5-2956-sp-vary-cell2-2-04.csv" 
## "CALAN201-trk6-3056-sp-vary-cell1-2-01.csv" 
## "CALAN202-trk3-3199-sp-vary135-cell3-2-04.csv"  
## "CALAN202-trk3-2809-sp-vary-cell1-2-06.csv" 
## "CALAN202-trk3-3327-sp-vary-cell2-2-02.csv" 
## "CALAN202-trk3-3395-sp-vary-cell1-2-02.csv" 

## ZF file names that could be used: 
## "ZF1506-trk1-5746-sp-vary-cell1-2-0A.csv" 
## "ZF1423-trk1-5581-sp-vary-cell3-2-04.csv" 
## "ZF1501-trk8-6095-sp-vary-cell1-2-01.csv" 
## "ZF1503-trk7-5268-sp-vary-cell2-2-02.csv" 
## "ZF1503-trk7-5330-sp-vary-cell1-2-01.csv" 
## "ZF1504-trk4-5045-sp-vary-cell3-2-05.csv" 
## "ZF1504-trk5-5253-sp-vary-cell2-2-03.csv" 
## "ZF1510-trk7-5020-sp-vary-cell3-2-05.csv" 
## "ZF1423-trk3-5772-sp-vary-cell4-2-0C.csv" 
## "ZF1425-trk5-5221-sp-vary-cell1-2-01.csv" 
## "ZF1425-trk5-5221-sp-vary-cell3-2-04.csv" 
## "ZF1503-trk7-5268-sp-vary-cell3-2-06.csv" 


fn <- 21 # put file list number here
title=list.files()[fn]     
hist1<-read.csv(list.files()[fn]) 

names(hist1) <- c('dir', 'speed', paste('bin', as.character(1:20), sep=""))
hist1[dim(hist1)[1]+1,]<-hist1[dim(hist1)[1]-1,]
rownames(hist1)<-1:dim(hist1)[1]


rowmeans<-rowMeans(hist1[,3:22])

baseline<-rowmeans[seq(1,length(rowmeans),by=2)]

percentchange <- (rowmeans[seq(2,length(rowmeans),by=2)] - mean(baseline))/mean(baseline) * 100
data<-cbind(percentchange, hist1[seq(2,length(rowmeans),by=2),1:2])
degrees <- acos(((sqrt((((data$speed)*144*0.02768)/2)^2+30^2))^2+(sqrt((((data$speed)*144*0.02768)/2)^2+30^2))^2-((data$speed)*144*0.02768)^2)/(2*(sqrt((((data$speed)*144*0.02768)/2)^2+30^2))*(sqrt((((data$speed)*144*0.02768)/2)^2+30^2))))*180/pi

degrees<-round(degrees,2)
data<-cbind(data,degrees)

means <- aggregate(data$percentchange, by=list(data$dir,degrees), 'mean')
names(means) <- c('dir', 'degrees', 'pchange') 

max_firing <- max(means[,3])
norm <- data$percentchange/max_firing
normmeans <- cbind(data, norm)
normmeansf <- aggregate(normmeans$norm, by=list(data$dir, data$degrees), 'mean')

myse <- function(x){
  sd(x)/sqrt(length(x))
}

se <- aggregate(normmeans$norm, by=list(normmeans$dir,normmeans$degrees), 'myse')
normmeansf$se <- se$x
names(normmeansf) <- c('dir', 'degrees', 'norm', 'SEM')

norm_df <- data.frame(dir = c(normmeansf[,1]), degrees = c(normmeansf[,2]), norm = c(normmeansf[,3]), sem = c(normmeansf[,4]))

pref<-norm_df[1,1]
anti<-pref+180

# #######log scale plot for small multiples#####
# plot(log10(norm_df$degrees[norm_df$dir==pref]), norm_df$norm[norm_df$dir==pref], type ="n",xlab="", ylab="", xlim=range(log10(norm_df$degrees)), ylim=range(-1,1.2), xaxt='n', yaxt='n')
# plot(log10(norm_df$degrees[norm_df$dir==pref]), norm_df$norm[norm_df$dir==pref], type ="n",xlab="", ylab="", xlim=range(log10(norm_df$degrees)), ylim=range(-2,1.5), main = c(fn, title))
# 
# arrows(log10(norm_df$degrees[norm_df$dir==pref]), norm_df$norm[norm_df$dir==pref]+norm_df$sem[norm_df$dir==pref], log10(norm_df$degrees[norm_df$dir==pref]), norm_df$norm[norm_df$dir==pref]-norm_df$sem[norm_df$dir==pref], angle=90,code=3, length=0.015)
# arrows(log10(norm_df$degrees[norm_df$dir==anti]), norm_df$norm[norm_df$dir==anti]+norm_df$sem[norm_df$dir==anti], log10(norm_df$degrees[norm_df$dir==anti]), norm_df$norm[norm_df$dir==anti]-norm_df$sem[norm_df$dir==anti], angle=90,code=3, length=0.015)
# 
# lines(log10(norm_df$degrees[norm_df$dir==pref]), norm_df$norm[norm_df$dir==pref], col="black", type="b", pch=0, lwd=3)
# lines(log10(norm_df$degrees[norm_df$dir==anti]), norm_df$norm[norm_df$dir==anti], col="black", type="b", pch=23, lwd=3)
# abline(0,0, col = "gray60", lty = 2, lwd = 2)



norm_speed <- as_tibble(norm_df)
norm_speed <- norm_speed %>%
  rename(response = norm, direction = dir)
norm_speed$direction <- as.factor(norm_speed$direction)

limits <- aes(ymax = response + sem, ymin=response - sem)
p <- ggplot(norm_speed, aes(x = degrees, y = response, color = direction, group = direction)) +
  scale_y_continuous(limits = c(-2,1.65)) + scale_x_log10() ## y axis was c(-1.7,1.5) 

p + geom_point(aes(shape = direction, size = 2)) +
  scale_shape_manual( values = c(15,23)) +
  geom_line() +
  geom_errorbar(limits)

# geom_errorbar(limits, width = 1)
