## Created by A. Gaede Oct 3, 2016
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

## ZF file names could be used: 
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


fn <- 44 ## put file number here
title = list.files()[fn]     
hist1 <- read.csv(list.files()[fn]) 

names(hist1) <- c('direction', 'speed', paste('bin', as.character(1:20), sep=""))

row_means <- rowMeans(hist1[,3:22])

baseline <- row_means[seq(1,length(row_means),by=2)]
baseline <- mean(baseline)

motion_firing <- row_means[seq(2,length(row_means), by = 2)]
firing_less_baseline <- motion_firing - baseline
max_firing <- max(abs(firing_less_baseline))

my_data <- cbind(firing_less_baseline, hist1[seq(2,length(row_means),by=2),1:2])
degrees <- acos(((sqrt((((my_data$speed)*144*0.02768)/2)^2+30^2))^2+(sqrt((((my_data$speed)*144*0.02768)/2)^2+30^2))^2-((my_data$speed)*144*0.02768)^2)/(2*(sqrt((((my_data$speed)*144*0.02768)/2)^2+30^2))*(sqrt((((my_data$speed)*144*0.02768)/2)^2+30^2))))*180/pi
degrees <- round(degrees,2)
my_data <- cbind(my_data,degrees)

my_data <- my_data %>% 
  mutate(norm_firing = firing_less_baseline/max_firing) 
  
norm_means <- aggregate(my_data$norm_firing, by=list(my_data$direction,degrees), 'mean')

names(norm_means) <- c('direction', 'degrees', 'norm_means') 

## Function for calculating SEM
myse <- function(x){
  sd(x)/sqrt(length(x))
}

se <- aggregate(my_data$norm_firing, by=list(my_data$direction,my_data$degrees), 'myse')
norm_means$se <- se$x
names(norm_means) <- c('direction', 'degrees', 'norm_means', 'SEM')

######### Graphing code:

norm_speed <- as_tibble(norm_means)
norm_speed$direction <- as.factor(norm_speed$direction)

limits <- aes(ymax = norm_means + SEM, ymin = norm_means - SEM)
p <- ggplot(norm_speed, aes(x = degrees, y = norm_means, color = direction, group = direction)) +
  scale_y_continuous(limits = c(-1.0, 1.0)) + scale_x_log10() 

p + geom_point(aes(shape = direction, size = 2)) +
  scale_shape_manual( values = c(15,23)) +
  geom_line() +
  geom_errorbar(limits)

## for quadrant plot sample curves
## ZF file numbers: 77, 56
## CALAN file number: 35

## "CALAN199-trk3-2889-sp-vary-cell3-2-01.csv"
## "ZF1503-trk7-5268-sp-vary-cell2-2-02.csv"
## "ZF1425-trk5-5221-sp-vary-cell1-2-01.csv"

# p2 <- ggplot(norm_speed, aes(x = degrees, y = norm_means, color = direction, group = direction)) +
#   scale_y_continuous(limits = c(-1.0, 1.0))
# 
# p2 + geom_point(aes(shape = direction, size = 2)) +
#   scale_shape_manual( values = c(15,23)) +
#   geom_line() +
#   geom_errorbar(limits, width = 1)

