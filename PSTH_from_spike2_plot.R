library(readr)
library(tidyverse)

dir <- INSERT PATH TO "PSTH_examples_figureS1"
setwd(dir)

## 0 degrees
psth_0_deg_100ms <- read_csv("psth_0_deg_100ms.csv")
names(psth_0_deg_100ms) <- c('time', 'spikes')
barplot(psth_0_deg_100ms$spikes, col = "black", ylim = c(0, 25), main =  "0 degrees")

## 45 degrees
psth_45_deg_100ms <- read_csv("psth_45_deg_100ms.csv")
names(psth_45_deg_100ms) <- c('time', 'spikes')
barplot(psth_45_deg_100ms$spikes, col = "black", ylim = c(0, 25), main =  "45 degrees")

## 90 degrees
psth_90_deg_100ms <- read_csv("psth_90_deg_100ms.csv")
names(psth_90_deg_100ms) <- c('time', 'spikes')
barplot(psth_90_deg_100ms$spikes, col = "black", ylim = c(0, 25), main =  "90 degrees")


## 135 degrees
psth_135_deg_100ms <- read_csv("psth_135_deg_100ms.csv")
names(psth_135_deg_100ms) <- c('time', 'spikes')
barplot(psth_135_deg_100ms$spikes, col = "black", ylim = c(0, 25), main =  "135 degrees")

## 180 degrees
psth_180_deg_100ms <- read_csv("psth_180_deg_100ms.csv")
names(psth_180_deg_100ms) <- c('time', 'spikes')
barplot(psth_180_deg_100ms$spikes, col = "black", ylim = c(0, 25), main =  "180 degrees")

## 225 degrees
psth_225_deg_100ms <- read_csv("psth_225_deg_100ms.csv")
names(psth_225_deg_100ms) <- c('time', 'spikes')
barplot(psth_225_deg_100ms$spikes, col = "black", ylim = c(0, 25), main =  "225 degrees")

## 270 degrees
psth_270_deg_100ms <- read_csv("psth_270_deg_100ms.csv")
names(psth_270_deg_100ms) <- c('time', 'spikes')
barplot(psth_270_deg_100ms$spikes, col = "black", ylim = c(0, 25), main =  "270 degrees")

## 315 degrees
psth_315_deg_100ms <- read_csv("psth_315_deg_100ms.csv")
names(psth_315_deg_100ms) <- c('time', 'spikes')
barplot(psth_315_deg_100ms$spikes, col = "black", ylim = c(0, 25), main =  "315 degrees")


