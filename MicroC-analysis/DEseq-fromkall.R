library("tximport")
library("readr")
library("tximportData")
library("rhdf5")
dir <- system.file("extdata", package="tximportData")
samples <- read.table(file.path(
