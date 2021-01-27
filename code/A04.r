# --------------------------------------------------------------------
#  _____.r
#  Author: James Robert White, PhD
#  Email: jwhite@respherabio.com
# --------------------------------------------------------------------
#  description of program
# --------------------------------------------------------------------
rm(list = ls())
library(RColorBrewer)
library(scales)
library(ggplot2)
require(gtools)
library(MASS)
# ------------------------------------------------------------------------
# set the working directory as the location of the script
commandArgsArray <- commandArgs(trailingOnly = F)
fileEntry        <- commandArgsArray[grep("--file",commandArgsArray)]
scriptPath       <- dirname(sub("--file=", "", fileEntry))
scriptPrefix     <- gsub(pattern=".r$", replacement="", x=tail(strsplit(fileEntry, "/")[[1]],n=1), ignore.case=TRUE)
scriptPrefix     <- sub("--file=", "", scriptPrefix)
setwd(scriptPath)
# ------------------------------------------------------------------------
# create output directory
analysisdir = paste("../analysis/", scriptPrefix, sep="")
unlink(analysisdir, recursive=TRUE)
dir.create(analysisdir)
# ------------------------------------------------------------------------
