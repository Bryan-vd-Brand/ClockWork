#!/usr/bin/env Rscript
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)
library(ggforce)
library(RColorBrewer)

#KnockoutReport is a tsv of with header
# Sample    SampleName	Reference	Unmodified	Reads_aligned	Knockout	Reason	Insertions	Deletions	Substitutions

koFile <- Sys.glob(file.path("KnockoutReport.tsv"))
koDF <- read.table(file = koFile, sep = "\t", header = 1)
koDF <- koDF %>% mutate(X = gsub('_','',substring(SampleName, first = 2,last=3))) %>% mutate(Y = substring(SampleName, first = 1, last = 1)) %>% mutate(X=factor(X,levels=unique(str_sort(X, numeric=TRUE)))) %>% mutate(Y=factor(Y,levels=unique(str_sort(Y, decreasing = TRUE, numeric=TRUE))))
meanVal <- mean(koDF$Reads_aligned)
koDF <- koDF %>% mutate(sizeScaler = (2 + ((Reads_aligned / meanVal) * 4)))

wellPlot <- ggplot(koDF, aes(x = X, y = Y, fill = Reason)) +
    geom_dotplot(binaxis = "y",stackdir = "center", dotsize = 2, stackgroups = TRUE) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ylab("WellRow") + xlab("WellColumn") +
    scale_color_manual(name="", values = c("yellowgreen","darkorange1","maroon","black")) +
    scale_fill_manual(name="", values = c("yellowgreen","darkorange1","maroon","black")) +
    facet_wrap(~Sample)

ggsave(wellPlot, file = "wellPlot.pdf", units = "cm", width = 50, height = 10, limitsize = FALSE)

readCountPlot <- ggplot(koDF, aes(x = X, y = Y, color = Reason)) +
    geom_point(size = koDF$sizeScaler) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    ylab("WellRow") + xlab("WellColumn") +
    scale_color_manual(name="", values = c("yellowgreen","darkorange1","maroon","black")) +
    scale_fill_manual(name="", values = c("yellowgreen","darkorange1","maroon","black")) +
    facet_wrap(~Reference)

ggsave(readCountPlot, file = "readCountPlot.pdf", units = "cm", width = 50, height = 10, limitsize = FALSE)
