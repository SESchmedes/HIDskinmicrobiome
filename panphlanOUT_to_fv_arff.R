#!/usr/bin/Rscript
##############################################
#Authored by Sarah E. Schmedes
#Date Released: 08/01/2017
#This script converts a panphlan output into a feature vector
# and .arff file for Weka
##############################################
library(foreign)
library(GetoptLong)

#Set flags
GetoptLong(
    "bodysite=s", "bodysite symbol, character, mandatory option",
    "clade=s", "clade name, character, mandatory option",
    "output=s", "output path, character, mandatory option"
    )
#Read in data
input.file <- sprintf("%s/%s/%s_%s_gene_presence_absence_noref", output, bodysite, bodysite, clade)
tbl <- read.table(input.file, header = T, sep="\t")

#If panphlan output, then reformat table and remove gene columns with no variance
rownames(tbl) <- tbl[,1]
tbl <- tbl[, -1]
tbl <- t(tbl)
tbl <- as.data.frame(tbl)
tbl$Sample <- as.factor(rownames(tbl))
tbl$Individual <- as.factor(substr(tbl$Sample, 1, 4))
alldevs <- apply(tbl,2,sd)
alldevs['Individual'] <- 1
alldevs['Sample'] <- 1
tbl2 <- subset(tbl, T, select=(alldevs > 1e-6))
tbl2$Individual <- as.factor(tbl2$Individual)
tbl2$Sample <- as.factor(tbl2$Sample)
rownames(tbl2) <- NULL
tbl2$Sample <- NULL
write.table(tbl2, file=sprintf("%s/%s/%s_%s_pangenePA_fv.txt", output, bodysite, bodysite, clade), sep="\t", quote=F, row.names=F, eol="\n")
write.arff(tbl2, file=sprintf("%s/%s/%s_%s_pangenePA_fv.arff", output, bodysite, bodysite, clade))
