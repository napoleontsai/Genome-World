#!/usr/bin/env Rscript
# Use API to access RNA modifications database
library(httr)
library(jsonlite)
library(stringr)
library(seqinr)
dir <- getwd()
setwd(dir)

# Parse modifications data
root_directory <- "www.genesilico.pl/modomics/api/"
modAPI <- paste(root_directory, "modifications", sep = "")
res <- GET(modAPI)
allMod <- fromJSON(rawToChar(res$content))
original_base <- list()
all_short_names <- list()
all_iso_mass <- list()
for (i in 1:length(allMod)){
  original_base[[i]] <- allMod[1:length(allMod)][[i]]$original_base
}
for (i in 1:length(allMod)){
  all_short_names[[i]] <- allMod[1:length(allMod)][[i]]$short_name
}
for (i in 1:length(allMod)){
  all_iso_mass[[i]] <- allMod[1:length(allMod)][[i]]$mass_monoiso
}
mod.df <- as.data.frame(cbind(unlist(original_base), unlist(all_short_names), unlist(all_iso_mass)))
colnames(mod.df) <- c("original_base","short_names","mass_monoiso")

# Perform mass calculations on the selected modification
# RNase digested precursor ions
InputMod <- as.character(read.table("FileDirect.txt")$V1[4])
InputMod_Mass <- as.numeric(as.character(mod.df[which(mod.df$short_names %in% InputMod),]$mass_monoiso)) - 2 + 63
InputOriginal <- as.character(mod.df[which(mod.df$short_names %in% InputMod),]$original_base)
rnase.digested.fragments <- read.csv(as.character(read.table("FileDirect.txt")$V1[2]), header = TRUE)
MatchedSequences <- str_count(as.character(unlist(rnase.digested.fragments$Digested_Sequence)), InputOriginal)

Original_Base_Detection <- function(unmodified_base){
  mass_Adenosine <- 329.0525*str_count(unmodified_base,"A")  
  mass_Guanosine <- 345.0474*str_count(unmodified_base,"G")
  mass_Cytidine <- 305.0413*str_count(unmodified_base,"C")
  mass_Uridine <- 306.0253*str_count(unmodified_base,"U")
  return(sum(mass_Adenosine, mass_Guanosine, mass_Cytidine, mass_Uridine))
}
MassOriginalBase <- Original_Base_Detection(InputOriginal)

potential.mod.bases <- sapply(MatchedSequences, seq_len)

modified_oligo_mass <- function(OligoMass, NumberOfModBases){
  ModOligoMass <- OligoMass - NumberOfModBases*(MassOriginalBase - InputMod_Mass)
  return(ModOligoMass)
}

modified.precursor.mass <- list()
for (i in 1:length(rnase.digested.fragments$Oligo_Mass_Da)){
  modified.precursor.mass[[i]] <- modified_oligo_mass(rnase.digested.fragments$Oligo_Mass_Da[i], potential.mod.bases[[i]])
}
names(modified.precursor.mass) <- rnase.digested.fragments$Digested_Sequence

DigestionType <- as.character(read.table("FileDirect.txt")$V1[2])
tagOne <- unlist(str_remove(unlist(str_split(DigestionType, "out_"))[2], ".csv"))
OutputNameOne <- paste("std_out_precursor", tagOne, "modmass.json", sep = "_")
write_json(modified.precursor.mass, OutputNameOne)

# CID fragmented product ions
cid.fragments <- read.csv(as.character(read.table("FileDirect.txt")$V1[3]), header = TRUE)
MatchedSequences <- str_count(as.character(unlist(cid.fragments[2])), InputOriginal)

potential.mod.bases <- sapply(MatchedSequences, seq_len)

modified.product.mass <- list()
for (i in 1:length(cid.fragments$Mass_Da)){
  modified.product.mass[[i]] <- modified_oligo_mass(cid.fragments$Mass_Da[i], potential.mod.bases[[i]])
}
names(modified.product.mass) <- unlist(cid.fragments[2])

ionType <- as.character(read.table("FileDirect.txt")$V1[3])
tagTwo <- unlist(str_remove(unlist(str_split(ionType, "out_"))[2], ".csv"))
OutputNameTwo <- paste("std_out_product", tagTwo, "modmass.json", sep = "_")
write_json(modified.product.mass, OutputNameTwo)
q()
