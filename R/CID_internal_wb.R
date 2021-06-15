#!/usr/bin/env Rscript
# Calculate all internal CID fragments
library(stringr)
library(seqinr)

dir <- print(getwd())
setwd(dir) # Automatically set work directory

rnase.digested.fragments <- read.csv(as.character(read.table("FileDirect.txt")$V1[2]),header = TRUE) # FileDirect.txt file path
internal.sequences <- substring(rnase.digested.fragments$Digested_Sequence[which(str_count(rnase.digested.fragments$Digested_Sequence)>3)],2,str_count(rnase.digested.fragments$Digested_Sequence[which(str_count(rnase.digested.fragments$Digested_Sequence)>3)])-1)
cid.frag.array <- list()
for (i in 1:length(internal.sequences)){
  cid.frag <- list()
  for (j in 1:str_count(internal.sequences[i])){
    cid.frag[[j]] <- substring(as.character(internal.sequences[i]),j,c(j+1:str_count(internal.sequences[i])))
  }
  cid.frag.array[[i]] <- cid.frag
}
cid.frag.array <- sapply(sapply(cid.frag.array,unlist),unique)
internal.frag <- unlist(cid.frag.array)[which(str_count(unlist(cid.frag.array))>1)]

# Calculate internal fragments mass
mass_A <- list()
mass_G <- list()
mass_C <- list()
mass_U <- list()
wb_ions_mass <- list()
for (i in 1:length(internal.frag)){
  mass_A[[i]] <- str_count(internal.frag[[i]],"A")*329.0525
  mass_G[[i]] <- str_count(internal.frag[[i]],"G")*345.0474
  mass_C[[i]] <- str_count(internal.frag[[i]],"C")*305.0413
  mass_U[[i]] <- str_count(internal.frag[[i]],"T")*306.0253
  wb_ions_mass[[i]] <- mass_A[[i]] + mass_C[[i]] + mass_G[[i]] + mass_U[[i]] + 78.9585 - 62.96359 + 2.01565
}
wb.ions.frag.mass <- cbind(internal.frag,unlist(wb_ions_mass))
colnames(wb.ions.frag.mass) <- c("w-b_fragments","Mass_Da")

DigestionType <- as.character(read.table("FileDirect.txt")$V1[2])
tag <- unlist(str_split(DigestionType, "out_"))[2]
OutputNameOne <- paste("std_out_wb_ions", tag, sep = "_")
write.csv(wa.ions.frag.mass, OutputNameOne) # Output spreadsheet
q()