
# Calculate all internal fragments
library(stringr)
library(seqinr)
rnase.digested.fragments <- read.csv(as.character(read.table("my path/FileDirect.txt")$V1[2]),header = TRUE) # FileDirect.txt file path
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
  mass_A[[i]] <- str_count(internal.frag[[i]],"A")*329.2
  mass_G[[i]] <- str_count(internal.frag[[i]],"G")*345.2
  mass_C[[i]] <- str_count(internal.frag[[i]],"C")*305.2
  mass_U[[i]] <- str_count(internal.frag[[i]],"T")*306.2
  wb_ions_mass[[i]] <- mass_A[[i]] + mass_C[[i]] + mass_G[[i]] + mass_U[[i]] + 79 - 63 + 2
}
wb.ions.frag.mass <- cbind(internal.frag,unlist(wb_ions_mass))
colnames(wb.ions.frag.mass) <- c("w-b_fragments","Mass_Da")
write.csv(wb.ions.frag.mass, 'output folder path/std_out_wb_ions.csv') # Output folder path
q()