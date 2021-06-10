# Calculate post enzyme-digested oligo sequences
library(stringr)
library(seqinr)
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Automatically set work directory
print(getwd())

FastaFile <- read.fasta(as.character(read.table("FileDirect.txt")$V1[1])) # FileDirect.txt file path
frag.one <- unlist(strsplit(str_flatten(FastaFile[[1]]),"c"))
digest.frag.one <- str_remove(paste(frag.one,"c")," ")
t_frag <- strsplit(digest.frag.one[grep("t",digest.frag.one)],"t")
t_frag_laststr <- list()
for(i in 1:length(t_frag)){
  t_frag_laststr[[i]] <- t_frag[[i]][length(t_frag[[i]])]
}
t_frag_laststr_rm <- list()
for(i in 1:length(t_frag)){
  t_frag_laststr_rm[[i]] <- t_frag[[i]][1:(length(t_frag[[i]])-1)]
}
t_paste <- list()
for(i in 1:length(t_frag_laststr_rm)){
  t_paste[[i]] <- str_remove(paste(t_frag_laststr_rm[[i]],"t")," ")
}
digest.frag.two <- list()
for(i in 1:length(t_paste)){
  digest.frag.two[[i]] <- unlist(c(t_paste[[i]],t_frag_laststr[[i]]))
}
digest.frag.one[grep("t",digest.frag.one)] = c(digest.frag.two)
digest.frag.two.df <- unlist(digest.frag.one)
if (FastaFile[[1]][length(FastaFile[[1]])] == "c"){
  digest.frag.two.df[length(digest.frag.two.df)] <- digest.frag.two.df[length(digest.frag.two.df)]
} else{
  digest.frag.two.df[length(digest.frag.two.df)] <- str_remove(digest.frag.two.df[length(digest.frag.two.df)],"c")
}
if(digest.frag.two.df[length(digest.frag.two.df)] == ""){
  digest.frag.two.df <- digest.frag.two.df[1:length(digest.frag.two.df)-1]
} else {
  digest.frag.two.df <- digest.frag.two.df[1:length(digest.frag.two.df)]
}
digest.frag.two.df <- as.data.frame(digest.frag.two.df)
RNaseA_digested_frag <- na.omit(apply(digest.frag.two.df,2,str_to_upper))

# Calculate mass of the digested fragments
temp_frag_split <- strsplit(RNaseA_digested_frag[1:length(RNaseA_digested_frag)],"")
mass_A <- str_count(temp_frag_split,"A")*329.2
mass_G <- str_count(temp_frag_split,"G")*345.2
mass_C <- str_count(temp_frag_split,"C")*305.2
mass_U <- str_count(temp_frag_split,"T")*306.2
oligo_mass <- list()
for(i in 1:length(temp_frag_split)){
  oligo_mass[i] <- mass_A[i] + mass_C[i] + mass_G[i] + mass_U[i] + 18
}
RNaseA_digested_frag <- cbind(RNaseA_digested_frag, unlist(oligo_mass))
colnames(RNaseA_digested_frag) <- c("Digested_Sequence","Oligo_Mass_Da")
write.csv(RNaseA_digested_frag,'std_out_RNaseA.csv') # Output spreadsheet
q()
