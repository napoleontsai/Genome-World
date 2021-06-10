
# Calculate all possible CID fragments
library(stringr)
library(seqinr)
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Automatically set work directory
print(getwd())

rnase.digested.fragments <- read.csv(as.character(read.table("FileDirect.txt")$V1[2]),header = TRUE) # FileDirect.txt file path
cid.frag.array <- list()
for (i in 1:nrow(rnase.digested.fragments)){
  cid.frag <- list()
  for (j in 1:str_length(rnase.digested.fragments$Digested_Sequence[i])){
    cid.frag[[j]] <- substring(as.character(rnase.digested.fragments$Digested_Sequence[i]),j,c(j:str_length(rnase.digested.fragments$Digested_Sequence[i])))
  }
  cid.frag.array[[i]] <- cid.frag
}
cid.frag.array <- sapply(cid.frag.array,unlist)

# Select a-w ions
a.ions <- list()
for (i in 1:length(cid.frag.array)){
  a.ions[[i]] <- cid.frag.array[[i]][1:max(str_count(cid.frag.array[[i]]))-1]
}
w.ions <- list()
for (i in 1:nrow(rnase.digested.fragments)){
  y.ions[[i]] <- substring(rnase.digested.fragments$Digested_Sequence[i],c(2:str_count(rnase.digested.fragments$Digested_Sequence[i])))
}
a.ions <- unlist(a.ions[which(str_count(rnase.digested.fragments$Digested_Sequence) > 1)])
w.ions <- unlist(w.ions[which(str_count(rnase.digested.fragments$Digested_Sequence) > 1)])

# Calculate a-w ion fragments mass
mass_A <- list()
mass_G <- list()
mass_C <- list()
mass_U <- list()
a_ions_mass <- list()
for (i in 1:length(a.ions)){
  mass_A[[i]] <- str_count(a.ions[[i]],"A")*329.2
  mass_G[[i]] <- str_count(a.ions[[i]],"G")*345.2
  mass_C[[i]] <- str_count(a.ions[[i]],"C")*305.2
  mass_U[[i]] <- str_count(a.ions[[i]],"T")*306.2
  a_ions_mass[[i]] <- mass_A[[i]] + mass_C[[i]] + mass_G[[i]] + mass_U[[i]] + 1 - 79
}
a.ions.frag.mass <- cbind(a.ions,unlist(a_ions_mass))

mass_A <- list()
mass_G <- list()
mass_C <- list()
mass_U <- list()
w_ions_mass <- list()
for (i in 1:length(w.ions)){
  mass_A[[i]] <- str_count(w.ions[[i]],"A")*329.2
  mass_G[[i]] <- str_count(w.ions[[i]],"G")*345.2
  mass_C[[i]] <- str_count(w.ions[[i]],"C")*305.2
  mass_U[[i]] <- str_count(w.ions[[i]],"T")*306.2
  w_ions_mass[[i]] <- mass_A[[i]] + mass_C[[i]] + mass_G[[i]] + mass_U[[i]] + 17 + 79
}
w.ions.frag.mass <- cbind(w.ions,unlist(w_ions_mass))
colnames(a.ions.frag.mass) <- c("a_fragments","Mass_Da")
colnames(w.ions.frag.mass) <- c("w_fragments","Mass_Da")
write.csv(a.ions.frag.mass, 'std_out_a_ions.csv') # Output spreadsheet
write.csv(w.ions.frag.mass, 'std_out_w_ions.csv') # Output spreadsheet
q()