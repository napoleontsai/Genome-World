
# Calculate all possible CID fragments
library(stringr)
library(seqinr)
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Automatically set working directory
print(getwd())

rnase.digested.fragments <- read.csv(as.character(read.table("my path/FileDirect.txt")$V1[2]),header = TRUE) # FileDirect.txt file path
cid.frag.array <- list()
for (i in 1:nrow(rnase.digested.fragments)){
  cid.frag <- list()
  for (j in 1:str_length(rnase.digested.fragments$Digested_Sequence[i])){
    cid.frag[[j]] <- substring(as.character(rnase.digested.fragments$Digested_Sequence[i]),j,c(j:str_length(rnase.digested.fragments$Digested_Sequence[i])))
  }
  cid.frag.array[[i]] <- cid.frag
}
cid.frag.array <- sapply(cid.frag.array,unlist)

# Select b-x ions
b.ions <- list()
for (i in 1:length(cid.frag.array)){
  b.ions[[i]] <- cid.frag.array[[i]][1:max(str_count(cid.frag.array[[i]]))-1]
}
x.ions <- list()
for (i in 1:nrow(rnase.digested.fragments)){
  x.ions[[i]] <- substring(rnase.digested.fragments$Digested_Sequence[i],c(2:str_count(rnase.digested.fragments$Digested_Sequence[i])))
}
b.ions <- unlist(b.ions[which(str_count(rnase.digested.fragments$Digested_Sequence) > 1)])
x.ions <- unlist(x.ions[which(str_count(rnase.digested.fragments$Digested_Sequence) > 1)])

# Calculate b-x ion fragments mass
mass_A <- list()
mass_G <- list()
mass_C <- list()
mass_U <- list()
b_ions_mass <- list()
for (i in 1:length(b.ions)){
  mass_A[[i]] <- str_count(b.ions[[i]],"A")*329.2
  mass_G[[i]] <- str_count(b.ions[[i]],"G")*345.2
  mass_C[[i]] <- str_count(b.ions[[i]],"C")*305.2
  mass_U[[i]] <- str_count(b.ions[[i]],"T")*306.2
  b_ions_mass[[i]] <- mass_A[[i]] + mass_C[[i]] + mass_G[[i]] + mass_U[[i]] + 1 - 63
}
b.ions.frag.mass <- cbind(b.ions,unlist(b_ions_mass))

mass_A <- list()
mass_G <- list()
mass_C <- list()
mass_U <- list()
x_ions_mass <- list()
for (i in 1:length(x.ions)){
  mass_A[[i]] <- str_count(x.ions[[i]],"A")*329.2
  mass_G[[i]] <- str_count(x.ions[[i]],"G")*345.2
  mass_C[[i]] <- str_count(x.ions[[i]],"C")*305.2
  mass_U[[i]] <- str_count(x.ions[[i]],"T")*306.2
  x_ions_mass[[i]] <- mass_A[[i]] + mass_C[[i]] + mass_G[[i]] + mass_U[[i]] + 17 + 63
}
x.ions.frag.mass <- cbind(x.ions,unlist(x_ions_mass))
colnames(b.ions.frag.mass) <- c("b_fragments","Mass_Da")
colnames(x.ions.frag.mass) <- c("x_fragments","Mass_Da")
write.csv(b.ions.frag.mass, 'std_out_b_ions.csv') # Output spreadsheet
write.csv(x.ions.frag.mass, 'std_out_x_ions.csv') # Output spreadsheet
q()
