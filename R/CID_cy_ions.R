
# Calculate all possible CID fragments
library(stringr)
library(seqinr)

dir <- print(getwd())
setwd(dir) # Automatically set work directory

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

# Select c-y ions
c.ions <- list()
for (i in 1:length(cid.frag.array)){
  c.ions[[i]] <- cid.frag.array[[i]][1:max(str_count(cid.frag.array[[i]]))-1]
}
y.ions <- list()
for (i in 1:nrow(rnase.digested.fragments)){
  y.ions[[i]] <- substring(rnase.digested.fragments$Digested_Sequence[i],c(2:str_count(rnase.digested.fragments$Digested_Sequence[i])))
}
c.ions <- unlist(c.ions[which(str_count(rnase.digested.fragments$Digested_Sequence) > 1)])
y.ions <- unlist(y.ions[which(str_count(rnase.digested.fragments$Digested_Sequence) > 1)])

# Calculate c-y ion fragments mass
mass_A <- list()
mass_G <- list()
mass_C <- list()
mass_U <- list()
c_ions_mass <- list()
for (i in 1:length(c.ions)){
  mass_A[[i]] <- str_count(c.ions[[i]],"A")*329.2
  mass_G[[i]] <- str_count(c.ions[[i]],"G")*345.2
  mass_C[[i]] <- str_count(c.ions[[i]],"C")*305.2
  mass_U[[i]] <- str_count(c.ions[[i]],"T")*306.2
  c_ions_mass[[i]] <- mass_A[[i]] + mass_C[[i]] + mass_G[[i]] + mass_U[[i]] + 1
}
c.ions.frag.mass <- cbind(c.ions,unlist(c_ions_mass))

mass_A <- list()
mass_G <- list()
mass_C <- list()
mass_U <- list()
y_ions_mass <- list()
for (i in 1:length(c.ions)){
  mass_A[[i]] <- str_count(y.ions[[i]],"A")*329.2
  mass_G[[i]] <- str_count(y.ions[[i]],"G")*345.2
  mass_C[[i]] <- str_count(y.ions[[i]],"C")*305.2
  mass_U[[i]] <- str_count(y.ions[[i]],"T")*306.2
  y_ions_mass[[i]] <- mass_A[[i]] + mass_C[[i]] + mass_G[[i]] + mass_U[[i]] + 17
}
y.ions.frag.mass <- cbind(y.ions,unlist(y_ions_mass))
colnames(c.ions.frag.mass) <- c("c_fragments","Mass_Da")
colnames(y.ions.frag.mass) <- c("y_fragments","Mass_Da")
write.csv(c.ions.frag.mass, 'std_out_c_ions.csv') # Output spreadsheet
write.csv(y.ions.frag.mass, 'std_out_y_ions.csv') # Output spreadsheet
q()
