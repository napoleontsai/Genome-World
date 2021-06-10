
# Calculate post enzyme-digested oligo sequences
library(stringr)
library(seqinr)
library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Automatically set work directory
print(getwd())

FastaFile <- read.fasta(as.character(read.table("FileDirect.txt")$V1[1])) # FindDirect.txt file path
frag.one <- unlist(strsplit(str_flatten(FastaFile[[1]]),"g"))
digest.frag.one <- str_remove(paste(frag.one,"g")," ")
if (FastaFile[[1]][length(FastaFile[[1]])] == "g"){
  digest.frag.one[length(digest.frag.one)] <- digest.frag.one[length(digest.frag.one)]
} else{
  digest.frag.one[length(digest.frag.one)] <- str_remove(digest.frag.one[length(digest.frag.one)],"g")
}
if(digest.frag.one[length(digest.frag.one)] == ""){
  digest.frag.one <- digest.frag.one[1:length(digest.frag.one)-1]
} else {
  digest.frag.one <- digest.frag.one[1:length(digest.frag.one)]
}
digest.frag.one <- as.data.frame(digest.frag.one)
digest.frag.one <- na.omit(apply(digest.frag.one,2,str_to_upper))

# Calculate mass of the digested fragments
temp_frag_split <- strsplit(digest.frag.one[1:length(digest.frag.one)],"")
mass_A <- str_count(temp_frag_split,"A")*329.2
mass_G <- str_count(temp_frag_split,"G")*345.2
mass_C <- str_count(temp_frag_split,"C")*305.2
mass_U <- str_count(temp_frag_split,"T")*306.2
oligo_mass <- list()
for(i in 1:length(temp_frag_split)){
  oligo_mass[i] <- mass_A[i] + mass_C[i] + mass_G[i] + mass_U[i] + 18
}
digest.frag.one <- cbind(digest.frag.one, unlist(oligo_mass))
colnames(digest.frag.one) <- c("Digested_Sequence","Oligo_Mass_Da")
write.csv(digest.frag.one,'std_out_RNaseT1.csv') # Output spreadsheet
q()
