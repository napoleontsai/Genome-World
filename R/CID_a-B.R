
# Load CID fragment sequences and mass table
library(stringr)
library(seqinr)

dir <- print(getwd())
setwd(dir) # Automatically set work directory


cid.fragments <- read.csv(as.character(read.table("FileDirect.txt")$V1[3]),header = TRUE)[,2:3] # FileDirect.txt file path

# Calculate base loss mass
temp.cid.frag <- cid.fragments[which(str_count(cid.fragments[,1])>1),]
temp.seq <- list()
for (i in 1:nrow(temp.cid.frag)){
  temp.seq[i] <- substring(as.character(temp.cid.frag[,1][i]),2,str_count(as.character(temp.cid.frag[,1][i])))
}
temp.seq <- unlist(temp.seq)
loss_A <- list()
loss_G <- list()
loss_C <- list()
loss_U <- list()
for (i in 1:length(temp.seq)){
  loss_A[i] <- (str_count(temp.seq[i],"A")/str_count(temp.seq[i],"A"))*(-135.13)
}
for (i in 1:length(temp.seq)){
  loss_G[i] <- (str_count(temp.seq[i],"G")/str_count(temp.seq[i],"G"))*(-151.13)
}
for (i in 1:length(temp.seq)){
  loss_C[i] <- (str_count(temp.seq[i],"C")/str_count(temp.seq[i],"C"))*(-111.10)
}
for (i in 1:length(temp.seq)){
  loss_U[i] <- (str_count(temp.seq[i],"T")/str_count(temp.seq[i],"T"))*(-112.09)
}
loss_mass_A <- list()
for (i in 1:nrow(temp.cid.frag)){
  loss_mass_A[i] <- temp.cid.frag$Mass_Da[i] + unlist(loss_A[i])
}
loss_mass_G <- list()
for (i in 1:nrow(temp.cid.frag)){
  loss_mass_G[i] <- temp.cid.frag$Mass_Da[i] + unlist(loss_G[i])
}
loss_mass_C <- list()
for (i in 1:nrow(temp.cid.frag)){
  loss_mass_C[i] <- temp.cid.frag$Mass_Da[i] + unlist(loss_C[i])
}
loss_mass_U <- list()
for (i in 1:nrow(temp.cid.frag)){
  loss_mass_U[i] <- temp.cid.frag$Mass_Da[i] + unlist(loss_U[i])
}
loss.base <- cbind(temp.cid.frag,unlist(loss_mass_A),unlist(loss_mass_G),unlist(loss_mass_C),unlist(loss_mass_U))
colnames(loss.base) <- c("Sequence","Mass_Da","Loss_A","Loss_G","Loss_C","Loss_U")
write.csv(loss.base, 'std_out_loss_base.csv') # Output spreadsheet
q()
