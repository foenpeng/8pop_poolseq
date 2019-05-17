library(data.table)
library(poolfstat)
setwd("D:/Foen Peng/PoolSeq_fish_6_pop")

load("./processed_data/all_pool_pairs_new_fst_pooldata_Indel")
load("./processed_data/all_pool_pairs_new_fst_Indel")


reads_only<-as.data.table(cbind(pooldata@snp.info,pooldata@readcoverage,pooldata@refallele.readcount))
reads_Fst<-cbind(reads_only,Fst_matrix$PairwiseSnpFST)
# change variable names
#Fst<- fread("./processed_data/Fst_8pop.txt", header = T)
setnames(reads_Fst, c("V1","V2","V3","V4","V5",
                             "V6","V7","V8","V9","V10","V11","V12","V13",
                             "V14","V15","V16","V17","V18","V19","V20","V21"),
                     c("LG", "Pos", "ref","Base_A", "Base_a",
                             "boot_Nreads","echo_Nreads","fred_Nreads","gos_Nreads","law_Nreads","pach_Nreads","rob_Nreads","say_Nreads",
                             "boot_N_A","echo_N_A","fred_N_A","gos_N_A","law_N_A","pach_N_A","rob_N_A","say_N_A"))

fwrite(reads_Fst,"../OneDrive - University of Connecticut/Poolseq_8pops/processed_data/Fst_8pops_WithIndel.csv")
write.table(as.data.frame(Fst_matrix$PairwiseFSTmatrix),sep=",","../OneDrive - University of Connecticut/Poolseq_8pops/processed_data/Pop_Fst_8pops_WithIndel.csv",row.names = T, col.names = T)
