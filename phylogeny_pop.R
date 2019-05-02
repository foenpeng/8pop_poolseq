library(data.table)
load("D:/Foen Peng/PoolSeq_fish_6_pop/processed_data/all_pool_pairs_new_fst_pooldata")
pop_Fst<-as.data.table(Fst_matrix[[1]])
fwrite(pop_Fst,"D:/Foen Peng/OneDrive - University of Connecticut/Poolseq_8pops/result/pop_pairwise_Fst.csv")


library(ape)
dis_matrx<-Fst_matrix[[1]]
tree<-nj(dis_matrx)
root(tree,outgroup = "say",resolve.root = TRUE)
png("D:/Foen Peng/OneDrive - University of Connecticut/Poolseq_8pops/result/NJ_tr_pop.png",width = 2000, height = 1600, res = 300)
plot(tree)
edgelabels(round(tree$edge.length,2), bg="white", col="black", font=1,cex=0.8)
dev.off()

