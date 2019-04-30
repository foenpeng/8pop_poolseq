# In the preprocessing steps for Fst_8pop.txt, coverage requirement was set to 10-5000 for all pops, and SNP contain indels were dropped

#################################################################
#  1.  Load packages
#################################################################

library(data.table)

setwd("D:/Foen Peng/OneDrive - University of Connecticut/Poolseq_8pops")

#################################################################
#  2.  Load data and organize as needed
#################################################################

# change variable names
#Fst<- fread("./processed_data/Fst_8pop.txt", header = T)
#setnames(Fst, c("V1","V2","V3","V4","V5",
#                             "V6","V7","V8","V9","V10","V11","V12","V13",
#                             "V14","V15","V16","V17","V18","V19","V20","V21"),
#                     c("LG", "Pos", "ref","Base_A", "Base_a",
#                             "boot_Nreads","echo_Nreads","fred_Nreads","gos_Nreads","law_Nreads","pach_Nreads","rob_Nreads","say_Nreads",
#                             "boot_N_A","echo_N_A","fred_N_A","gos_N_A","law_N_A","pach_N_A","rob_N_A","say_N_A"))
#fwrite(Fst,"Fst_8pops.csv")

Fst<- fread("./processed_data/Fst_8pops.csv", header = T)
Fst_sample<- fread("./processed_data/Fst_8pops.csv", header = T, nrows=10)

# change LG name, 1000 means mitochondial DNA, NA means unmapped
Fst[substr(LG, 1, 2) == "ch", LGn := as.numeric(as.roman(substr(LG, 4, nchar(LG))))]
# replace all the negative and na values with 0
Fst[,(22:49):=lapply(.SD, function(x) replace(x, which(x < 0 | is.na(x)), 0)),.SDcols=22:49]
# replace all the 1 with 0.99
Fst[,(22:49):=lapply(.SD, function(x) replace(x, which(x ==1), 0.99)),.SDcols=22:49]

pop_list<-c("boot","echo","fred","gos","law","pach","rob","say")

pairwise_data<-function(pop1,pop2,outgroup="say",data=Fst){
  pop1_N_A<-paste0(pop1,"_N_A")
  pop2_N_A<-paste0(pop2,"_N_A")
  outgroup_N_A<-paste0(outgroup,"_N_A")
  
  pop1_Nreads<-paste0(pop1,"_Nreads")
  pop2_Nreads<-paste0(pop2,"_Nreads")
  outgroup_Nreads<-paste0(outgroup,"_Nreads")
  
  pop1_FreqA<-paste0(pop1,"_FreqA")
  pop2_FreqA<-paste0(pop2,"_FreqA")
  
  pop1_pop2_Fst<-paste0(pop1,"_vs_",pop2)
  pop1_outgroup_Fst<-paste0(pop1,"_vs_",outgroup)
  pop2_outgroup_Fst<-paste0(pop2,"_vs_",outgroup)
  
  #extract infor from big table
  col_keep<-append(colnames(data)[1:5],c("LGn",pop1_N_A,pop1_Nreads,pop2_N_A,pop2_Nreads,outgroup_N_A,outgroup_Nreads,pop1_pop2_Fst,pop1_outgroup_Fst,pop2_outgroup_Fst))
  return(data[,..col_keep])
}

# calculate allele frequency
pop1_pop2_comp[,pop1_FreqA := get(pop1_N_A)/get(pop1_Nreads)]
pop1_pop2_comp[,pop2_FreqA := pop2_N_A/pop2_Nreads]

# calculate populatin branch statistics
pop1_pop2_comp[,paste0("pbs.",pop1) := ((-log(1- pop1_pop2_Fst)) + (- log(1-pop1_outgroup_Fst)) - (- log(1-pop2_outgroup_Fst )))/2]
pop1_pop2_comp[,paste0("pbs.",pop2) := ((-log(1- pop1_pop2_Fst)) + (- log(1-pop2_outgroup_Fst)) - (- log(1-pop1_outgroup_Fst )))/2]
pop1_pop2_comp[,paste0("pbs.",outgroup) := ((-log(1- pop1_outgroup_Fst)) + (- log(1-pop2_outgroup_Fst)) - (- log(1-pop1_pop2_Fst)))/2]



