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


# read only relevant data
# pop1 should be the one with first letter in front of pop2, e.g. "g" > "r", for the sake of mathing pairwise comparison in gos_vs_rob
pop1<-"boot"
pop2<-"gos"
outgroup="say"


col_keep<-append(c("LGn", "Pos", "Base_A", "Base_a"), 
                 c(paste0(pop1,"_N_A"),paste0(pop2,"_N_A"),paste0(outgroup,"_N_A"),
                 paste0(pop1,"_Nreads"),paste0(pop2,"_Nreads"),paste0(outgroup,"_N_A"),
                 paste0(pop1,"_vs_",pop2),paste0(pop1,"_vs_",outgroup),paste0(pop2,"_vs_",outgroup)))
pop1_pop2_Fst<-Fst[,..col_keep]
setnames(pop1_pop2_Fst,-(1:4),c("pop1_N_A","pop2_N_A","outgroup_N_A",
                                "pop1_Nreads","pop2_Nreads","outgroup_Nreads",
                                "pop1_pop2_Fst","pop1_outgroup_Fst","pop2_outgroup_Fst"))

#################################################################
#  3.  Calculate PBS and other statistics
#################################################################

{
  # calculate allele frequency
  pop1_pop2_Fst[,pop1_FreqA := pop1_N_A/pop1_Nreads]
  pop1_pop2_Fst[,pop2_FreqA := pop2_N_A/pop2_Nreads]
  pop1_pop2_Fst[,AFD_pop1_pop2 := abs(pop1_FreqA-pop2_FreqA)]
  
  # calculate populatin branch statistics
  pop1_pop2_Fst[,pbs.1 := ((-log(1- pop1_pop2_Fst)) + (- log(1-pop1_outgroup_Fst)) - (- log(1-pop2_outgroup_Fst )))/2]
  pop1_pop2_Fst[,pbs.2 := ((-log(1- pop1_pop2_Fst)) + (- log(1-pop2_outgroup_Fst)) - (- log(1-pop1_outgroup_Fst )))/2]
  pop1_pop2_Fst[,pbs.o := ((-log(1- pop1_outgroup_Fst)) + (- log(1-pop2_outgroup_Fst)) - (- log(1-pop1_pop2_Fst)))/2]
}
#################################################################
#  5. sliding window
#################################################################
{
  windowsize <- 50000
  stepsize<-10000
  
  chr <- pop1_pop2_Fst[, .(chr_length = max(Pos)), by = LGn]
  #chr[,cumulative_chrlengths := cumsum(chr_length)]
  #chr[,cumulative_chrSTART := shift(cumulative_chrlengths,1,0,"lag")]
  chr_slide<-chr[rep(1:.N,(chr_length)%/%stepsize)][,.(window_start=(0:.N)*stepsize, window_end=((0:.N)*stepsize+windowsize),chr_length=chr_length), by = LGn]
  chr_slide<-chr_slide[window_end<chr_length+stepsize]
  pop1_pop2_Fst[,Pos_join:=Pos]
  pop1_pop2_slide<-pop1_pop2_Fst[chr_slide, 
                    on=c("LGn","Pos_join>window_start","Pos_join<window_end"), 
                    allow.cartesian=T][,.("PBS.nSNPs" = .N,
                                          "PBS.1" = mean(pbs.1,na.rm=T), 
                                          "PBS.2" = mean(pbs.2,na.rm=T),
                                          "PBS.o" = mean(pbs.o,na.rm=T),
                                          "Fst_1_2" = mean(pop1_pop2_Fst,na.rm=T), 
                                          "AFD_1_2"= mean(AFD_pop1_pop2,na.rm=T),
                                          "window_end" = Pos_join.1[1]),
                                       by=.(LGn,Pos_join)]
  setnames(pop1_pop2_slide, "Pos_join", "window_start")
  pop1_pop2_slide[,Pos:=(window_end+window_start)/2]
}
#################################################################
#  6.  Plotting results
#################################################################

##### Plot binned genomic map, modified from Dan's plot function from file Graphics function.R

plotgene_bin <- function(dat, pop, focalLG, specifybps = F, minpos, maxpos, Genestart, Geneend, Genename, statstoplot){
    npops <- length(statstoplot)
    par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
    if(specifybps == F){
      minpos <- 0
      maxpos <- dat[LGn==focalLG,max(Pos)]
    }
    options(scipen=5)
    col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
    for(i in 1:length(statstoplot)){
      
      plot(dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos,.(Pos/1000000,eval(parse(text = statstoplot[i])))],pch = 16, cex = 0.6, axes = F, xlab = paste("Chromosome",focalLG), ylab = statstoplot[i],col=col[i])
      abline(v=dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos & eval(parse(text = paste(statstoplot[i],".focal",sep=""))),Pos/1000000],col=rgb(red=0,green=0,blue=0,alpha=100,maxColorValue = 255))
      ticks<-seq(0,max(dat[LGn %in% c(focalLG),Pos])/1e6,1)
      axis(1, at=seq(0,maxpos/1e6,0.5), labels = F)
      axis(1, at=seq(0,maxpos/1e6,1), labels = seq(0,maxpos/1e6,1))
      axis(2)
    }
    mtext(paste("Chromosome",focalLG, "1:",pop[1],"  2:",pop[2]),side=1,line=1, outer = TRUE)
  }
  
# to plot the data
{
  dat_plot <- pop1_pop2_slide
  statstoplot = c("PBS.1", "PBS.2", "PBS.o","Fst_1_2","AFD_1_2")
  
  cutoff_calculation<-function(var, data=dat_plot, threshold=0.99){
    var.cutoff <- quantile(data[,get(var)],threshold,na.rm=T)
    data[,paste0(var,".focal"):= ifelse(get(var)> var.cutoff, T, F)]
  }
  
  lapply(statstoplot,cutoff_calculation)
  
  
  for (chromosome in 2) {
    png(filename=sprintf("./result/chr_map_50k_slide/chr%s.%s_%s.divergence_0.99.png", chromosome,pop1,pop2),width = 3000, height = 2000,res=300)
    plotgene_bin(dat_plot,pop=c(pop1,pop2), focalLG=chromosome, statstoplot=statstoplot)
    dev.off()
  } 
}