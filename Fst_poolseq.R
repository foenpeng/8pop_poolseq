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

Fst<- fread("./processed_data/Fst_8pops_WithIndel.csv", header = T)
Fst_sample<- fread("./processed_data/Fst_8pops_WithIndel.csv", header = T, nrows=10)
# read population Fst as data frame, becasue data table does not support row names
Fst_genome<-read.csv("./processed_data/Pop_Fst_8pops_WithIndel.csv", header = T)
#row.names(Fst_genome)<-colnames(Fst_genome)
Fst_genome<-as.matrix(Fst_genome)

# change LG name, 1000 means mitochondial DNA, NA means unmapped
Fst[substr(LG, 1, 2) == "ch", LGn := as.numeric(as.roman(substr(LG, 4, nchar(LG))))]
Fst[is.na(LGn),LGn:=100]
# replace all the negative and na values with 0
Fst[,(22:49):=lapply(.SD, function(x) replace(x, which(x < 0 | is.na(x)), 0)),.SDcols=22:49]
# replace all the 1 with 0.99
Fst[,(22:49):=lapply(.SD, function(x) replace(x, which(x ==1), 0.99)),.SDcols=22:49]

pop_list<-c("boot","echo","fred","gos","law","pach","rob","say")

# import gene information
## function to convert linkage group names from string to numeric number
{
  convert_LG_name<-function(dataset, LG_column){
    dataset[substr(get(LG_column), 1, 2) == "MT", LGn := 0]
    dataset[substr(get(LG_column), 1, 2) == "gr", LGn := as.numeric(as.roman(substr(get(LG_column), 6, nchar(get(LG_column)))))]
    dataset[substr(get(LG_column), 1, 2) == "sc", LGn := as.numeric(substr(get(LG_column), 10, nchar(get(LG_column))))]
  }
  gene_loc<-unique(fread("../Analysis_Expression/Data files RAW/GeneLocations.csv",sep=","))
  gene_name<-fread("../Analysis_Expression/Data files RAW/GeneID_to_GeneName.csv")
  gene_LG<-fread("../Analysis_Expression/Data files RAW/GeneID to LinkageGroup.csv",header = TRUE)
  setnames(gene_name,c("Ensembl Gene ID","Associated Gene Name"),c("GeneID","gene.name"))
  convert_LG_name(gene_LG, "LinkageGroup")
  # LG file has a space at the begining of GeneID.
  gene_LG[,GeneID:=substring(V1,2)]
  
  setkey(gene_loc,GeneID)
  setkey(gene_LG, GeneID)
  setkey(gene_name,GeneID)
  
  # merge three tables to gene_info
  gene_info<-gene_LG[gene_loc, on=c("GeneID")]
  gene_info<-gene_info[gene_name]
  gene_info<-gene_info[,!c("V1", "LinkageGroup")]
  setkey(gene_info,LGn,Start,Stop)
  
  # seperate the full length from the exons
  gene_info[,full:=FALSE]
  gene_info[,c("gene_start","gene_end"):=list(min(Start),max(Stop)),by=GeneID]
  gene_info[Start==gene_start&Stop==gene_end,full:=TRUE]
  gene_info[,c("gene_start","gene_end"):=NULL]

}


#################################################################
#  2.  Pairwise analysis
#################################################################

# read only relevant data
# pop1 should be the one with first letter in front of pop2, e.g. "g" > "r", for the sake of mathing pairwise comparison in gos_vs_rob
pop1<-"gos"
pop2<-"rob"
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
#  1.  Calculate PBS and other statistics


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
#  2. sliding window

{
  windowsize <- 10000
  stepsize<-5000
  
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
#  3.  Plotting results


##### Plot binned genomic map, modified from Dan's plot function from file Graphics function.R

plotgene_bin <- function(dat, pop, focalLG, specifybps = F, minpos, maxpos, Genestart, Geneend, Genename, statstoplot, y_range){
    npops <- length(statstoplot)
    par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
    if(specifybps == F){
      minpos <- 0
      maxpos <- dat[LGn==focalLG,max(Pos)]
    }
    options(scipen=5)
    col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
    for(i in 1:length(statstoplot)){
      y_lim = unlist(y_range[,get(statstoplot[i])])
      plot(dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos,.(Pos/1000000,eval(parse(text = statstoplot[i])))],pch = 16, cex = 0.6, axes = F, xlab = paste("Chromosome",focalLG), ylab = statstoplot[i],col=col[i], ylim=y_lim)
      abline(v=dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos & eval(parse(text = paste(statstoplot[i],".focal",sep=""))),Pos/1000000],col=rgb(red=0,green=0,blue=0,alpha=50,maxColorValue = 255))
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
  y_range = dat_plot[,..statstoplot][,lapply(.SD,function(x) list(min(x,na.rm = T),max(x,na.rm = T)))]
  
  for (chromosome in 1:21) {
    png(filename=sprintf("./result/chr_map_10k_slide/chr%s.%s_%s.divergence_0.99.png", chromosome,pop1,pop2),width = 3000, height = 2000,res=300)
    plotgene_bin(dat_plot,pop=c(pop1,pop2), focalLG=chromosome, statstoplot=statstoplot,y_range=y_range)
    dev.off()
  } 
}




#################################################################
#  3.  All comparison analysis
#################################################################

#fib_list<-c("boot","echo","fred","gos","law","pach","rob")
#pop_h_list<-c("boot","fred","law","rob")
#pop_l_list<-c("echo","gos","pach")

fib_list<-c("boot","echo","gos","law","rob")
pop_h_list<-c("law","rob")
pop_l_list<-c("boot","echo","gos")


generate_pop_list<-function(list){
  all_comb<-t(combn(list,2))
  return(paste0(all_comb[,1],"_vs_",all_comb[,2]))
}

hh_list<-generate_pop_list(pop_h_list)
ll_list<-generate_pop_list(pop_l_list)
all_list<-generate_pop_list(fib_list)
hl_list<-all_list[!(all_list %in% hh_list | all_list %in% ll_list)]

setcolorder(Fst,c(hh_list,hl_list,ll_list))

# create sliding window frames
windowsize <- 1000
stepsize<-500
setkeyv(Fst,c("LGn","Pos"))
chr <- Fst[, .(chr_snp=.N,chr_length = max(Pos)), by = LGn]
chr[,Cumulative_Chrlengths:=cumsum(chr_length)]
chr[,Cumulative_ChrSTART := c(0, Cumulative_Chrlengths[-length(Cumulative_Chrlengths)])]
chr_slide<-chr[rep(1:.N,(chr_snp)%/%stepsize)][,.(window_start=(0:.N)*stepsize, window_end=((0:.N)*stepsize+windowsize),chr_snp=chr_snp,chr_length, cum_chr_start=Cumulative_ChrSTART), by = LGn]
chr_slide<-chr_slide[window_end<chr_snp+stepsize]

#dXY and Kirkpatrick's method
dxy<-data.table(LGn=Fst[,LGn],Pos=Fst[,Pos])

dxy_cal<-function(comparison, dat){
  pops <- unlist(strsplit(comparison, split="_vs_"))
  pop1_freq<-paste0(pops[1],"_freq")
  pop2_freq<-paste0(pops[2],"_freq")
  dxy <- dat[,get(pop1_freq)*(1-get(pop2_freq))+get(pop2_freq)*(1-get(pop1_freq))]
  return(dxy)
}

for(i in pop_list){
  dxy[,paste0(i,"_freq"):= Fst[,get(paste0(i,"_N_A"))/get(paste0(i,"_Nreads"))]]
}
for(i in all_list){
  dxy[,paste0(i,"_dxy"):= dxy[,dxy_cal(i,dxy)]]
}


dxy[,Pos_join:= rowid(LGn)]
Ra_slide_temp<-dxy[chr_slide, 
                      on=c("LGn","Pos_join>window_start","Pos_join<window_end"), 
                      allow.cartesian=T]
Ra_slide_part1<-Ra_slide_temp[,.("PBS.nSNPs" = .N,
                        "Pos_start" = min(Pos),
                        "Pos_end" = max(Pos),
                        "Pos" = (max(Pos)+min(Pos))/2, 
                        "Pos_cum" = (max(Pos)+min(Pos))/2 + max(cum_chr_start)),
                     by=.(LGn,Pos_join)]
Ra_slide_part2<-Ra_slide_temp[,lapply(.SD,mean),.SDcol=paste0(all_list,"_dxy"), by=.(LGn,Pos_join)]
Ra_slide<-Ra_slide_part1[Ra_slide_part2,on=c("LGn","Pos_join")]
rm(Ra_slide_temp,Ra_slide_part1,Ra_slide_part2)

sum_dist<-function(y) sum(unlist(lapply(y,function(x) -log(1-x))),na.rm=T)


pbs_cal<-function(data,dist_fun,hl=hl_list,hh=hh_list,ll=ll_list){
  snp_num<-nrow(data)
  data[,D_hl:=sum_dist(.SD),.SD=hl_list,by=seq_len(snp_num)]
  data[,D_hh:=sum_dist(.SD),.SD=hh_list,by=seq_len(snp_num)]
  data[,D_ll:=sum_dist(.SD),.SD=ll_list,by=seq_len(snp_num)]
  data[,pbs:=(D_hl-(D_hh*3/5+D_ll*6/2))/12,by=seq_len(snp_num)]
}



# Tau_Foen: calculate the branch lenth between h and L pops: (D_hl - D_hh*N_ll/(N_hh-1) -D_ll*N_hh/(N_ll-1))/N_hl
sum_dist<-function(y) sum(unlist(lapply(y,function(x) -log(1-x))),na.rm=T)


pbs_cal<-function(data,dist_fun,hl=hl_list,hh=hh_list,ll=ll_list){
  snp_num<-nrow(data)
  data[,D_hl:=sum_dist(.SD),.SD=hl_list,by=seq_len(snp_num)]
  data[,D_hh:=sum_dist(.SD),.SD=hh_list,by=seq_len(snp_num)]
  data[,D_ll:=sum_dist(.SD),.SD=ll_list,by=seq_len(snp_num)]
  data[,pbs:=(D_hl-(D_hh*3/5+D_ll*6/2))/12,by=seq_len(snp_num)]
}


# read popgen files, which is very large
if(file.exists("./result/Tau_Foen_RobLaw_BootEchoGos.csv")){
  Tau_foen_data<-fread("./result/Tau_Foen_RobLaw_BootEchoGos.csv",header = T)
  Fst<-Fst[Tau_foen_data,on=c("LGn","Pos")]
}else{
pbs_cal(data=Fst,dist_fun = sum_dist)}

# calculate average Distance
Fst[,D_hl:=D_hl/length(hl_list)]
Fst[,D_hh:=D_hh/length(hh_list)]
Fst[,D_ll:=D_ll/length(ll_list)]


#fwrite(Fst[,.(LGn,Pos,pbs)],"./result/Tau_Foen_RobLaw_BootEchoGos.csv")

# Tau_Foen: correct for phylogeny
phy_correction<-function(comparison, value, matrix){
  pops <- unlist(strsplit(comparison, split="_vs_"))
  correct <- value-matrix[pops[1],pops[2]]
  return(correct)
}


Fst_phy_correct<-data.table(LGn=Fst[,LGn],Pos=Fst[,Pos])
for(i in all_list){
  Fst_phy_correct[,(i):=phy_correction(i,Fst[,get(i)],Fst_genome)]
}
pbs_cal(data=Fst_phy_correct,dist_fun = sum_dist)



###### creating sliding window by SNP numbers
slide_data<-Fst
slide_data[,Pos_join:= rowid(LGn)]
Tau_slide<-slide_data[chr_slide, 
               on=c("LGn","Pos_join>window_start","Pos_join<window_end"), 
               allow.cartesian=T][,.("PBS.nSNPs" = .N,
                                     "Pos_start" = min(Pos),
                                     "Pos_end" = max(Pos),
                                     "Pos" = (max(Pos)+min(Pos))/2, 
                                     "Pos_cum" = (max(Pos)+min(Pos))/2 + max(cum_chr_start),
                                     "pbs_slide" = mean(pbs,na.rm=T)),
                                  by=.(LGn,Pos_join)]

# function to plot region view of every SNP
{
  plot_region <- function(dat, focalLG, minpos, maxpos, gene.name, gene.id, exon, statstoplot, y_range){
    npops <- length(statstoplot)
    par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
    
    options(scipen=5)
    col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
    for(i in 1:length(statstoplot)){
      y_lim = unlist(y_range[,get(statstoplot[i])])
      plot(dat[,.(Pos,eval(parse(text = statstoplot[i])))], pch = 16, cex = 0.5, axes = F, ylim=y_lim*1.2, ylab = statstoplot[i],col=col[i])
      rect(exon[,Start], y_lim[1],  exon[,Stop],y_lim[2], border = rgb(1,0,0,0) , lwd = 1, col = rgb(1,0,0,0.1))
      text(exon[,(Start+Stop)/2],y_lim[2]/1.2,labels = paste(substr(gene.id,14,nchar(gene.id)),gene.name,sep=" "),cex=0.7, pos=1,srt=90)
      axis(1, at=seq(minpos,maxpos,0.025E6), labels = F)
      axis(1, at=seq(minpos,maxpos,0.05E6), labels = seq(minpos/1e6,maxpos/1e6,0.05))
      axis(2)
    }
    mtext(paste("LG",focalLG,":", minpos,"-",maxpos),side=1,line=1, outer = TRUE)
  }
  
  region_to_plot<-data.table(
    LGn = 7,
    start = 13.1e6,
    end = 13.7e6
  )
  statstoplot = c("pbs","gos_vs_rob","D_hl","D_ll","D_hh")
  dat_plot = Fst
  y_range = dat_plot[,..statstoplot][,lapply(.SD,function(x) list(min(x),max(x)))]
  for (i in 1:region_to_plot[,.N]) {
    LGn_plot<-region_to_plot[i,LGn]
    region_start_plot<-region_to_plot[i,start]
    region_end_plot<-region_to_plot[i,end]
    
    #png(filename=sprintf("./result/Region_zoomin/LG%s_%s-%s_RobLaw_BootEchoGos.png",LGn_plot,region_start_plot/1e6,region_end_plot/1e6),width = 3500, height = 2000,res=300)
    
    gene_info_plot<-gene_info[Start>region_start_plot & Stop<region_end_plot & LGn==LGn_plot,]
    gene_name_plot<-gene_info_plot[full==TRUE,gene.name]
    gene_id_plot<-gene_info_plot[full==TRUE,GeneID]
    exon_plot<-gene_info_plot[full==TRUE,.(Start, Stop)]
    SNP_plot<-dat_plot[LGn==LGn_plot & Pos >= region_start_plot & Pos <= region_end_plot,.SD,keyby=Pos]
    plot_region(SNP_plot,focalLG=LGn_plot, minpos=region_start_plot, maxpos=region_end_plot, gene.name=gene_name_plot,gene.id=gene_id_plot, exon=exon_plot,
                statstoplot = statstoplot, y_range=y_range)
    #dev.off()
  } 
}

# function to plot per chromosome map
{
  plot_bin_LG <- function(dat, pop, focalLG, specifybps = F, minpos, maxpos, Genestart, Geneend, Genename, statstoplot){
  npops <- length(statstoplot)
  par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
  if(specifybps == F){
    minpos <- 0
    maxpos <- dat[LGn==focalLG,max(Pos)]
  }
  options(scipen=5)
  col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
  for(i in 1:length(statstoplot)){
    y_lim=c(-0.7,0.7)
    plot(dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos,.(Pos/1000000,eval(parse(text = statstoplot[i])))],pch = 16, cex = 0.6, axes = F,ylim=y_lim, xlab = paste("Chromosome",focalLG), ylab = statstoplot[i],col=col[i])
    abline(v=dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos & eval(parse(text = paste(statstoplot[i],".focal",sep=""))),Pos/1000000],col=rgb(red=0,green=0,blue=0,alpha=100,maxColorValue = 255))
    ticks<-seq(0,max(dat[LGn %in% c(focalLG),Pos])/1e6,1)
    axis(1, at=seq(0,maxpos/1e6,0.5), labels = F)
    axis(1, at=seq(0,maxpos/1e6,1), labels = seq(0,maxpos/1e6,1))
    axis(2)
  }
  mtext(paste("Chromosome",focalLG),side=1,line=1, outer = TRUE)
  }
}

# function to plot genome map
{
  plot_bin_genome <- function(dat, statstoplot){

    chrcol <- ifelse(dat[,eval(parse(text = paste(statstoplot,".focal",sep="")))], 1, "grey60")
    y_lim<-c(dat[,min(eval(parse(text = statstoplot)))*1.2],dat[,max(eval(parse(text = statstoplot)))*1.2])
    plot(dat[,Pos_cum],dat[,eval(parse(text = statstoplot))], col=chrcol,axes=F, xlab = "Linkage Group", ylab = statstoplot,pch=16,cex=0.6, ylim = y_lim)
    #text(dat[sig==T,Cumulative_position],dat[sig==T,abs_Z], labels = dat[sig==T,snp.id],cex=0.7, pos=2)
    abline(h=dat[eval(parse(text = paste(statstoplot,".focal",sep="")))==FALSE, max(eval(parse(text = statstoplot)))], lty=5, col='black')
    abline(v=dat[eval(parse(text = paste(statstoplot,".focal",sep="")))==TRUE, Pos_cum], lty=1, col=rgb(0,0,0,alpha=0.3))
    Chr.mid <- dat[,(min(Pos_cum)+max(Pos_cum))/2,by=LGn][,V1]
    axis(2)
    axis(1, at= c(dat[,min(Pos_cum),by=LGn][,V1],dat[,max(Pos_cum)]), labels = c(1:21,"Un",NA))
    rect(dat[,min(Pos_cum),by=LGn][,V1],y_lim[1],dat[,max(Pos_cum),by=LGn][,V1],y_lim[2],border = NA, lwd = 0.1, col= rep(c(NA,rgb(0,0,0,alpha=0.1)),12) )

  }
}






{
  dat_plot <- Tau_slide
  statstoplot = c("pbs_slide")
  
  cutoff_calculation<-function(var, data=dat_plot, threshold=0.99){
    var.cutoff <- quantile(data[,get(var)],threshold,na.rm=T)
    data[,paste0(var,".focal"):= ifelse(get(var)> var.cutoff, T, F)]
  }
  
  lapply(statstoplot,cutoff_calculation)
  
  png(sprintf("./result/%s_RobLaw_BootEchoGos_1000SNP.png",statstoplot), width = 5000, height = 2000,res=350)  
  plot_bin_genome(dat=dat_plot,statstoplot)
  dev.off()
  
  for (chromosome in 1:21) {
    png(filename=sprintf("./result/chr_1000SNP_RobLaw_BootEchoGos/chr%s.divergence_0.99.png", chromosome),width = 4000, height = 2000,res=300)
    plot_bin_LG(dat_plot,pop=c("",""), focalLG=chromosome, statstoplot=statstoplot)
    dev.off()
  } 
}







##################################################################################################################################
# plot for evolution poster

{
  plot_bin_LG <- function(dat, pop, focalLG, specifybps = F, minpos, maxpos, Genestart, Geneend, Genename, statstoplot){
    npops <- length(statstoplot)
    par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
    if(specifybps == F){
      minpos <- 0
      maxpos <- dat[LGn==focalLG,max(Pos)]
    }
    options(scipen=5)
    col<-c("darkblue","darkred","darkgreen","darkorange","darkmagenta")
    for(i in 1:length(statstoplot)){
      y_lim=c(-0.7,0.7)
      plot(dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos,.(Pos/1000000,eval(parse(text = statstoplot[i])))],pch = 20,cex = 0.5,axes = F,ylim=y_lim, xlab = paste("Chromosome",focalLG), ylab = "Populatin Distance",col="black", cex.lab=1.5)
      abline(v=dat[LGn %in% c(focalLG) & Pos > minpos & Pos < maxpos & eval(parse(text = paste(statstoplot[i],".focal",sep=""))),Pos/1000000],col=rgb(red=0,green=0,blue=0,alpha=100,maxColorValue = 255))
      ticks<-seq(0,max(dat[LGn %in% c(focalLG),Pos])/1e6,1)
      axis(1, at=seq(0,maxpos/1e6,0.5), labels = F)
      axis(1, at=seq(0,maxpos/1e6,1), labels = seq(0,maxpos/1e6,1),cex.axis=1.5)
      axis(2,cex.axis=1.2)
    }
    mtext(paste("Chromosome",focalLG),side=1,line=1, outer = TRUE, cex=1.5)
  }



{
  dat_plot <- Tau_slide
  statstoplot = c("pbs_slide")
  
  cutoff_calculation<-function(var, data=dat_plot, threshold=0.99){
    var.cutoff <- quantile(data[,get(var)],threshold,na.rm=T)
    data[,paste0(var,".focal"):= ifelse(get(var)> var.cutoff, T, F)]
  }
  
  lapply(statstoplot,cutoff_calculation)
  
  
  for (chromosome in 12) {
    pdf(file="./result/Evol2019_poster_pbs_LG12.pdf", width = 15, height = 4)

    plot_bin_LG(dat_plot,pop=c("",""), focalLG=chromosome, statstoplot=statstoplot)
    dev.off()
  } 
}
}




{
  plot_region <- function(dat, focalLG, minpos, maxpos, gene.name, gene.id, exon, statstoplot, y_range, y_lab){
    npops <- length(statstoplot)
    par(mfrow = c(npops,1),mar = c(1.5,4,0.5,0.5), mgp = c(2, 0.75, 0), oma = c(3,0,2,0))
    
    options(scipen=5)
    col<-c("black","black","darkgreen","darkorange","darkmagenta")

    for(i in 1:length(statstoplot)){
      y_lim = unlist(y_range[,get(statstoplot[i])])
      plot(dat[,.(Pos,eval(parse(text = statstoplot[i])))], pch = 20, cex = 0.5,cex.lab = 1.5, axes = F, ylim=y_lim, ylab = y_lab[i],col=col[i])
      rect(exon[,Start], y_lim[1],  exon[,Stop],y_lim[2], border = NA , lwd = 1, col = rgb(0,0,0,0.1))
      #text(exon[,(Start+Stop)/2],y_lim[2]/1.2,labels = paste(substr(gene.id,14,nchar(gene.id)),gene.name,sep=" "),cex=0.7, pos=1,srt=90)
      axis(1, at=seq(minpos,maxpos,0.025E6), labels = F)
      axis(1, at=seq(minpos,maxpos,0.05E6), labels = seq(minpos/1e6,maxpos/1e6,0.05), cex.axis=1.5)
      axis(2, cex.axis=1.2)
    }
    mtext(paste("LG",focalLG,":", minpos,"-",maxpos),side=1,line=1, outer = TRUE,cex=1.5)
  }
  
  region_to_plot<-data.table(
    LGn = 12,
    start = 12.75e6,
    end = 12.9e6
  )
  statstoplot = c("gos_vs_rob")
  dat_plot = Fst
  y_range = dat_plot[,..statstoplot][,lapply(.SD,function(x) list(min(x),max(x)))]
  for (i in 1:region_to_plot[,.N]) {
    LGn_plot<-region_to_plot[i,LGn]
    region_start_plot<-region_to_plot[i,start]
    region_end_plot<-region_to_plot[i,end]
    
    #png(filename=sprintf("./result/Region_zoomin/LG%s_%s-%s_RobLaw_BootEchoGos.png",LGn_plot,region_start_plot/1e6,region_end_plot/1e6),width = 3500, height = 2000,res=300)
    pdf(file="./result/Evol2019_poster_pbs_zoominLG12.75_FstOnly.pdf", width = 15, height = 4)
    gene_info_plot<-gene_info[Start>region_start_plot & Stop<region_end_plot & LGn==LGn_plot,]
    gene_name_plot<-gene_info_plot[full==TRUE,gene.name]
    gene_id_plot<-gene_info_plot[full==TRUE,GeneID]
    exon_plot<-gene_info_plot[full==TRUE,.(Start, Stop)]
    SNP_plot<-dat_plot[LGn==LGn_plot & Pos >= region_start_plot & Pos <= region_end_plot,.SD,keyby=Pos]
    plot_region(SNP_plot,focalLG=LGn_plot, minpos=region_start_plot, maxpos=region_end_plot, gene.name=gene_name_plot,gene.id=gene_id_plot, exon=exon_plot,
                statstoplot = statstoplot, y_range=y_range, y_lab<-c("Fst (Rob-Gos)"))
    dev.off()
  } 
}

# plot coverage
{
  gos_cov<- fread("./processed_data/Read_coverage/Coverage_Gos.csv", header = T)
  rob_cov<- fread("./processed_data/Read_coverage/Coverage_Rob.csv", header = T)
  gos_cov_12.75<-gos_cov[Position>=12.75e6 & Position <=12.9e6]
  rob_cov_12.75<-rob_cov[Position>=12.75e6 & Position <=12.9e6]
  rm(gos_cov,rob_cov)
  
  
  slide_data<-rob_cov_12.75
  slide_data[,LGn:=12]
  slide_data[,Pos_join:= rowid(LGn)]
  cov_slide<-slide_data[chr_slide, 
                        on=c("LGn","Pos_join>window_start","Pos_join<window_end"), 
                        allow.cartesian=T][,.("Pos" = (max(Position)+min(Position))/2, 
                                              "Pos_cum" = (max(Position)+min(Position))/2 + max(cum_chr_start),
                                              "Coverage" = mean(Coverage,na.rm=T)),
                                           by=.(LGn,Pos_join)]


  region_to_plot<-data.table(
    LGn = 12,
    start = 12.75e6,
    end = 12.9e6
  )
  dat_plot=cov_slide
  statstoplot = c("Coverage")
  y_range = data.table(Coverage=c(0,300))
  for (i in 1:region_to_plot[,.N]) {
    LGn_plot<-region_to_plot[i,LGn]
    region_start_plot<-region_to_plot[i,start]
    region_end_plot<-region_to_plot[i,end]
    
    #png(filename=sprintf("./result/Region_zoomin/LG%s_%s-%s_RobLaw_BootEchoGos.png",LGn_plot,region_start_plot/1e6,region_end_plot/1e6),width = 3500, height = 2000,res=300)
    pdf(file="./result/Evol2019_poster_pbs_zoominLG12.75_rob_cov.pdf", width = 15, height = 4)
    gene_info_plot<-gene_info[Start>region_start_plot & Stop<region_end_plot & LGn==LGn_plot,]
    gene_name_plot<-gene_info_plot[full==TRUE,gene.name]
    gene_id_plot<-gene_info_plot[full==TRUE,GeneID]
    exon_plot<-gene_info_plot[full==TRUE,.(Start, Stop)]
    SNP_plot<-dat_plot[LGn==LGn_plot & Pos >= region_start_plot & Pos <= region_end_plot,.SD,keyby=Pos]
    plot_region(SNP_plot,focalLG=LGn_plot, minpos=region_start_plot, maxpos=region_end_plot, gene.name=gene_name_plot,gene.id=gene_id_plot, exon=exon_plot,
                statstoplot = statstoplot, y_range=y_range, y_lab="Coverage")
    dev.off()
  } 
}