library(data.table)
library(ggplot2)
library(gggenes) 

setwd("/Users/pengfoen/Documents/Research/My_Paper/2021_Jesse_QTL/CYP3a48_figure")
  genes<-fread( "./LG12.75_GFF3.txt", header = F, sep="\t") 
  genes<-genes[2:9,]
  genes[,c("V2","V3","V6","V8"):=NULL]
  setnames(genes,c("V1","V4","V5","V7","V9"),c("molecule","start","end","strand","gene"))
  genes[,orientation:=ifelse(start<end,1,-1)]
  genes[,gene:=rep("",8)]
  
  gene_view<-ggplot(genes, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
    geom_gene_arrow(arrowhead_height = unit(3, "mm"), arrowhead_width = unit(1, "mm")) +
    geom_gene_label(align = "left") +
    facet_wrap(~ molecule, scales = "free", ncol = 1) +
    scale_fill_brewer(palette = "Set3") +
    theme_genes()
  ggsave(file="gene_view.svg", plot=gene_view, width=10, height=4)
  
  gos_cov<- fread("./gos_cov_12.75.csv", header = T)
  rob_cov<- fread("./rob_cov_12.75.csv", header = T)
  say_cov<- fread("./say_cov_12.75.csv", header = T)
  cov<-gos_cov[rob_cov,on=("Position")][,.("Pos"=Position,
                                                             "rob_cov"=i.Coverage,
                                                             "gos_cov"=Coverage)]
  cov<-say_cov[cov,on=c("Position"="Pos")]
  setnames(cov,"Coverage","say_cov")
  rm(gos_cov,rob_cov,say_cov)
  region_start=12750000
  region_end=12850000
  step=500
  cov_slide=cov[,.("Say_Coverage"=mean(say_cov),"Rob_Coverage"=mean(rob_cov),"Gos_Coverage"=mean(gos_cov),"Pos"=Position[1]+step/2) , keyby = .(gr=cut(Position, breaks=seq(region_start+step, region_end, by=step)))]
 
  cov_view<-ggplot(cov_slide, aes(x=Pos)) + 
    geom_line(aes(y = Rob_Coverage,colour="Rob_Coverage")) + 
    geom_line(aes(y = Gos_Coverage,colour="Gos_Coverage")) +
    geom_line(aes(y = Say_Coverage,colour="Say_Coverage")) +
    scale_colour_manual("", 
                        values = c("Rob_Coverage"="green", "Gos_Coverage"="blue", "Say_Coverage"="black")) +
    scale_y_continuous("Coverage", limits = c(0,300)) +
    scale_x_continuous("Region in LG12 (Mb)",breaks=seq(region_start,region_end,0.01E6), labels = seq(region_start/1e6,region_end/1e6,0.01))+
    theme_classic() 
  ggsave(file="cov_view.svg", plot=cov_view, width=10, height=4)
  

  
