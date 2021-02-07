library(data.table)

setwd("/Users/pengfoen/OneDrive - University of Connecticut/Poolseq_8pops")

# load and process Amanda's data on 2018
fib_dat_2018<-fread("./processed_data/Fibrosis_wild_Hund_Wild_Parasites_Data_2018.csv")
fib_dat_2018<-fib_dat_2018[,.("lake"=Lake,
                              "sex"=Sex_Fish,
                              "infected"=as.numeric(as.logical(Infected)),
                              "fibrosis"=as.numeric(as.logical(as.numeric(substr(Fibrosis_Score,1,1)))),
                              "worm_count"=WormCount,
                              "worm_mass_total"=Mass_All_worms)]
fib_dat_2018[,lake:=paste0(lake,"_2018")]


# load and process Will's data on 2016
fib_dat_2016<-fread("./processed_data/Fibrosis_wild_Shimcurateddata.csv")
pop_list<-c("Boot","Echo","Frederick ","Gosling","Lawson","Pachena","Roberts","Rosselle")
fib_dat_2016<-subset(fib_dat_2016, lake%in%pop_list, c(3,6,7,9,41,46))
fib_dat_2016[,lake:=paste0(lake,"_2016")]
setnames(fib_dat_2016,c("Nworms","wormmass.total"), c("worm_count","worm_mass_total"))

# combine datasets and calculate fibrosis rate, 
# assumbing fish with fibrosis but no worm was exposed to worm previously, fibrosis rate is calculated as (fishes with fibrosis/(fishes with fibrosis+fish has worms but no fibrosis))
fib_both<-rbind(fib_dat_2016,fib_dat_2018)
fib_both<-na.omit(fib_both, cols=c("infected","fibrosis"))  
fib_both[,exposed:=0]
fib_both[!(infected==0 & fibrosis==0),exposed:=1]
fib_summary<-fib_both[,.(num_fishes=.N,
                         num_fibrosis=sum(fibrosis),
                         num_exposed=sum(exposed),
                         num_infected=sum(infected),
                         fibrosis_rate=sum(fibrosis)/sum(exposed)),
                      by=lake]
setkey(fib_summary,fibrosis_rate)

# plot the results
png(filename="./result/Fibrosis_rate_wild_pop.png",width = 4000, height = 2000,res=300)
xx<-barplot(fib_summary[,fibrosis_rate*100],names.arg=fib_summary[,lake],ylim=c(0,120),ylab = "Fibrosis rate",cex.names=0.9)
text(x = xx, y = fib_summary[,fibrosis_rate*100], label = paste0("N=",fib_summary[,num_fishes]," ",fib_summary[,round(fibrosis_rate,2)*100],"%"), pos = 3, cex = 0.9, col = "red")
dev.off()

fib_dat_2016[,.(num_fishes=.N,
             worm_mass_total=mean(worm_mass_total,na.rm = T),
             worm_mass_avg=mean(worm_mass_total/worm_count,na.rm = T),
             num_infected=sum(infected)),by=lake]

fib_both[,.(num_fishes=.N,
                  worm_mass_total=mean(worm_mass_total,na.rm = T),
                  worm_mass_avg=mean(worm_mass_total/worm_count,na.rm = T),
                  num_infected=sum(infected)),by=lake]

png(filename="./result/Avg_worm_mass_wild_pop.png",width = 4000, height = 2000,res=300)
tt<-barplot(fib_summary[,worm_mass_avg],names.arg=fib_summary[,lake],ylim=c(0,120),ylab = "Average Worm Mass",cex.names=0.9)
text(x = tt, y = fib_summary[,worm_mass_avg], label = paste0("N=",fib_summary[,num_fishes]," ",fib_summary[,round(worm_mass_avg,2)],"g"), pos = 3, cex = 0.9, col = "red")
dev.off()

