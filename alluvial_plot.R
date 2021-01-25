library(RColorBrewer)
library(ggalluvial)
library(reshape2)
library(dplyr)
library(data.table)


# Read metadata
# Metadata needs columns: file_name, donor, cell_type, time_point
## file_name: contains tab-delimited .txt file name
## donor: unique name for samles donor
## cell_type: cell type
## time_point: time point when sample was harvested. Must be in numeric format.

## example:
## file_name	sample_id	..filter..	group	time_point	cell_type
## ch45_14.txt	ch45_14	conv:MiXcr,ds:24957,ncfilter:remove	Kadyrov	0.5	nCD8
## ch45_15.txt	ch45_15	conv:MiXcr,ds:24957,ncfilter:remove	Kadyrov	1	nCD8
## ch45_16.txt	ch45_16	conv:MiXcr,ds:24957,ncfilter:remove	Kadyrov	0.5	Treg
## ch45_17.txt	ch45_17	conv:MiXcr,ds:24957,ncfilter:remove	Kadyrov	1	Treg
## ch45_18.txt	ch45_18	conv:MiXcr,ds:24957,ncfilter:remove	Logacheva	1	mature_nCD4

args = commandArgs(trailingOnly = TRUE)
fn = args[1]
param = args[2]
metadata = read.table(fn, header = T)

mode = param

#metadata = read.table("/home/eshnayder/TCR/TCR_2_points/not_downsample/top1000nc/metadata.txt", header = T)

## Grouping metadata
groups = as.data.frame(metadata %>% group_by(donor, cell_type) %>% summarise(n = n()))
filtered = subset(groups, n > 1)
filtered = merge(filtered, data, by = "donor", all = F)

## isolation names, cell_types and time points
cell_types = unique(groups$cell_type)
name = unique(filtered$donor)
time_points = as.numeric(unique(filtered$time_point))

## Common clonotypes table function
# name - from "donor" column
# cell - from "cell_type" column
# mode - include or not V-segment. mode = "CDR3_V" - include V-segment

clonotypes_tables = function(name, cell, mode){
  sub = subset(metadata, metadata$donor == as.character(name))
  sub =  subset(sub, sub$cell_type == as.character(cell))
  if (nrow(sub) > 1){
    sub = sub[order(sub$time_point), ]
  
    if (file.exists(paste(getwd(), sub$file_name[1], sep="/")) & file.exists(paste(getwd(), sub$file_name[2], sep="/"))){
      point1 = read.table(sub$file_name[1], header = T)
      #cat("check1")
      point2 = read.table(sub$file_name[2], header = T)
      
    
      if (nrow(point1) >= 1000 & nrow(point2) >= 1000){
        #common.clonotype = merge(point1[1:1000,], point2[1:1000,], by = c("cdr3aa", "v"))
        if (mode == "CDR3_V"){
          common.clonotype = merge(point1, point2, by = c("cdr3aa", "v"))
          
          #common.clonotype = inner_join(point1[,c(1,2,4,5)], point2[,c(1,2,4,5)], by = c("cdr3aa", "v"))
          common.clonotype$clon = paste(common.clonotype$cdr3aa, common.clonotype$v, sep = "~")
        }
        else{
          common.clonotype = merge(point1, point2, by = c("cdr3aa"))
          #common.clonotype = inner_join(point1[,c(1,2,4)], point2[,c(1,2,4)], by = c("cdr3aa"))
          common.clonotype$clon = common.clonotype$cdr3aa
        }
        common.clonotype = common.clonotype[,c("clon", "freq.x", "freq.y")]
        colnames(common.clonotype) = c("Clonotypes", "0.5", "1")
        dim(common.clonotype)
        common.clonotype = common.clonotype[order(common.clonotype[2], decreasing = T), ]
        common.clonotype
      }
      #return(common.clonotype)
    }
  else{
    print(paste(name, cell, "hasn't two time points", sep = " ")) 
   }
  }
}

## Alluvium plot function
flow_chart = function(commonclon, name, cell_type){
  melted = melt(commonclon, id = "Clonotypes")
  n = nrow(commonclon)
  p = ggplot(melted, aes(x=variable, y=value, stratum =Clonotypes ,
                     alluvium = Clonotypes, fill = Clonotypes))+ 
    scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(n))+
    geom_flow()+ ggtitle(paste("Clonotypes number:", nrow(commonclon)))+
    geom_stratum()+theme_classic()+ 
      xlab("time point")+ylab("ratio")+
      theme(panel.grid.major.y = element_line( size=.1, color="black" , linetype = "dotted"))
  ggsave(paste0(name, "_",cell_type,".png"),plot = p, width = 14, height = 7)
}


## Check time points' format and make plots
if (any(is.na(time_points))){
  print("Time points aren't in numeric format")
}else{
  for (i in 1:length(name)){
    for (j in 1:length(cell_types)){
      tab = clonotypes_tables(name[i], cell_types[j], mode)
      cat(name[i])
      if (is.data.frame(tab) ){
        tab$Clonotypes[which(duplicated(tab$Clonotypes))] = paste(tab$Clonotypes[which(duplicated(tab$Clonotypes))], which(duplicated(tab$Clonotypes)), sep = "_")
        flow_chart(tab, name[i], cell_types[j])
      }
    }
  }
}

