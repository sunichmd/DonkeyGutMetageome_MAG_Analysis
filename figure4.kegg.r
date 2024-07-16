## kegg
all_pro = read.csv("kegg/MAG_gene_kegg_counts1.txt",header = F,sep="\t")
all_protein = sum(all_pro[,2])
all_protein
all_cazy_protein = sum(as.numeric(all_pro[,3]))
all_cazy_protein
### 扇形图(unann,ann)------------

ratio <- c(all_cazy_protein/all_protein*100,(100-all_cazy_protein/all_protein*100))
options(digits=3)
disease <- c(paste0("KEGG annotation\n",all_cazy_protein,"(",round(all_cazy_protein/all_protein*100,2),"%)"),
             paste0("Others\n",all_protein-all_cazy_protein,"(",round(100-all_cazy_protein/all_protein*100,2),"%)"))

pdf("KEGG/Anno_kegg_percent_pie.pdf")
pie(ratio, labels=disease,
    radius = 0.8,clockwise=45,
    main = "All proteins",border = F,
    col = c("#BEAED4","#FEE9b2"))
dev.off()

tax_count = read.csv("kegg/all_ko_L1.txt",header=F,sep="\t")
colnames(tax_count) = c("Var1","Freq")
p = ggplot(data=tax_count, mapping=aes(x=Freq,y=Var1,fill=Var1))+
  geom_bar(stat = "identity") +   
  geom_text(aes(label = Freq),size = 3,vjust = 0.5,hjust=-0.2) +
  theme_classic() + #去除边框 
  xlab("KEGG Level1 Number") + ylab(NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,max(tax_count$Freq)+max(tax_count$Freq)/20)) +
  theme(legend.position="none")
p
ggsave("kegg/Anno_kegg_level1_count_boxplot.pdf",p,width = 15,height = 5)  

p = tax_count %>%
  mutate(p = Freq / sum(Freq)) %>%
  mutate(labels = scales::label_percent(accuracy = 0.01)(Freq / sum(Freq))) %>%
  ggplot(aes("",p, fill = Var1)) +
  geom_col() +
  geom_text(aes(label = labels), position = position_stack(0.3),size=2.5) +
  coord_polar("y") +
  scale_fill_manual(values = c('#D6E7A3', "lightblue", "mistyrose","#FEE9b2", 
                               "lightcyan", "lavender", "cornsilk", "#F4CAE4"),
                    name="KEGG Level1") +
  theme_void()
p
ggsave("kegg/Anno_kegg_level1_count_pie.pdf",p,width = 9,height = 7)

tax_count = read.csv("kegg/all_ko_L2.txt",header=F,sep="\t")
colnames(tax_count) = c("Var1","Freq")
p = ggplot(data=tax_count, mapping=aes(x=Freq,y=Var1,fill=Var1))+
  geom_bar(stat = "identity") +   
  geom_text(aes(label = Freq),size = 3,vjust = 0.5,hjust=-0.2) +
  theme_classic() + #去除边框 
  xlab("KEGG Level2 Number") + ylab(NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,max(tax_count$Freq)+max(tax_count$Freq)/20)) +
  theme(legend.position="none")
p
ggsave("kegg/Anno_kegg_level2_count_boxplot.pdf",p,width = 15,height = 6)

### 热图-----------
library(tidyr)# 使用的gather & spread
library(reshape2) # 使用的函数 melt & dcast 
library(RColorBrewer)

#sample_ann = read.table("物种注释和beta多样性/tax/taxonomy.txt",header = T,row.names = 1,sep="\t")

for(tax in c("L1","L2","L3")){
  n=20
  sample_count = read.csv(paste0("kegg/mag_ko_",tax,"_count.txt"),sep="\t",header = F,row.names = NULL)
  
  gd1_wide<-spread(sample_count,V2,V3)
  rownames(gd1_wide) = gd1_wide[,1]
  
  sample_ann = read.table("物种注释和beta多样性/tax/taxonomy.txt",header = T,row.names = 1,sep="\t")
  
  row_all = sample_ann[rownames(gd1_wide),c("Phylum","Domain"),drop=F]
  row_all = row_all[order(row_all$Domain,row_all$Phylum),]
  
  gd1_wide[is.na(gd1_wide)] <- 0
  #View(gd1_wide[rownames(row_all),-1])
  library(RColorBrewer)
  pheatmap::pheatmap(log2(gd1_wide[rownames(row_all),-1]+1),scale = "column",
                     show_rownames = F,,angle_col = 90,
                     annotation_row = row_all,annotation_names_row = F,
                     cluster_rows = F,cluster_cols = F,
                     color = colorRampPalette(c("navy","white","firebrick3"))(100),
                     border_color = "NA",filename = paste0("kegg/mag_ko_",tax,"_heatmap.pdf"),
                     width = n,height = n
  )
  pheatmap::pheatmap(log2(gd1_wide[rownames(row_all),-1]+1),scale = "column",
                     show_rownames = F,,angle_col = 90,
                     annotation_row = row_all,annotation_names_row = F,
                     cluster_rows = T,cluster_cols = T,
                     color = colorRampPalette(c("navy","white","firebrick3"))(100),
                     border_color = "NA",filename = paste0("kegg/mag_ko_",tax,"_heatmap_cluster.pdf"),
                     width = n,height = n
  )
}
dev.off()

#### kegg统计----------
tax_count = read.csv("result2/result2/kegg/all_ko_L1.txt",header=F,sep="\t")
library(ggplot2)
colnames(tax_count) = c("Var1","Freq")
p = ggplot(data=tax_count, mapping=aes(x=Freq,y=Var1,fill=Var1))+
  geom_bar(stat = "identity") +   
  geom_text(aes(label = Freq),size = 3,vjust = 0.5,hjust=-0.2) +
  theme_classic() + #去除边框 
  xlab("KEGG Level1 Number") + ylab(NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,max(tax_count$Freq)+max(tax_count$Freq)/20)) +
  theme(legend.position="none")
p
ggsave("result2/result2/kegg/Anno_kegg_level1_count_boxplot.pdf",p,width = 15,height = 5)  

tax_count = read.csv("result2/result2/kegg/all_ko_L2.txt",header=F,sep="\t")
library(ggplot2)
colnames(tax_count) = c("Var1","Freq")
p = ggplot(data=tax_count, mapping=aes(x=Freq,y=Var1,fill=Var1))+
  geom_bar(stat = "identity") +   
  geom_text(aes(label = Freq),size = 3,vjust = 0.5,hjust=-0.2) +
  theme_classic() + #去除边框 
  xlab("KEGG Level2 Number") + ylab(NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,max(tax_count$Freq)+max(tax_count$Freq)/20)) +
  theme(legend.position="none")
p
ggsave("result2/result2/kegg/Anno_kegg_level2_count_boxplot.pdf",p,width = 15,height = 8)  

tax_count = read.csv("result2/result2/kegg/all_ko_L3.txt",header=F,sep="\t")
colnames(tax_count) = c("Var1","Freq")
p = ggplot(data=tax_count, mapping=aes(x=Freq,y=Var1,fill=Var1))+
  geom_bar(stat = "identity") +   
  geom_text(aes(label = Freq),size = 3,vjust = 0.5,hjust=-0.2) +
  theme_classic() + #去除边框 
  xlab("KEGG Level3 Number") + ylab(NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,max(tax_count$Freq)+max(tax_count$Freq)/20)) +
  theme(legend.position="none")
p
ggsave("result2/result2/kegg/Anno_kegg_level3_count_boxplot.pdf",p,width = 15,height = 45)
 
## Prevotella kegg 统计
tax_count = read.csv("result2/result2/kegg/Prevotella_ko_L1.txt",header=F,sep="\t")
library(ggplot2)
colnames(tax_count) = c("Var1","Freq")
p = ggplot(data=tax_count, mapping=aes(x=Freq,y=Var1,fill=Var1))+
  geom_bar(stat = "identity") +   
  geom_text(aes(label = Freq),size = 3,vjust = 0.5,hjust=-0.2) +
  theme_classic() + #去除边框 
  xlab("KEGG Level1 Number") + ylab(NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,max(tax_count$Freq)+max(tax_count$Freq)/20)) +
  theme(legend.position="none")
p
ggsave("result2/result2/kegg/Prevotella_Anno_kegg_level1_count_boxplot.pdf",p,width = 15,height = 5)  

tax_count = read.csv("result2/result2/kegg/Prevotella_ko_L2.txt",header=F,sep="\t")
library(ggplot2)
colnames(tax_count) = c("Var1","Freq")
p = ggplot(data=tax_count, mapping=aes(x=Freq,y=Var1,fill=Var1))+
  geom_bar(stat = "identity") +   
  geom_text(aes(label = Freq),size = 3,vjust = 0.5,hjust=-0.2) +
  theme_classic() + #去除边框 
  xlab("KEGG Level2 Number") + ylab(NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,max(tax_count$Freq)+max(tax_count$Freq)/20)) +
  theme(legend.position="none")
p
ggsave("result2/result2/kegg/Prevotella_Anno_kegg_level2_count_boxplot.pdf",p,width = 15,height = 8)  

tax_count = read.csv("result2/result2/kegg/Prevotella_ko_L3.txt",header=F,sep="\t")
colnames(tax_count) = c("Var1","Freq")
p = ggplot(data=tax_count, mapping=aes(x=Freq,y=Var1,fill=Var1))+
  geom_bar(stat = "identity") +   
  geom_text(aes(label = Freq),size = 3,vjust = 0.5,hjust=-0.2) +
  theme_classic() + #去除边框 
  xlab("KEGG Level3 Number") + ylab(NULL) +
  scale_x_continuous(expand = c(0,0),limits = c(0,max(tax_count$Freq)+max(tax_count$Freq)/20)) +
  theme(legend.position="none")
p
ggsave("result2/result2/kegg/Prevotella_Anno_kegg_level3_count_boxplot.pdf",p,width = 15,height = 30) 

