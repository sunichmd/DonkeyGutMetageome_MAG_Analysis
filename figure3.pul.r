## PUL------
all_pro = read.table("pul/sample_pro_pul_counts.txt",header = F,sep="\t")
all_protein = sum(all_pro[,2])
all_protein
all_cazy_protein = sum(all_pro[,3])
### 扇形图(unann,ann)-------

ratio <- c(all_cazy_protein/all_protein*100,(100-all_cazy_protein/all_protein*100))
options(digits=3)
disease <- c(paste0("PULs annotation\n",all_cazy_protein,"(",round(all_cazy_protein/all_protein*100,2),"%)"),
             paste0("Others\n",all_protein-all_cazy_protein,"(",round(100-all_cazy_protein/all_protein*100,2),"%)"))

colors <-c('', '#53A85F', '#F1BB72', '#F3B1A0', 
           '#D6E7A3', '#57C3F3', '#476D87',
           '#E59CC4', '#AB3282', '#23452F', '#BD956A')
colors = ["#8AC6D1","#BBDED6","#FAE3D9","#FFB6B9"]
pdf("pul/Anno_pul_percent_pie.pdf")
pie(ratio, labels=disease,
    radius = 1.0,clockwise=45,
    main = "All proteins",border = F,
    col = c('#BBDED6',"#FEE9b2"))
dev.off()

### 箱线图（identity分布）----------
all_type = read.table("pul/all_pul_diamond.txt",header = F,sep="\t")
table(all_type$V3 > 80)
library(ggplot2)
p = ggplot(data = all_type,aes(x="PUL",y=V3)) + 
  geom_boxplot(fill="#BBDED6") +
  theme_classic() +
  xlab(NULL) + ylab("percent identitiles")
ggsave("pul/Anno_pul_type_identity_boxplot.pdf",p,width = 5,height = 3)  

### 统计门-种水平的pul个数---------
for(tax in c("species","genus","family","order","class","phylum")){
  tax_count = read.table(paste0("pul/",tax,"_pul_count.txt"),header=T,sep="\t")
  p = ggplot(data=tax_count, mapping=aes(y=tax_count[,1],x=Count,fill=tax_count[,1]))+
    geom_bar(stat = "identity") +   
    geom_text(aes(label = Count),size = 3,vjust = 0.3,hjust=-0.2) +
    theme_classic() + #去除边框 
    xlab("PUL Count") + ylab(NULL) +
    scale_x_continuous(expand = c(0,0),limits = c(0,max(tax_count$Count)+max(tax_count$Count)/30)) +
    theme(legend.position="none")
  ggsave(paste0("pul/",tax,"_pul_count_barplot.pdf"),p,width = 20,height =10)
}

for(tax in c("species","genus")){
  tax_count = read.table(paste0("pul/",tax,"_pul_count.txt"),header=T,sep="\t")
  p = ggplot(data=tax_count, mapping=aes(y=tax_count[,1],x=Count,fill=tax_count[,1]))+
    geom_bar(stat = "identity") +   
    geom_text(aes(label = Count),size = 3,vjust = 0.3,hjust=-0.2) +
    theme_classic() + #去除边框 
    xlab("PUL Count") + ylab(NULL) +
    scale_x_continuous(expand = c(0,0),limits = c(0,max(tax_count$Count)+max(tax_count$Count)/30)) +
    theme(legend.position="none")
  ggsave(paste0("pul/",tax,"_pul_count_barplot.pdf"),p,width = 20,height =30)
}
