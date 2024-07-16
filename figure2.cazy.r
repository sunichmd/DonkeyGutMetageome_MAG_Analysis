## cazy------
all_pro = read.table("sample_pro_cazy_counts.txt",header = F,sep="\t")
all_protein = sum(all_pro[,2])
all_protein
all_cazy_protein = sum(all_pro[,3])
### 扇形图(unann,ann)------------

ratio <- c(all_cazy_protein/all_protein*100,(100-all_cazy_protein/all_protein*100))
options(digits=3)
disease <- c(paste0("CAZys annotation\n",all_cazy_protein,"(",round(all_cazy_protein/all_protein*100,2),"%)"),
             paste0("Others\n",all_protein-all_cazy_protein,"(",round(100-all_cazy_protein/all_protein*100,2),"%)"))

colors <-c('', '#53A85F', '#F1BB72', '#F3B1A0', 
           '#D6E7A3', '#57C3F3', '#476D87',
           '#E59CC4', '#AB3282', '#23452F', '#BD956A')
colors = ["#8AC6D1","#BBDED6","#FAE3D9","#FFB6B9"]
pdf("dbcan2/Anno_cazy_percent_pie.pdf")
pie(ratio, labels=disease,
    radius = 1.0,clockwise=45,
    main = "All proteins",border = F,
    col = c('#F3B1A0',"#FEE9b2"))
dev.off()

### 统计CBM、CE、GH、GT、PL、AA --------
awk 'BEGIN{OFS=FS="\t"}{
if(FNR==1) print "gene\ttype\tidenti";\
if($2~/CBM/) print $1,"CBM",$3;\
if($2~/CE/) print $1,"CE",$3;\
if($2~/GH/) print $1,"GH",$3;\
if($2~/GT/) print $1,"GT",$3;\
if($2~/PL/) print $1,"PL",$3;\
if($2~/AA/) print $1,"AA",$3}' all_cazy_diamond.txt > cazy_type_count.txt


all_cazy = read.table("dbcan2/all_cazy_diamond.txt",header = F,sep="\t")
all_type = read.table("dbcan2/cazy_type_count.txt",header=T,sep="\t")

for(i in 1:nrow(all_cazy)){
  if(grepl("CBM",all_cazy$V2[i])){
    tmp_type = data.frame(gene = all_cazy[i,1],
                          type = "CBM",
                          identi = all_cazy[i,3])
    all_type = rbind(all_type,tmp_type)
  }
  if(grepl("CE",all_cazy$V2[i])){
    tmp_type = data.frame(gene = all_cazy[i,1],
                          type = "CE",
                          identi = all_cazy[i,3])
    all_type = rbind(all_type,tmp_type)
  }
  if(grepl("GH",all_cazy$V2[i])){
    tmp_type = data.frame(gene = all_cazy[i,1],
                          type = "GH",
                          identi = all_cazy[i,3])
    all_type = rbind(all_type,tmp_type)
  }
  if(grepl("GT",all_cazy$V2[i])){
    tmp_type = data.frame(gene = all_cazy[i,1],
                          type = "GT",
                          identi = all_cazy[i,3])
    all_type = rbind(all_type,tmp_type)
  }
  if(grepl("PL",all_cazy$V2[i])){
    tmp_type = data.frame(gene = all_cazy[i,1],
                          type = "PL",
                          identi = all_cazy[i,3])
    all_type = rbind(all_type,tmp_type)
  }
  if(grepl("AA",all_cazy$V2[i])){
    tmp_type = data.frame(gene = all_cazy[i,1],
                          type = "AA",
                          identi = all_cazy[i,3])
    all_type = rbind(all_type,tmp_type)
  }
}

### 饼图 Cazy类型组成--------
awk 'BEGIN{OFS=FS="\t"}FNR>1{a[$2]+=1}\
END{print "type\tCount";for(i in a) print i,a[i]}' cazy_type_count.txt >cazy_type_count1.txt

type_count = read.table("dbcan2/cazy_type_count1.txt",header=T,sep = "\t")
ratio <- type_count$Count/sum(type_count$Count)
options(digits=3)
disease <- paste0(type_count$type,"\n",round(ratio*100,2),"%")

colors <-c('', '#53A85F', '#F1BB72', '#F3B1A0', 
           '#D6E7A3', '#57C3F3', '#476D87',
           '#E59CC4', '#AB3282', '#23452F', '#BD956A')
colors = ["#8AC6D1","#BBDED6","#FAE3D9","#FFB6B9"]

library(dplyr)
p = type_count %>%
    mutate(p = Count / sum(Count)) %>%
    mutate(labels = scales::label_percent(accuracy = 0.01)(Count / sum(Count))) %>%
    ggplot(aes("",p, fill = type)) +
    geom_col() +
    geom_text(aes(label = labels), position = position_stack(0.3),size=3.2) +
    coord_polar("y") +
    scale_fill_manual(values = c('#D6E7A3', "lightblue", "mistyrose", 
                                 "lightcyan", "lavender", "cornsilk")) +
    theme_void()
p
ggsave("dbcan2/Anno_cazy_type_percent_pie.pdf",p,width = 5,height = 5)
### 箱线图（identity分布）----------
table(all_type$identi > 85)
library(ggplot2)
p = ggplot(data = all_type,aes(x=type,y=identi,fill=type)) + 
  geom_boxplot() +
  theme_classic() +
  xlab(NULL) + ylab("percent identitiles") +
  scale_fill_manual(values = c('#D6E7A3', "lightblue", "mistyrose", 
                               "lightcyan", "lavender", "cornsilk"))
ggsave("dbcan2/Anno_cazy_type_identity_boxplot.pdf",p,width = 5,height = 3)  