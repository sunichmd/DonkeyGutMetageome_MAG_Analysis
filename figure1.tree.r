# gtdbtk建树-----
## gtdbtk infer --msa_file temp/gtdb_classify/align/tax.bac120.user_msa.fasta.gz --out_dir temp/gtdb_infer
## time gtdbtk infer --msa_file temp/gtdb_classify/align/tax.ar53.user_msa.fasta.gz --out_dir temp/gtdb_infer
## bac ---------
library(ggplot2)
library(ggtree)
tree=read.tree("MAGs进化树/filter.unrooted.tree")
data=fortify(tree)

map=read.csv("MAGs进化树/annotation.txt", header=T, sep="\t")
map$Phylum = gsub("p__","",map$Phylum)
# 圆形（线型）

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

colr = c(
  "#FF0000",  # 红色
  "#00FF00",  # 绿色
  "#0000FF",  # 蓝色
  "#FFFF00",  # 黄色
  "#800080",  # 紫色
  "#FFA500",  # 橙色
  "#FFC0CB",  # 粉红色
  "#00FFFF",  # 青色
  "#FFD700",  # 金色
  "#C0C0C0",  # 银色
  "#000000",  # 黑色
  "#FFFFFF",  # 白色
  "#A52A2A",  # 棕色
  "#808080",  # 灰色
  "#8B0000",  # 深红色
  "#006400",  # 深绿色
  "#00008B",  # 深蓝色
  "#FFD700",  # 深黄色
  "#800080",  # 深紫色
  "#FF8C00",  # 深橙色
  "#FF1493",  # 深粉红色
  "#008B8B"   # 深青色
)

custom_colors <- c("#1f77b4", "#aec7e8","#98df8a","#ff7f0e","#d62728","#FFD999", "#9467bd", "#8c564b",
                   "#e377c2", "#7f7f7f", "#2ca02c","#bcbd22", "#17becf", "#ff9896",
                   "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
                   "#ffbb78",  "#A52A2A", "#FFC0CB")
gra=ggtree(tree,layout="circular", size=0.1) %<+% map +
  # 树型、线粗细、末端颜色 + 注释信息
  geom_tiplab(aes(label=NA, col=Phylum), hjust=-2, align=TRUE,linetype = "solid", linesize=1.6,size=1)+
  # 注释、颜色、高度、对其、虚点大小
  theme(legend.title=element_text(face="bold"), legend.position="none",legend.text=element_text(size=rel(0.5)))+
  # 图例位置、文字大小
  xlim(NA, max(data$x)*1.3) +  scale_colour_manual(values = custom_colors)
gra

# ggtree(tree,layout="circular", size=0.1) %<+% map +
#   #geom_highlight(node=842,fill="red")+
#   #geom_nodelab(aes(subset=!isTip,label=node),hjust=-.3,color="blue")
#   # 树型、线粗细、末端颜色 + 注释信息
#   geom_tiplab(aes(label=NA, col=Phylum), hjust=-2, align=TRUE,linetype = "solid", linesize=1.6,size=1)+
#   # 注释、颜色、高度、对其、虚点大小
#   theme(legend.title=element_text(face="bold"), legend.position="none",legend.text=element_text(size=rel(0.5)))+
#   # 图例位置、文字大小
#   xlim(NA, max(data$x)*1.3) +  scale_colour_manual(values = custom_colors)

pdf("tree_circular_line.pdf",width=8,height = 8)
gra
dev.off()

# facet_plot(gra,"1",data = map ,
#            geom = geom_point,mapping = aes(x = map$Phylum))

gra_node = gra$data[1:nrow(map),]
rownames(gra_node) = gra_node$node
gra_colo = rbind(ggplot_build(gra)$data[[3]],ggplot_build(gra)$data[[5]])
rownames(gra_colo) = gra_colo$node
gra_node$colr = gra_colo[rownames(gra_node),"colour"]
table(gra_node$colr)               

gra_node_colr = gra_node[,c("Phylum","colr")]
gra_node_colr = gra_node_colr[!duplicated(gra_node_colr),]

gra_node_num = data.frame(table(gra_node$Phylum))
colnames(gra_node_num)[1] = "Phylum"
gra_node_num = merge(gra_node_num,gra_node_colr,by="Phylum")

gra_node_num$group="Bacteria"
#gra_node_num$group[2:4] <- "b"
p1 = ggplot(data=gra_node_num, mapping=aes(x=Phylum,y=Freq,fill=Phylum))+
  geom_bar(stat = "identity") +   
  scale_fill_manual(values = gra_node_num$colr) +
  coord_flip() + theme(legend.position="none") +
  geom_text(aes(label = Freq),position=position_dodge(width = 0.9),size = 2.5,vjust = 0.4) +
  theme(axis.ticks = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#去除背景
        panel.border = element_blank()) + #去除边框 
  xlab(NULL) + ylab("Number") + facet_grid(group~., switch = "y") +
  theme(strip.placement='outside') +
  theme(strip.background.y = element_rect(fill = "#007896")) +
  scale_y_continuous(expand = c(0,7))
ggsave("MAGs_num.pdf",p1,width = 5,height = 5)
library(ggpubr)
p2 = ggarrange(gra,p1)
ggsave("MAGs_tree.pdf",p2,width = 8)

## 古菌------------
custom_colors <- c("#1f77b4", "#FFC0CB")

library(ggplot2)
library(ggtree)
tree=read.tree("MAGs进化树/filter_ar.unrooted.tree")
data=fortify(tree)

map=read.csv("MAGs进化树/ar_annotation.txt", header=T, sep="\t")
map$Phylum = gsub("p__","",map$Phylum)

gra=ggtree(tree, size=0.1) %<+% map +
  # 树型、线粗细、末端颜色 + 注释信息
  geom_tiplab(aes(label=Phylum, col=Phylum), hjust=0, align=F,linetype = "solid", linesize=1,size=3) +
  geom_point(aes(color=Phylum), size=3) +
  # 注释、颜色、高度、对其、虚点大小
  theme(legend.title=element_text(face="bold"), legend.position="right")+
  # 图例位置、文字大小
  #xlim(NA, max(data$x)*1.3) +
  scale_colour_manual(values = custom_colors)
gra

# ggtree(tree,layout="circular", size=0.1) %<+% map +
#   #geom_highlight(node=842,fill="red")+
#   #geom_nodelab(aes(subset=!isTip,label=node),hjust=-.3,color="blue")
#   # 树型、线粗细、末端颜色 + 注释信息
#   geom_tiplab(aes(label=NA, col=Phylum), hjust=-2, align=TRUE,linetype = "solid", linesize=1.6,size=1)+
#   # 注释、颜色、高度、对其、虚点大小
#   theme(legend.title=element_text(face="bold"), legend.position="none",legend.text=element_text(size=rel(0.5)))+
#   # 图例位置、文字大小
#   xlim(NA, max(data$x)*1.3) +  scale_colour_manual(values = custom_colors)

pdf("tree_circular_line.pdf",width=8,height = 8)
gra
dev.off()

# facet_plot(gra,"1",data = map ,
#            geom = geom_point,mapping = aes(x = map$Phylum))

gra_node = gra$data[1:nrow(map),]
rownames(gra_node) = gra_node$node
gra_colo = rbind(ggplot_build(gra)$data[[3]])
rownames(gra_colo) = gra_colo$node
gra_node$colr = gra_colo[rownames(gra_node),"colour"]
table(gra_node$colr)               

gra_node_colr = gra_node[,c("Phylum","colr")]
gra_node_colr = gra_node_colr[!duplicated(gra_node_colr),]

gra_node_num = data.frame(table(gra_node$Phylum))
colnames(gra_node_num)[1] = "Phylum"
gra_node_num = merge(gra_node_num,gra_node_colr,by="Phylum")

gra_node_num$group="Archaea"
#gra_node_num$group[2:4] <- "b"
p1 = ggplot(data=gra_node_num, mapping=aes(x=Phylum,y=Freq,fill=Phylum))+
  geom_bar(stat = "identity") +   
  scale_fill_manual(values = gra_node_num$colr) +
  coord_flip() + theme(legend.position="none") +
  geom_text(aes(label = Freq),position=position_dodge(width = 0.9),size = 2.5,vjust = 0.4) +
  theme(
    panel.grid.major =element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),#去除背景
    panel.border = element_blank(),
    axis.ticks = element_blank()) + #去除边框 
  xlab(NULL) + ylab("Number") + facet_grid(group~., switch = "y") +
  theme(strip.placement='outside') +
  theme(strip.background.y = element_rect(fill = "#A52A2A"))
p1
ggsave("MAGs_num.pdf",p1,width = 5,height = 5)
library(ggpubr)
p2 = ggarrange(gra,p1)
ggsave("MAGs_tree.pdf",p2,width = 8)

# iqtree 使用 ModelFinder 来自动选择最佳的氨基酸替代模型----
# iqtree2 -s ar.fa -m MFP  -nt 16
# iqtree2 -s bac.fa -m MFP -bb 1000 -nt 16
## bac ---------
library(ggplot2)
library(ggtree)
tree=read.tree("MCI1/tree/bac.fa.iqtree.treefile")
data=fortify(tree)

map=read.csv("MAGs进化树/annotation.txt", header=T, sep="\t")
map$Phylum = gsub("p__","",map$Phylum)
# 圆形（线型）

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

colr = c(
  "#FF0000",  # 红色
  "#00FF00",  # 绿色
  "#0000FF",  # 蓝色
  "#FFFF00",  # 黄色
  "#800080",  # 紫色
  "#FFA500",  # 橙色
  "#FFC0CB",  # 粉红色
  "#00FFFF",  # 青色
  "#FFD700",  # 金色
  "#C0C0C0",  # 银色
  "#000000",  # 黑色
  "#FFFFFF",  # 白色
  "#A52A2A",  # 棕色
  "#808080",  # 灰色
  "#8B0000",  # 深红色
  "#006400",  # 深绿色
  "#00008B",  # 深蓝色
  "#FFD700",  # 深黄色
  "#800080",  # 深紫色
  "#FF8C00",  # 深橙色
  "#FF1493",  # 深粉红色
  "#008B8B"   # 深青色
)

custom_colors <- c("#1f77b4", "#aec7e8","#98df8a","#ff7f0e","#d62728","#FFD999", "#9467bd", "#8c564b",
                   "#e377c2", "#7f7f7f", "#2ca02c","#bcbd22", "#17becf", "#ff9896",
                   "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
                   "#ffbb78",  "#A52A2A", "#FFC0CB")
gra=ggtree(tree,layout="circular", size=0.1) %<+% map +
  # 树型、线粗细、末端颜色 + 注释信息
  geom_tiplab(aes(label=NA, col=Phylum), hjust=-2, align=TRUE,linetype = "solid", linesize=1.6,size=1)+
  # 注释、颜色、高度、对其、虚点大小
  theme(legend.title=element_text(face="bold"), legend.position="none",legend.text=element_text(size=rel(0.5)))+
  # 图例位置、文字大小
  xlim(NA, max(data$x)*1.3) +  scale_colour_manual(values = custom_colors)
gra

# ggtree(tree,layout="circular", size=0.1) %<+% map +
#   #geom_highlight(node=842,fill="red")+
#   #geom_nodelab(aes(subset=!isTip,label=node),hjust=-.3,color="blue")
#   # 树型、线粗细、末端颜色 + 注释信息
#   geom_tiplab(aes(label=NA, col=Phylum), hjust=-2, align=TRUE,linetype = "solid", linesize=1.6,size=1)+
#   # 注释、颜色、高度、对其、虚点大小
#   theme(legend.title=element_text(face="bold"), legend.position="none",legend.text=element_text(size=rel(0.5)))+
#   # 图例位置、文字大小
#   xlim(NA, max(data$x)*1.3) +  scale_colour_manual(values = custom_colors)

pdf("MCI1/tree/tree_circular_line.pdf",width=8,height = 8)
gra
dev.off()

# facet_plot(gra,"1",data = map ,
#            geom = geom_point,mapping = aes(x = map$Phylum))

gra_node = gra$data[1:nrow(map),]
rownames(gra_node) = gra_node$node
gra_colo = rbind(ggplot_build(gra)$data[[3]],ggplot_build(gra)$data[[5]])
rownames(gra_colo) = gra_colo$node
gra_node$colr = gra_colo[rownames(gra_node),"colour"]
table(gra_node$colr)               

gra_node_colr = gra_node[,c("Phylum","colr")]
gra_node_colr = gra_node_colr[!duplicated(gra_node_colr),]

gra_node_num = data.frame(table(gra_node$Phylum))
colnames(gra_node_num)[1] = "Phylum"
gra_node_num = merge(gra_node_num,gra_node_colr,by="Phylum")

gra_node_num$group="Bacteria"
#gra_node_num$group[2:4] <- "b"
p1 = ggplot(data=gra_node_num, mapping=aes(x=Phylum,y=Freq,fill=Phylum))+
  geom_bar(stat = "identity") +   
  scale_fill_manual(values = gra_node_num$colr) +
  coord_flip() + theme(legend.position="none") +
  geom_text(aes(label = Freq),position=position_dodge(width = 0.9),size = 2.5,vjust = 0.4) +
  theme(axis.ticks = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#去除背景
        panel.border = element_blank()) + #去除边框 
  xlab(NULL) + ylab("Number") + facet_grid(group~., switch = "y") +
  theme(strip.placement='outside') +
  theme(strip.background.y = element_rect(fill = "#007896")) +
  scale_y_continuous(expand = c(0,7))
ggsave("MCI1/tree/MAGs_num.pdf",p1,width = 5,height = 5)
library(ggpubr)
p2 = ggarrange(gra,p1)
ggsave("MCI1/tree/Bac_MAGs_iqtree_tree.pdf",p2,width = 11,height =7)

## 古菌------------
custom_colors <- c("#1f77b4", "#FFC0CB")

library(ggplot2)
library(ggtree)
tree=tree=read.tree("MCI1/tree/ar.fa.iqtree.treefile")
data=fortify(tree)

map=read.csv("MAGs进化树/ar_annotation.txt", header=T, sep="\t")
map$Phylum = gsub("p__","",map$Phylum)

gra=ggtree(tree, size=0.1) %<+% map +
  # 树型、线粗细、末端颜色 + 注释信息
  geom_tiplab(aes(label=Phylum, col=Phylum), hjust=0, align=F,linetype = "solid", linesize=1,size=3) +
  geom_point(aes(color=Phylum), size=3) +
  # 注释、颜色、高度、对其、虚点大小
  theme(legend.title=element_text(face="bold"), legend.position="right")+
  # 图例位置、文字大小
  #xlim(NA, max(data$x)*1.3) +
  scale_colour_manual(values = custom_colors)
gra

# ggtree(tree,layout="circular", size=0.1) %<+% map +
#   #geom_highlight(node=842,fill="red")+
#   #geom_nodelab(aes(subset=!isTip,label=node),hjust=-.3,color="blue")
#   # 树型、线粗细、末端颜色 + 注释信息
#   geom_tiplab(aes(label=NA, col=Phylum), hjust=-2, align=TRUE,linetype = "solid", linesize=1.6,size=1)+
#   # 注释、颜色、高度、对其、虚点大小
#   theme(legend.title=element_text(face="bold"), legend.position="none",legend.text=element_text(size=rel(0.5)))+
#   # 图例位置、文字大小
#   xlim(NA, max(data$x)*1.3) +  scale_colour_manual(values = custom_colors)

pdf("MCI1/tree/tree_circular_line.pdf",width=8,height = 8)
gra
dev.off()

# facet_plot(gra,"1",data = map ,
#            geom = geom_point,mapping = aes(x = map$Phylum))

gra_node = gra$data[1:nrow(map),]
rownames(gra_node) = gra_node$node
gra_colo = rbind(ggplot_build(gra)$data[[3]])
rownames(gra_colo) = gra_colo$node
gra_node$colr = gra_colo[rownames(gra_node),"colour"]
table(gra_node$colr)               

gra_node_colr = gra_node[,c("Phylum","colr")]
gra_node_colr = gra_node_colr[!duplicated(gra_node_colr),]

gra_node_num = data.frame(table(gra_node$Phylum))
colnames(gra_node_num)[1] = "Phylum"
gra_node_num = merge(gra_node_num,gra_node_colr,by="Phylum")

gra_node_num$group="Archaea"
#gra_node_num$group[2:4] <- "b"
p1 = ggplot(data=gra_node_num, mapping=aes(x=Phylum,y=Freq,fill=Phylum))+
  geom_bar(stat = "identity") +   
  scale_fill_manual(values = gra_node_num$colr) +
  coord_flip() + theme(legend.position="none") +
  geom_text(aes(label = Freq),position=position_dodge(width = 0.9),size = 2.5,vjust = 0.4) +
  theme(
    panel.grid.major =element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),#去除背景
    panel.border = element_blank(),
    axis.ticks = element_blank()) + #去除边框 
  xlab(NULL) + ylab("Number") + facet_grid(group~., switch = "y") +
  theme(strip.placement='outside') +
  theme(strip.background.y = element_rect(fill = "#A52A2A"))
p1
ggsave("MAGs_num.pdf",p1,width = 5,height = 5)
library(ggpubr)
p2 = ggarrange(gra,p1)
ggsave("MCI1/tree/Ar_MAGs_iqtree_tree.pdf",p2,width = 12,height = 5)


# fasttree 构建最大似然树----
# FastTree -lg < bac.fa >bac_fastree.tree
# FastTree -lg < ar.fa >ar_fastree.tree
## bac ---------
library(ggplot2)
library(ggtree)
tree=read.tree("MCI1/tree/bac_fastree.tree")
data=fortify(tree)

map=read.csv("MAGs进化树/annotation.txt", header=T, sep="\t")
map$Phylum = gsub("p__","",map$Phylum)
# 圆形（线型）

library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

colr = c(
  "#FF0000",  # 红色
  "#00FF00",  # 绿色
  "#0000FF",  # 蓝色
  "#FFFF00",  # 黄色
  "#800080",  # 紫色
  "#FFA500",  # 橙色
  "#FFC0CB",  # 粉红色
  "#00FFFF",  # 青色
  "#FFD700",  # 金色
  "#C0C0C0",  # 银色
  "#000000",  # 黑色
  "#FFFFFF",  # 白色
  "#A52A2A",  # 棕色
  "#808080",  # 灰色
  "#8B0000",  # 深红色
  "#006400",  # 深绿色
  "#00008B",  # 深蓝色
  "#FFD700",  # 深黄色
  "#800080",  # 深紫色
  "#FF8C00",  # 深橙色
  "#FF1493",  # 深粉红色
  "#008B8B"   # 深青色
)

custom_colors <- c("#1f77b4", "#aec7e8","#98df8a","#ff7f0e","#d62728","#FFD999", "#9467bd", "#8c564b",
                   "#e377c2", "#7f7f7f", "#2ca02c","#bcbd22", "#17becf", "#ff9896",
                   "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5",
                   "#ffbb78",  "#A52A2A", "#FFC0CB")
gra=ggtree(tree,layout="circular", size=0.1) %<+% map +
  # 树型、线粗细、末端颜色 + 注释信息
  geom_tiplab(aes(label=NA, col=Phylum), hjust=-2, align=TRUE,linetype = "solid", linesize=1.6,size=1)+
  # 注释、颜色、高度、对其、虚点大小
  theme(legend.title=element_text(face="bold"), legend.position="none",legend.text=element_text(size=rel(0.5)))+
  # 图例位置、文字大小
  xlim(NA, max(data$x)*1.3) +  scale_colour_manual(values = custom_colors)
gra

# ggtree(tree,layout="circular", size=0.1) %<+% map +
#   #geom_highlight(node=842,fill="red")+
#   #geom_nodelab(aes(subset=!isTip,label=node),hjust=-.3,color="blue")
#   # 树型、线粗细、末端颜色 + 注释信息
#   geom_tiplab(aes(label=NA, col=Phylum), hjust=-2, align=TRUE,linetype = "solid", linesize=1.6,size=1)+
#   # 注释、颜色、高度、对其、虚点大小
#   theme(legend.title=element_text(face="bold"), legend.position="none",legend.text=element_text(size=rel(0.5)))+
#   # 图例位置、文字大小
#   xlim(NA, max(data$x)*1.3) +  scale_colour_manual(values = custom_colors)

pdf("MCI1/tree/tree_circular_line.pdf",width=8,height = 8)
gra
dev.off()

# facet_plot(gra,"1",data = map ,
#            geom = geom_point,mapping = aes(x = map$Phylum))

gra_node = gra$data[1:nrow(map),]
rownames(gra_node) = gra_node$node
gra_colo = rbind(ggplot_build(gra)$data[[3]],ggplot_build(gra)$data[[5]])
rownames(gra_colo) = gra_colo$node
gra_node$colr = gra_colo[rownames(gra_node),"colour"]
table(gra_node$colr)               

gra_node_colr = gra_node[,c("Phylum","colr")]
gra_node_colr = gra_node_colr[!duplicated(gra_node_colr),]

gra_node_num = data.frame(table(gra_node$Phylum))
colnames(gra_node_num)[1] = "Phylum"
gra_node_num = merge(gra_node_num,gra_node_colr,by="Phylum")

gra_node_num$group="Bacteria"
#gra_node_num$group[2:4] <- "b"
p1 = ggplot(data=gra_node_num, mapping=aes(x=Phylum,y=Freq,fill=Phylum))+
  geom_bar(stat = "identity") +   
  scale_fill_manual(values = gra_node_num$colr) +
  coord_flip() + theme(legend.position="none") +
  geom_text(aes(label = Freq),position=position_dodge(width = 0.9),size = 2.5,vjust = 0.4) +
  theme(axis.ticks = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#去除背景
        panel.border = element_blank()) + #去除边框 
  xlab(NULL) + ylab("Number") + facet_grid(group~., switch = "y") +
  theme(strip.placement='outside') +
  theme(strip.background.y = element_rect(fill = "#007896")) +
  scale_y_continuous(expand = c(0,7))
ggsave("MCI1/tree/MAGs_num.pdf",p1,width = 5,height = 5)
library(ggpubr)
p2 = ggarrange(gra,p1)
ggsave("MCI1/tree/Bac_MAGs_tree.pdf",p2,width = 11,height =7)

## 古菌------------
custom_colors <- c("#1f77b4", "#FFC0CB")

library(ggplot2)
library(ggtree)
tree=read.tree("MCI1/tree/ar_fastree.tree")
data=fortify(tree)

map=read.csv("MAGs进化树/ar_annotation.txt", header=T, sep="\t")
map$Phylum = gsub("p__","",map$Phylum)

gra=ggtree(tree, size=0.1) %<+% map +
  # 树型、线粗细、末端颜色 + 注释信息
  geom_tiplab(aes(label=Phylum, col=Phylum), hjust=0, align=F,linetype = "solid", linesize=1,size=3) +
  geom_point(aes(color=Phylum), size=3) +
  # 注释、颜色、高度、对其、虚点大小
  theme(legend.title=element_text(face="bold"), legend.position="right")+
  # 图例位置、文字大小
  #xlim(NA, max(data$x)*1.3) +
  scale_colour_manual(values = custom_colors)
gra

# ggtree(tree,layout="circular", size=0.1) %<+% map +
#   #geom_highlight(node=842,fill="red")+
#   #geom_nodelab(aes(subset=!isTip,label=node),hjust=-.3,color="blue")
#   # 树型、线粗细、末端颜色 + 注释信息
#   geom_tiplab(aes(label=NA, col=Phylum), hjust=-2, align=TRUE,linetype = "solid", linesize=1.6,size=1)+
#   # 注释、颜色、高度、对其、虚点大小
#   theme(legend.title=element_text(face="bold"), legend.position="none",legend.text=element_text(size=rel(0.5)))+
#   # 图例位置、文字大小
#   xlim(NA, max(data$x)*1.3) +  scale_colour_manual(values = custom_colors)

pdf("MCI1/tree/tree_circular_line.pdf",width=8,height = 8)
gra
dev.off()

# facet_plot(gra,"1",data = map ,
#            geom = geom_point,mapping = aes(x = map$Phylum))

gra_node = gra$data[1:nrow(map),]
rownames(gra_node) = gra_node$node
gra_colo = rbind(ggplot_build(gra)$data[[3]])
rownames(gra_colo) = gra_colo$node
gra_node$colr = gra_colo[rownames(gra_node),"colour"]
table(gra_node$colr)               

gra_node_colr = gra_node[,c("Phylum","colr")]
gra_node_colr = gra_node_colr[!duplicated(gra_node_colr),]

gra_node_num = data.frame(table(gra_node$Phylum))
colnames(gra_node_num)[1] = "Phylum"
gra_node_num = merge(gra_node_num,gra_node_colr,by="Phylum")

gra_node_num$group="Archaea"
#gra_node_num$group[2:4] <- "b"
p1 = ggplot(data=gra_node_num, mapping=aes(x=Phylum,y=Freq,fill=Phylum))+
  geom_bar(stat = "identity") +   
  scale_fill_manual(values = gra_node_num$colr) +
  coord_flip() + theme(legend.position="none") +
  geom_text(aes(label = Freq),position=position_dodge(width = 0.9),size = 2.5,vjust = 0.4) +
  theme(
    panel.grid.major =element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),#去除背景
    panel.border = element_blank(),
    axis.ticks = element_blank()) + #去除边框 
  xlab(NULL) + ylab("Number") + facet_grid(group~., switch = "y") +
  theme(strip.placement='outside') +
  theme(strip.background.y = element_rect(fill = "#A52A2A"))
p1
ggsave("MAGs_num.pdf",p1,width = 5,height = 5)
library(ggpubr)
p2 = ggarrange(gra,p1)
ggsave("MCI1/tree/Ar_MAGs_tree.pdf",p2,width = 12,height = 5)
