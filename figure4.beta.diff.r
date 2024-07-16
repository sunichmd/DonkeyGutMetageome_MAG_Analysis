## Beta多样性
# Beta多样性----------

## all---------
mkdir -p 物种注释和beta多样性/beta/
usearch -beta_div ra_mag.txt \
-filename_prefix 物种注释和beta多样性/beta/

for i in bray_curtis euclidean jaccard bray_curtis_binary
do
Rscript ${db}/script/beta_pcoa.R \
--input 物种注释和beta多样性/beta/${i}.txt --design metadata.txt \
--group Group --label F --width 89 --height 59 \
--output 物种注释和beta多样性/beta/${i}.pcoa.pdf
Rscript ${db}/script/beta_cpcoa.R \
--input 物种注释和beta多样性/beta/${i}.txt --design metadata.txt \
--group Group --label F --width 89 --height 59 \
--output 物种注释和beta多样性/beta/${i}.cpcoa.pdf
done

i=bray_curtis
Rscript ${db}/script/beta_pcoa.R \
--input ${i}.txt --design metadata.txt \
--group Group --label F --width 89 --height 59 \
--output ${i}.pcoa.pdf

## pair--------
mkdir -p 物种注释和beta多样性/beta/pair3_beta/

for i in bray_curtis euclidean jaccard bray_curtis_binary
do
Rscript ${db}/script/beta_pcoa.R \
--input 物种注释和beta多样性/beta/${i}.txt --design metadata1.txt \
--group Group --label  T --width 89 --height 59 \
--output 物种注释和beta多样性/beta/pair3_beta/${i}.pcoa.pdf
Rscript ${db}/script/beta_cpcoa.R \
--input 物种注释和beta多样性/beta/${i}.txt --design metadata1.txt \
--group Group --label T --width 89 --height 59 \
--output 物种注释和beta多样性/beta/pair3_beta/${i}.cpcoa.pdf
done

## 差异分析
## 两组差异分析-----------
### 两组的所有样本的wilcox test分析-------------

library(vegan)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(stringr)
library(ggrepel)

metadata = read.table("metadata.txt",header=T,sep='\t')
rownames(metadata)=metadata$Sample

#metadata$pair = metadata$Person_name
metadata$group = metadata$Group

pairdata = read.table("MAGs_ra.txt",header=T,sep='\t',row.names = 1,comment.char = "")

idx = metadata$Sample %in% colnames(pairdata)
metadata = metadata[idx,,drop=F]
pairdata = pairdata[, metadata$Sample]

# 标准化为百分比
if(F){
  norm = t(t(pairdata)/colSums(pairdata,na=T)*100)
  # 按丰度筛选标准化特征表和原始值
  idx = rowMeans(norm) > 0
  norm = norm[idx, ]
  colSums(norm)
  pairdata = pairdata[idx, ]
  
  norm1_rm = t(t(pairdata1)/colSums(pairdata1,na=T)*100)
  # 按丰度筛选标准化特征表和原始值
  idx = rowMeans(norm1_rm) > 0
  norm1_rm = norm1_rm[idx, ]
  colSums(norm1_rm)
  pairdata1 = pairdata1[idx, ]
}else{
  norm = pairdata
  norm1_rm = pairdata
}

dim(norm)
print("Your are using unpair Wilcoxon test!")
unique(metadata$Group)
group_listL = c("C","VC","DC")

for(m in 1:2){
  for(n in (m+1):3){
    group_list_compare = paste0(group_listL[m],"-",group_listL[n])
    print(group_list_compare)
    group_list = c(group_listL[m],group_listL[n])
    
    idx = metadata$group %in% group_list[1]
    GroupA = norm[,rownames(metadata[idx,,drop=F])]
    idx = metadata$group %in% group_list[2]
    GroupB = norm[,rownames(metadata[idx,,drop=F])]
    nrDAF = data.frame(list=rownames(norm), row.names =rownames(norm) )
    for(i in 1:nrow(nrDAF)){
      #print(i)
      # analydata = data.frame(Subject = metadata$Group,
      #                        Group = metadata$Group,
      #                        Value = norm[i,metadata$Sample])
      # colnames(analydata)[3]="Value"
      # 对每行Feature进行秩合检验
      FC = (mean(as.numeric(GroupA[i,]))+0.0000001)/(mean(as.numeric(GroupB[i,]))+0.0000001)
      nrDAF[i,2]=log2(FC)
      nrDAF[i,3]=log2(max(as.numeric(c(GroupA[i,],GroupB[i,])))*10000)
      nrDAF[i,4]= suppressWarnings(wilcox.test(as.numeric(GroupA[i,]),as.numeric(GroupB[i,]))$p.value)
    }
    nrDAF=nrDAF[,-1]
    colnames(nrDAF)=c("logFC", "logCPM", "PValue")
    nrDAF$FDR = p.adjust(nrDAF$PValue, method="fdr", dim(nrDAF)[1])
    nrDAF$logCPM=round(nrDAF$logCPM,3)
    pvalue = 0.05
    fdr = 1
    nrDAF$level = ifelse(nrDAF$logFC>0 & nrDAF$PValue<pvalue & nrDAF$FDR<fdr, "Enriched",
                         ifelse(nrDAF$logFC<0 & nrDAF$PValue<pvalue & nrDAF$FDR<fdr, "Depleted",
                                "NotSig"))
    # nrDAF$level=factor(nrDAF$level,levels = c("Enriched","Depleted","NotSig"))
    print(table(nrDAF$level))

    # Add MeanA and MeanB in percentage
    # calculate groupA mean
    A_list = subset(metadata, group %in% group_list[1])
    A_norm = norm[, rownames(A_list)]
    A_mean = as.data.frame(rowMeans(A_norm))
    colnames(A_mean)=group_list[1]
    # calculate groupB mean
    B_list = subset(metadata, group %in% group_list[2])
    B_norm = norm[, rownames(B_list)]
    B_mean = as.data.frame(rowMeans(B_norm))
    colnames(B_mean)=group_list[2]
    
    # merge and reorder
    Mean = round(cbind(A_mean, B_mean, 
                       A_norm, B_norm), 3) #
    Mean = Mean[rownames(nrDAF),]
    # 修正列名
    colnames(nrDAF)[1] = "log2FC"
    colnames(nrDAF)[2] = "log2CPM"
    output=cbind(nrDAF,Mean)
    write.table(output,paste0("compare/",
                              group_list_compare,"-compare-unpair.txt"),
                sep="\t",row.names = T,col.names = NA,quote=F)
  }
}

#### kegg富集分析------------

library(stringr)
kegg_file = read.table("kegg/all_mag_keggpathway.txt",sep="\t",header = T)
length(unique(kegg_file$KEGG.Pathway))

kegg_term_file <- read.table("C:/EasyMicrobiome/kegg/khier1.tsv",sep="\t",header = T)

kegg_term2name = data.frame(ko_term = kegg_term_file$pathwayid,
                            name = kegg_term_file$pathway)

for(group_list_compare in c("C-DC","C-VC","VC-DC")){
  group=str_split(group_list_compare,"-",simplify = T)
  compare_res = read.table(paste0("compare/",
                                  group_list_compare,"-compare-unpair.txt"),sep="\t",header = T,row.names = 1)
  # enrich
  diff_deg = MAG=rownames(compare_res[compare_res$level=="Enriched",])
  print(paste0(group_list_compare,": Enriched MAGs in ",group[1]," ",length(diff_deg)))
  kegg_term2gene = kegg_file[kegg_file$MAG %in% diff_deg,]
  colnames(kegg_term2gene) = c("gene","ko_term")
  kegg_term2gene = kegg_term2gene[,c("ko_term","gene")]
  
  library(clusterProfiler)
  # KEGG
  # 使用enricher函数进行KEGG的富集分析
  kegg_enrich <- enricher(gene=diff_deg,pvalueCutoff = 1,pAdjustMethod = "none",
                          TERM2GENE = kegg_term2gene,TERM2NAME = kegg_term2name,
                          minGSSize = 0,qvalueCutoff = 1)
  kegg_res = as.data.frame(kegg_enrich)[,c(1,2,8,9)]
  kegg_res_ann = merge(kegg_res,kegg_term_file,all.x = T,by.x="Description",by.y="pathway")
  write.csv(kegg_res_ann,paste0("compare/",group_list_compare,"-KEGG_enrichment_Enriched_",group[1],".csv"),row.names = F)
  # detel
  diff_deg = MAG=rownames(compare_res[compare_res$level=="Depleted",])
  print(paste0(group_list_compare,": Enriched MAGs in ",group[2]," ",length(diff_deg)))
  kegg_term2gene = kegg_file[kegg_file$MAG %in% diff_deg,]
  colnames(kegg_term2gene) = c("gene","ko_term")
  kegg_term2gene = kegg_term2gene[,c("ko_term","gene")]
  
  library(clusterProfiler)
  # KEGG
  # 使用enricher函数进行KEGG的富集分析
  kegg_enrich <- enricher(gene=diff_deg,pvalueCutoff = 1,pAdjustMethod = "none",
                          TERM2GENE = kegg_term2gene,TERM2NAME = kegg_term2name,
                          minGSSize = 0,qvalueCutoff = 1)
  kegg_res = as.data.frame(kegg_enrich)[,c(1,2,8,9)]
  kegg_res_ann1 = merge(kegg_res,kegg_term_file,all.x = T,by.x="Description",by.y="pathway")
  write.csv(kegg_res_ann1,paste0("compare/",group_list_compare,"-KEGG_enrichment_Enriched_",group[2],".csv"),row.names = F)
  #View(kegg_res_ann1[kegg_res_ann1$Count>9,])
  library(ggvenn)
  plot_dat = list(kegg_res_ann$Description,kegg_res_ann1$Description)
  names(plot_dat) = group
  p = ggvenn(plot_dat,
             stroke_color = "white",
             fill_color = c("#E41A1C","#1E90FF"),
             set_name_color =c("#E41A1C","#1E90FF"))
  ggsave(filename = paste0("compare/",group_list_compare,"-KEGG_enrichment-venn",".pdf"),p)
  
  g1_uni = kegg_res_ann[kegg_res_ann$ID %in% setdiff(kegg_res_ann$ID,kegg_res_ann1$ID),]
  write.csv(g1_uni,paste0("compare/",group_list_compare,"-KEGG_enrichment_Enriched_",group[1],"-Unique.csv"),row.names = F)
  
  g2_uni = kegg_res_ann1[kegg_res_ann1$ID %in% setdiff(kegg_res_ann1$ID,kegg_res_ann$ID),]
  write.csv(g2_uni,paste0("compare/",group_list_compare,"-KEGG_enrichment_Enriched_",group[2],"-Unique.csv"),row.names = F)
}
```
#### Unique 通路---------------
groupL = c("C","VC","DC")
g1 = read.csv(paste0("compare/C-VC-KEGG_enrichment_Enriched_",groupL[1],"-Unique.csv"))
g11 = read.csv(paste0("compare/C-DC-KEGG_enrichment_Enriched_",groupL[1],"-Unique.csv"))
c_all_un = rbind(g1,g11)

g2 = read.csv(paste0("compare/C-VC-KEGG_enrichment_Enriched_",groupL[2],"-Unique.csv"))
g22 = read.csv(paste0("compare/VC-DC-KEGG_enrichment_Enriched_",groupL[2],"-Unique.csv"))
vc_all_un = rbind(g2,g22)

g3 = read.csv(paste0("compare/C-DC-KEGG_enrichment_Enriched_",groupL[3],"-Unique.csv"))
g33 = read.csv(paste0("compare/VC-DC-KEGG_enrichment_Enriched_",groupL[3],"-Unique.csv"))
dc_all_un = rbind(g3,g33)

library(ggvenn)
plot_dat = list(unique(c_all_un$Description),unique(vc_all_un$Description),unique(dc_all_un$Description))
names(plot_dat) = groupL
p = ggvenn(plot_dat,
           stroke_color = "white",
           fill_color = c("#E41A1C","#1E90FF",'#3cb44b'),
           set_name_color =c("#E41A1C","#1E90FF",'#3cb44b'))
p
ggsave(filename = paste0("compare/Unique-KEGG_enrichment-venn.pdf"),p)

c_id = setdiff(c_all_un$Description,c(intersect(c_all_un$Description,vc_all_un$Description),intersect(c_all_un$Description,dc_all_un$Description)))

file_res = c_all_un[c_all_un$Description%in% c_id,]
View(file_res[,c(3,6,7,1,4)])
write.csv(file_res,paste0("compare/Unique-C-all-KEGG_enrichment_Enriched.csv"),row.names = F)

vc_id = setdiff(vc_all_un$Description,c(intersect(c_all_un$Description,vc_all_un$Description),intersect(vc_all_un$Description,dc_all_un$Description)))
length(vc_id)
file_res = vc_all_un[vc_all_un$Description%in% vc_id,]
View(file_res[,c(3,6,7,1,4)])
write.csv(file_res,paste0("compare/Unique-VC-all-KEGG_enrichment_Enriched.csv"),row.names = F)


dc_id = setdiff(dc_all_un$Description,c(intersect(vc_all_un$Description,dc_all_un$Description),intersect(c_all_un$Description,dc_all_un$Description)))
length(dc_id)
file_res = dc_all_un[dc_all_un$Description%in% dc_id,]
View(file_res[,c(3,6,7,1,4)])
write.csv(file_res,paste0("compare/Unique-DC-all-KEGG_enrichment_Enriched.csv"),row.names = F)


### 两组的配对样本的t.text ----------
group_listL = c("C","VC","DC")
metadata1 = read.table("metadata1.txt",sep="\t",header = T)
metadata1$group = metadata1$Group

for(m in 1:2){
  for(n in (m+1):3){
    group_list_compare = paste0(group_listL[m],"-",group_listL[n])
    #print(group_list_compare)
    group_list = c(group_listL[m],group_listL[n])
    
    metadata2 = metadata1[order(metadata1$pair),]
    rownames(metadata2) = metadata2$Sample
    metadata2 = metadata2[metadata2$Group %in% group_list,]
    print(paste(group_list_compare,"pair sample number:",nrow(metadata2)))
    
    norm1 = norm1_rm[,metadata2$Sample]
    idx = metadata2$group %in% group_list[1]
    GroupA = norm1[,rownames(metadata2[idx,,drop=F])]
    idx = metadata2$group %in% group_list[2]
    GroupB = norm1[,rownames(metadata2[idx,,drop=F])]
    
    
    nrDAF = data.frame(list=rownames(norm1), row.names =rownames(norm1) )
    for(i in 1:nrow(nrDAF)){
      #print(i)
      # 对每行Feature进行秩合检验
      FC = (mean(as.numeric(GroupA[i,]))+0.0000001)/(mean(as.numeric(GroupB[i,]))+0.0000001)
      if(max(as.numeric(c(GroupA[i,],GroupB[i,])))==0){
        next
      }
      nrDAF[i,2]=log2(FC)
      nrDAF[i,3]=log2(max(as.numeric(c(GroupA[i,],GroupB[i,])))*10000)
      nrDAF[i,4]= suppressWarnings(t.test(as.numeric(GroupA[i,]),
                                          as.numeric(GroupB[i,]),paired = T)$p.value)
    }
    nrDAF=nrDAF[,-1]
    colnames(nrDAF)=c("logFC", "logCPM", "PValue")
    nrDAF$FDR = p.adjust(nrDAF$PValue, method="fdr", dim(nrDAF)[1])
    nrDAF$logCPM=round(nrDAF$logCPM,3)
    pvalue = 0.05
    fdr = 1
    nrDAF$level = ifelse(nrDAF$logFC>0 & nrDAF$PValue<pvalue & nrDAF$FDR<fdr, "Enriched",
                         ifelse(nrDAF$logFC<0 & nrDAF$PValue<pvalue & nrDAF$FDR<fdr, "Depleted",
                                "NotSig"))
    # nrDAF$level=factor(nrDAF$level,levels = c("Enriched","Depleted","NotSig"))
    print(table(nrDAF$level))

    # Add MeanA and MeanB in percentage
    # calculate groupA mean
    A_list = subset(metadata2, group %in% group_list[1])
    A_norm = norm1[, rownames(A_list)]
    A_mean = as.data.frame(rowMeans(A_norm))
    colnames(A_mean)=group_list[1]
    # calculate groupB mean
    B_list = subset(metadata2, group %in% group_list[2])
    B_norm = norm1[, rownames(B_list)]
    B_mean = as.data.frame(rowMeans(B_norm))
    colnames(B_mean)=group_list[2]
    
    # merge and reorder
    Mean = round(cbind(A_mean, B_mean, 
                       A_norm, B_norm), 3) #
    Mean = Mean[rownames(nrDAF),]
    # 修正列名
    colnames(nrDAF)[1] = "log2FC"
    colnames(nrDAF)[2] = "log2CPM"
    output=cbind(nrDAF,Mean)
    write.table(output,paste0("compare/pair/",
                              group_list_compare,"-compare-pair-ttest.txt"),
                sep="\t",row.names = T,col.names = NA,quote=F)
  }
}

#### kegg富集分析----------
library(stringr)
kegg_file = read.table("kegg/all_mag_keggpathway.txt",sep="\t",header = T)
length(unique(kegg_file$KEGG.Pathway))

kegg_term_file <- read.table("C:/EasyMicrobiome/kegg/khier1.tsv",sep="\t",header = T)

kegg_term2name = data.frame(ko_term = kegg_term_file$pathwayid,
                            name = kegg_term_file$pathway)

for(group_list_compare in c("C-DC","C-VC","VC-DC")){
  group=str_split(group_list_compare,"-",simplify = T)
  compare_res = read.table(paste0("compare/pair/",
                                  group_list_compare,"-compare-pair-ttest.txt"),sep="\t",header = T,row.names = 1)
  # enrich
  diff_deg = MAG=rownames(compare_res[compare_res$level=="Enriched",])
  print(paste0(group_list_compare,": Enriched MAGs in ",group[1]," ",length(diff_deg)))
  kegg_term2gene = kegg_file[kegg_file$MAG %in% diff_deg,]
  colnames(kegg_term2gene) = c("gene","ko_term")
  kegg_term2gene = kegg_term2gene[,c("ko_term","gene")]
  
  library(clusterProfiler)
  # KEGG
  # 使用enricher函数进行KEGG的富集分析
  kegg_enrich <- enricher(gene=diff_deg,pvalueCutoff = 1,pAdjustMethod = "none",
                          TERM2GENE = kegg_term2gene,TERM2NAME = kegg_term2name,
                          minGSSize = 0,qvalueCutoff = 1)
  kegg_res = as.data.frame(kegg_enrich)[,c(1,2,8,9)]
  kegg_res_ann = merge(kegg_res,kegg_term_file,all.x = T,by.x="Description",by.y="pathway")
  write.csv(kegg_res_ann,paste0("compare/pair/",group_list_compare,"-KEGG_enrichment_Enriched_",group[1],".csv"),row.names = F)
  # detel
  diff_deg = MAG=rownames(compare_res[compare_res$level=="Depleted",])
  print(paste0(group_list_compare,": Enriched MAGs in ",group[2]," ",length(diff_deg)))
  kegg_term2gene = kegg_file[kegg_file$MAG %in% diff_deg,]
  colnames(kegg_term2gene) = c("gene","ko_term")
  kegg_term2gene = kegg_term2gene[,c("ko_term","gene")]
  
  library(clusterProfiler)
  # KEGG
  # 使用enricher函数进行KEGG的富集分析
  kegg_enrich <- enricher(gene=diff_deg,pvalueCutoff = 1,pAdjustMethod = "none",
                          TERM2GENE = kegg_term2gene,TERM2NAME = kegg_term2name,
                          minGSSize = 0,qvalueCutoff = 1)
  kegg_res = as.data.frame(kegg_enrich)[,c(1,2,8,9)]
  kegg_res_ann1 = merge(kegg_res,kegg_term_file,all.x = T,by.x="Description",by.y="pathway")
  write.csv(kegg_res_ann1,paste0("compare/pair/",group_list_compare,"-KEGG_enrichment_Enriched_",group[2],".csv"),row.names = F)
  #View(kegg_res_ann1[kegg_res_ann1$Count>9,])
  library(ggvenn)
  plot_dat = list(kegg_res_ann$Description,kegg_res_ann1$Description)
  names(plot_dat) = group
  p = ggvenn(plot_dat,
             stroke_color = "white",
             fill_color = c("#E41A1C","#1E90FF"),
             set_name_color =c("#E41A1C","#1E90FF"))
  ggsave(filename = paste0("compare/pair/",group_list_compare,"-KEGG_enrichment-venn",".pdf"),p)
  
  g1_uni = kegg_res_ann[kegg_res_ann$ID %in% setdiff(kegg_res_ann$ID,kegg_res_ann1$ID),]
  write.csv(g1_uni,paste0("compare/pair/",group_list_compare,"-KEGG_enrichment_Enriched_",group[1],"-Unique.csv"),row.names = F)
  
  g2_uni = kegg_res_ann1[kegg_res_ann1$ID %in% setdiff(kegg_res_ann1$ID,kegg_res_ann$ID),]
  write.csv(g2_uni,paste0("compare/pair/",group_list_compare,"-KEGG_enrichment_Enriched_",group[2],"-Unique.csv"),row.names = F)
}
#### Unique 通路------------
groupL = c("C","VC","DC")
g1 = read.csv(paste0("compare/pair/C-VC-KEGG_enrichment_Enriched_",groupL[1],"-Unique.csv"))
g11 = read.csv(paste0("compare/pair/C-DC-KEGG_enrichment_Enriched_",groupL[1],"-Unique.csv"))
c_all_un = rbind(g1,g11)

g2 = read.csv(paste0("compare/pair/C-VC-KEGG_enrichment_Enriched_",groupL[2],"-Unique.csv"))
g22 = read.csv(paste0("compare/pair/VC-DC-KEGG_enrichment_Enriched_",groupL[2],"-Unique.csv"))
vc_all_un = rbind(g2,g22)

g3 = read.csv(paste0("compare/pair/C-DC-KEGG_enrichment_Enriched_",groupL[3],"-Unique.csv"))
g33 = read.csv(paste0("compare/pair/VC-DC-KEGG_enrichment_Enriched_",groupL[3],"-Unique.csv"))
dc_all_un = rbind(g3,g33)

library(ggvenn)
plot_dat = list(unique(c_all_un$Description),unique(vc_all_un$Description),unique(dc_all_un$Description))
names(plot_dat) = groupL
p = ggvenn(plot_dat,
           stroke_color = "white",
           fill_color = c("#E41A1C","#1E90FF",'#3cb44b'),
           set_name_color =c("#E41A1C","#1E90FF",'#3cb44b'))
p
ggsave(filename = paste0("compare/pair/Unique-KEGG_enrichment-venn.pdf"),p)

c_id = setdiff(c_all_un$Description,c(intersect(c_all_un$Description,vc_all_un$Description),intersect(c_all_un$Description,dc_all_un$Description)))
length(c_id)
file_res = c_all_un[c_all_un$Description%in% c_id,]
View(file_res[,c(3,6,7,1,4)])
write.csv(file_res,paste0("compare/pair/Unique-C-all-KEGG_enrichment_Enriched.csv"),row.names = F)

vc_id = setdiff(vc_all_un$Description,c(intersect(c_all_un$Description,vc_all_un$Description),intersect(vc_all_un$Description,dc_all_un$Description)))
length(vc_id)
file_res = vc_all_un[vc_all_un$Description%in% vc_id,]
View(file_res[,c(3,6,7,1,4)])
write.csv(file_res,paste0("compare/pair/Unique-VC-all-KEGG_enrichment_Enriched.csv"),row.names = F)


dc_id = setdiff(dc_all_un$Description,c(intersect(vc_all_un$Description,dc_all_un$Description),intersect(c_all_un$Description,dc_all_un$Description)))
length(dc_id)
file_res = dc_all_un[dc_all_un$Description%in% dc_id,]
View(file_res[,c(3,6,7,1,4)])
write.csv(file_res,paste0("compare/pair/Unique-DC-all-KEGG_enrichment_Enriched.csv"),row.names = F)

