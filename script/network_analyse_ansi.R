#!/usr/bin/env Rscript
# ?????
package_list <- c("reshape2","ggplot2","devtools","bindrcpp", "VennDiagram",
                  "ggthemes","agricolae","dplyr","igraph", "psych","sqldf","ggrepel",
                  "digest","AnnotationDbi", "impute", "GO.db", "preprocessCore","WGCNA","multtest")

for(p in package_list){
  suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
}

# ??????????????
zi.pi<-function(nodes_bulk, z.bulk, modularity_class, degree){
	
	z.bulk[abs(z.bulk)>0] <- 1
	module<-which(colnames(nodes_bulk)==modularity_class)
	module.max<-max(nodes_bulk[,module])
	degree<-which(colnames(nodes_bulk)==degree)
	
	#??????????????????
	bulk.module<-list(NA)
	length(bulk.module)<-module.max
	
	for(i in 1:max(nodes_bulk[,module])){
		bulk.module[[i]]<-z.bulk[which(nodes_bulk[,module]==i),which(nodes_bulk[,module]==i)]
		
		bulk.module[[i]] <- as.data.frame(bulk.module[[i]])
		rownames(bulk.module[[i]])<-rownames(z.bulk)[which(nodes_bulk[,module]==i)]
		colnames(bulk.module[[i]])<-colnames(z.bulk)[which(nodes_bulk[,module]==i)]
	}
	# ???????????????0-1????
	bulk.module
	
	# ???????????????within-module degree z
	
	## ????list
	z_bulk<-list(NA)
	length(z_bulk)<-module.max
	
	for(i in 1:length(z_bulk)){
		z_bulk[[i]]<-bulk.module[[i]][,1]
		z_bulk[[i]]<-as.data.frame(z_bulk[[i]])
		colnames(z_bulk[[i]])<-"z"
		rownames(z_bulk[[i]])<-rownames(bulk.module[[i]])
	}
	
	z_bulk
	
	## ????z???
	for(i in 1:max(nodes_bulk[,module])){
		# ??????????????????1??????????0??z score = 0
		if(length(bulk.module[[i]])==1){
			z_bulk[[i]][,1]<-0
		}else if(sum(bulk.module[[i]])==0){
			z_bulk[[i]][,1]<-0
		}else{
			# z???????????????????????????????????????????????scale????????????
			k <- rowSums(bulk.module[[i]])
			mean <- mean(k)
			sd <- sd(k)
			sd
			if (sd==0){
				z_bulk[[i]][,1]<-0
			}else{
				z_bulk[[i]][,1]<-(k-mean)/sd
			}
		}
	}
	
	## z????
	for(i in 2:max(nodes_bulk[,module])) {
		z_bulk[[i]]<-rbind(z_bulk[[i-1]],z_bulk[[i]])
	}
	z_bulk<-z_bulk[[module.max]]
	z_bulk
	table(z_bulk>1.5)
	
	# ???????????????????????????????????????????????????
	bulk.module1<-list(NA)
	length(bulk.module1)<-module.max
	
	for(i in 1:max(nodes_bulk[,module])){
		bulk.module1[[i]] <- z.bulk[,which(nodes_bulk[,module]==i)]
		bulk.module1[[i]]<- as.data.frame(bulk.module1[[i]])
		rownames(bulk.module1[[i]])<-rownames(z.bulk)
		colnames(bulk.module1[[i]])<-colnames(z.bulk)[which(nodes_bulk[,module]==i)]
	}
	bulk.module1
	
	# ?????????????? among-module connectivity c
	
	## ????list
	c_bulk<-list(NA)
	length(c_bulk)<-module.max
	
	for(i in 1:length(c_bulk)){
		c_bulk[[i]] <- z.bulk[,1]
		c_bulk[[i]]<-as.matrix(c_bulk[[i]])
		colnames(c_bulk[[i]])<-"c"
		rownames(c_bulk[[i]])<-rownames(z.bulk)
		c_bulk[[i]][,1]<-NA
	}
	c_bulk
	
	## ?????????????????????????????
	for(i in 1:max(nodes_bulk[,module])){
		c_bulk[[i]]<- rowSums(bulk.module1[[i]])
		c_bulk[[i]]<-as.matrix(c_bulk[[i]])
		#c_bulk[[i]][which(nodes_bulk$modularity == i),] = c_bulk[[i]][which(nodes_bulk$modularity == i),] -1
		c_bulk[[i]]<-c_bulk[[i]]*c_bulk[[i]]
		colnames(c_bulk[[i]])<-"c"
		rownames(c_bulk[[i]])<-rownames(z.bulk)
	}
	c_bulk[[1]]
	
	## ?????????????????????????????
	for(i in 2:max(nodes_bulk[,module])){
		c_bulk[[i]] <- c_bulk[[i]] + c_bulk[[i-1]]
	}
	c_bulk<-c_bulk[[module.max]]
	
	# 1-????????????????????????/????????????
	# ?????????????????????????????????????????????????y????Pi??????1????????????????????????????????????????????y????Pi??0??
	# ????????????0????????????????????????????????????????????????????????????????????????????????????????????
	for(i in 1:length(c_bulk)){
		if(nodes_bulk$degree[i]==0){
			c_bulk[i] <- 0
		}else{
			c_bulk[i] <- 1-(c_bulk[i]/(nodes_bulk$degree[i]*nodes_bulk$degree[i]))
		}
	}
	
	colnames(c_bulk)<-"c"
	table(c_bulk<0.62)
	
	#z,c????
	z_c_bulk<-c_bulk
	z_c_bulk<-as.data.frame(z_c_bulk)
	z_c_bulk$z<-z_bulk[match(rownames(c_bulk),rownames(z_bulk)),]
	z_c_bulk<-z_c_bulk[,c(2,1)]
	names(z_c_bulk)[1:2]<-c('within_module_connectivities','among_module_connectivities')
	
	z_c_bulk$nodes_id<-rownames(z_c_bulk)
	nodes_bulk$nodes_id<-rownames(nodes_bulk)
	z_c_bulk<-merge(z_c_bulk,nodes_bulk,by='nodes_id')
	z_c_bulk
	
}

# ???????????????????????
matrix2igraph <- function(matr,matr1,r.threshold,p.threshold,cor_method,p_adjMethod){
	#matr=otu
	#matr=otu1
	if(is.na(matr1[1])){
		occor <- corAndPvalue(matr,method = c(cor_method))
	}else{
		occor <- corAndPvalue(matr,matr1,method = c(cor_method))
	}
	##R?
	occor.r <- occor$cor
	
	# multiple test the p values
	if(p_adjMethod != "None"){
		mtadj <- mt.rawp2adjp(unlist(occor$p),proc=p_adjMethod)
		adpcor <- mtadj$adjp[order(mtadj$index),2]
		occor.p <- matrix(adpcor,dim(matr)[2])
		rownames(occor.p) = rownames(occor.r)
		colnames(occor.p) = colnames(occor.r)
	}else{
		occor.p <- occor$p
		table(occor$p<0.5)
	}
	
	write.table(cbind(data.frame(cor=rownames(occor.r)),occor.r), paste0(otu_sample_file, opts$output, "_", i, "_", cor_method,"_all_correlation_result.txt"),col.names = T, row.names = F, quote = F, sep ="\t")
	write.table(cbind(data.frame(cor_pvlaue=rownames(occor.p)),occor.p), paste0(otu_sample_file, opts$output, "_", i, "_", cor_method, "_", p_adjMethod, "_all_correlationPvlaue_result.txt"),col.names = T, row.names = F, quote = F, sep ="\t")
	# ?????????????????????????????????????R????????????????????0
	occor.r[is.na(occor.r)] <- 0
	occor.r[occor$p > p.threshold|abs(occor$cor) < r.threshold] = 0
	
	if(is.na(matr1[1])){
		diag(occor.r) <- 0
	}

	# ????igraph????
	# diag=FALSE ?????????????
	if(is.na(matr1[1])){
		igraph <- graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
	}else{
		expand.matrix <- function(A){
			m <- nrow(A)
			n <- ncol(A)
			B <- matrix(0,nrow = m, ncol = m)
			rownames(B)=rownames(A)
			C <- matrix(0,nrow = n, ncol = n)
			rownames(C)=colnames(A)
			cbind(rbind(B,t(A)),rbind(A,C))
		}
		occor.r = expand.matrix(occor.r)
		igraph <- graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=TRUE,diag=FALSE)
	}
	
	cor_table <<- occor.r
	write.table(cbind(data.frame(cor_pvlaue=rownames(cor_table)),cor_table), paste0(otu_sample_file, opts$output, "_", i, "_", cor_method, "_", p_adjMethod, "_",r.threshold,"_correlationPvlaue_result.txt"),col.names = T, row.names = F, quote = F, sep ="\t")
	print(paste("Points:", length(V(igraph)),"edges:",length(E(igraph))))
	# NOTE:????????weighted=NULL,???????????????????????????t???????????????????????????????????
	# ????????????????????
	# occor.r[occor.r!=0] <- 1
	# igraph <- graph_from_adjacency_matrix(occor.r,mode="undirected",weighted=NULL,diag=FALSE)
	
	# ?????????????????????????????
	# remove isolated nodes?????????????otu???????????otu ???????????????????
	bad.vs <- V(igraph)[degree(igraph) < 1]
	#bad.vs = NULL
	igraph <- delete.vertices(igraph, bad.vs)
	print(paste("After reomoving the single points, points:", length(V(igraph)),"edges:",length(E(igraph))))
	return(igraph)
}

node_pro <- function(igraph,outdir){  
	# ?????
	igraph.degree<-igraph::degree(igraph)
	# ???????????
	igraph.cen.degree<-centralization.degree(igraph)$res
	# ?????????????
	igraph.betweenness<-centralization.betweenness(igraph)$res
	# ?????????
	igraph.closeness<-centralization.closeness(igraph)$res
	
	igraph.node.pro <- cbind(igraph.degree,igraph.closeness,igraph.betweenness,igraph.cen.degree)
	colnames(igraph.node.pro)<-c("igraph.degree","igraph.closeness","igraph.betweenness","igraph.cen.degree")
	igraph.node.pro
}

net_pro<-function(igraph,outdir){
	# network property
	# ?????? The size of the graph (number of edges)
	num.edges <- length(E(igraph)) # length(curve_multiple(igraph))
	num.edges
	# ???????? Order (number of vertices) of a graph
	num.vertices <- length(V(igraph))# length(diversity(igraph, weights = NULL, vids = 	V(igraph)))
	num.vertices
	# ??????(connectance) ??????????????????????????????????????????????????????????????????????????????????????????????
	connectance <- edge_density(igraph,loops=FALSE)# ? graph.density;loops?????TRUE,?????????????self loops??A--A??B--B???????
	connectance
	# ?????(Average degree)
	average.degree <- mean(igraph::degree(igraph))# ?????2M/N,????M ??N ???????????????????????
	average.degree
	# ???????????(Average path length)
	average.path.length <- average.path.length(igraph) # ?mean_distance(igraph) # mean_distance calculates the average path length in a graph
	average.path.length
	# ???(Diameter)
	diameter <- diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL)
	diameter
	# ?????? edge connectivity / group adhesion
	edge.connectivity <- edge_connectivity(igraph)
	edge.connectivity
	# ??????(Clustering coefficient)??????????????????????????????????????????????????????????????????????????????????C????????????????????????????????
	clustering.coefficient <- transitivity(igraph) 
	clustering.coefficient
	no.clusters <- no.clusters(igraph)
	no.clusters
	# ????????(Degree centralization)
	centralization.degree <- centralization.degree(igraph)$centralization
	centralization.degree
	# ??????????(Betweenness centralization)
	centralization.betweenness <- centralization.betweenness(igraph)$centralization 
	centralization.betweenness
	# ??????????(Closeness centralization)
	centralization.closeness <- centralization.closeness(igraph)$centralization
	centralization.closeness
	
	num.pos.edges<-sum(E(igraph)$weight>0)# number of postive correlation
	num.neg.edges<-sum(E(igraph)$weight<0)# number of negative correlation
	
	igraph.network.pro <- rbind(num.edges,num.pos.edges,num.neg.edges,num.vertices,connectance,average.degree,average.path.length,diameter,edge.connectivity,clustering.coefficient,no.clusters,centralization.degree,centralization.betweenness,centralization.closeness)
	rownames(igraph.network.pro)<-c("num.edges","num.pos.edges","num.neg.edges","num.vertices","connectance","average.degree","average.path.length","diameter","edge.connectivity","clustering.coefficient","no.clusters","centralization.degree","centralization.betweenness","centralization.closeness")
	colnames(igraph.network.pro)<- "value"
	igraph.network.pro
}


site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# ???????????????????????????????
if (!suppressWarnings(suppressMessages(require("optparse", character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))) {
  install.packages(p, repos=site)
  require("optparse",character.only=T)
}

# ????
if (TRUE){
  option_list = list(
    make_option(c("-i", "--input"), type="character", default="result/L/otutab1.txt",
                help="Input file; ????????????? [default %default]"),
    make_option(c("-I", "--input1"), type="character", default="",
    						help="Input another file with the same column of input; ????????????? [default %default]"),
    make_option(c("-d", "--design"), type="character", default="result/metadata_L_C.txt",
                help="Design file; ?????????? [default %default]"),
    make_option(c("-o", "--output"), type="character", default=".A",
    						help="Output prefix; ??????.txt??/pdf? [default %default]"),
    make_option(c("-l", "--node_label"), type="character", default="",
                help="Node label of all input file.If not input file, the network will show the input file row name; ????????????? [default %default]"),
    make_option(c("-L", "--node_label_col"), type="numeric", default=6,
                help="Node label file cols.If not input file, the network will show the input file row name.Default 1-7 response kingdom-species; ????????????? [default %default]"),
    make_option(c("-n", "--groupColname"), type="character", default="Treat",
                help="Group column name; ???????? [default %default]"),
    make_option(c("-s", "--normlized"), type="logical", default=T,
    						help="Whether normlize the data; ????? [default %default]"),
    make_option(c("-t", "--filter_abundance"), type="numeric", default=0.0001,
    						help="Filter low abundance node; ???????? [default %default]"),
    make_option(c("-g", "--group"), type="character", default="",
                help="Input the group name.If not input a group name, the script will calculate all group network; If input multi groups,please spearate them by space; ?????????????? [default %default]"),
    make_option(c("-p", "--pvalue"), type="numeric", default=0.05,
                help="Threshold of P-value, ????????? [default %default]"),
    make_option(c("-r", "--cor"), type="numeric", default=0.8,
                help="Threshold of Correlation, ????????? [default %default]"),
    make_option(c("-m", "--cor_method"), type="character", default="spearman",
                help="Choose a method to calculate the correlation of input file??we support<pearson, kendall, spearman>; ???????????? [default %default]"),
    make_option(c("-a", "--p_adjMethod"), type="character", default="Bonferroni",
                help="Choose a method to adjust the pvalue of correlation??we support<Bonferroni, Holm, Hochberg, SidakSS, SidakSD, BH, BY, ABH, TSBH,None>; ????P????? [default %default]"),
    make_option(c("-M", "--calculateModule"), type="logical", default=F,
                help="Whether calculate the module of network; ??????????? [default %default]"),
    make_option(c("-N", "--calculateCoreNode"), type="logical", default=T,
    						help="Whether calculate the core node of network; ?????????????? [default %default]"),
    make_option(c("-S", "--nodeSize"), type="character", default="12,20",
                help="The range of Node size; ?????????? [default %default]"),
    make_option(c("-w", "--width"), type="numeric", default=14,
                help="Figure width; ???? [default %default]"),
    make_option(c("-e", "--height"), type="numeric", default=7,
                help="Figure heidth; ???? [default %default]")
    )
  opts = parse_args(OptionParser(option_list=option_list))
  # usage = paste0("usage: %prog -i ",opts$input," -d ",opts$design," -n ",opts$groupColname," -l ",opts$node_label," -L ",opts$node_label_col,
  # 							 " -p ",opts$pvalue," -r ",opts$cor," -m ",opts$cor_method," -a ",opts$p_adjMethod," -s ",opts$normlized," -t ",opts$filter_abundance,
  # 							 " -M ",opts$calculateModule," -N ",opts$calculateCoreNode," -w ",opts$width," -e ",opts$height," -o ",opts$output)
  # print(usage)
}

otu_sample_file = opts$input
otu_sample_file1 = opts$input1
design_file = opts$design
r.threshold = opts$cor
p.threshold = opts$pvalue
groupColname = opts$groupColname
width = opts$width
height = opts$height
normlized = opts$normlized
calculateCoreNode = opts$calculateCoreNode

library(stringr)
sizeMin=as.numeric(str_split(opts$nodeSize,",",simplify = T)[1])
sizeMax=as.numeric(str_split(opts$nodeSize,",",simplify = T)[2])

print(paste("The threshold of correlation is", r.threshold))
print(paste("The threshold of Pvalue is",opts$p_adjMethod, p.threshold))

# ?????????????????????????????
maxtaxnum=1000
# ????vertex????

# ??????
otu <- read.table(otu_sample_file, header=T, row.names= 1, sep="\t", comment.char="") 
design <- read.table(design_file, header=T, row.names= 1, sep="\t")
idx = rownames(design) %in% colnames(otu)
design = design[idx, ]
otu = otu[, rownames(design)]
if(normlized){
	norm = t(t(otu)/colSums(otu,na=T))
}else{
	norm = as.matrix(otu)
}

#print(colnames(norm))

if(opts$group==""){
	use_group = unique(design[,groupColname])
}else{
	use_group = strsplit(opts$group," ")[[1]]
}

otu_raw = norm

if(file.exists(otu_sample_file1)){
	otu1 <- read.table(otu_sample_file1, header=T, row.names= 1, sep="\t", comment.char="") 
	idx = rownames(design) %in% colnames(otu1)
	design = design[idx, ]
	otu1 = otu1[, rownames(design)]
	#print(head(otu1))
	if(normlized){
		norm1 = t(t(otu1)/colSums(otu1,na=T))
	}else{
		norm1 = as.matrix(otu1)
	}
	otu_raw1 = norm1
}
#print(head(norm1))
plot_stat = c()
for( i in use_group){
	print(i)
	otu = otu_raw
	cx = grep(i, design[,groupColname])
	
	otu<-t(otu)
	otu = otu[cx,]
	#print(head(otu))
	# ?????????
	if(TRUE){
		otu <- otu[,colSums(otu)/sum(otu) > opts$filter_abundance]
	}
	dim(otu)
	# # ?????????????????????????taxonomy.txt?????Phylum
	# ?????????????????????????taxonomy.txt?????Phylum
	if(opts$node_label!=""){
		gcol = opts$node_label_col
		glab = opts$node_label_col
		
		if(class(try(read.table(opts$node_label, header=F, row.names= 1, sep="\t", comment.char=""),silent=T)) == "try-error"){
			otu_tax <- data.frame(x=c(rep(NA,length(colnames(otu)))),row.names=colnames(otu))
		}else{
			otu_tax <- read.table(opts$node_label, header=F, row.names= 1, sep="\t", comment.char="")
			otu_tax <- as.data.frame(otu_tax[colnames(otu),])
		}
		
	}else{
		gcol = 1
		glab = 1
		otu_tax <- as.data.frame(colnames(otu),colnames(otu))
	}
	#removeCols = which(colSums(otu)==0)
	# if(length(removeCols)!=0){
	# 	otu_abundance <- colSums(otu[,-removeCols])
	# 	otu_pro <- cbind(otu_tax[-removeCols,],otu_abundance)
	# 	igraph <- matrix2igraph(otu[, -removeCols], r.threshold, p.threshold, opts$cor_method, opts$p_adjMethod)
	# }else{
	otu_abundance <- colSums(otu)
	otu_pro <- cbind(otu_tax,otu_abundance)
	rownames(otu_pro) = names(otu_abundance)
	#print(head(otu_pro))
	if(file.exists(otu_sample_file1)){
		otu1 = otu_raw1
		otu1 <- t(otu1)
		otu1 = otu1[cx,]
		if(TRUE){
			otu1 <- otu1[,colSums(otu1)/sum(otu1) > opts$filter_abundance]
		}
		
		# # ?????????????????????????taxonomy.txt?????Phylum
		# ?????????????????????????taxonomy.txt?????Phylum
		if(opts$node_label!=""){
			gcol = opts$node_label_col
			glab = opts$node_label_col
			
			if(class(try(read.table(opts$node_label, header=F, row.names= 1, sep="\t", comment.char=""),silent=T)) == "try-error"){
				otu_tax1 <- data.frame(x=c(rep(NA,length(colnames(otu1)))),row.names=colnames(otu1))
			}else{
				otu_tax1 <- read.table(opts$node_label, header=F, row.names= 1, sep="\t", comment.char="")
				otu_tax1 <- as.data.frame(otu_tax1[colnames(otu1),])
			}
			
		}else{
			gcol = 1
			glab = 1
			otu_tax1 <- as.data.frame(colnames(otu1),colnames(otu1))
		}
		#print(head(otu_raw1))
		#removeCols = which(colSums(otu)==0)
		# if(length(removeCols)!=0){
		# 	otu_abundance <- colSums(otu[,-removeCols])
		# 	otu_pro <- cbind(otu_tax[-removeCols,],otu_abundance)
		# 	igraph <- matrix2igraph(otu[, -removeCols], r.threshold, p.threshold, opts$cor_method, opts$p_adjMethod)
		# }else{
		otu_abundance <- colSums(otu1)
		otu_pro1 <- cbind(otu_tax1,otu_abundance)
		rownames(otu_pro1) = names(otu_abundance)
		
		
		if(T){
			otu_pro1 = otu_pro1[!is.na(otu_pro1[,1]),]
			#print(otu1)
			otu1 = otu1[,rownames(otu_pro1)]
		}
		
		igraph <- matrix2igraph(otu,otu1,r.threshold, p.threshold, opts$cor_method, opts$p_adjMethod)
		# yunxu liang ge leixing de wuzhong zunei xiangguan fenxi 
		if(F){
			igraph <- matrix2igraph(cbind(otu,otu1), NA , r.threshold, p.threshold, opts$cor_method, opts$p_adjMethod)
		}
	}else{
			
		igraph <- matrix2igraph(otu, NA , r.threshold, p.threshold, opts$cor_method, opts$p_adjMethod)
		
	}
		igraph.weight <- E(igraph)$weight
		if(T){
		# ????????????????????????????????????
		E(igraph)$weight <-  2**((E(igraph)$weight - 
																min(E(igraph)$weight)) / diff(range(E(igraph)$weight)))
		min(E(igraph)$weight);max(E(igraph)$weight)
		}
	#}
		if(file.exists(otu_sample_file1)){
			otu_pro = rbind(otu_pro,otu_pro1)
		}
	# ??????????
	if(calculateCoreNode){
		V(igraph)$degree <- degree(igraph)
		V(igraph)$modularity <- membership(cluster_fast_greedy(igraph,weights = NULL))
		
		#print(table(V(igraph)$name))
		nodes_list <- data.frame(
			nodes_id = V(igraph)$name, 
			degree = V(igraph)$degree, 
			modularity = V(igraph)$modularity,
			row.names = V(igraph)$name 
		)

		#nodes_list = nodes_list[colnames(occor.r),]
		zi_pi <- zi.pi(nodes_list, cor_table[nodes_list$nodes_id,nodes_list$nodes_id], degree = 'degree', modularity_class = 'modularity')
		zi_pi <- na.omit(zi_pi)   #NA ??????????????? 0 ????
		
		within_module_connectivities_thre = 2.5
		among_module_connectivities_thre = 0.62
		
		table(zi_pi$within_module_connectivities>within_module_connectivities_thre)
		table(zi_pi$among_module_connectivities<among_module_connectivities_thre)
		
		
		zi_pi[which(zi_pi$within_module_connectivities < within_module_connectivities_thre &
									zi_pi$among_module_connectivities < among_module_connectivities_thre),'type'] <- 'Peripherals'
		zi_pi[which(zi_pi$within_module_connectivities < within_module_connectivities_thre &
									zi_pi$among_module_connectivities > among_module_connectivities_thre),'type'] <- 'Connectors'
		zi_pi[which(zi_pi$within_module_connectivities > within_module_connectivities_thre &
									zi_pi$among_module_connectivities < among_module_connectivities_thre),'type'] <- 'Provincial hubs '
		zi_pi[which(zi_pi$within_module_connectivities > within_module_connectivities_thre &
									zi_pi$among_module_connectivities > among_module_connectivities_thre),'type'] <- 'Kinless hubs'
		
		rownames(zi_pi) = zi_pi$nodes_id
		zi_pi_plot = merge(zi_pi,otu_pro[,glab,drop=F],all.x=T,by="row.names")
		write.table(zi_pi_plot, paste0(otu_sample_file, opts$output, ".", i,'_z_c_score_result.txt'), sep = '\t', row.names = FALSE, quote = FALSE)
		
		if(T){
			zi_pi_plot = merge(zi_pi,otu_pro[,glab,drop=F],all.x=T,by="row.names")
			zi_pi_plot[,ncol(zi_pi_plot)] <- ifelse(zi_pi_plot$type!="Peripherals",zi_pi_plot[,ncol(zi_pi_plot)],NA)
			
			colnames(zi_pi_plot)[ncol(zi_pi_plot)] = "label"
			colnames(zi_pi_plot)
			p = ggplot(zi_pi_plot, aes(among_module_connectivities, within_module_connectivities,label=label)) +
				geom_point(aes(color = type),alpha = 0.5, size = 2) +
				scale_color_manual(values = c('gray','red','blue','purple'),
													 limits = c('Peripherals', 'Connectors', 'Provincial hubs ', 'Kinless hubs'))+
				theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
							panel.background = element_blank(), legend.key = element_blank()) +
				labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
				geom_vline(xintercept = among_module_connectivities_thre) +
				geom_hline(yintercept = within_module_connectivities_thre)+
				geom_text_repel(point.padding=NA,max.overlaps = 1000)
			#geom_text(aes(among_module_connectivities, within_module_connectivities, label=label))
			p
			ggsave(paste0(otu_sample_file, opts$output, ".", i,'_net_core_annotation.pdf'),p)
			
		}
		if(T){
			zi_pi_plot = merge(zi_pi,otu_pro[,glab,drop=F],all.x=T,by="row.names")
			zi_pi_plot[,ncol(zi_pi_plot)] <- ifelse(zi_pi_plot$type!="Peripherals",zi_pi_plot$nodes_id,NA)
			
			length(table(zi_pi_plot[,ncol(zi_pi_plot)]))
			colnames(zi_pi_plot)[ncol(zi_pi_plot)] = "label"
			colnames(zi_pi_plot)
			p = ggplot(zi_pi_plot, aes(among_module_connectivities, within_module_connectivities,label=label)) +
				geom_point(aes(color = type),alpha = 0.5, size = 2) +
				scale_color_manual(values = c('gray','red','blue','purple'),
													 limits = c('Peripherals', 'Connectors', 'Provincial hubs ', 'Kinless hubs'))+
				theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'),
							panel.background = element_blank(), legend.key = element_blank()) +
				labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
				geom_vline(xintercept = among_module_connectivities_thre) +
				geom_hline(yintercept = within_module_connectivities_thre)+
				geom_text_repel(point.padding=NA,max.overlaps = 1000)
			#geom_text(aes(among_module_connectivities, within_module_connectivities, label=label))
			p
			ggsave(paste0(otu_sample_file, opts$output, ".", i,'_net_core_otu.pdf'),p)
			}
	}
	
	#print(V(igraph))
	netpro_result <- net_pro(igraph)
	nodepro_result <- node_pro(igraph)
	cat(paste(i, "Pos Cor:", sum(igraph.weight>0),"\n"))# number of postive correlation
	cat(paste(i, "Neg Cor:", sum(igraph.weight<0),"\n"))# number of negative correlation
	
	if(file.exists(otu_sample_file1)){
		value = c(sum(igraph.weight>0),sum(igraph.weight<0),length(grep("Bac_",V(igraph)$name)),length(grep("Fung_",V(igraph)$name)))
		type=c("pos","neg","bac","fung")
	}else{
		value = c(sum(igraph.weight>0),sum(igraph.weight<0),length(V(igraph)))
		type=c("pos","neg","all")
	}
	
	plot_stat = rbind(plot_stat,
					  data.frame(
						   value = value,
						   group = rep(i,length(value)),
						   type = type,
						   variable = c(rep("edge",2),rep("node",length(value)-2)))
					  )
	
	# set edge color??postive correlation ?څ?red, negative correlation?څ?blue
	E.color <- igraph.weight
	E.color <- ifelse(E.color>0, "red",ifelse(E.color<0, "blue","grey"))
	E(igraph)$color <- as.character(E.color)
	
	igraph.otu <- merge(otu_pro,
											data.frame(name = V(igraph)$name,row.names = V(igraph)$name),
											all.y=T,by="row.names") # ?????OTU????
	
	rownames(igraph.otu) = igraph.otu$name
	igraph.otu = igraph.otu[,c(-1,-ncol(igraph.otu))]
	#print(igraph.otu)
	## ????????????????????????????????sinzeMax-sizeMin??
	a= as.numeric(igraph.otu$otu_abundance)
	a=a[!is.na(a)]
	igraph.size <- (sizeMax-sizeMin)/(max(a)-min(a))*(a-min(a)) + sizeMin
	V(igraph)$size <- igraph.size
	V(igraph)$size
	# ???y?????
	#cat(head(as.matrix(igraph.otu)))
	
	if(opts$node_label!=""){
	  igraph.otu[,gcol][is.na(igraph.otu[,gcol])] = "Unassigned"
	  node.col<-igraph.otu[,gcol]
	  node.col<-as.character(node.col)
  	fre.tax <- names(sort(table(node.col),decreasing =T)[1:maxtaxnum])
  	node.col[!(node.col %in% fre.tax)]<-"Rare_groups"
  	node.col<-as.factor(node.col)
  	otu.tax.levels<-levels(node.col)
  	levels(node.col) <- rainbow(length(otu.tax.levels)) # ??????levles?????????????????????????????????
  	V(igraph)$color <- as.character(node.col)
  	node.label <- as.character(igraph.otu[,glab])
  	V(igraph)$label <- node.label
	}else{
	  top_n=1
	  all_degree<- sort(degree(igraph), decr=T)
	  all_degree = all_degree[V(igraph)$name]
	  node.col<-ifelse(all_degree>=quantile(all_degree, prob=1-top_n/100),"Core","Not Core")
	  
	  node.col<-as.factor(node.col)
	  levels(node.col) = c("purple","lightblue")
	  V(igraph)$color <- as.character(node.col)
	  
	  node.label <- V(igraph)$name
	  V(igraph)$label <- node.label
	}
	
	
	#print(node.label)
	
	write_graph(igraph, paste0(otu_sample_file, opts$output, ".", i,"_co-occurrence_network.graphml"),"graphml")
	if(opts$calculateModule){
		## ???????????
		# ????? modularity
		fc <- cluster_fast_greedy(igraph,weights = NULL)# cluster_walktrap 	cluster_edge_betweenness, cluster_fast_greedy, cluster_spinglass
		modularity <- modularity(igraph,membership(fc))
		# ??????????????
		comps <- membership(fc)
		colbar <- rainbow(max(comps))
		V(igraph)$color <- colbar[comps]
		module = "_module"
		nodepro_result = cbind(nodepro_result, data.frame(igraph.module=comps[rownames(nodepro_result)]))
		module_igraph = data.frame()
		for(j in (1:max(comps))){
			sub_graph = induced_subgraph(igraph,comps==j)
			sub_graph_num.vertices <- length(V(sub_graph))
			sub_graph_edge <- length(E(sub_graph))
			sub_graph_connectance <- edge_density(sub_graph,loops=FALSE)
			sub_graph_average.degree <- mean(igraph::degree(sub_graph))
			sub_graph_average.path.length <- average.path.length(sub_graph)
			sub_graph_diameter <- diameter(sub_graph, directed = FALSE, unconnected = TRUE, weights = NULL)
			sub_graph_edge.connectivity <- edge_connectivity(sub_graph)
			sub_graph_clustering.coefficient <- transitivity(sub_graph)
			sub_graph_no.clusters <- no.clusters(sub_graph)
			sub_graph_centralization.degree <- centralization.degree(sub_graph)$centralization
			sub_graph_centralization.betweenness <- centralization.betweenness(sub_graph)$centralization 
			sub_graph_centralization.closeness <- centralization.closeness(sub_graph)$centralization
			sub_graph_num.pos.edges <- sum(E(sub_graph)$weight>0)# number of postive correlation
			sub_graph_num.neg.edges <- sum(E(sub_graph)$weight<0)
			
			module_igraph = rbind(module_igraph, c(paste0("Modules",j),sub_graph_edge,sub_graph_num.pos.edges,sub_graph_num.neg.edges,sub_graph_num.vertices,sub_graph_connectance,sub_graph_average.degree,sub_graph_average.path.length,sub_graph_diameter,sub_graph_edge.connectivity,sub_graph_clustering.coefficient,sub_graph_no.clusters,sub_graph_centralization.degree,sub_graph_centralization.betweenness,sub_graph_centralization.closeness))
		}
		colnames(module_igraph) = c("modules","num.edges","num.pos.edges","num.neg.edges","num.vertices","connectance","average.degree","average.path.length","diameter","edge.connectivity","clustering.coefficient","no.clusters","centralization.degree","centralization.betweenness","centralization.closeness")
		write.table(module_igraph, paste0(otu_sample_file, opts$output,"_", i , module, "_Attributes.txt"),col.names = T, row.names = F, quote = F, sep ="\t")
		
	}else{
		module = ""
	}
	
	nodepro_result = cbind(nodepro_result, data_frame(node_label=node.label),data.frame(type=V(igraph)$color))
	nodepro_result = cbind(nodepro_result, data_frame(size=igraph.size))
	#print(nodepro_result)
	# ???????igraph??weight??????????????layout??????????
	#E(igraph)$weight <- NA
	e <- get.edgelist(igraph,names=F)
	
	library(qgraph)
	
	l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(igraph),
	                                       area=8*(vcount(igraph)^2),
	                                       repulse.rad=(vcount(igraph)^2.7))
	# ????????
	pdf(file = paste0(otu_sample_file, opts$output, ".", i, module ,"_co-occurrence_network.pdf"),
	    width = width, height = height)
	
	# ?????????????vertex.frame.color?????????
	set.seed(321)
	#coords<-layout_with_fr(igraph,niter=9999,grid="nogrid")
	if(length(E(igraph)) > 500){
	  edge_width = 0.01
	  plot(igraph,main=paste(toupper(i), "Co-occurrence network"),
	       vertex.label = NA,
	       vertex.frame.color="black",
	       vertex.label.cex=0.5,
	       edge.lty=1,
	       edge.curved=F,
	       edge.width = edge_width,
	       margin=c(0,0,0,0),vertex.shape="circle",
	       layout=l)
	  #V(igraph)$label = NA
	}else{
	  edge_width = 1
	  plot(igraph,main=paste(toupper(i), "Co-occurrence network"),
	       #vertex.label = NA,
	       vertex.frame.color="black",
	       vertex.label.cex=0.5,
	       edge.lty=1,
	       edge.curved=F,
	       edge.width = edge_width,
	       margin=c(0,0,0,0),vertex.shape="circle",
	       layout=l)
	}
	
	if(opts$node_label!=""){
	  if(length(unique(node.label))==1 && is.na(unique(node.label))){
	    cat(paste("Warning: The nodes of network can not be noted!\n"))
	  }else{
	    legend(-2,1, legend=otu.tax.levels,col=levels(node.col), pch=16,cex=0.8,bty="n")
	  }
	}
	if(opts$calculateModule){
		legend(1.2,0.75,legend=c(paste0("Modules",c(1:max(comps)))),col=colbar,pch=16,cex=1,bty="n")
	}
	legend(1.2,1, legend=c("pos.cor","neg.cor"),col=c("red","blue"),lty=1,lwd=2, bty="n")
	# ????????????????????vertex.label.dist????????vertex.label.color????????vertex.label.degree
	
	dev.off()
	
	# ???????????
	write_graph(igraph, paste0(otu_sample_file, opts$output, "_", i, "_igraph_edgelist.txt"), format="edgelist")
	E(igraph)$weight <- igraph.weight 
	write_graph(igraph, paste0(otu_sample_file, opts$output, "_", i, "_igraph_col.txt"), format="ncol")#?????????
	write.table(cbind(data.frame(network_attributes=rownames(netpro_result)), netpro_result), paste0(otu_sample_file, opts$output, "_", i, "_igraph_networkAttributes.txt"),col.names = T, row.names = F, quote = F, sep ="\t")
	write.table(cbind(data.frame(node_attributes=rownames(nodepro_result)), nodepro_result), paste0(otu_sample_file, opts$output, "_", i, module, "_igraph_nodeAttributes.txt"),col.names = T, row.names = F, quote = F, sep ="\t")
	
}
library(ggplot2)

p = ggplot(data = plot_stat,aes(x=variable,y=value,fill=type))+
	geom_bar(stat="identity",position = position_dodge())+
	theme_bw() + geom_text(aes(label=value),
						   position=position_dodge(.9),
						   size=5,vjust=-0.5,check_overlap = T) +
	facet_grid(.~group)
	
ggsave(filename=paste0(otu_sample_file, opts$output, "stat.pdf"),p)
