


#####################################################################################
##########################    require packages     ##################################
#####################################################################################

library(Seurat)
library(DUBStepR)
library(SingleCellExperiment)
library(geneBasisR)
library(Seurat)
library(HighlyRegionalGenes)
library(FEAST)
library(rhdf5)
library(ranger)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(reshape2)
library(tidyverse)
library(viridis)
library(paletteer)
library(gridExtra)
library(fpc)

#####################################################################################
##########################   clustering methods   ###################################
#####################################################################################

## Seurat clustering with PCA
Cluster=function(data, gene_list){
  SeuObj = CreateSeuratObject(counts=data)
  SeuObj <- NormalizeData(SeuObj,verbose = FALSE)
  all.gene = rownames(SeuObj)
  SeuObj <- ScaleData(SeuObj, features = all.gene, verbose = FALSE)
  if(ncol(data)<50){
    SeuObj <- RunPCA(SeuObj, features = gene_list,verbose = FALSE,npcs=20)
  }else{
    SeuObj <- RunPCA(SeuObj, features = gene_list,verbose = FALSE)
  }
  if(length(gene_list)<10){
    SeuObj <- FindNeighbors(SeuObj, verbose = FALSE, dims = 1:(length(gene_list)-1))
  }else{
    SeuObj <- FindNeighbors(SeuObj, verbose = FALSE)
  }
  SeuObj <- FindClusters(SeuObj, verbose = FALSE)
  cluster_ident=as.numeric(Idents(SeuObj))
  return(cluster_ident)
}

## Seurat clustering without PCA
Cluster=function(data, gene_list){
  SeuObj = CreateSeuratObject(counts=data)
  SeuObj <- NormalizeData(SeuObj,verbose = FALSE)
  SeuObj <- FindNeighbors(SeuObj, verbose = FALSE, dims = NULL, features = gene_list)
  SeuObj <- FindClusters(SeuObj, verbose = FALSE)
  cluster_ident=as.numeric(Idents(SeuObj))
  return(cluster_ident)
}

#####################################################################################
############## gene selection with other feature selection strategies  ##############
#####################################################################################

HVGfun=function(data, nfeatures=2000, method="vst"){
  pbmc = CreateSeuratObject(counts=data)
  pbmc <- NormalizeData(object = pbmc)
  pbmc <- FindVariableFeatures(pbmc, selection.method = method, nfeatures = nfeatures)
  return(VariableFeatures(pbmc))
}

FEASTfun=function(data, k, n_feature=4000){
  # n_feature=4000
  Y = process_Y(data, thre = 2)
  # ixs = FEAST(Y, k=k)
  ixs  = FEAST_fast(Y, k=k)
  return(rownames(Y)[ixs][1:n_feature])
}

geneBasisRfun=function(data){
  data = as.data.frame(data)
  sce=SingleCellExperiment(assays = list(counts = data, logcounts = log2(data+1)))
  sce = retain_informative_genes(sce)
  return(gene_search(sce, n_genes_total = 50)$gene)
}

DUBStepRfun=function(data, pcs=20, min.cells = 0.05*ncol(data)){
  res <- DUBStepR(input.data = data, num.pcs = pcs, min.cells = min.cells)
  return(res$optimal.feature.genes)
}


HRGfun=function(data){
  pbmc = CreateSeuratObject(counts=data)
  pbmc <- NormalizeData(pbmc,verbose = FALSE)
  all.genes = rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes, verbose = FALSE)
  if(ncol(data)<50){
    pbmc <- RunPCA(pbmc, features = all.genes,verbose = FALSE,npcs=20)
  }else{
    pbmc <- RunPCA(pbmc, features = all.genes,verbose = FALSE)
  }
  pbmc=FindRegionalGenes(pbmc,dims = 1:10,nfeatures = 2000,overlap_stop = 0.95)
  gene_num = HRG_elbowplot(pbmc)
  return(RegionalGenes(pbmc,nfeatures = gene_num))
}


datalist = c("")
medlist = c("DUBStepR","FEAST","geneBasisR","HVG","HRG")
runtimeMat = matrix(nrow = length(medlist),ncol = length(datalist))
colnames(runtimeMat) = datalist
rownames(runtimeMat) = medlist
for (d in datalist) {
  print(d)
  load("")
  for (m in medlist) {
    print(m)
    if(m=="DUBStepR"){
      sgenes = DUBStepRfun(data[rownames(x),])
    }else if(m=="FEAST"){
      sgenes = FEASTfun(data,k = length(table(label)),n_feature = 2000)
    }else if(m=="geneBasisR"){
      sgenes = geneBasisRfun(data)
    }else if(m=="HRG"){
      sgenes = HRGfun(data)
    }else{
      sgenes = HVGfun(data)
    }
    runningtime <- proc.time() - timestart
    runtimeMat[m, as.character(d)] <- runningtime[[3]]
  }
}

#####################################################################################
# 1. Comparison of clustering performance with other feature selection strategies  ##
#####################################################################################

# for Rdata scRNAseq data sets
datalist = c()
medlist = c("CellBRF", "DUBStepR","Feats","FEAST","geneBasisR","HVG", "HRG")
resmat1 = matrix(nrow = length(datalist), ncol = length(medlist))
colnames(resmat1) = medlist
rownames(resmat1) = datalist
resmat2 = resmat1
for (d in datalist) {
  print(d)
  load("")
  for (m in medlist) {
    print(m)
    if(m=="CellBRF"){
      sgenes = read.table(paste0("./",d,"_CellBRF_gs_res.txt"))$V1
    }else if(m=="HVG"){
      sgenes <- HVGfun(data)
    }else {
      sgenes = read.table(paste0("./",m,"_",d,".txt"))$V1
    }
    print(length(sgenes))
    rownames(data) = gsub("[_]","-",rownames(data))
    sgenes = gsub("[_]","-",sgenes)
    clusters <- Cluster(data = data, gene_list = sgenes)
    resmat1[d,m] <- evalcluster(label, clusters)[[1]]
    resmat2[d,m] <- evalcluster(label, clusters)[[3]]
  }
}

# for h5 data
datalist = list.files("./h5data/")
datalist = c()
medlist = c("CellBRF", "DUBStepR","Feats","FEAST","geneBasisR","HVG", "HRG")
resmat1 = matrix(nrow = length(datalist), ncol = length(medlist))
colnames(resmat1) = medlist
rownames(resmat1) = datalist
resmat2 = resmat1
for (d in datalist) {
  print(d)
  x = h5read(paste0("./h5data/",datalist[i]),"X")
  y = h5read(paste0("./h5data/",datalist[i]),"Y")
  data = x
  label = y
  rownames(data) = paste0("gene",1:nrow(data))
  colnames(data) = paste0("cell",1:ncol(data))
  for (m in medlist) {
    print(m)
    if(m=="CellBRF"){
      sgenes = read.table(paste0("./",d,"_CellBRF_gs_res.txt"))$V1
    }else if(m=="HVG"){
      sgenes <- HVGfun(data)
    }else {
      sgenes = read.table(paste0("./",m,"_",d,".txt"))$V1
    }
    print(length(sgenes))
    rownames(data) = gsub("[_]","-",rownames(data))
    sgenes = gsub("[_]","-",sgenes)
    clusters <- Cluster(data = data, gene_list = sgenes)
    resmat1[d,m] <- evalcluster(label, clusters)[[1]]
    resmat2[d,m] <- evalcluster(label, clusters)[[3]]
  }
}
  
# plot

plotdf = resmat1
data_melt<-melt(plotdf)
rankmat <- matrix(nrow = 0,ncol = ncol(plotdf))
for (i in 1:nrow(plotdf)) {
  rankmat <- rbind(rankmat, rank(-plotdf[i,], ties.method='min'))
}
data_melt <- cbind(data_melt, melt(rankmat)[,3])
data_melt <- cbind(data_melt, rep(rownames(plotdf),7))
names(data_melt) = c('Method', 'ARI', 'Ranking', 'Dataset')
data_melt %>% 
  # Add a column called 'type': do we want to highlight the group or not?
  mutate(type=ifelse(Method=="CellBRF","Highlighted","Normal")) %>%
  # Build the boxplot. In the 'fill' argument, give this column
  ggplot(aes(x=Method, y=ARI, fill=type, alpha=type)) + 
  geom_boxplot() +
  # stat_summary(fun.data=MinMeanSEMMax, geom="boxplot") +
  scale_fill_manual(values=c("#69b3a2", "grey")) +
  scale_alpha_manual(values=c(1,0.1)) +
  xlab("") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 13, family = "serif", face = 'bold'),
        axis.text.y = element_text(size = 13, family = "serif"),
        axis.title.y = element_text(size = 13, family = "serif", face = 'bold'),
        legend.position = "none")
ggsave(path = "./", filename = "Fig2C.png",width = 8, height = 6, device='png', dpi=300)

data_melt %>% 
  ggplot(aes(x=Method, y=Ranking, fill=Dataset, color=Dataset)) + 
  geom_dotplot(binaxis='y', stackdir='center',  binwidth = 0.18, stackgroups = TRUE, binpositions="all") +
  # xlab("Feature Selection Method") +
  xlab("") +
  ylab("Rank") +
  scale_y_continuous(breaks = c(1:7))+
  theme(panel.grid.major = element_line(colour="gray"),
        panel.grid.minor = element_line(colour="gray"),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 13, family = "serif", face = 'bold'),
        axis.text.y = element_text(size = 13, family = "serif"),
        axis.title.y = element_text(size = 13, family = "serif", face = 'bold'),
        legend.position = "top",
        legend.key = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 13, family = "serif"))
ggsave(path = "./", filename = "Fig2B.png",width = 8, height = 6, device='png', dpi=300)


#####################################################################################
################################   neighbors comparison ############################
#####################################################################################

medlist = c("CellBRF", "DUBStepR","Feats","FEAST","geneBasisR", "HVG", "HRG")
datalist = c("Buettner","Chu_celltime", "Chung", "Darmanis", "Deng", "Engel", "Goolam", "Kim", "Koh",
             "Kolodz", "Kumar", "Leng", "Li", "Maria2", "Pollen", "Robert", "Ting", "Treutlein",
             "Usoskin", "Yan", "Yeo", "Zhou")
stadf = as.data.frame(matrix(nrow = 0, ncol = 4))
for (i in datalist) {
  print(i)
  tmp = read.table(paste0("./neighbor_auc/",i,".txt"))
  colnames(tmp) = medlist
  rownames(tmp) = c(1:15)
  tmp = tmp[-c(2,4,6,7,8,9,11,12,13,14),]
  tmp = cbind(melt(tmp), rep(c(1,3,5,10,15), length(medlist)), rep(i, nrow(tmp)))
  colnames(tmp) = c("Method", "AUC", "K_neighbor", "Dataset")
  stadf = rbind(stadf, tmp)
}

klist = c(1,3,5,10,15)
meanmat = matrix(0,nrow = 5,ncol = 7)
rownames(meanmat) = klist
colnames(meanmat) = medlist
sdmat = meanmat
for (i in klist) {
  tmp = stadf[stadf$K_neighbor==i, ]
  for (m in medlist) {
    meanmat[as.character(i), m] = mean(tmp[tmp$Method==m,]$AUC)
    sdmat[as.character(i), m] = sd(tmp[tmp$Method==m,]$AUC)
  }
}

stadf$K_neighbor = as.character(stadf$K_neighbor)

ggplot(stadf,aes(x = K_neighbor, y=AUC, fill = Method))+
  geom_boxplot() +
  ylab("K-nearest neighbor consistency") +
  xlab("Number of nearest neighbors")+
  theme_bw() +
  scale_x_discrete(limits = as.character(c(1,3,5,10,15))) +
  scale_fill_manual(values = paletteer_d("ggthemes::Tableau_10")) +
  theme(axis.title = element_text(size = 15, family = "serif", face = 'bold'),
        legend.title = element_text(size = 15, family = "serif", face = 'bold'),
        strip.text = element_text(size = 15, family = "serif", face = 'bold'),
        panel.grid =element_blank(),
        legend.position = "top",
        legend.key = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 15, family = "serif"))

ggsave(path = "./", filename = "Fig3A.png",width = 8, height = 6, device='png', dpi=300)


#### small dataset silhouette index

medlist = c("CellBRF", "DUBStepR","Feats","FEAST","geneBasisR", "HVG", "HRG")
datalist = c("Buettner","Chu_celltime", "Chung", "Darmanis", "Deng", "Engel", "Goolam", "Kim", "Koh",
             "Kolodz", "Kumar", "Leng", "Li", "Maria2", "Pollen", "Robert", "Ting", "Treutlein",
             "Usoskin", "Yan", "Yeo", "Zhou")
resmat = matrix(nrow = length(datalist), ncol = length(medlist))
colnames(resmat) = medlist
rownames(resmat) = datalist
for (d in datalist) {
  print(d)
  load(paste0("./",d,".Rdata"))
  
  l = rep(0, ncol(data))
  for (t in names(table(label))){
    l[which(label==t)] = which(names(table(label))==t)
  }
  for (m in medlist) {
    print(m)
    if(m=="CellBRF"){
      sgenes = read.table(paste0("./",d,"_CellBRF_gs_res.txt"))$V1
    }else if(m=="HVG"){
      sgenes <- HVGfun(data)
    }else {
      sgenes = read.table(paste0("./",m,"_",d,".txt"))$V1
    }
    
    if(length(intersect(rownames(data), sgenes))==0){
      rownames(data) = gsub("[_]","-",rownames(data))
      sgenes = gsub("[_]","-",sgenes)
    }  
    tmp = data[intersect(rownames(data), sgenes),]
    res <- cluster.stats(dist(t(tmp)), l)
    resmat[d,m] <- res$avg.silwidth
  }
  
}

save(resmat, file = "./si.Rdata")

load(file = "./si.Rdata")
colnames(resmat) = c("CellBRF", colnames(resmat)[-1])
data_melt <- melt(resmat)
data_melt <- data_melt[,-1]
names(data_melt) = c('Method', 'SI')
data_melt %>% 
  mutate(type=ifelse(Method=="CellBRF","Highlighted","Normal")) %>%
  ggplot(aes(x=Method, y=SI, fill=type, alpha=type)) + 
  geom_boxplot() +
  scale_fill_manual(values=c("#69b3a2", "grey")) +
  scale_alpha_manual(values=c(1,0.1)) +
  xlab("") +
  ylab("Silhouette Index") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 15, family = "serif", face = 'bold'),
        axis.text.y = element_text(size = 15, family = "serif"),
        axis.title.y = element_text(size = 15, family = "serif", face = 'bold'),
        legend.position = "none")
ggsave(path = "./", filename = "Fig3B.png",width = 8, height = 6, device='png', dpi=300)



#####################################################################################
################################   CellBRF steps analysis ###########################
#####################################################################################

# 1. data balancing
plotdf = read.table("clipboard", header = T)
data_melt<-melt(plotdf[,1:4])
names(data_melt) = c('Criteria', 'Value')
data_melt$Criteria[which(data_melt$Criteria=="NMI.1")] = "NMI"
data_melt$Criteria[which(data_melt$Criteria=="ARI.1")] = "ARI"
data_melt = cbind(data_melt, rep(rownames(plotdf), 4))
data_melt = cbind(data_melt, c(rep("CellBRF", 24),rep("Unbalanced", 24)))
names(data_melt) = c('Criteria', 'Value', 'Dataset', 'Method')
data_melt$Dataset = factor(data_melt$Dataset, levels = rownames(plotdf))

ggplot(data_melt, aes(x=Dataset, y=Value, col = Method, shape = Criteria)) + 
  geom_point(size=8, alpha = 0.5)+
  scale_color_manual(values=paletteer_d("ggsci::default_nejm")) +
  scale_shape_manual(values = c(15, 17)) +
  theme(panel.grid = element_line(color = "gray"),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 13, family = "serif", face = 'bold'),
        axis.text.y = element_text(size = 13, family = "serif"),
        axis.title.y = element_text(size = 13, family = "serif", face = 'bold'),
        axis.title.x = element_text(size = 13, family = "serif", face = 'bold'),
        legend.position = "top",
        legend.key = element_rect(fill = "white", color = NA),
        legend.title = element_text(size = 13, family = "serif"),
        legend.text = element_text(size = 13))

ggplot(data_melt, aes(x=Dataset, y=Value, col = Criteria)) + 
  geom_point(size=8, alpha = 0.5)+
  scale_color_manual(values=paletteer_d("ggsci::default_nejm")) +
  theme(panel.grid = element_line(color = "gray"),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 15, family = "serif", face = 'bold', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15, family = "serif"),
        axis.title.y = element_text(size = 15, family = "serif", face = 'bold'),
        axis.title.x = element_text(size = 15, family = "serif", face = 'bold'),
        legend.position = "top",
        legend.key = element_rect(fill = "white", color = NA),
        legend.title = element_text(size = 15, family = "serif"),
        legend.text = element_text(size = 15))
ggsave(path = "./", filename = "Fig4B.png",width = 8, height = 6, device='png', dpi=300)

df = melt(plotdf[,5])
df = cbind(df, rownames(plotdf))
colnames(df) = c("SE","Dataset")
df$Dataset = factor(df$Dataset, levels = rownames(plotdf))
ggplot(data=df, aes(x=Dataset, y=SE, group=1)) +
  geom_point(size=5, pch=18)+
  ylab("Balance Entropy (BE)") +
  theme(panel.grid = element_line(color = "gray"),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 15, family = "serif", face = 'bold', angle = 45, hjust = 1),
        axis.text.y = element_text(size = 15, family = "serif"),
        axis.title.y = element_text(size = 15, family = "serif", face = 'bold'),
        axis.title.x = element_text(size = 15, family = "serif", face = 'bold'))
ggsave(path = "./", filename = "Fig4A.png",width = 8, height = 6, device='png', dpi=300)


# 2. random forest importance

datalist<-list.files("./total/")
datalist = as.character(sapply(datalist,function(x){strsplit(x,"[.]")[[1]][1]}))
datalist = datalist[-c(1:7,16,17,26,31)]
medlist = c("FEAST","vst","dispersion", "HRG")
for (d in datalist) {
  print(d)
  load(paste0("./total/",d,".Rdata"))
  gn = rownames(data)
  cn = colnames(data)
  data = apply(data, 2, as.numeric)
  rownames(data) = gn
  colnames(data) = cn
  for (m in medlist) {
    print(m)
    if(m=="FEAST"){
      sgenes = FEASTfun(data, k=length(table(as.character(label))))
    }else if(m=="vst"){
      sgenes = HVGfun(data, n_feature=4000, method="vst")
    }else if(m=="dispersion"){
      sgenes = HVGfun(data, n_feature=4000, method="dispersion")
    }else{
      sgenes = HRGfun(data, n_feature=4000)
    }
    save(sgenes, file = paste0("./Re_gimp_test_res/",d,"_",m,".Rdata"))
  }
}

datalist = list.files("./h5data/")
datalist = as.character(sapply(datalist,function(x){strsplit(x,"[.]")[[1]][1]}))
medlist = c("FEAST","vst","dispersion", "HRG")
for (d in datalist) {
  print(d)
  x = h5read(paste0("./h5data/",d,".h5"),"X")
  y = h5read(paste0("./h5data/",d,".h5"),"Y")
  data = x
  label = y
  rownames(data) = paste0("gene",1:nrow(data))
  colnames(data) = paste0("cell",1:ncol(data))
  for (m in medlist) {
    print(m)
    if(m=="FEAST"){
      sgenes = FEASTfun(data, k=length(table(as.character(label))))
    }else if(m=="vst"){
      sgenes = HVGfun(data, n_feature=4000, method="vst")
    }else if(m=="dispersion"){
      sgenes = HVGfun(data, n_feature=4000, method="dispersion")
    }else{
      sgenes = HRGfun(data, n_feature=4000)
    }
    save(sgenes, file = paste0("./Re_gimp_test_res/",d,"_",m,".Rdata"))
  }
}

gslist = c(20, 30, 50, 100, 200, 400, 700, 1200, 2000, 4000)
medlist = c("CellBRF", "FEAST","vst","dispersion","HRG")
datalist<-list.files("./total/")
datalist = as.character(sapply(datalist,function(x){strsplit(x,"[.]")[[1]][1]}))
datalist = datalist[-c(1:7,16,17,26,31)]
for (d in datalist) {
  resmat1 = matrix(nrow = length(medlist), ncol = length(gslist))
  colnames(resmat1) = gslist
  rownames(resmat1) = medlist
  resmat2 = resmat1
  print(d)
  load(paste0("./",d,".Rdata"))
  for (m in medlist) {
    print(m)
    if(m=="CellBRF"){
      g_imp = read.table(paste0("./",d,"_CellBRF_gene_imp_res.txt"))$V1
      g_name = read.table(paste0("./",d,"_CellBRF_genenames.txt"))$V1
      names(g_imp) = g_name
      sgenes = names(sort(g_imp, decreasing = T))
    }else{
      load(file = paste0("./",d,"_",m,".Rdata"))
    }
    rownames(data) = gsub("[_]","-",rownames(data))
    sgenes = gsub("[_]","-",sgenes)
    
    for (g in gslist) {
      print(g)
      clusters1 <- Cluster(data = data, gene_list = sgenes[1:g])
      resmat1[m,as.character(g)] <- evalcluster(label, clusters1)[[1]]
      resmat2[m,as.character(g)] <- evalcluster(label, clusters1)[[3]]
    }
  }
}

datalist = list.files("./h5data/")
datalist = as.character(sapply(datalist,function(x){strsplit(x,"[.]")[[1]][1]}))
for (d in datalist) {
  print(d)
  resmat1 = matrix(nrow = length(medlist), ncol = length(gslist))
  colnames(resmat1) = gslist
  rownames(resmat1) = medlist
  resmat2 = resmat1
  x = h5read(paste0("./h5data/",d,".h5"),"X")
  y = h5read(paste0("./h5data/",d,".h5"),"Y")
  data = x
  label = y
  rownames(data) = paste0("gene",1:nrow(data))
  colnames(data) = paste0("cell",1:ncol(data))
  
  label = as.character(label)
  l = rep(0, ncol(data))
  for (t in names(table(label))){
    l[which(label==t)] = which(names(table(label))==t)
  }
  
  for (m in medlist) {
    print(m)
    if(m=="CellBRF"){
      g_imp = read.table(paste0("./",d,"_CellBRF_gene_imp_res.txt"))$V1
      g_name = read.table(paste0("./",d,"_CellBRF_genenames.txt"))$V1
      names(g_imp) = g_name
      sgenes = names(sort(g_imp, decreasing = T))
      sgenes = rownames(data)[as.numeric(sgenes)+1]
    }else{
      load(file = paste0("./Re_gimp_test_res/",d,"_",m,".Rdata"))
    }
    rownames(data) = gsub("[_]","-",rownames(data))
    sgenes = gsub("[_]","-",sgenes)
    
    for (g in gslist) {
      print(g)
      clusters1 <- Cluster(data = data, gene_list = sgenes[1:g])
      resmat1[m,as.character(g)] <- evalcluster(label, clusters1)[[1]]
      resmat2[m,as.character(g)] <- evalcluster(label, clusters1)[[3]]
    }
  }
}

datalist<-list.files("./")
datalist = unique(as.character(sapply(datalist,function(x){strsplit(x,"_A")[[1]][1]})))
datalist = unique(as.character(sapply(datalist,function(x){strsplit(x,"_N")[[1]][1]})))
ARIsat = data.frame()
for (d in datalist) {
  print(d)
  ARIres = read.csv(paste0("./",d,"_ARI.csv"))
  rownames(ARIres) = c("CellBRF", "FEAST", "HVGvst", "HVGdisp", "HRG")
  ARIres = ARIres[,-1]
  colnames(ARIres) = c(20, 30, 50, 100, 200, 400, 700, 1200, 2000, 4000)
  tmp1 = melt(ARIres)
  tmp1 = cbind(rep(rownames(ARIres),length(colnames(ARIres))), tmp1)
  colnames(tmp1) = c("Method","gsn","ARI")
  
  ARIsat = rbind(ARIsat, tmp1)
}
gslist = c(20, 30, 50, 100, 200, 400, 700, 1200, 2000, 4000)
medlist = c("CellBRF", "FEAST", "HVGvst", "HVGdisp", "HRG")
ARIfinal = data.frame()
for (i in medlist) {
  for (j in gslist) {
    idx = intersect(which(ARIsat$Method==i), which(ARIsat$gsn==j))
    ARIfinal = rbind(ARIfinal, c(i,j,mean(ARIsat$ARI[idx])))
  }
}
colnames(ARIfinal) = c("Method","gsn","ARI")

ARIfinal = read.table("clipboard", header = T)
ARIfinal$ARI = as.numeric(ARIfinal$ARI)
ARIfinal$gsn = factor(ARIfinal$gsn, levels = gslist)

# mean(gsnum)
ggplot(data=ARIfinal, aes(x=gsn, y=ARI, color=Method, group=Method)) +
  geom_line(linetype = "dashed")+
  geom_point(size=3, pch=17)+
  scale_color_manual(values = paletteer_d("ggsci::default_nejm")) +
  xlab("Number of selected features") +
  ylab("ARI") +
  theme(panel.grid = element_line(color = "gray"),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 15, family = "serif"),
        axis.text.y = element_text(size = 15, family = "serif"),
        axis.title.y = element_text(size = 15, family = "serif", face = 'bold'),
        axis.title.x = element_text(size = 15, family = "serif", face = 'bold'),
        legend.position = "top",
        legend.key = element_rect(fill = "white", color = NA),
        legend.text = element_text(size = 15, family = "serif"))
ggsave(path = "./", filename = "Fig5A.png",width = 10, height = 5, device='png', dpi=300)

# 3. gene set optimal
datalist<-list.files("./")
datalist = unique(as.character(sapply(datalist,function(x){strsplit(x,"_A")[[1]][1]})))
datalist = unique(as.character(sapply(datalist,function(x){strsplit(x,"_N")[[1]][1]})))
ARImat = data.frame()
gn = c()
for (d in datalist) {
  print(d)
  sgenes = read.table(paste0("./",d,"_CellBRF_gs_res.txt"))$V1
  gn = c(gn, length(sgenes))
  ARIres = read.csv(paste0("./",d,"_ARI.csv"))
  rownames(ARIres) = c("CellBRF", "FEAST", "HVGvst", "HVGdisp", "HRG")
  ARIres = ARIres[,-1]
  colnames(ARIres) = c(20, 30, 50, 100, 200, 400, 700, 1200, 2000, 4000)
  ARImat = rbind(ARImat, ARIres[1,])
}
rownames(ARImat) = datalist

ARImat = read.table("clipboard",header = T)
MinMeanSEMMax <- function(x) {
  v <- c(min(x), quantile(x, 0.25), mean(x), quantile(x, 0.75), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  v
}
colnames(ARImat) = c(c(20, 30, 50, 100, 200, 400, 700, 1200, 2000, 4000),"3σ_rule based")
df = melt(ARImat)
colnames(df) = c("gsn", "ARI")
df$gsn = factor(df$gsn, levels = c("3σ_rule based", c(20, 30, 50, 100, 200, 400, 700, 1200, 2000, 4000)))
df %>% 
  mutate(type=ifelse(gsn=="3σ_rule based","Highlighted","Normal")) %>%
  ggplot(aes(x=gsn, y=ARI, fill=type, alpha=type)) + 
  # geom_boxplot() +
  stat_summary(fun.data=MinMeanSEMMax, geom="boxplot") +
  # ggplot(aes(x=gsn, y=Rank, fill=Dataset, color=Dataset)) + 
  # geom_dotplot(binaxis='y', stackdir='center',  binwidth = 0.2, stackgroups = TRUE, binpositions="all") +
  xlab("Number of selected features") +
  ylab("ARI") +
  scale_fill_manual(values=c("#69b3a2", "grey")) +
  scale_alpha_manual(values=c(1,0.1)) +
  # scale_y_continuous(breaks = c(1:13))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA),
        axis.text.x = element_text(size = 15, family = "serif"),
        axis.text.y = element_text(size = 15, family = "serif"),
        axis.title.y = element_text(size = 15, family = "serif", face = 'bold'),
        axis.title.x = element_text(size = 15, family = "serif", face = 'bold'),
        legend.position = "none")
ggsave(path = "./", filename = "Fig5B.png",width = 10, height = 5, device='png', dpi=300)


#####################################################################################
################################   time course data analysis ########################
#####################################################################################

medlist = c("CellBRF", "DUBStepR","Feats","FEAST","geneBasisR", "HVG", "HRG")
load(paste0("./Chu_celltime.Rdata"))
data=log2(data+1)
sgenes = read.table(paste0("./Chu_celltime_CellBRF_gs_res.txt"))$V1
# sgenes = read.table(paste0("./Feats_Chu_celltime.txt"))$V1
pbmc <- RunTSNE(t(as.matrix(data[intersect(sgenes, rownames(data)),])))
tsne <- pbmc@cell.embeddings
tsne <- as.data.frame(tsne)
tsne <- cbind(tsne, label)
colnames(tsne) <- c("tSNE_1", "tSNE_2", "type")
ggplot(tsne, aes(x = tSNE_1, y = tSNE_2, color = type)) +
  geom_point(size=1) +
  scale_color_manual(values = paletteer_d("ggthemes::Superfishel_Stone")) +
  ylab("tSNE_2") +
  xlab("tSNE_1")+
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(size = 15, family = "serif"),
        axis.text.y = element_text(size = 15, family = "serif"),
        axis.title.y = element_text(size = 15, family = "serif", face = 'bold'),
        axis.title.x = element_text(size = 15, family = "serif", face = 'bold'),
        legend.text = element_text(size = 15, family = "serif"),
        panel.grid =element_blank())
ggsave(path = "./", filename = "Fig6_CBRF.png",width = 6, height = 6, device='png', dpi=300)

load(paste0("./Chu_celltime.Rdata"))
load(paste0("./Chu_DEgenes.Rdata"))
allrefgene = intersect(refgenelist[,1],rownames(data))
mostpatternlist = names(table(refgenelist[,2])[which(table(refgenelist[,2])>100)])
centralpatternlist = names(table(refgenelist[,2])[which(table(refgenelist[,2])>10 & table(refgenelist[,2])<=100)])
rarepatternlist = names(table(refgenelist[,2])[which(table(refgenelist[,2])<=10)])

mostlist = c()
for (p in mostpatternlist) {
  mostlist = c(mostlist, refgenelist[which(refgenelist[,2]==p),1])
}
centrallist = c()
for (p in centralpatternlist) {
  centrallist = c(centrallist, refgenelist[which(refgenelist[,2]==p),1])
}
rarelist = c()
for (p in rarepatternlist) {
  rarelist = c(rarelist, refgenelist[which(refgenelist[,2]==p),1])
}

overlapRatioMat = matrix(nrow = length(medlist), ncol = 3)
colnames(overlapRatioMat) = c("Central","Rare", "Most")
rownames(overlapRatioMat) = medlist
for (m in medlist) {
  if(m=="CellBRF"){
    sgenes = read.table(paste0("./Chu_celltime_CellBRF_gs_res.txt"))$V1
  }else {
    sgenes = read.table(paste0("./",m,"_Chu_celltime.txt"))$V1
  }
  overlapRatioMat[m,"Most"] = length(intersect(sgenes,mostlist))/length(sgenes)
  overlapRatioMat[m,"Central"] = length(intersect(sgenes,centrallist))/length(sgenes)
  overlapRatioMat[m,"Rare"] = length(intersect(sgenes,rarelist))/length(sgenes)
}
set.seed(2022)
df = melt(overlapRatioMat)
colnames(df) = c("Method","Group","Number")
df$Group = factor(df$Group, levels = c("Most","Central", "Rare"))
ggplot(df, aes(x = Method, y = Number, group = Group)) + 
  geom_col(aes(fill=Group))  +
  scale_fill_manual(values = paletteer_d("ggthemes::Superfishel_Stone"),
                    labels = c("Most","Rare","Extremely Rare")) +
  ylim(0,0.9)+
  ylab("Overlap") +
  theme_bw() +
  theme(
    text=element_text(size=15,  family="serif"),
    axis.text = element_text(size = 15, family = "serif", face = 'bold'),
    axis.title.y = element_text(size = 15, family = "serif", face = 'bold'),
    legend.text = element_text(size = 15, family = "serif"),
    legend.title = element_text(size = 15, family = "serif", face = 'bold'),
    axis.text.x = element_text(size = 15, family = "serif", face = 'bold', angle = 45, hjust = 1),
    axis.title.x = element_blank(),
    strip.text = element_text(size=10),
    panel.grid =element_blank(),
    legend.position = c(.99, .99),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6))
ggsave(path = "./", filename = "Fig6_overlap.png",width = 6, height = 6, device='png', dpi=300)

load(paste0("./Chu_celltime.Rdata"))
plotdf = matrix(nrow = 5, ncol = 2)
clusters1 = Cluster(data = data, gene_list = rownames(data))
plotdf[1,1] = evalcluster(label, clusters1)[[1]]
plotdf[1,2] = evalcluster(label, clusters1)[[3]]
clusters1 = Cluster(data = data, gene_list = allrefgene)
plotdf[2,1] = evalcluster(label, clusters1)[[1]]
plotdf[2,2] = evalcluster(label, clusters1)[[3]]
clusters1 = Cluster(data = data, gene_list = mostlist)
plotdf[3,1] = evalcluster(label, clusters1)[[1]]
plotdf[3,2] = evalcluster(label, clusters1)[[3]]
clusters1 = Cluster(data = data, gene_list = centrallist)
plotdf[4,1] = evalcluster(label, clusters1)[[1]]
plotdf[4,2] = evalcluster(label, clusters1)[[3]]
clusters1 = Cluster(data = data, gene_list = rarelist)
plotdf[5,1] = evalcluster(label, clusters1)[[1]]
plotdf[5,2] = evalcluster(label, clusters1)[[3]]

colnames(plotdf) = c("NMI", "ARI")
rownames(plotdf) = c("ALL", "All Pattern", "Main Pattern", "Rare Pattern", "Extremely Rare Pattern")
plotdf = read.table("clipboard", sep = "\t", header = T)
colnames(plotdf) = c("Genes","NMI", "ARI")
plotdf = plotdf[-1,]
plotdf$Genes = c("ALL", "Main", "Rare", "Extremely Rare")
df = melt(plotdf)
colnames(df) = c("Genes", "Criteria", "Value")
df$Genes = factor(df$Genes, levels = c("ALL", "Main", "Rare", "Extremely Rare"))
ggplot(df, aes(y = Value, x = Genes, fill=Criteria)) + 
  geom_bar(stat = "identity", position = position_dodge())  +
  scale_fill_manual(values=paletteer_d("ggthemes::Superfishel_Stone")) +
  xlab("Gene expression pattern") +
  coord_cartesian(ylim=c(0.6, 0.9)) +
  theme_bw() +
  theme(text=element_text(size=10,  family="serif",face = 'bold'),
        panel.grid =element_blank(),
        axis.title.y = element_blank(),
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        axis.text.x = element_text(size = 15, family = "serif"),
        axis.text.y = element_text(size = 15, family = "serif"),
        axis.title.x = element_text(size = 15, family = "serif", face = 'bold'),
        legend.text = element_text(size = 15, family = "serif"),
        legend.title = element_text(size = 15, family = "serif"))

ggsave(path = "./", filename = "Fig6_clustering.png",width = 6, height = 6, device='png', dpi=300)

medlist = c("CellBRF", "DUBStepR","Feats","FEAST","geneBasisR", "HVG", "HRG")
load(paste0("./Chu_celltime.Rdata"))
rownames(data) = gsub("[_]","-",rownames(data))
SeuObj = CreateSeuratObject(counts=data)
SeuObj <- NormalizeData(SeuObj,verbose = FALSE)
all.gene = rownames(SeuObj)
SeuObj <- ScaleData(SeuObj, features = all.gene, verbose = FALSE)

resmat1 = matrix(nrow = length(medlist), ncol = 5)
rownames(resmat1) = medlist
colnames(resmat1) = as.character(c(1,3,5,10,15))
resmat2 = matrix(nrow = length(medlist), ncol = 1)
rownames(resmat2) = medlist
colnames(resmat2) = c("silhouette")
l = rep(0, ncol(data))
for (t in names(table(label))){
  l[which(label==t)] = which(names(table(label))==t)
}

testidx = sample(ncol(data), floor(ncol(data)*0.2), replace = F)
trainidx = setdiff(1:ncol(data), testidx)

for (m in medlist) {
  if(m=="CellBRF"){
    sgenes = read.table(paste0("./Chu_celltime_CellBRF_gs_res.txt"))$V1
  }else{
    sgenes = read.table(paste0("./",m,"_Chu_celltime.txt"))$V1
  }
  sgenes = gsub("[_]","-",sgenes)
  pbmc <- RunTSNE(object = SeuObj,seed.use=2022, features = sgenes, dims=NULL)
  umap <- pbmc@reductions$tsne@cell.embeddings
  umap <- as.data.frame(umap)
  umap <- cbind(umap, label)
  
  train_cl <- umap[trainidx,]
  test_cl <- umap[testidx,]
  train_scale <- scale(train_cl[, 1:2])
  test_scale <- scale(test_cl[, 1:2])
  
  for (k in c(1,3,5,10,15)) {
    classifier_knn <- knn(train = train_scale,
                          test = test_scale,
                          cl = train_cl$label,
                          k = k)
    cm <- table(test_cl$label, classifier_knn)
    misClassError <- mean(classifier_knn != test_cl$label)
    print(paste('Accuracy =', 1-misClassError))
    resmat1[m,as.character(k)] = 1-misClassError
  }
  
  res <- cluster.stats(dist(umap), l)
  resmat2[m,1] <- res$avg.silwidth
}

###############################################################################################################

load(paste0("./Puram.Rdata"))
gn = rownames(data)
cn = colnames(data)
data = apply(data, 2, as.numeric)
rownames(data) = gn
colnames(data) = cn
data=log2(data+1)
sgenes = read.table(paste0("./Puram_CellBRF_gs_res.txt"))$V1
# sgenes = read.table(paste0("./Feats_Puram.txt"))$V1
pbmc <- RunTSNE(t(as.matrix(data[sgenes,])))
umap <- pbmc@cell.embeddings
umap <- as.data.frame(umap)
umap <- cbind(umap, label)
colnames(umap) <- c("UMAP_1", "UMAP_2", "type")
ggplot(umap, aes(x = UMAP_1, y = UMAP_2, color = type)) +
  geom_point(size=0.5) +
  scale_color_manual(values = paletteer_d("ggthemes::Superfishel_Stone")) +
  ylab("tSNE_2") +
  xlab("tSNE_1")+
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme(legend.title = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(size = 15, family = "serif"),
        axis.text.y = element_text(size = 15, family = "serif"),
        axis.title.y = element_text(size = 15, family = "serif", face = 'bold'),
        axis.title.x = element_text(size = 15, family = "serif", face = 'bold'),
        legend.text = element_text(size = 15, family = "serif"),
        panel.grid =element_blank())
rm(pbmc)
ggsave(path = "./", filename = "Fig7_CBRF.png",width = 6, height = 6, device='png', dpi=300)

refmat = read.table("./PanglaoDB_markers_27_Mar_2020.tsv",
                    sep = "\t", quote="")
colnames(refmat) = refmat[1,]
refmat = refmat[-1,]

type = c("B cells", "Dendritic cells", "Endothelial cells","Fibroblasts","Macrophages","Mast cells","Myocytes","T cells")
refgenes = c()
for (t in type) {
  print(t)
  tmp = intersect(rownames(data),refmat[which(refmat$`cell type`==t),"official gene symbol"])
  # tmp = intersect(derefgenes, tmp)
  print(length(tmp))
  refgenes = c(refgenes, list(tmp))
}
type = c("B cell", "Dendritic", "Endothelial","Fibroblast","Macrophage","Mast","Myocyte","T cell")
names(refgenes) = type

overlapRatioMat = matrix(nrow = length(medlist), ncol = length(refgenes))
colnames(overlapRatioMat) = type
rownames(overlapRatioMat) = medlist
for (m in medlist) {
  if(m=="CellBRF"){
    sgenes = read.table(paste0("./Puram_CellBRF_gs_res.txt"))$V1
  }else {
    sgenes = read.table(paste0("./",m,"_Puram.txt"))$V1
  }
  for (i in 1:length(refgenes)) {
    overlapRatioMat[m,names(refgenes)[i]] = length(intersect(sgenes,refgenes[[i]]))/length(intersect(sgenes, rownames(data)))
  }
}

df = melt(overlapRatioMat)
colnames(df) = c("Method","Group","Number")
ggplot(df, aes(x = Method, y = Number, group = Group)) +
  geom_col(aes(fill=Group))  +
  # ggplot(df, aes(x = Method, y = Number)) +
  # geom_bar(stat="identity", fill="#6F99ADFF", width=0.75)  +
  scale_fill_manual(values = paletteer_d("ggthemes::Superfishel_Stone")) +
  ylab("Overlap") +
  xlab("Methods")+
  theme_bw() +
  theme(text=element_text(size=15,  family="serif"),
        axis.text = element_text(size = 15, family = "serif", face = 'bold'),
        axis.title.y = element_text(size = 15, family = "serif", face = 'bold'),
        legend.text = element_text(size = 15, family = "serif"),
        axis.text.x = element_text(size = 15, family = "serif", face = 'bold', angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "right",
        panel.grid =element_blank())
ggsave(path = "D:/CSU/cell_forest/NC/plot_ismb/", filename = "Fig7_overlap.png",width = 6, height = 6, device='png', dpi=300)

sgenes = read.table(paste0("./Puram_CellBRF_gs_res.txt"))$V1

df = apply(data[sgenes,which(label=="Fibroblast")], 1, as.numeric)
sdlist = apply(df, 2, sd)
col = as.character(paletteer_c("ggthemes::Blue", 12))
heatmap(t(df[,order(sdlist, decreasing = T)[1:30]]), scale = "none", 
        col =  col, labCol = NA, keep.dendro = F)

df = apply(data[sgenes,which(label=="Fibroblast")], 1, as.numeric)
sdlist = apply(df, 2, sd)
df = t(df[,order(sdlist, decreasing = T)[1:30]])
colnames(df) = which(label=="Fibroblast")
hc<-hclust(dist(df))
rowInd<-hc$order
hc<-hclust(dist(t(df)))
colInd<-hc$order
df.m<-df[rowInd,colInd]
df.m<-apply(df.m,1,function(x){(x-min(x))/(max(x)-min(x))})
df.m<-t(df.m)

subdf.m<-matrix(ncol = 0, nrow = nrow(df.m))
rownames(subdf.m)<-rownames(df.m)
for (i in seq(1,ncol(df.m),50)) {
  subdf.m<-cbind(subdf.m, df.m[,i])
}
colnames(subdf.m)<-colnames(df.m)[seq(1,ncol(df.m),50)]
df.m<-subdf.m

coln<-colnames(df.m)
rown<-rownames(df.m)
colnames(df.m)<-1:ncol(df.m)
rownames(df.m)<-1:nrow(df.m)
df.m<-melt(df.m)
base_size<-12
ggplot(df.m, aes(Var2, Var1)) +
  geom_tile(aes(fill = value), 
            color = "black",
            lwd = 0.5,
            linetype = 1) + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  scale_x_continuous(expand = c(0, 0),labels=coln,breaks=1:length(coln)) + 
  scale_y_continuous(expand = c(0, 0),labels=rown,breaks=1:length(rown)) +
  theme(text=element_text(size=15,  family="serif"),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_text(size = 15, family = "serif", face = 'bold'),
        # legend.title = element_blank())
        legend.position="none")
ggsave(path = "./", filename = "Fig7_heatmap.png",width =6, height = 6, device='png', dpi=300)

SeuObj = CreateSeuratObject(counts=data)
SeuObj <- NormalizeData(SeuObj,verbose = FALSE)
all.gene = rownames(SeuObj)
SeuObj <- ScaleData(SeuObj, features = all.gene, verbose = FALSE)

medlist = c("CellBRF", "DUBStepR","Feats","FEAST","geneBasisR", "HVG", "HRG")
resmat1 = matrix(nrow = length(medlist), ncol = 5)
rownames(resmat1) = medlist
colnames(resmat1) = as.character(c(1,3,5,10,15))
resmat2 = matrix(nrow = length(medlist), ncol = 1)
rownames(resmat2) = medlist
colnames(resmat2) = c("silhouette")
l = rep(0, ncol(data))
for (t in names(table(label))){
  l[which(label==t)] = which(names(table(label))==t)
}

testidx = sample(ncol(data), floor(ncol(data)*0.2), replace = F)
trainidx = setdiff(1:ncol(data), testidx)

for (m in medlist) {
  if(m=="CellBRF"){
    sgenes = read.table(paste0("./Puram_ClustFS_gs_res.txt"))$V1
  }else{
    sgenes = read.table(paste0("./",m,"_Puram.txt"))$V1
  }
  sgenes = gsub("[_]","-",sgenes)
  pbmc <- RunTSNE(object = SeuObj,seed.use=2022, features = sgenes, dims=NULL)
  umap <- pbmc@reductions$tsne@cell.embeddings
  umap <- as.data.frame(umap)
  umap <- cbind(umap, label)
  
  train_cl <- umap[trainidx,]
  test_cl <- umap[testidx,]
  train_scale <- scale(train_cl[, 1:2])
  test_scale <- scale(test_cl[, 1:2])
  
  for (k in c(1,3,5,10,15)) {
    classifier_knn <- knn(train = train_scale,
                          test = test_scale,
                          cl = train_cl$label,
                          k = k)
    cm <- table(test_cl$label, classifier_knn)
    misClassError <- mean(classifier_knn != test_cl$label)
    print(paste('Accuracy =', 1-misClassError))
    resmat1[m,as.character(k)] = 1-misClassError
  }
  
  res <- cluster.stats(dist(umap), l)
  resmat2[m,1] <- res$avg.silwidth
}

###################  rare cell analysis  ####################

counts <- ReadMtx(mtx = "./jurkat/matrix.mtx",
                  cells = "./jurkat/barcodes.tsv",
                  features = "./jurkat/genes.tsv")
annotation <- read.table("clipboard")

names = rownames(annotation)[c(which(annotation$SNV=="species1" | annotation$SNV=="species2"))]
annotation = annotation[names,1]
names(annotation) = names
counts = counts[,names]

annotation[which(annotation=="species1")] <- "293T"
annotation[which(annotation=="species2")] <- "jurkat"

save(counts, annotation, file = "./jurkat/jurkat_raw.Rdata")

rare_w = c(0.005, 0.01, 0.015, 0.02, 0.025, 0.05)
set.seed(2022)
for (r in 1:length(rare_w)) {
  id = c()
  id = c(id, sample(which(annotation=="jurkat"), floor(rare_w[r]*length(which(annotation=="jurkat")))))
  id = c(id, which(annotation=="293T"))
  data = as.matrix(counts[,id])
  label = annotation[id]
  save(data, label, file = paste0("./jurkat_sub_datasets/jurkat_sub",r,".Rdata"))
}

datalist = list.files("./jurkat_sub_datasets/")
for (i in 1:length(datalist)) {
  # print(i)
  dn = strsplit(datalist[i],"[.]")[[1]][1]
  load(paste0("./jurkat_sub_datasets/",datalist[i]))
  write.csv(data, file = paste0("./",dn,".csv"))
  write.csv(label, file = paste0("./",dn,"_label.csv"))
}

datalist = list.files("./jurkat_sub_datasets/")
medlist = c("DUBStepR","FEAST","geneBasisR","HVG", "HRG")
for (d in datalist) {
  print(d)
  d = strsplit(d,"[.]")[[1]][1]
  load(paste0("./jurkat_sub_datasets/",d,".Rdata"))
  gn = rownames(data)
  cn = colnames(data)
  data = apply(data, 2, as.numeric)
  rownames(data) = gn
  colnames(data) = cn
  for (m in medlist) {
    print(m)
    if(m=="DUBStepR"){
      x=log2(data+1)
      kvar = which(apply(x,1,sd)>0)
      x=x[kvar,]
      sgenes = DUBStepRfun(data[rownames(x),])
    }else if(m=="FEAST"){
      sgenes = FEASTfun(data,k = length(table(label)),n_feature = 2000)
    }else if(m=="geneBasisR"){
      sgenes = geneBasisRfun(data)
    }else if(m=="HRG"){
      sgenes = HRGfun(data)
    }else{
      sgenes = HVGfun(data)
    }
    write.table(sgenes,
                file = paste0("./",m,"_", d,".txt"),
                row.names = FALSE, col.names = FALSE, quote=FALSE)
  }
}


datalist = list.files("./jurkat_sub_datasets/")
medlist = c("CellBRF", "DUBStepR","Feats","FEAST","geneBasisR", "HVG", "HRG", "Raw")
for (i in 2:length(datalist)) {
  print(i)
  d = strsplit(datalist[i],"[.]")[[1]][1]
  load(file = paste0("./jurkat_sub_res/",d,"_res.Rdata"))
  
  rownames(stamat)[8] = "original FiRE"
  stamat[which(is.na(stamat))] = 0
  stamat = cbind(rep(0, 8), stamat)
  colnames(stamat)[1] = "0"
  pd = melt(stamat)
  colnames(pd) = c("Method", "selected_N", "rare_N")
  
  ggplot(pd, aes(x = selected_N, y = rare_N, color = Method, fill=Method)) +
    geom_line(size = 2, alpha=0.8) +
    scale_color_manual(values = paletteer_d("ggsci::default_nejm")) +
    ylab("Numbers of rare (Jurkat) cells") +
    xlab("Number of selected rare cells")+
    theme_bw() +
    guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(text=element_text(size=15,  family="serif",face = 'bold'),
          strip.text = element_text(size=15, family="serif",face = 'bold'),
          panel.grid =element_blank(),
          legend.position = c(.01, .99),
          legend.justification = c("left", "top"),
          legend.box.just = "left",
          legend.margin = margin(0, 0, 0, 0),
          legend.key.width=unit(0.1,"cm"),
          legend.key.height=unit(0.1,"cm"))
  ggsave(path = "./", filename = paste0("Fig8_",i,".png"),width = 10, height = 6, device='png', dpi=300)
}


load("./jurkat_sub5.Rdata")
sgenes = read.table(paste0("./jurkat_sub5_CellBRF_gs_res.txt"))$V1
pbmc = CreateSeuratObject(counts=data)
pbmc <- NormalizeData(pbmc,verbose = FALSE)
df = pbmc@assays$RNA@data[sgenes,order(label)]
df_lab = label[order(label)]

col <- paletteer_c("ggthemes::Blue", 256)
col_col <- as.character(df_lab)
col_col[which(col_col=="293T")] = "#F59885"
col_col[which(col_col=="jurkat")] = "#D44D44"
gn = rownames(df)
cn = colnames(df)
df = apply(df, 2, as.numeric)
rownames(df) = gn
colnames(df) = cn

png("./Fig8_heatmap.png",width=9.5,height=6,units="in",res=300, family = "serif")
heatmap(df, scale = "none", col =  col, ColSideColors = col_col, Colv = NA, labCol = NA,
        xlab = NA, ylab = NA, margins = c(1,6), cexRow = 0.3)
dev.off()


df<-melt(df)
base_size<-12
ggplot(df, aes(Var2, Var1)) +
  geom_tile(aes(fill = value), 
            color = "black",
            lwd = 0.5,
            linetype = 1) + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_grey(base_size = base_size) + labs(x = "", y = "") + 
  theme(text=element_text(family="serif"),
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        # legend.title = element_blank())
        legend.position="none")

###################  run time analysis  ####################

load("./68KPBMCsdata.Rdata")
set.seed(2022)

nsamples = c(100,200,500,1000,2000,5000,10000,20000)
for (n in nsamples) {
  type_num = ceiling(table(label)/sum(table(label))*n)
  idx = c()
  for (t in names(table(label))) {
    idx = c(idx, sample(which(label==t),type_num[t]))
  }
  subdata = data[idx,]
  sublab = label[idx]
  save(subdata,sublab,file = paste0("./time_test_datalist/",n,".Rdata"))
  write.csv(subdata,file = paste0("./time_test_datalist/",n,".csv"))
  write.csv(sublab,file = paste0("./time_test_datalist/",n,"_label.csv"))
}

datalist = c("100","200","500","1000","2000","5000","10000","20000")
# medlist = c("DUBStepR","FEAST","geneBasisR","HVG")
medlist = c("HRG")
runtimeMat = matrix(nrow = length(medlist),ncol = length(datalist))
colnames(runtimeMat) = datalist
rownames(runtimeMat) = medlist
for (d in datalist) {
  load(paste0("./time_test_datalist/",d,".Rdata"))
  subdata = t(subdata)
  
  for (m in medlist) {
    print(m)
    timestart <- proc.time()
    if(m=="DUBStepR"){
      res <- DUBStepR(input.data = subdata)
    }else if(m=="FEAST"){
      Y = process_Y(subdata, thre = 2)
      ixs  = FEAST_fast(Y, k=length(unique(as.character(sublab))))
    }else if(m=="geneBasisR"){
      data = as.data.frame(as.matrix(subdata))
      sce=SingleCellExperiment(assays = list(counts = data, logcounts = log2(data+1)))
      sce = retain_informative_genes(sce)
      # run gene selection
      genes = gene_search(sce, n_genes_total = 50)
    }else{
      genes = HRGfun(subdata)
    }
    runningtime <- proc.time() - timestart
    runtimeMat[m, as.character(d)] <- runningtime[[3]]
  }
}

save(runtimeMat, file = "./HRG_runtime.Rdata")

# plot
df = read.table("clipboard")
# df = log2(df+1)

png("./Fig9.png",width=8,height=6,units="in",res=300)
plot(as.numeric(df[1,])~c(1:8) , type="b" , bty="l" , 
     xlab="Number of cells" , ylab="Runtime (secs)" , col="#BC3C29", 
     lwd=3 , pch=17 , ylim=c(0,800), family = "serif", xaxt = "n")
axis(1, at=1:8, labels=c("100", "200", "500", "1000", "2000", "5000", "10000", "20000"))
lines(as.numeric(df[2,])~c(1:8), col="#0072B5" , lwd=3 , pch=15 , type="b" )
lines(as.numeric(df[3,])~c(1:8), col="#E18727" , lwd=3 , pch=16 , type="b" )
lines(as.numeric(df[4,])~c(1:8), col="#20854E" , lwd=3 , pch=18 , type="b" )
lines(as.numeric(df[5,])~c(1:8), col="#7876B1" , lwd=3 , pch=0 , type="b" )
lines(as.numeric(df[6,])~c(1:8), col="#6F99AD" , lwd=3 , pch=1 , type="b" )
lines(as.numeric(df[7,])~c(1:8), col="#B07AA1" , lwd=3 , pch=2 , type="b" )
legend("topleft", 
       legend = c("CellBRF", "DUBStepR","Feats","FEAST","geneBasisR", "HVG", "HRG"),
       col = c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD", "#B07AA1"), 
       pch = c(17,15,16,18, 0, 1, 2), 
       bty = "n",
       pt.cex = 1, 
       cex = 0.8, 
       text.col = "black", 
       horiz = F , 
       inset = c(0, 0))
dev.off()


