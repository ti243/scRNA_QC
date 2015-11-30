#Author: Tomislav Ilicic
#Date: 2/07/15
#Organisation: WTSI
#Description: SCRIPT TO EXTRACT BIOLOGICAL AND TECHNICAL FEATURES FORM HT-SEQ COUNT DATA. 

library(getopt)
library(topGO)
library(DESeq)
library(biomaRt)
library(ggplot2)
library(e1071)
library(mvoutlier)
library(grid)

#EXTRACT FEATURES FOR DATASET
extract_features<-function(read_metrics, counts, genes, file_name, output_dir, common_features_input, GO_terms_input, extra_genes_input, organism) {
  info ("Extract features")
  info(paste0("Gene expression: ", input))
  info(paste0("Features output dir: ", output_dir))
  info(paste0("GO Terms to use: ", GO_terms_input))
  info(paste0("Common Features: ", common_features_input))
  info(paste0("Extra genes Features: ", extra_genes_input))
  
  
  output_dir= paste(output_dir, file_name ,sep="/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  #LOAD ALL DATA
  GO_terms<-read.table(GO_terms_input, header=FALSE, stringsAsFactors=FALSE)
  extra_genes=read.table(extra_genes_input, header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
  common_features=unlist(read.table(common_features_input, header=FALSE, stringsAsFactors=FALSE, fill=TRUE))
  extra_genes=apply(extra_genes, 1, function(x){
    return(subset(x, x != ""))
  })
  
  info(paste0("DATA Loaded."))
  
  
  cell_names=colnames(counts)
  cell_names<-cell_names[-1]
  rownames(counts)<-genes
  
  #GENERATE ALL FEATURES       
  features_all<-feature_generation(read_metrics, counts, GO_terms, extra_genes, organism)
  info(paste0("FEATURES generated."))
  sds<-apply(features_all, 2, sd)
  #REMOVE 0-VARIANCE VALUE FEATURES
  features_all<-features_all[,sds != 0]
  types<-c("all", "common")
  
  features_common=features_all[,which(colnames(features_all) %in% common_features)]
  f_all=paste0(output_dir, "/", file_name, ".", types[1], ".features")
  f_common=paste0(output_dir, "/", file_name, ".", types[2], ".features")
  write.table(features_common, f_common)
  write.table(features_all, f_all)
  
  info(paste0("Features saved: ", f_all))
  info(paste0("Features saved: ", f_common))
  
  return(list(features_all, features_common))
  
}


#GENERATE ALL FEATURES
feature_generation<-function(read_metrics, counts_nm, GO_terms, extra_genes, organism) {
  
  features<-list()
  
  #REMOVE ALL 0 GENES
  counts_nm=data.frame(counts_nm)
  
  genes_mean<-rowMeans(counts_nm)
  genes_zero<-which(genes_mean == 0)
  genes_mean<-genes_mean[-genes_zero]
  counts_nm_mean<-counts_nm[-genes_zero,]/genes_mean
  
  ########################################
  ####TECHINCAL FEATURES##################
  ########################################
  
  #ONLY CONSIDER READS MAPPES TO GENES (EXCLUDING ERCCS)
  number_mapped_reads_prop=(read_metrics$mapped-read_metrics$ercc)/read_metrics$total
  read_metrics_depending_mapped_reads=read_metrics[, -c(1:3)]
  read_metrics_depending_mapped_reads_prop=read_metrics_depending_mapped_reads/read_metrics$mapped
  
  
  #HOPE THAT THIS WORKS ALSO WHEN REGRESSION NORMALIZATION HAS BEEN APPLIED
  #AS SOME REGRESSION METHODS ALSO PUSH 0 TO ANOTHER VALUE
  detected_genes<-apply(counts_nm, 2, function (x) {
    return(length(which(x > 0)))
  })
  
  counts_nm_mean_log<-log(counts_nm_mean + 0.001)
  genes_var<-apply(counts_nm_mean_log,1,var)
  genes_means_log<-log(genes_mean)
  
  #SELECT ONLY HIGHLY VARIABLES GENES WITH STRONG EXPRESSION BASED ON LAST QUANTILE
  i<-which(genes_var > quantile(genes_var)[4] & genes_means_log > quantile(genes_means_log)[4])
  counts_nm_mean_log_high_var_mean<-counts_nm_mean_log[i,]
  
  #VARIANCE ACROSS HIGHLY EXPRESSED GENES
  transcriptome_variance<-apply(counts_nm_mean_log_high_var_mean, 2, var)
  
  #NUMBER OF HIGHLY EXPRESSED GENES
  m<-quantile(counts_nm_mean[i,2])[4]
  number_of_highly_expressed_variable_genes<-apply(counts_nm_mean[i,], 2, function(x) {
    return(length(which(x > m)))
  })
  
  #NUMBER OF LOW TO HIGH EXPRESED AND VARIABLE GENES PER INTERVAL
  num_of_high_var_exp_genes_interval<-apply(counts_nm_mean_log_high_var_mean, 2, function(x) {
    
    return(hist(x, breaks=c(-1000, -4, -2, 0, 2, 1000), plot=FALSE)$counts)
  })
  num_of_high_var_exp_genes_interval<-t(num_of_high_var_exp_genes_interval)
  colnames(num_of_high_var_exp_genes_interval)<-paste0("num_of_high_var_exp_genes_interval_", 1:ncol(num_of_high_var_exp_genes_interval))
  
  cell_to_mean_corr<-cor(counts_nm, rowMeans(counts_nm), method="spearman")
  colnames(cell_to_mean_corr)=c("cell_to_mean_corr")
  
  
  ########################################
  ####BIOLOGICAL FEATURES##################
  ########################################
  
  #ASSUME IT IS MOUSE
  mouse_GO_BP <- annFUN.org("BP", mapping = "org.Mm.eg.db", ID = "ensembl")
  mouse_GO_CC <- annFUN.org("CC", mapping = "org.Mm.eg.db", ID = "ensembl")
  
  if(organism == "human") {
    mouse_GO_CC <- annFUN.org("CC", mapping = "org.Hs.eg.db", ID = "ensembl")
    mouse_GO_BP <- annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "ensembl")  
  }

  mouse_GO=c(mouse_GO_BP, mouse_GO_CC)
  
  #PROPORTION OF MAPPED READS MAPPED TO THE GO TERM
  go_prop=sapply(unlist(GO_terms), function (go_id) {
    prop=sum_prop(counts_nm, unlist(mouse_GO[go_id])) 
    return(prop)
  }, simplify=FALSE)
  go_prop=do.call(cbind, go_prop)
  go_names=Term(unlist(GO_terms))
  colnames(go_prop)=go_names  
  
  #PROPORTION OF MAPPED READS MAPPED TO SPECIFIC GENES
  extra_genes_prop=sapply(extra_genes, function (extra_g) {
    #print(extra_g[1])
    prop=sum_prop(counts_nm, extra_g[-1]) 
    return(prop)
  }, simplify=FALSE)
  extra_genes_prop=do.call(cbind, extra_genes_prop)
  extra_genes_colnames=sapply(extra_genes, function (extra_g) {
    return(extra_g[1])
  }, simplify=FALSE)
  colnames(extra_genes_prop)=unlist(extra_genes_colnames)
  
  total_reads=read_metrics$total
  features<-cbind(total_reads, number_mapped_reads_prop, read_metrics_depending_mapped_reads_prop, 
                  detected_genes, transcriptome_variance, num_of_high_var_exp_genes_interval, cell_to_mean_corr, number_of_highly_expressed_variable_genes, 
                  go_prop, extra_genes_prop)
  
  rownames(features)=colnames(counts_nm)
  return(features)
}

sum_prop=function(counts, genes_interest) {
  genes_interest_i<-which(rownames(counts) %in% unlist(genes_interest))
  genes_interest_counts_prop<-colSums(counts[genes_interest_i,])
  return(genes_interest_counts_prop)
}

#ASSES QUALITY OF A CELL USING FEATURES AND LABELS FROM 960mES CELLS
asses_cell_quality_SVM<-function(training_set_features, training_set_labels, test_set_features, ensemble_param, output_dir) {
  
  data_set<-data.frame(l=training_set_labels[,2], unlist(as.matrix(training_set_features)))
  formula = l ~ .
  data_set$l <-factor(data_set$l)
  
  test_data<-data.frame(as.matrix(test_set_features))
  
  final_results<-sapply(1:nrow(ensemble_param), function(x) {
    parameters<-ensemble_param[x,]
    weights <- table(data_set$l) 
    weights[1]<-3
    weights[2] <- 1
    kernel<-"radial"
    
    model<-svm(formula, data=data_set, gamma=parameters[1], cost=parameters[2], kernel=kernel, class.weights=weights)
    
    pred_test<- predict(model, test_data)
    svm_test<-as.numeric(levels(pred_test))[pred_test]
    return(svm_test)
  }, simplify=FALSE)
  final_results<-do.call(cbind, final_results)
  final<-vote(final_results)
  return(as.matrix(final))
  
}

#ASSES CELL QUALITY USING PCA AND OUTLIER DETECTION ON ALL DIMENSIONS
asses_cell_quality_PCA<-function(test_set_features, output_dir, prefix) {
  
  pca<-prcomp(test_set_features, scale=TRUE, center=TRUE)
  pca_var_explained=summary(pca)
  pcout_c<-pcout(pca$x[, 1:2])
  
  dim=min(10, ncol(test_set_features))
  uni_2<-uni.plot(pca$x[,1:dim])
  
  low_qual_i=which(pcout_c$wfinal01==0)
  #low_qual_i = which(uni_2$outliers == TRUE)
  
  mtdna=t.test(test_set_features[,6][low_qual_i], test_set_features[,6][-low_qual_i], alternative = "greater")$p.value
  mapped_prop=t.test(test_set_features[,2][low_qual_i], test_set_features[,2][-low_qual_i], alternative = "less")$p.value
  
  
  if ((mapped_prop) > 0.5) {
    types=rep(0, nrow(test_set_features))
    types[low_qual_i] = 1
    
  } else{
    types=rep(1, nrow(test_set_features))
    types[low_qual_i] = 0
  }
  
  data_frame<-data.frame(type=as.character(types), pca$x)  
  col=c("0" = "red","1" = "darkgreen")
  
  
  plot<-ggplot(data_frame, aes(x=PC1, y=PC2)) + geom_point(aes(colour=type)) + theme_bw()  + theme(axis.line=element_blank(),
                                                                                                   axis.text.x=element_text(),axis.text.y=element_text(),
                                                                                                   axis.ticks.length = unit(0, "mm"),
                                                                                                   axis.title.x=element_blank(),
                                                                                                   axis.title.y=element_blank(),
                                                                                                   panel.background=element_blank(),
                                                                                                   panel.border= element_rect(fill=NA,color="black", 
                                                                                                                              linetype="solid"),
                                                                                                   panel.grid.major=element_blank(),
                                                                                                   panel.grid.minor=element_blank(),
                                                                                                   plot.background=element_blank()) + scale_colour_manual(values=col)
  
  top5_PC1_i=order(abs(pca$rotation[,1]), decreasing = TRUE)
  top5_f_pc1=names(pca$rotation[(top5_PC1_i[1:3]),1])
  names=colnames(test_set_features)
  test_set_features_top5_pc1_i=which(names %in% top5_f_pc1)
  size=0.5
  text_size=12
  border_size=1
  plotsPC1=sapply(test_set_features_top5_pc1_i, function(f) {
    
    
    feature=test_set_features[,f]
    df=data.frame(counts=log(feature+0.0001), type=types)
    plot<-ggplot(df, aes(x=type)) + geom_boxplot(aes(colour=factor(type), y = counts), trim=TRUE, alpha=0.3, adjust=1, size=size, outlier.size = 0) +ggtitle(names[f]) + theme_bw()  + theme(axis.line=element_blank(),
                                                                                                                                                                                             axis.text.x=element_blank(),axis.text.y=element_blank(),
                                                                                                                                                                                             axis.ticks.length = unit(0, "mm"),
                                                                                                                                                                                             axis.title.x=element_blank(),
                                                                                                                                                                                             axis.title.y=element_blank(),
                                                                                                                                                                                             legend.position="none",
                                                                                                                                                                                             plot.title=element_text(size=text_size),
                                                                                                                                                                                             panel.background=element_blank(),
                                                                                                                                                                                             panel.border= element_rect(fill=NA,color="black", size=border_size, 
                                                                                                                                                                                                                        linetype="solid"),
                                                                                                                                                                                             panel.grid.major=element_blank(),
                                                                                                                                                                                             panel.grid.minor=element_blank(),
                                                                                                                                                                                             plot.background=element_blank()) + scale_color_manual(values=col)
    return(plot)
  }, simplify=FALSE)
  
  
  top5_PC2_i=order(abs(pca$rotation[,2]), decreasing = TRUE)
  top5_f_pc2=names(pca$rotation[(top5_PC2_i[1:3]),2])
  names=colnames(test_set_features)
  test_set_features_top5_pc2_i=which(names %in% top5_f_pc2)
  
  plotPC2=sapply(test_set_features_top5_pc2_i, function(f) {
    
    feature=test_set_features[,f]
    df=data.frame(counts=log(feature+0.0001), type=types)
    plot<-ggplot(df, aes(x=type)) + geom_boxplot(aes(colour=factor(type), y = counts), trim=TRUE, alpha=0.3, adjust=1, size=size, outlier.size = 0) +ggtitle(names[f]) + theme_bw()  + theme(axis.line=element_blank(),
                                                                                                                                                                                             axis.text.x=element_blank(),axis.text.y=element_blank(),
                                                                                                                                                                                             axis.ticks.length = unit(0, "mm"),
                                                                                                                                                                                             axis.title.x=element_blank(),
                                                                                                                                                                                             axis.title.y=element_blank(),
                                                                                                                                                                                             legend.position="none",
                                                                                                                                                                                             panel.background=element_blank(),
                                                                                                                                                                                             plot.title=element_text(size=text_size),
                                                                                                                                                                                             panel.border= element_rect(fill=NA,color="black", size=border_size, 
                                                                                                                                                                                                                        linetype="solid"),
                                                                                                                                                                                             panel.grid.major=element_blank(),
                                                                                                                                                                                             panel.grid.minor=element_blank(),
                                                                                                                                                                                             plot.background=element_blank()) + scale_color_manual(values=col)
    return(plot)
  }, simplify=FALSE)
  
  l=matrix(c(2,2,3,3,4,4, rep(1,5), 5,rep(1,5), 5,rep(1,5), 6,rep(1,5), 6,rep(1,5), 7),nrow=6)
  l=cbind(l, c(rep(1,5), 7))
  
  t=matrix(c(2,2,3,3,4,4,10, rep(1,5), 5,5,rep(1,5), 5,5,rep(1,5), 6,6,rep(1,5), 6,6,rep(1,5), 7,7),nrow=7)
  t=cbind(t, c(rep(1,5), 7,7))
  
  o<-paste0(output_dir, "/",prefix, "_pca_quality_check.pdf")
  pdf(o, width=10, height=7)
  multiplot(plot,  plotlist=c(plotsPC1,(plotPC2)),layout=l)
  dev.off()
  annot = cbind(rownames(test_set_features), types)
  return(annot)
}


multiplot <- function(..., plotlist=NULL, file, cols=6, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#VOTING SCHEME FOR MODEL ENSEMBLE
vote<-function(predictions) {
  freq<-apply(predictions, 1, function (x) {
    f<-table(x)
    return(names(f)[which.max(f)])
  })
  
  freq<-as.numeric(freq)
  return(freq)
}

#CALCULATE ACCURACY
calculate_accuracy<-function(pred, lab) {
  table(lab[,2], pred)
  success_i<-which(lab[,2]==pred)  
  bad_original_i<-which(lab[,2]==0 & lab[,3]!="S" & lab[,3]!="G")
  bad_all_i<-which(lab[,2]==0)
  good_i<-which(lab[,2]==1)
  
  
  sucess_bad_all_i<-intersect(success_i, bad_all_i)
  sucess_good_i<-intersect(success_i, good_i)
  TP<-length(sucess_bad_all_i)
  TN<-length(sucess_good_i)
  FN<-length(bad_all_i)-TP
  FP<-length(good_i)-TN
  
  specificity = TN/(TN+FP)
  sensitivity = TP/(TP+FN)
  
  r = (c(sensitivity, specificity))
  names(r) = c("Sensitivity", "Specificity")
  return(r)
}

info<-function(string) {
  print(paste0("[INFO]:", string))
}

normalise_by_factor<-function(counts, factor) { 
  t(t(counts)/factor)
}


get_file_name_ext<-function(path) {
  split<-strsplit(path, "/")
  split<-split[[1]]
  return (split[length(split)]) 
}


get_file_name_no_ext<-function(path) {
  file_name<-get_file_name_ext(path)
  file_name<-gsub("\\..*", "", file_name)
  return (file_name) 
}
spec = matrix(
  c("input"  , "i", 1, "character",
    "input_labels"  , "l", 1, "character",
    "output", "o", 1, "character"), 
  byrow=TRUE, ncol=4);
opt = getopt(spec);

input=opt$input
input_labels=opt$input_labels
output=opt$output



################################################################################
###############################TRAINING DATA RAW COUNTS####################################
################################################################################


output_dir="/Users/ti1/Google\ Drive/projects/quality_control/data/not_annotated_data/"
GO_terms_input="/Users/ti1/Google\ Drive/projects/quality_control/tables/biological_features.txt"
extra_genes_input="/Users/ti1/Google\ Drive/projects/quality_control/tables/extra_genes.txt"
common_features_input="/Users/ti1/Google Drive/projects/quality_control/tables/common_features.txt"


input="/Users/ti1/Google\ Drive/projects/quality_control/data/original_labels/ola_mES/ola_mES.sorted.txt"
data<-read.table(input, header=TRUE)
genes=data[, 1]
data=data[,-1]
rownames(data)=genes
data=data[,order(colnames(data))]

input="/Users/ti1/Google\ Drive/projects/quality_control/data/original_labels/ola_mES/ola_mES.stats"
read_metrics=read.table(input, sep=",", header=TRUE)
read_metrics=read_metrics[which(read_metrics[,1] %in% colnames(data)),]
read_metrics=read_metrics[order((read_metrics[,1])),]

file_name<-get_file_name_no_ext(input)

qual_i<-grep("_", rownames(data))
qual_s<-as.character(rownames(data)[qual_i])
data_qual<-data[qual_i,]
data[qual_i,]=data[qual_i,]

ercc_i<-grep("ERCC", rownames(data))
ercc_s<-as.character(rownames(data)[ercc_i])
counts_ercc<-data[ercc_i,]
ercc=colSums(counts_ercc)
ht_seq_features= t(data_qual[c(1, 2),])


read_metrics=read_metrics[,c(2,3,5,6)]
read_metrics = cbind(read_metrics, ht_seq_features, ercc)
counts<-data[-c(ercc_i,qual_i),]
genes=rownames(counts)

counts_nm<-data.frame(normalise_by_factor(counts, colSums(counts)))

output_dir="/Users/ti1/Google Drive/projects/quality_control/data/23_11_15"
file_name=paste0(file_name, "_raw_counts")
common_feature_names=c("Mapped reads %", "Non-exon reads %", "Multi-mapped reads %", "Transcriptome variance", "Cytoplasm %", "mtDNA %", "Mitochondria %")
all_feature_names= c("Total reads", "Mapped reads %", "Multi-mapped reads %", "Non-exon reads%", "Ambigious genes %", "ERCC ratio", "# detected genes", "Transcriptome variance",
                     "Highly variable Int.1", "Highly variable Int.2", "Highly variable Int.3", "Highly variable Int.4", "Highly variable Int.5", 
                     "#Highly expressed", "Apoptosis %", "Metablosim %", "Ribosomes %", "Membrane %", "Cytoplasm %", "Extracellular region %", "mtDNA", "Mitochondria upreg. %", "Mitochondria downreg. %", "Actb %", "Gadph%")

features=(extract_features(read_metrics, counts_nm, genes, file_name,output_dir , common_features_input, GO_terms_input, extra_genes_input, "human"))

quality_PCA = asses_cell_quality_PCA(features[[1]], paste0(output_dir, "/", file_name), paste0(file_name, "_all_features"))
quality_PCA = asses_cell_quality_PCA(features[[2]], paste0(output_dir, "/", file_name), paste0(file_name, "_common_features"))

labels_old=read.table("/Users/ti1/Google Drive/projects/quality_control/data/2_7/ola_mES/re-labeling/labels_and_names_detailed_ola_mES_outliers_uni2_i.txt", header=TRUE)
labels_original=read.table("/Users/ti1/Google\ Drive/projects/quality_control/data/original_labels/ola_mES/labels_and_names_detailed_ola_mES.txt", header=TRUE)

#missed by new
labels_old=labels_old[order(labels_old[,1]),]
quality_PCA=quality_PCA[order(quality_PCA[,1]),]
labels_original=labels_original[order(labels_original[,1]),]


differnt_than_origin=quality_PCA[-which(quality_PCA[,2]==as.character(labels_original[,2])),]
differnt_than_origin_good_to_bad=differnt_than_origin[which(differnt_than_origin[,2]==0),1]
labels_new = labels_original
labels_new[,3]=as.character(labels_original[,3])
labels_new[which(labels_original[,1] %in% differnt_than_origin_good_to_bad), 3]="S"
labels_new[which(labels_original[,1] %in% differnt_than_origin_good_to_bad), 2]=0

table(quality_PCA[which(quality_PCA[,2]==as.character(labels_original[,2])),2])[1]/78
table(quality_PCA[which(quality_PCA[,2]==as.character(labels_original[,2])),2])[2]/882

table(labels_new[,3])
table(labels_original[,3])

################################################################################
###############################TRAINING DATA CUFFLINKS SUPPLEMENT####################################
################################################################################

output_dir="/Users/ti1/Google\ Drive/projects/quality_control/data/not_annotated_data/"
GO_terms_input="/Users/ti1/Google\ Drive/projects/quality_control/tables/biological_features.txt"
extra_genes_input="/Users/ti1/Google\ Drive/projects/quality_control/tables/extra_genes.txt"
common_features_input="/Users/ti1/Google Drive/projects/quality_control/tables/common_features.txt"

# MERGE ALL FPKM COUNTS AND TRANSFORM THEM TO TPM
# fpkmToTpm <- function(fpkm)
# {
#   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
# }
# 
# 
# genes_expr<-list.files(path = "/Users/ti1/Google\ Drive/projects/quality_control/data/original_labels/ola_mES/", pattern = ".*_cufflinks-2_2_1_Linux_x86_64.counts$", all.files = FALSE,
#                        full.names = TRUE, recursive = TRUE)
# 
# data_complete=sapply(genes_expr, function(x) {
#   
#   data<-read.table(x, header=FALSE, stringsAsFactors = F)
#   return(data[, -1])
# }, simplify=FALSE)
# 
# genes<-read.table(genes_expr[1], header=TRUE)[,1]
# data=do.call(cbind, data_complete)
# names=data[1, ]
# data=data[-1, ]
# data <- apply(data,2, as.numeric)
# colnames(data)=names
# 
# genes_dup=which(duplicated(genes)==TRUE)
# isoforms=unique(genes[genes_dup])
# 
# genes_dup=sapply(isoforms, function(x) {
#   
#   i=which(genes==x)
#   return(i)  
# }, simplify=FALSE)
# genes_dup=unlist(genes_dup)
# 
# sums=sapply(isoforms, function(x) {
#   
#   i=which(genes==x)
#   return(colSums(data[i, ]))
#   
# }, simplify=FALSE)
# 
# duplicates_sums=do.call(rbind, sums)
# rownames(duplicates_sums)=isoforms
# 
# data=data[-genes_dup, ]
# genes=genes[-genes_dup]
# rownames(data)=genes
# 
# which(rownames(data) %in% rownames(duplicates_sums))
# data_com=rbind(data, duplicates_sums)
# genes=rownames(data_com)
# data_com=fpkmToTpm(data_com) 
# data_com=data_com[,order(colnames(data_com))]

#write.table(data_com, "/Users/ti1/Google\ Drive/projects/quality_control/data/original_labels/ola_mES/ola_mES_cufflinks.counts")
data=read.table("/Users/ti1/Google\ Drive/projects/quality_control/data/original_labels/ola_mES/ola_mES_cufflinks.counts", header=TRUE)
labels_original=read.table("/Users/ti1/Google\ Drive/projects/quality_control/data/original_labels/ola_mES/labels_and_names_detailed_ola_mES.txt", header=TRUE)

#data=data[,which(colnames(data) %in% labels_original[,1])]

input="/Users/ti1/Google\ Drive/projects/quality_control/data/original_labels/ola_mES/ola_mES.stats"
read_metrics=read.table(input, sep=",", header=TRUE)
read_metrics=read_metrics[which(read_metrics[,1] %in% colnames(data)),]
read_metrics=read_metrics[order((read_metrics[,1])),]

file_name<-get_file_name_no_ext(input)

ercc_i<-grep("ERCC", rownames(data))
ercc_s<-as.character(rownames(data)[ercc_i])
counts_ercc<-data[ercc_i,]
ercc=colSums(counts_ercc)

read_metrics=read_metrics[,c(2,3,5,6)]
read_metrics = cbind(read_metrics, ercc)
counts_nm=data
genes=rownames(counts_nm)

features=(extract_features(read_metrics, counts_nm, genes, paste0(file_name, "_cufflinks"), "/Users/ti1/Google Drive/projects/quality_control/data/23_11_15", common_features_input, GO_terms_input, extra_genes_input))
quality_PCA = asses_cell_quality_PCA(features[[1]], paste0(output_dir, "/", file_name))
#quality_PCA = asses_cell_quality_PCA(features[[2]], paste0(output_dir, "/", file_name))



#missed by new
quality_PCA=quality_PCA[order(quality_PCA[,1]),]
labels_original=labels_original[order(labels_original[,1]),]

differnt_than_origin=quality_PCA[-which(quality_PCA[,2]==as.character(labels_original[,2])),]
differnt_than_origin_good_to_bad=differnt_than_origin[which(differnt_than_origin[,2]==0),1]
labels_new = labels_original
labels_new[,3]=as.character(labels_original[,3])
labels_new[which(labels_original[,1] %in% differnt_than_origin_good_to_bad), 3]="S"
labels_new[which(labels_original[,1] %in% differnt_than_origin_good_to_bad), 2]=0

table(quality_PCA[which(quality_PCA[,2]==as.character(labels_original[,2])),2])[1]/78
table(quality_PCA[which(quality_PCA[,2]==as.character(labels_original[,2])),2])[2]/882

table(labels_new[,3])
table(labels_original[,3])


################################################################################
###############################HUMAN SUPPLEMENT DATA############################
################################################################################
genes_expr<-list.files(path = "/Users/ti1/Google\ Drive/projects/quality_control/data/not_annotated_data/", pattern = ".*\\.counts.txt$", all.files = FALSE,
                       full.names = TRUE, recursive = TRUE)
stats<-list.files(path = "/Users/ti1/Google\ Drive/projects/quality_control/data/not_annotated_data/", pattern = ".*\\.stats$", all.files = FALSE,
                  full.names = TRUE, recursive = TRUE)
extra_genes_human_input="/Users/ti1/Google\ Drive/projects/quality_control/data/not_annotated_data/extra_genes_humans.txt"
output_dir="/Users/ti1/Google\ Drive/projects/quality_control/data/not_annotated_data/"


features_data_sets=sapply(1:length(genes_expr), function(x){
  
  data<-read.table(genes_expr[[x]], header=TRUE)
  data=data[,order(colnames(data))]
  
  read_metrics=read.table(stats[[x]], sep=",", header=TRUE)
  read_metrics=read_metrics[order((read_metrics[,1])),]
  
  file_name<-get_file_name_no_ext(genes_expr[[x]])
  
  qual_i<-grep("_", rownames(data))
  qual_s<-as.character(rownames(data)[qual_i])
  data_qual<-data[qual_i,]
  data[qual_i,]=data[qual_i,]
  
  ercc_i<-grep("ERCC", rownames(data))
  ercc_s<-as.character(rownames(data)[ercc_i])
  counts_ercc<-data[ercc_i,]
  ercc=colSums(counts_ercc)
  ht_seq_features= t(data_qual[c(1, 2),])
  
  read_metrics=read_metrics[,c(2,3,5,6)]
  read_metrics = cbind(read_metrics, ht_seq_features, ercc)
  counts<-data[-c(ercc_i,qual_i),]
  genes=rownames(counts)
  counts_nm<-data.frame(normalise_by_factor(counts, colSums(counts)))
  
  
  features=(extract_features(read_metrics, counts_nm, genes, file_name, output_dir, common_features_input, GO_terms_input, extra_genes_human_input, "human"))
  all=features[[1]]
  common=features[[2]]
  
  colnames(all)=all_feature_names
  colnames(common)=common_feature_names
  quality_PCA = asses_cell_quality_PCA(all, paste0(output_dir, "/",file_name), paste0(file_name, "_all"))
  
  return(features)
}, simplify=FALSE)


combined=rbind(features_data_sets[[1]][[2]], features_data_sets[[2]][[2]])


training_set_common_features=read.table("/Users/ti1/Google Drive/projects/quality_control/data/23_11_15/ola_mES_raw_counts/ola_mES_raw_counts.common.features", header=TRUE)
training_set_all_features=read.table("/Users/ti1/Google Drive/projects/quality_control/data/23_11_15/ola_mES_raw_counts/ola_mES_raw_counts.all.features", header=TRUE)
training_set_labels=read.table("/Users/ti1/Google Drive/projects/quality_control/data/30_12/ola_mES/pre_clustering/ola_mES_no_corr_features/labels_and_names_detailed_ola_mES_pca.txt", header=TRUE)
paramaters_core<-read.table("/Users/ti1/Google\ Drive/projects/quality_control/data/20_11/ola_mES/model_ensemble_param_generic_pca_no_corr_features/ola_mES_generic.ensemble_parameters_boost.txt")



quality_SVM = asses_cell_quality_SVM(training_set_common_features, training_set_labels, features_data_sets[[1]][[2]], paramaters_core, output_dir)

