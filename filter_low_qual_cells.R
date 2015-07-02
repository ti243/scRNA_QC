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


#EXTRACT FEATURES FOR DATASET
extract_features<-function(input, output_dir, common_features_input, GO_terms_input, extra_genes_input) {
  info ("Extract features")
  info(paste0("Gene expression: ", input))
  info(paste0("Features output dir: ", output_dir))
  info(paste0("GO Terms to use: ", GO_terms_input))
  info(paste0("Common Features: ", common_features_input))
  info(paste0("Extra genes Features: ", extra_genes_input))
  
  
  file_name<-get_file_name_no_ext(input)
  output_dir= paste(output_dir, file_name ,sep="/")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  #LOAD ALL DATA
  data<-read.table(input, header=TRUE)
  GO_terms<-read.table(GO_terms_input, header=FALSE, stringsAsFactors=FALSE)
  extra_genes=read.table(extra_genes_input, header=FALSE, stringsAsFactors=FALSE, fill=TRUE)
  common_features=unlist(read.table(common_features_input, header=FALSE, stringsAsFactors=FALSE, fill=TRUE))
  extra_genes=apply(extra_genes, 1, function(x){
    return(subset(x, x != ""))
  })
  
  info(paste0("DATA Loaded."))
  

  cell_names=colnames(data)
  cell_names<-cell_names[-1]
  
  genes<-data[,1]
  counts<-data[,-1]
  rownames(counts)<-genes
  
  #GENERATE ALL FEATURES       
  features_all<-feature_generation(counts, FALSE, GO_terms, extra_genes)
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
feature_generation<-function(counts, use_ercc, GO_terms, extra_genes) {
  
  features<-list()
  
  qual_i<-grep("_", rownames(counts))
  qual_s<-as.character(rownames(counts)[qual_i])
  counts_qual<-counts[qual_i,]
  
  #ERCCs (IF AVAILABLE)
  ercc_i<-grep("ERCC", rownames(counts))
  ercc_s<-as.character(rownames(counts)[ercc_i])
  counts_ercc<-counts[ercc_i,]
  
  #NORMALIZE CELLS BY TOTAL NUMBER OF GENERATED READS (INCLUDING ERCC, GENES, UNMAPPED ETC)
  number_total_reads<-colSums(counts)
  counts<-counts[-c(ercc_i,qual_i),]
  counts_nm<-normalise_by_factor(counts, number_total_reads)
  
  #UNMAPPED, MAPPED TO INTRONS ETC.
  counts_ht_seq_features_prop<-t(counts_qual)/number_total_reads
  
  #REMOVE ALL 0 GENES
  genes_mean<-rowMeans(counts_nm)
  genes_zero<-which(genes_mean == 0)
  genes_mean<-genes_mean[-genes_zero]
  counts_nm_mean<-counts_nm[-genes_zero,]/genes_mean
  
  ########################################
  ####TECHINCAL FEATURES##################
  ########################################
  
  #ONLY CONSIDER READS MAPPES TO GENES (EXCLUDING ERCCS)
  number_mapped_reads<-colSums(counts)
  number_mapped_reads_prop<-number_mapped_reads/number_total_reads
  
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
    
    return(hist(x, breaks=c(-1000, -4, -2, 0, 2, 4 ,1000), plot=FALSE)$counts)
  })
  num_of_high_var_exp_genes_interval<-t(num_of_high_var_exp_genes_interval)
  colnames(num_of_high_var_exp_genes_interval)<-paste0("num_of_high_var_exp_genes_interval_", 1:ncol(num_of_high_var_exp_genes_interval))
    
  if (use_ercc == TRUE) {
    ercc_ratio<-colSums(counts_ercc)/colSums(counts)
  } else {
    ercc_ratio<-rep(0,ncol(counts_ercc))
  }
  cell_to_mean_corr<-cor(counts, rowMeans(counts), method="spearman")
  colnames(cell_to_mean_corr)=c("cell_to_mean_corr")
  ########################################
  ####BIOLOGICAL FEATURES##################
  ########################################
  mouse_GO_BP <- annFUN.org("BP", mapping = "org.Mm.eg.db", ID = "ensembl")
  mouse_GO_CC <- annFUN.org("CC", mapping = "org.Mm.eg.db", ID = "ensembl")
  mouse_GO=c(mouse_GO_BP, mouse_GO_CC)
  
  #PROPORTION OF MAPPED READS MAPPED TO THE GO TERM
  go_prop=sapply(unlist(GO_terms), function (go_id) {
    prop=calulate_prop(counts, unlist(mouse_GO[go_id]), number_mapped_reads) 
    return(prop)
  }, simplify=FALSE)
  go_prop=do.call(cbind, go_prop)
  go_names=Term(unlist(GO_terms))
  colnames(go_prop)=go_names  

  #PROPORTION OF MAPPED READS MAPPED TO SPECIFIC GENES
  extra_genes_prop=sapply(extra_genes, function (extra_g) {
    prop=calulate_prop(counts, extra_g[-1], number_mapped_reads) 
    return(prop)
  }, simplify=FALSE)
  extra_genes_prop=do.call(cbind, extra_genes_prop)
  extra_genes_colnames=sapply(extra_genes, function (extra_g) {
    return(extra_g[1])
  }, simplify=FALSE)
  colnames(extra_genes_prop)=unlist(extra_genes_colnames)
  

  features<-cbind(number_total_reads,number_mapped_reads_prop, counts_ht_seq_features_prop, 
                  detected_genes, transcriptome_variance, num_of_high_var_exp_genes_interval, ercc_ratio, cell_to_mean_corr, number_of_highly_expressed_variable_genes, 
                  go_prop, extra_genes_prop)
  
  return(features)
}

calulate_prop=function(counts, genes_interest, number_mapped_reads) {
  genes_interest_i<-which(rownames(counts) %in% unlist(genes_interest))
  genes_interest_counts_prop<-colSums(counts[genes_interest_i,])/number_mapped_reads
  return(genes_interest_counts_prop)
}

#ASSES QUALITY OF A CELL USING FEATURES AND LABELS FROM 960mES CELLS
asses_cell_quality<-function(training_set_features, training_set_labels, test_set_features, ensemble_param, output_dir) {
  
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


output_dir="/Users/ti1/Google\ Drive/projects/quality_control/data/2_7"
GO_terms_input="/Users/ti1/Google\ Drive/projects/quality_control/tables/biological_features.txt"
extra_genes_input="/Users/ti1/Google\ Drive/projects/quality_control/tables/extra_genes.txt"
common_features_input="/Users/ti1/Google Drive/projects/quality_control/tables/common_features.txt"

#TRAINING DATA
input="/Users/ti1/Google\ Drive/projects/quality_control/data/original_labels/ola_mES/ola_mES.sorted.txt"
extract_features(input, output_dir, common_features_input, GO_terms_input, extra_genes_input)

sorted_gene_expr<-list.files(path = "/Users/ti1/Google Drive/projects/quality_control/data/original_labels//", pattern = ".*\\.sorted.txt$", all.files = FALSE,
                      full.names = TRUE, recursive = TRUE,

features_data_sets=sapply(sorted_gene_expr[7:length(sorted_gene_expr)], function(x){
  features=extract_features(x, output_dir, common_features_input, GO_terms_input, extra_genes_input)
})
                           
training_set_all_features=read.table("/Users/ti1/Google Drive/projects/quality_control/data/2_7/ola_mES/ola_mES.all.features", header=TRUE)
training_set_common_features=read.table("/Users/ti1/Google Drive/projects/quality_control/data/2_7/ola_mES//ola_mES.common.features", header=TRUE)
training_set_labels=read.table("/Users/ti1/Google Drive/projects/quality_control/data/30_12/ola_mES/pre_clustering/ola_mES_no_corr_features/labels_and_names_detailed_ola_mES_pca.txt", header=TRUE)
  
paramaters_core<-read.table("/Users/ti1/Google\ Drive/projects/quality_control/data/20_11/ola_mES/model_ensemble_param_generic_pca_no_corr_features/ola_mES_generic.ensemble_parameters_boost.txt")

new_labels<-asses_cell_quality(training_set_common_features, training_set_labels, kedar_features_common, paramaters_core, output_dir)
accuracy = calculate_accuracy(final, kedar_labels)


