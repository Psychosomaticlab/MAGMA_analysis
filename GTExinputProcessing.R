library(conflicted)
library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)

library(dplyr)
library(tidyr)

setwd("path for the working directory")
file="./GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
data <- read_tsv(file,skip=2, show_col_types = FALSE) %>% dplyr::rename(gene_id = Name)

data <- data %>% dplyr::select(!matches(c("Bladder","Cervix - Ectocervix","Cervix - Endocervix",
                                   "Cells - Cultured fibroblasts","Cells - EBV-transformed lymphocytes",
                                   "Fallopian Tube","Kidney - Cortex","Kidney - Medulla","Testis")))

colnames(data) <- gsub(" - | ", "_", colnames(data))

#examples for grouping GTEx data
folder<-"./GTEX/brain-cortex"
folder<-"./GTEX/brain-other"


filelist <- list.files(path = folder, pattern = "\\.gct$", full.names = TRUE)
data_list <- lapply(filelist, function(fi) {
  gct_obj <- read_tsv(fi,skip=2, show_col_types = FALSE) %>% dplyr::rename(gene_id = Name)
  return(gct_obj) 
})
merged_data <-merged_data %>% dplyr::select(-contains("id."))
row_medians <- apply(merged_data[, 3:ncol(merged_data)], MARGIN = 1, FUN = median, na.rm = TRUE)
merged_data$Median <- row_medians


braincortex <- merged_data %>% dplyr::select(gene_id, Description, Median) %>% dplyr::rename(Brain_cortex = Median)
brainother <- merged_data %>% dplyr::select(gene_id, Description, Median) %>% dplyr::rename(Brain_other = Median)

data <- data %>% dplyr::select(-contains(c("Adipose","Brain","Artery","Colon","Esophagus","Heart","Skin")))

merged_data <- adipose %>%
  full_join(artery, by = c("gene_id", "Description")) %>%
  full_join(braincortex, by = c("gene_id", "Description")) %>%
  full_join(brainother, by = c("gene_id", "Description"))%>%
  full_join(colon, by = c("gene_id", "Description"))%>%
  full_join(data, by = c("gene_id", "Description"))%>%
  full_join(esophagus, by = c("gene_id", "Description"))%>%
  full_join(heart, by = c("gene_id", "Description")) %>%
  full_join(skin, by = c("gene_id", "Description"))

merged_data <- merged_data %>% add_count(gene_id) %>% dplyr::filter(n==1) %>%dplyr::select(-n) 
merged_data <- merged_data %>%dplyr::mutate(gene_id=gsub("\\..+","",gene_id))

entrez_ensembl <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL)
entrez_ensembl_unique_genes_entrez <- entrez_ensembl %>% count(gene_id) %>% dplyr::filter(n==1)
entrez_ensembl_unique_genes_ens <- entrez_ensembl %>% count(ensembl_id) %>% dplyr::filter(n==1)
entrez_ensembl <- dplyr::filter(entrez_ensembl,gene_id%in%entrez_ensembl_unique_genes_entrez$gene_id & ensembl_id %in% entrez_ensembl_unique_genes_ens$ensembl_id)
colnames(entrez_ensembl)[1] <- "Gene"
colnames(entrez_ensembl)[2] <- "gene_id"
entrez_ensembl$Gene <- as.character(entrez_ensembl$Gene)

save(merged_data, entrez_ensembl, file = "./magma.Rdata")


##for gene binary analysis
exp<- merged_data%>%tidyr::gather(key = Tissue,value=Expr,-gene_id,-Description) %>%as_tibble()
exp1<-inner_join(exp,entrez_ensembl,by = "gene_id")

not_expressed <- exp1 %>% group_by(Gene) %>% dplyr::filter(sum(Expr) == 0) %>%  dplyr::select(Gene) %>% unique()
exp1 <- dplyr::filter(exp1,!Gene%in%not_expressed$Gene)

calculate_tau <- function(Expr_values) {
  normalized_expr <- Expr_values / max(Expr_values)
  tau_value <- sum(1 - normalized_expr) / (length(Expr_values) - 1)
  return(tau_value)
}
calculate_entropy <- function(Expr_values) {
  e <- Expr_values[Expr_values > 0]
  if (length(e) == 0) return(0)
  e <- e / sum(e)
  entropy <- -sum(e * log2(e))  
  return(entropy)
}

exp1 <- exp1 %>% group_by(Gene) %>% mutate(specificity=Expr/sum(Expr)) %>% mutate(max_spe = max(specificity, na.rm = TRUE)) %>%  mutate(tau=calculate_tau(Expr)) %>%  mutate(entropy = calculate_entropy(Expr)) %>%  mutate(specentropy = 1 - (entropy / log2(n()))) %>%ungroup()

n_genes <- length(unique(exp1$Gene))
n_genes_to_keep <- (n_genes * 0.1) %>% round()
topN <- exp1 %>% dplyr::filter(Expr>1) %>% group_by(Tissue) %>% top_n(.,n_genes_to_keep,specificity)

gtexset<-dplyr::select(topN,Tissue,Gene)
write.table(gtexset, file = "./gtexset_1_10_genes.txt", sep = "\t", col.names = FALSE,row.names = FALSE, quote = FALSE)



#for gene property analysis
expp<-merged_data
expp$Average <- rowMeans(expp[, 3:ncol(expp)], na.rm = TRUE)
expp[, 3:ncol(expp)] <- log2(expp[, 3:ncol(expp)])
expp <- expp[is.finite(expp$Average), ]

exp2 <- inner_join(expp,entrez_ensembl,by = "gene_id")
exp2 <- exp2 %>%  dplyr::mutate(dplyr::across(everything(), ~ifelse(. == -Inf, 0, .))) %>%
                  dplyr::select(Gene, everything()) %>%
                  dplyr::select(-c(2, 3))

write.table(exp2, file ="./GTEX_V8_26_log2.txt", sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)


