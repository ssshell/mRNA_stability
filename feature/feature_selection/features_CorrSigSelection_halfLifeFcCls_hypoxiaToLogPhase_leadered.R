
############## feature selection by Spearman correlation and 
############## significance difference of features among target classes 
############## within transcript region and feature type
############## adding correlation with target classes as one of the criteria

library(tidyverse)
library(GGally)
library(RColorBrewer)
library(scales)


####### turn off scientific format
options(scipen = 999)


####### feature table with all features
###### use log phase abundance
path = '../FeatureSet_complete_logPhase_leadered/'
files = list.files(path, pattern = '.txt')
setwd('../FeatureSet_complete_logPhase_leadered/')
table_list = lapply(files, read.table, header = TRUE)
features_all = do.call(cbind, table_list)
setwd('../feature_selection/')


####### get leadered genes features
gene_leadered <- read.table("../../index/msmeg_CombinedAnnotation_CDS_leadered.bed", col.names = c("gene", "start", "end", "strand"))
features_leadered <- features_all %>% 
        rownames_to_column('gene') %>% 
        inner_join(gene_leadered, by = "gene") %>% 
        select(-start, -end, -strand) %>% 
        drop_na()


####### remove zero variance feature
############## test session
#### alternative way to test result
# removeZeroVar <- function(df){
#         df[, !sapply(df, function(x) min(x) == max(x))]
# }
# test <- removeZeroVar(features_leadered)

variances <- apply(features_leadered[-c(1)], 2, var)
variances[which(variances == 0)]
features_leadered_rmZeroVar <- features_leadered %>% 
        select(-AGAAAGGAGGT_5p25ntUpStream)


####### visualizing feature correlation
# features_leadered_rmZeroVar %>% 
#         select(-gene) %>% 
#         ggcorr(size = 0, 
#                legend.position = "None",
#                label = FALSE, 
#                method = c("all.obs", "spearman"), 
#                low = "#053061", mid = "#F7F7F7", high = "#67001F")


####### Spearman correlation (rs) among features
features_leadered_corr_p1 <- cor(features_leadered_rmZeroVar[-c(1)], method = "spearman")

features_leadered_corr_rs <- features_leadered_corr_p1 %>% 
        as.table() %>% as.data.frame() %>% 
        rename(feature1 = Var1, feature2 = Var2, rs = Freq) %>% 
        mutate(rs_abs = abs(rs)) %>% 
        filter(rs_abs >= 0.6 & feature1 != feature2)

features_leadered_corr_rs_uniq <- features_leadered_corr_rs[!duplicated(data.frame(t(apply(features_leadered_corr_rs[1:2], 1, sort)))), ] %>% 
        select(-rs)


####### get feature type annotation
###### Nucleotide features
path = '../FeatureSet_nucleotide_leadered/'
files = list.files(path, pattern = '.txt')
setwd('../FeatureSet_nucleotide_leadered/')
table_list = lapply(files, read.table, header = TRUE)
features_nucleotide = do.call(cbind, table_list)
setwd('../feature_selection')
features_nucleotide_id <- data.frame(colnames(features_nucleotide)) %>% 
        rename(feature = colnames.features_nucleotide.) %>% 
        mutate(type = "nucleotide")

###### Codon features
path = '../FeatureSet_codon_leadered/'
files = list.files(path, pattern = '.txt')
setwd('../FeatureSet_codon_leadered/')
table_list = lapply(files, read.table, header = TRUE)
features_codon = do.call(cbind, table_list)
setwd('../feature_selection')
features_codon_id <- data.frame(colnames(features_codon)) %>% 
        rename(feature = colnames.features_codon.) %>% 
        mutate(type = "codon")

###### Secondary structure features
path = '../FeatureSet_secondaryStructure_leadered/'
files = list.files(path, pattern = '.txt')
setwd('../FeatureSet_secondaryStructure_leadered/')
table_list = lapply(files, read.table, header = TRUE)
features_secondaryStructure = do.call(cbind, table_list)
setwd('../feature_selection')
features_secondaryStructure_id <- data.frame(colnames(features_secondaryStructure)) %>% 
        rename(feature = colnames.features_secondaryStructure.) %>% 
        mutate(type = "secondaryStructure")

###### Ribosome features
path = '../FeatureSet_ribosome_leadered/'
files = list.files(path, pattern = '.txt')
setwd('../FeatureSet_ribosome_leadered/')
table_list = lapply(files, read.table, header = TRUE)
features_ribosome = do.call(cbind, table_list)
setwd('../feature_selection')
features_ribosome_id <- data.frame(colnames(features_ribosome)) %>% 
        rename(feature = colnames.features_ribosome.) %>% 
        mutate(type = "ribosome")

###### Other properties
path = '../FeatureSet_otherProperties_logPhase_leadered/'
files = list.files(path, pattern = '.txt')
setwd('../FeatureSet_otherProperties_logPhase_leadered/')
table_list = lapply(files, read.table, header = TRUE)
features_otherProperties_type = do.call(cbind, table_list)
setwd('../feature_selection')
features_otherProperties_type_id <- data.frame(colnames(features_otherProperties_type)) %>% 
        rename(feature = colnames.features_otherProperties_type.) %>% 
        mutate(type = "otherProperties")


####### get transcript region annotation
###### feature table with 5'UTR features
path = '../FeatureSet_5pUTR/'
files = list.files(path, pattern = '.txt')
setwd('../FeatureSet_5pUTR/')
table_list = lapply(files, read.table, header = TRUE)
features_5pUTR = do.call(cbind, table_list)
setwd('../feature_selection')
features_5pUTR_id <- data.frame(colnames(features_5pUTR)) %>% 
        rename(feature = colnames.features_5pUTR.) %>% 
        mutate(region = "fpr_UTR")

###### feature table with 5' transcript features
path = '../FeatureSet_5pTranscript/'
files = list.files(path, pattern = '.txt')
setwd('../FeatureSet_5pTranscript/')
table_list = lapply(files, read.table, header = TRUE)
features_5pTrans = do.call(cbind, table_list)
setwd('../feature_selection')
features_5pTrans_id <- data.frame(colnames(features_5pTrans)) %>% 
        rename(feature = colnames.features_5pTrans.) %>% 
        mutate(region = "fpr_transcript")

###### feature table with CDS features
path = '../FeatureSet_CDS_logPhase/'
files = list.files(path, pattern = '.txt')
setwd('../FeatureSet_CDS_logPhase/')
table_list = lapply(files, read.table, header = TRUE)
features_CDS = do.call(cbind, table_list)
setwd('../feature_selection')
features_CDS_id <- data.frame(colnames(features_CDS)) %>% 
        rename(feature = colnames.features_CDS.) %>% 
        mutate(region = "CDS")

###### feature table with 3'UTR features
path = '../FeatureSet_3pUTR/'
files = list.files(path, pattern = '.txt')
setwd('../FeatureSet_3pUTR/')
table_list = lapply(files, read.table, header = TRUE)
features_3pUTR = do.call(cbind, table_list)
setwd('../feature_selection')
features_3pUTR_id <- data.frame(colnames(features_3pUTR)) %>% 
        rename(feature = colnames.features_3pUTR.) %>% 
        mutate(region = "tpr_UTR")


####### combine feature tables to include transcript region & type annotation
feature_id_region <- bind_rows(features_5pUTR_id, features_5pTrans_id,
                               features_CDS_id, features_3pUTR_id)

feature_id_type <- bind_rows(features_nucleotide_id, features_codon_id, 
                             features_secondaryStructure_id, features_ribosome_id,
                             features_otherProperties_type_id)


####### add transcript region & type to the Spearman correlation table
features_leadered_corr_f1 <- features_leadered_corr_rs_uniq %>% 
        select(-feature2) %>% 
        rename(feature = feature1) %>% 
        inner_join(feature_id_region, by = "feature") %>% 
        inner_join(feature_id_type, by = "feature")
features_leadered_corr_f2 <- features_leadered_corr_rs_uniq %>% 
        select(-feature1) %>% 
        rename(feature = feature2) %>% 
        inner_join(feature_id_region, by = "feature") %>% 
        inner_join(feature_id_type, by = "feature")


####### statistical tests on significant difference of half-life in target classes
half_life_class <- read.table('../../degradation_classes/half_life/CDS_halfLife_completeClass.txt', header = TRUE)

half_life_class_hypoxiaToLogPhase <- half_life_class %>% 
        select(gene, HalfLife_FCcls_hypoxiaToLogPhase) %>% 
        drop_na()

fc_hypoxiaToLogPhase <- features_leadered_rmZeroVar %>% 
        inner_join(half_life_class_hypoxiaToLogPhase, by = 'gene')

###### Kruskal-Wallis 
kruskal_pval_fc_hypoxiaToLogPhase <- lapply(fc_hypoxiaToLogPhase[, -c(1, ncol(fc_hypoxiaToLogPhase))], function(x) format.pval(kruskal.test(x ~ fc_hypoxiaToLogPhase$HalfLife_FCcls_hypoxiaToLogPhase)$p.value))
kruskal_pval_adjust_fc_hypoxiaToLogPhase <- as.data.frame(do.call(rbind, kruskal_pval_fc_hypoxiaToLogPhase)) %>% 
        rename(kruskal_pval_fc_hypoxiaToLogPhase = V1) %>% 
        mutate(across(contains("kruskal_pval_fc_hypoxiaToLogPhase"), ~ p.adjust(.x, method = "fdr", n = length(.x)), .names="{.col}_fdr")) %>% 
        mutate(kruskal_padjust_log10_fc_hypoxiaToLogPhase = log10(as.numeric(kruskal_pval_fc_hypoxiaToLogPhase))) %>% 
        rownames_to_column("feature")

###### Kolmogorowv-Smirnov
ks_pval_fc_hypoxiaToLogPhase <- list()
for (i in 2:(ncol(fc_hypoxiaToLogPhase) - 1)) {
        feature_i <- fc_hypoxiaToLogPhase[, c(i, ncol(fc_hypoxiaToLogPhase))]
        c1 <- feature_i[which(feature_i$HalfLife_FCcls_hypoxiaToLogPhase == "Large"), ]
        c2 <- feature_i[which(feature_i$HalfLife_FCcls_hypoxiaToLogPhase == "Med-large"), ]
        c3 <- feature_i[which(feature_i$HalfLife_FCcls_hypoxiaToLogPhase == "Med-small"), ]
        c4 <- feature_i[which(feature_i$HalfLife_FCcls_hypoxiaToLogPhase == "Small"), ]
        ks_pval_i = c(colnames(fc_hypoxiaToLogPhase)[i], format.pval(ks.test(c1[, 1], c2[, 1])$p.value),
                      format.pval(ks.test(c1[, 1], c3[, 1])$p.value),
                      format.pval(ks.test(c1[, 1], c4[, 1])$p.value),
                      format.pval(ks.test(c2[, 1], c3[, 1])$p.value),
                      format.pval(ks.test(c2[, 1], c4[, 1])$p.value),
                      format.pval(ks.test(c3[, 1], c4[, 1])$p.value))
        ks_pval_fc_hypoxiaToLogPhase[[i]] <- ks_pval_i
}

ks_pval_adjust_fc_hypoxiaToLogPhase <- as.data.frame(do.call("rbind", ks_pval_fc_hypoxiaToLogPhase)) %>% 
        rename(feature = V1, ks12 = V2, ks13 = V3, ks14 = V4, 
               ks23 = V5, ks24 = V6, ks34 = V7) %>%
        mutate(across(contains("ks12"), ~ p.adjust(.x, method = "fdr", n = length(.x)), .names="{.col}_fdr")) %>% 
        mutate(across(contains("ks13"), ~ p.adjust(.x, method = "fdr", n = length(.x)), .names="{.col}_fdr")) %>% 
        mutate(across(contains("ks14"), ~ p.adjust(.x, method = "fdr", n = length(.x)), .names="{.col}_fdr")) %>% 
        mutate(across(contains("ks23"), ~ p.adjust(.x, method = "fdr", n = length(.x)), .names="{.col}_fdr")) %>% 
        mutate(across(contains("ks24"), ~ p.adjust(.x, method = "fdr", n = length(.x)), .names="{.col}_fdr")) %>% 
        mutate(across(contains("ks34"), ~ p.adjust(.x, method = "fdr", n = length(.x)), .names="{.col}_fdr")) %>% 
        rowwise() %>% 
        mutate(ks_pval_fdr_fc_hypoxiaToLogPhase = min(as.numeric(ks12_fdr), as.numeric(ks13_fdr), as.numeric(ks14_fdr),
                                             as.numeric(ks23_fdr), as.numeric(ks24_fdr), as.numeric(ks34_fdr))) %>% 
        mutate(ks_padjust_log10_fc_hypoxiaToLogPhase = log10(ks_pval_fdr_fc_hypoxiaToLogPhase)) %>% 
        select(feature, ks_pval_fdr_fc_hypoxiaToLogPhase, ks_padjust_log10_fc_hypoxiaToLogPhase)

###### Combining two statistical test results 
fc_hypoxiaToLogPhase_classSig <- kruskal_pval_adjust_fc_hypoxiaToLogPhase %>% 
        inner_join(ks_pval_adjust_fc_hypoxiaToLogPhase, by = "feature") %>% 
        rowwise() %>% 
        mutate(avg_fdrPval_fc_hypoxiaToLogPhase = mean(c(kruskal_pval_fc_hypoxiaToLogPhase_fdr, ks_pval_fdr_fc_hypoxiaToLogPhase))) %>% 
        select(feature, kruskal_pval_fc_hypoxiaToLogPhase_fdr, ks_pval_fdr_fc_hypoxiaToLogPhase, avg_fdrPval_fc_hypoxiaToLogPhase)


####### Kendall's tau correlation (tau) between features and target class
fc_hypoxiaToLogPhase_ordinalCls <- fc_hypoxiaToLogPhase %>% 
        mutate(HalfLife_FCcls_hypoxiaToLogPhase_ordinal = case_when(HalfLife_FCcls_hypoxiaToLogPhase == "Small" ~ 1,
                                                         HalfLife_FCcls_hypoxiaToLogPhase == "Med-small" ~ 2,
                                                         HalfLife_FCcls_hypoxiaToLogPhase == "Med-large" ~ 3,
                                                         HalfLife_FCcls_hypoxiaToLogPhase == "Large" ~ 4)) %>% 
        select(-c(gene, HalfLife_FCcls_hypoxiaToLogPhase))
fc_hypoxiaToLogPhase_ordinalCls_features <- select(fc_hypoxiaToLogPhase_ordinalCls, -HalfLife_FCcls_hypoxiaToLogPhase_ordinal)

features_leadered_corrTarget_tau_p1 <- lapply(fc_hypoxiaToLogPhase_ordinalCls_features, function(x) cor(x, fc_hypoxiaToLogPhase_ordinalCls$HalfLife_FCcls_hypoxiaToLogPhase_ordinal, method = "kendall"))
features_leadered_corrTarget_tau <- as.data.frame(do.call(rbind, features_leadered_corrTarget_tau_p1)) %>% 
        rename(kendall_fc_hypoxiaToLogPhase = V1) %>% 
        mutate(kendall_fc_hypoxiaToLogPhase_abs = abs(kendall_fc_hypoxiaToLogPhase)) %>% 
        rownames_to_column("feature")


####### Combining all feature selection metrics
features_leadered_corr_f1_classSig_fc_hypoxiaToLogPhase <- features_leadered_corr_f1 %>% 
        inner_join(fc_hypoxiaToLogPhase_classSig, by = 'feature') %>% 
        inner_join(features_leadered_corrTarget_tau, by = 'feature') %>% 
        select(-kendall_fc_hypoxiaToLogPhase) %>% 
        rename(feature1 = feature, rs_abs_feature1 = rs_abs, region_feature1 = region, type_feature1 = type,
               kruskal_pval_fc_hypoxiaToLogPhase_fdr_feature1 = kruskal_pval_fc_hypoxiaToLogPhase_fdr,
               ks_pval_fdr_fc_hypoxiaToLogPhase_feature1 = ks_pval_fdr_fc_hypoxiaToLogPhase,
               avg_fdrPval_fc_hypoxiaToLogPhase_feature1 = avg_fdrPval_fc_hypoxiaToLogPhase,
               kendall_fc_hypoxiaToLogPhase_abs_feature1 = kendall_fc_hypoxiaToLogPhase_abs)

features_leadered_corr_f2_classSig_fc_hypoxiaToLogPhase <- features_leadered_corr_f2 %>% 
        inner_join(fc_hypoxiaToLogPhase_classSig, by = 'feature') %>%
        inner_join(features_leadered_corrTarget_tau, by = 'feature') %>% 
        select(-kendall_fc_hypoxiaToLogPhase) %>%
        rename(feature2 = feature, rs_abs_feature2 = rs_abs, region_feature2 = region, type_feature2 = type,
               kruskal_pval_fc_hypoxiaToLogPhase_fdr_feature2 = kruskal_pval_fc_hypoxiaToLogPhase_fdr,
               ks_pval_fdr_fc_hypoxiaToLogPhase_feature2 = ks_pval_fdr_fc_hypoxiaToLogPhase,
               avg_fdrPval_fc_hypoxiaToLogPhase_feature2 = avg_fdrPval_fc_hypoxiaToLogPhase,
               kendall_fc_hypoxiaToLogPhase_abs_feature2 = kendall_fc_hypoxiaToLogPhase_abs)

features_leadered_corrSig_fc_hypoxiaToLogPhase <- bind_cols(features_leadered_corr_f1_classSig_fc_hypoxiaToLogPhase, features_leadered_corr_f2_classSig_fc_hypoxiaToLogPhase) %>%
        select(-rs_abs_feature2) %>%
        rename(rs_abs = rs_abs_feature1)

write.table(features_leadered_corrSig_fc_hypoxiaToLogPhase, "./corrSig_halfLifeFcCls_hypoxiaToLogPhase_leadered.txt", quote = FALSE, row.names = FALSE, sep = "\t")

