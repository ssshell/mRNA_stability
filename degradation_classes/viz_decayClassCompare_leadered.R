
####### compare half-life and pattern degradation class for leadered genes
###### UpSetR: https://github.com/hms-dbmi/UpSetR/blob/master/README.md

library(tidyverse)
library(UpSetR)

####### complete half-life class
hl_CDS_complete_class <- read.table('./half_life/CDS_halfLife_completeClass.txt', header = TRUE)

####### get leadered CDS annotations
gene_leadered <- read.table("../index/msmeg_CombinedAnnotation_CDS_leadered.bed", col.names = c("gene", "start", "end", "strand"))

####### log phase 
###### get half-life class
hl_CDS_leadered_logPhase <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_cls_logPhase) %>% 
        drop_na() %>% 
        mutate(decay_class = sub("$", "_HalfLife", HalfLife_cls_logPhase)) %>% 
        select(gene, decay_class) %>% 
        inner_join(gene_leadered, by = "gene") %>% 
        select(-start, -end, -strand)

###### get decay pattern class
pattern_CDS_leadered_logPhase <- read.table("./pattern/CDS_pattern_logPhaseClass.txt",
                                   header = TRUE) %>% 
        # drop_na() %>% 
        mutate(decay_class = sub("$", "_Pattern", Pattern_cls_logPhase)) %>% 
        select(gene, decay_class) %>% 
        inner_join(gene_leadered, by = "gene") %>% 
        select(-start, -end, -strand)

###### get genes have both half-life and decay pattern class
hlAndPattern_logPhase <- inner_join(hl_CDS_leadered_logPhase, pattern_CDS_leadered_logPhase, by = "gene",
                                    suffix = c("_HalfLife", "_Pattern"))

hl_CDS_leadered_logPhase_common <- hl_CDS_leadered_logPhase[(hl_CDS_leadered_logPhase$gene %in% hlAndPattern_logPhase$gene), ]
pattern_CDS_leadered_logPhase_common <- pattern_CDS_leadered_logPhase[(pattern_CDS_leadered_logPhase$gene %in% hlAndPattern_logPhase$gene), ]

decay_complete_logPhase <- bind_rows(hl_CDS_leadered_logPhase_common, pattern_CDS_leadered_logPhase_common) %>% 
        mutate(decay_class = factor(decay_class))

###### get decay class list
class_msmeg_logPhase <- sort(unique(decay_complete_logPhase$decay_class))
class_msmeg_logPhase_list = list()
for (class_i in class_msmeg_logPhase) {
        class_temp_i <- decay_complete_logPhase %>% 
                filter(decay_class == class_i)
        class_msmeg_logPhase_list[[class_i]] <- class_temp_i$gene
}

###### compare decay class
############# test session for sets order
upset(fromList(class_msmeg_logPhase_list), nsets = 8, order.by = "freq",
      # sets = c("Fast_HalfLife", "Med-fast_HalfLife",
      #          "Med-slow_HalfLife", "Slow_HalfLife",
      #          "Fast_Pattern", "Med-fast_Pattern",
      #          "Med-slow_Pattern", "Slow_Pattern"),
      sets = c("Slow_HalfLife", "Slow_Pattern",
               "Med-slow_HalfLife", "Med-slow_Pattern",
               "Med-fast_HalfLife", "Med-fast_Pattern",
               "Fast_HalfLife", "Fast_Pattern"),
      keep.order = TRUE)
############# end of test session

upset_logPhase <- upset(fromList(class_msmeg_logPhase_list), nsets = 8, order.by = "freq",
                        point.size = 3,
                        line.size = 1,
                        text.scale = c(1.2, 1, 1.2, 1, 1, 1.5),
                        sets = c("Slow_HalfLife", "Slow_Pattern",
                                 "Med-slow_HalfLife", "Med-slow_Pattern",
                                 "Med-fast_HalfLife", "Med-fast_Pattern",
                                 "Fast_HalfLife", "Fast_Pattern"),
                        mainbar.y.label = "Intersected gene counts of decay class in log phase", sets.x.label = "Gene counts of decay class",
                        keep.order = TRUE)
png("./viz_toCompareDecayClass_byHalfLifeAndPattern/CDS_leadered_decayClass_commonHalfLifeAndPattern_logPhase.png", 
    res = 600, width = 5000, height = 3000)
print(upset_logPhase)
dev.off()

####### hypoxia
###### get half-life class
hl_CDS_leadered_hypoxia <- hl_CDS_complete_class %>% 
        select(gene, HalfLife_cls_hypoxia) %>% 
        drop_na() %>%
        mutate(decay_class = sub("$", "_HalfLife", HalfLife_cls_hypoxia)) %>% 
        select(gene, decay_class) %>% 
        inner_join(gene_leadered, by = "gene") %>% 
        select(-start, -end, -strand)

###### get decay class
pattern_CDS_leadered_hypoxia <- read.table("./pattern/CDS_pattern_hypoxiaClass.txt",
                                  header = TRUE) %>% 
        # drop_na() %>%
        mutate(decay_class = sub("$", "_Pattern", Pattern_cls_hypoxia)) %>% 
        select(gene, decay_class) %>% 
        inner_join(gene_leadered, by = "gene") %>% 
        select(-start, -end, -strand)

###### get genes have both half-life and decay pattern class
hlAndPattern_hypoxia <- inner_join(hl_CDS_leadered_hypoxia, pattern_CDS_leadered_hypoxia, by = "gene",
                                   suffix = c("_HalfLife", "_Pattern"))

hl_CDS_leadered_hypoxia_common <- hl_CDS_leadered_hypoxia[(hl_CDS_leadered_hypoxia$gene %in% hlAndPattern_hypoxia$gene), ]
pattern_CDS_leadered_hypoxia_common <- pattern_CDS_leadered_hypoxia[(pattern_CDS_leadered_hypoxia$gene %in% hlAndPattern_hypoxia$gene), ]

decay_complete_hypoxia <- bind_rows(hl_CDS_leadered_hypoxia_common, pattern_CDS_leadered_hypoxia_common) %>% 
        mutate(decay_class = factor(decay_class))

###### get decay class list
class_msmeg_hypoxia <- sort(unique(decay_complete_hypoxia$decay_class))
class_msmeg_hypoxia_list = list()
for (class_i in class_msmeg_hypoxia) {
        class_temp_i <- decay_complete_hypoxia %>% 
                filter(decay_class == class_i)
        class_msmeg_hypoxia_list[[class_i]] <- class_temp_i$gene
}

###### compare decay class
upset_hypoxia <- upset(fromList(class_msmeg_hypoxia_list), nsets = 8, order.by = "freq",
                       point.size = 3,
                       line.size = 1,
                       text.scale = c(1.2, 1, 1.2, 1, 1, 1.5),
                       sets = c("Slow_HalfLife", "Slow_Pattern",
                                "Med-slow_HalfLife", "Med-slow_Pattern",
                                "Med-fast_HalfLife", "Med-fast_Pattern",
                                "Fast_HalfLife", "Fast_Pattern"),
                       mainbar.y.label = "Intersected gene counts of decay class in hypoxia", sets.x.label = "Gene counts of decay class",
                       keep.order = TRUE)
png("./viz_toCompareDecayClass_byHalfLifeAndPattern/CDS_leadered_decayClass_commonHalfLifeAndPattern_hypoxia.png", 
    res = 600, width = 5000, height = 3000)
print(upset_hypoxia)
dev.off()
