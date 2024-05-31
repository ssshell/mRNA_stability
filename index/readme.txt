
#######
get tss_max_coverage.txt from ~/mRNA_stability/index_msmeg/CarlaAnnotations/
get annotation files from ~/mRNA_stability/index_msmeg/
#######

####### get leadered, leaderless & leader_undefined CDS
sed 1d msmeg_5pUTR_length.txt | awk '$2!="NA" && $2!=0' | cut -f1 > msmeg_CombinedAnnotation_CDS_leadered.txt
awk 'NR==FNR {id[$1]; next} $1 in id' msmeg_CombinedAnnotation_CDS_leadered.txt msmeg_CombinedAnnotation_CDS.bed > msmeg_CombinedAnnotation_CDS_leadered.bed

sed 1d msmeg_5pUTR_length.txt | awk '$2==0' | cut -f1 > msmeg_CombinedAnnotation_CDS_leaderless.txt
awk 'NR==FNR {id[$1]; next} $1 in id' msmeg_CombinedAnnotation_CDS_leaderless.txt msmeg_CombinedAnnotation_CDS.bed > msmeg_CombinedAnnotation_CDS_leaderless.bed

sed 1d msmeg_5pUTR_length.txt | awk '$2=="NA"' | cut -f1 > msmeg_CombinedAnnotation_CDS_leaderUndefined.txt
awk 'NR==FNR {id[$1]; next} $1 in id' msmeg_CombinedAnnotation_CDS_leaderUndefined.txt msmeg_CombinedAnnotation_CDS.bed > msmeg_CombinedAnnotation_CDS_leaderUndefined.bed

cat msmeg_CombinedAnnotation_CDS_leadered.bed msmeg_CombinedAnnotation_CDS_leaderUndefined.bed > msmeg_CombinedAnnotation_CDS_leaderedAndleaderUndefined.bed

####### get operon predictions
###### ~/index_msmeg/CarlaAnnotations/msmeg_predictedOperon.csv

####### add Ms1 to annotation
###### to get coverage for half-life with Ms1
##### Ms1: 6242368-6242625. +
msmeg_CombinedAnnotation.bed -> msmeg_CombinedAnnotation_Ms1Added.bed
msmeg_CDS_length.txt -> msmeg_Ms1Added_length.txt

