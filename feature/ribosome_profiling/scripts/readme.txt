
####### versions
sratoolkit.3.0.0
FastQC v0.11.9
Trimmomatic v0.39
Bowtie2 v2.4.5
samtools v1.16.1
StringTie v2.2.1

####### seq_processing
###### this folder contains all the sequencing libraries processing steps
##### download SRR files and convert to fastq files
#### ./TotalRibo
prefetch -v SRR8668649
prefetch -v SRR8668650
prefetch -v SRR8668657
prefetch -v SRR8668658

fasterq-dump SRR8668649
fasterq-dump SRR8668650
fasterq-dump SRR8668657
fasterq-dump SRR8668658

#### ./TotalmRNA
prefetch -v SRR8668651
prefetch -v SRR8668652
prefetch -v SRR8668659
prefetch -v SRR8668660

fasterq-dump SRR8668651
fasterq-dump SRR8668652
fasterq-dump SRR8668659
fasterq-dump SRR8668660

##### preprocessing of ribo seq
./fastqc_raw_ribo.sh
./trim_ribo.sh
./fastqc_trim_ribo.sh

##### preprocessing of mRNA seq
./fastqc_raw_mrna.sh

##### alignment of ribo seq
./bowtie2_ribo_tRNArRNA.sh
./bowtie2_ribo_genome.sh
./sam2bam_ribo.sh
./sort_bam_ribo.sh
./filter_ribo.sh

##### alignment of mRNA seq
./bowtie2_mRNA_tRNArRNA.sh
./bowtie2_mRNA_genome.sh
./sam2bam_mRNA.sh
./sort_bam_mRNA.sh
./filter_mRNA.sh


####### coverage
###### this folder contains all the quantification of coverage
##### get annotations gtf
python3 gene_gtf.py ../seq_processing/index/msmeg_CombinedAnnotation_CDS.bed ../seq_processing/index/msmeg_CombinedAnnotation_CDS_leaderless.txt CDS > ../seq_processing/index/msmeg_CombinedAnnotation_CDS.gtf 
python3 gene_gtf.py ../seq_processing/index/msmeg_CombinedAnnotation_CDS.bed ../seq_processing/index/msmeg_CombinedAnnotation_CDS_leaderless.txt CDS_20up > ../seq_processing/index/msmeg_CombinedAnnotation_CDS_20up.gtf
python3 gene_gtf.py ../seq_processing/index/msmeg_CombinedAnnotation_CDS.bed ../seq_processing/index/msmeg_CombinedAnnotation_CDS_leaderless.txt 5p_end > ../seq_processing/index/msmeg_CombinedAnnotation_CDS_5pend.gtf
python3 gene_gtf.py ../seq_processing/index/msmeg_CombinedAnnotation_CDS.bed ../seq_processing/index/msmeg_CombinedAnnotation_CDS_leaderless.txt 5p_end_excl > ../seq_processing/index/msmeg_CombinedAnnotation_CDS_5pend_excl.gtf
python3 gene_gtf.py ../seq_processing/index/msmeg_CombinedAnnotation_CDS_splitByPortions_last20percent.bed ../seq_processing/index/msmeg_CombinedAnnotation_CDS_leaderless.txt CDS > ../seq_processing/index/msmeg_CombinedAnnotation_CDS_splitByPortions_last20percent.gtf

#### adding regional annotation for ribosome profile
python3 gene_gtf.py ../seq_processing/index/msmeg_CombinedAnnotation_CDS_oneThird_5p.bed ../seq_processing/index/msmeg_CombinedAnnotation_CDS_leaderless.txt CDS > ../seq_processing/index/msmeg_CombinedAnnotation_CDS_oneThird_5p.gtf
python3 gene_gtf.py ../seq_processing/index/msmeg_CombinedAnnotation_CDS_oneThird_mid.bed ../seq_processing/index/msmeg_CombinedAnnotation_CDS_leaderless.txt CDS > ../seq_processing/index/msmeg_CombinedAnnotation_CDS_oneThird_mid.gtf
python3 gene_gtf.py ../seq_processing/index/msmeg_CombinedAnnotation_CDS_oneThird_3p.bed ../seq_processing/index/msmeg_CombinedAnnotation_CDS_leaderless.txt CDS > ../seq_processing/index/msmeg_CombinedAnnotation_CDS_oneThird_3p.gtf

##### get TPM value of each replicate
#### ribo seq coverage
./coverage_cds_ribo.sh
./coverage_cds_20up_ribo.sh
./coverage_5pend_ribo.sh
./coverage_5pend_excl_ribo.sh
./coverage_3pLast20percent_ribo.sh
./coverage_cds_oneThird5p_ribo.sh
./coverage_cds_oneThirdMid_ribo.sh
./coverage_cds_oneThird3p_ribo.sh

#### mRNA seq coverage
./coverage_cds_mRNA.sh
./coverage_cds_20up_mRNA.sh
./coverage_5pend_mRNA.sh
./coverage_5pend_excl_mRNA.sh
./coverage_3pLast20percent_mRNA.sh
./coverage_cds_oneThird5p_mRNA.sh
./coverage_cds_oneThirdMid_mRNA.sh
./coverage_cds_oneThird3p_mRNA.sh

#### extract TPM
./coverage_tpm_ribo.sh
./coverage_tpm_mRNA.sh
./coverage_tpm_ribo_last20percent.sh
./coverage_tpm_mRNA_last20percent.sh
./coverage_tpm_ribo_cds_oneThird5p.sh
./coverage_tpm_ribo_cds_oneThirdMid.sh
./coverage_tpm_ribo_cds_oneThird3p.sh
./coverage_tpm_mRNA_oneThird5p.sh
./coverage_tpm_mRNA_oneThirdMid.sh
./coverage_tpm_mRNA_oneThird3p.sh

#### combine samples
./combined_tpm_ribo.sh
./combined_tpm_mRNA.sh
./combined_tpm_ribo_last20percent.sh
./combined_tpm_mRNA_last20percent.sh
./combined_tpm_ribo_oneThird5p.sh
./combined_tpm_ribo_oneThirdMid.sh
./combined_tpm_ribo_oneThird3p.sh
./combined_tpm_mRNA_oneThird5p.sh
./combined_tpm_mRNA_oneThirdMid.sh
./combined_tpm_mRNA_oneThird3p.sh

##### normalize coverage
python3 normalized_ribo_profile.py ../seq_processing/bowtie2_ribo_genome/coverage_tpm/combined_coverage_tpm_ribo.txt ../seq_processing/bowtie2_mRNA_genome/coverage_tpm/combined_coverage_tpm_mRNA.txt > ../ribo_normalizedBy_totalmRNA.txt

python3 normalized_ribo_profile_3pLast20percent.py ../seq_processing/bowtie2_ribo_genome/coverage_tpm/combined_coverage_tpm_ribo_3pLast20percent.txt ../seq_processing/bowtie2_mRNA_genome/coverage_tpm/combined_coverage_tpm_mRNA_3pLast20percent.txt > ../ribo_normalizedBy_totalmRNA_3pLast20percent.txt

python3 normalized_ribo_profile_oneThird.py ../seq_processing/bowtie2_ribo_genome/coverage_tpm/combined_coverage_tpm_ribo_oneThird5p.txt ../seq_processing/bowtie2_mRNA_genome/coverage_tpm/combined_coverage_tpm_mRNA_oneThird5p.txt > ../ribo_normalizedBy_totalmRNA_oneThird5p.txt

python3 normalized_ribo_profile_oneThird.py ../seq_processing/bowtie2_ribo_genome/coverage_tpm/combined_coverage_tpm_ribo_oneThirdMid.txt ../seq_processing/bowtie2_mRNA_genome/coverage_tpm/combined_coverage_tpm_mRNA_oneThirdMid.txt > ../ribo_normalizedBy_totalmRNA_oneThirdMid.txt

python3 normalized_ribo_profile_oneThird.py ../seq_processing/bowtie2_ribo_genome/coverage_tpm/combined_coverage_tpm_ribo_oneThird3p.txt ../seq_processing/bowtie2_mRNA_genome/coverage_tpm/combined_coverage_tpm_mRNA_oneThird3p.txt > ../ribo_normalizedBy_totalmRNA_oneThird3p.txt

paste ../ribo_normalizedBy_totalmRNA_oneThird5p.txt ../ribo_normalizedBy_totalmRNA_oneThirdMid.txt ../ribo_normalizedBy_totalmRNA_oneThird3p.txt | cut -f1-2,4,6 > ../ribo_normalizedBy_totalmRNA_oneThird.txt

