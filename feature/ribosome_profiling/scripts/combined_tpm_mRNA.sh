#!/bin/bash

cd ../seq_processing/bowtie2_mRNA_genome/coverage_tpm

paste SRR8668651_cds_tpm_mRNA.txt SRR8668652_cds_tpm_mRNA.txt SRR8668651_cds_20up_tpm_mRNA.txt SRR8668652_cds_20up_tpm_mRNA.txt SRR8668651_5pend_tpm_mRNA.txt SRR8668652_5pend_tpm_mRNA.txt SRR8668651_5pend_excl_tpm_mRNA.txt SRR8668652_5pend_excl_tpm_mRNA.txt SRR8668659_cds_tpm_mRNA.txt SRR8668660_cds_tpm_mRNA.txt SRR8668659_cds_20up_tpm_mRNA.txt SRR8668660_cds_20up_tpm_mRNA.txt SRR8668659_5pend_tpm_mRNA.txt SRR8668660_5pend_tpm_mRNA.txt SRR8668659_5pend_excl_tpm_mRNA.txt SRR8668660_5pend_excl_tpm_mRNA.txt | cut -f1-2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32 > combined_coverage_tpm_mRNA.txt
