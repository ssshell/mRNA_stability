#!/bin/bash

cd ../seq_processing/bowtie2_mRNA_genome/coverage_tpm

paste SRR8668651_3pLast20percent_tpm_mRNA.txt SRR8668652_3pLast20percent_tpm_mRNA.txt SRR8668659_3pLast20percent_tpm_mRNA.txt SRR8668660_3pLast20percent_tpm_mRNA.txt | cut -f1-2,4,6,8 > combined_coverage_tpm_mRNA_3pLast20percent.txt 
