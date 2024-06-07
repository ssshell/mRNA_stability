#!/bin/bash

cd ../seq_processing/bowtie2_mRNA_genome/coverage_tpm

paste SRR8668651_oneThird3p_tpm_mRNA.txt SRR8668652_oneThird3p_tpm_mRNA.txt SRR8668659_oneThird3p_tpm_mRNA.txt SRR8668660_oneThird3p_tpm_mRNA.txt | cut -f1-2,4,6,8 > combined_coverage_tpm_mRNA_oneThird3p.txt

