#!/bin/bash

cd ../seq_processing/bowtie2_mRNA_genome/coverage_tpm

paste SRR8668651_oneThirdMid_tpm_mRNA.txt SRR8668652_oneThirdMid_tpm_mRNA.txt SRR8668659_oneThirdMid_tpm_mRNA.txt SRR8668660_oneThirdMid_tpm_mRNA.txt | cut -f1-2,4,6,8 > combined_coverage_tpm_mRNA_oneThirdMid.txt

