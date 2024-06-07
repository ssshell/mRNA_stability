#!/bin/bash

cd ../seq_processing/bowtie2_ribo_genome/coverage_tpm

paste SRR8668649_oneThirdMid_tpm_ribo.txt SRR8668650_oneThirdMid_tpm_ribo.txt SRR8668657_oneThirdMid_tpm_ribo.txt SRR8668658_oneThirdMid_tpm_ribo.txt | cut -f1-2,4,6,8 > combined_coverage_tpm_ribo_oneThirdMid.txt

