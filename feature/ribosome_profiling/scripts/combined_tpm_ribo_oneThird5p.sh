#!/bin/bash

cd ../seq_processing/bowtie2_ribo_genome/coverage_tpm

paste SRR8668649_oneThird5p_tpm_ribo.txt SRR8668650_oneThird5p_tpm_ribo.txt SRR8668657_oneThird5p_tpm_ribo.txt SRR8668658_oneThird5p_tpm_ribo.txt | cut -f1-2,4,6,8 > combined_coverage_tpm_ribo_oneThird5p.txt

