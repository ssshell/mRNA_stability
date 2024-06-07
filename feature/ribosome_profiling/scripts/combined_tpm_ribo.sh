#!/bin/bash

cd ../seq_processing/bowtie2_ribo_genome/coverage_tpm

paste SRR8668649_cds_tpm_ribo.txt SRR8668650_cds_tpm_ribo.txt SRR8668649_cds_20up_tpm_ribo.txt SRR8668650_cds_20up_tpm_ribo.txt SRR8668649_5pend_tpm_ribo.txt SRR8668650_5pend_tpm_ribo.txt SRR8668649_5pend_excl_tpm_ribo.txt SRR8668650_5pend_excl_tpm_ribo.txt SRR8668657_cds_tpm_ribo.txt SRR8668658_cds_tpm_ribo.txt SRR8668657_cds_20up_tpm_ribo.txt SRR8668658_cds_20up_tpm_ribo.txt SRR8668657_5pend_tpm_ribo.txt SRR8668658_5pend_tpm_ribo.txt SRR8668657_5pend_excl_tpm_ribo.txt SRR8668658_5pend_excl_tpm_ribo.txt | cut -f1-2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32 > combined_coverage_tpm_ribo.txt
