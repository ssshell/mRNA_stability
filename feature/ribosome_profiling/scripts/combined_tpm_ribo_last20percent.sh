#!/bin/bash

cd ../seq_processing/bowtie2_ribo_genome/coverage_tpm

paste SRR8668649_3pLast20percent_tpm_ribo.txt SRR8668650_3pLast20percent_tpm_ribo.txt SRR8668657_3pLast20percent_tpm_ribo.txt SRR8668658_3pLast20percent_tpm_ribo.txt | cut -f1-2,4,6,8 > combined_coverage_tpm_ribo_3pLast20percent.txt

