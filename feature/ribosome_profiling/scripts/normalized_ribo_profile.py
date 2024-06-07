
import sys
import numpy as np


class LoadTPM:
    """get raw TPM of 4 regions: cds, cds_20up, cds_5pend, cds_5pend_excl.
    Each with 2 replicates"""

    def __init__(self, filename_TPM_ribo, filename_TPM_totalmRNA):
        self.filename_TPM_ribo = filename_TPM_ribo
        self.filename_TPM_totalmRNA = filename_TPM_totalmRNA

    def get_coverage(self, file):
        dict_coverage = {}
        for line in file:
            line = line.rstrip().split('\t')
            dict_coverage[line[0]] = []
            for i in range(1, len(line)):
                dict_coverage[line[0]].append(float(line[i]) + 1)
        return dict_coverage


class AvgRep:
    """Take the average of replicates for each strain"""

    def __init__(self, dict_ribo, dict_totalmRNA, num_strain, num_replicate, num_region):
        self.ribo_ = dict_ribo
        self.totalmRNA_ = dict_totalmRNA
        self.strain_ = num_strain
        self.num_replicate_ = num_replicate
        self.num_region_ = num_region

    def get_avg_coverage(self, dict_coverage):
        dict_coverage_avg_rep = {}
        total_sample = self.strain_ * self.num_replicate_ * self.num_region_
        for k_gene, v_coverage in dict_coverage.items():
            dict_coverage_avg_rep[k_gene] = []
            for i in range(0, total_sample, self.num_replicate_):
                dict_coverage_avg_rep[k_gene].append(np.mean(v_coverage[i : i + self.num_replicate_]))
        return dict_coverage_avg_rep


def main():

    f_ribo = open(sys.argv[1], "r")
    f_totalmRNA = open(sys.argv[2], "r")
    coverage = LoadTPM(f_ribo, f_totalmRNA)
    coverage_ribo = coverage.get_coverage(coverage.filename_TPM_ribo)
    # print(coverage_ribo["MSMEG_0001"])
    coverage_totalmRNA = coverage.get_coverage(coverage.filename_TPM_totalmRNA)
    # print(coverage_totalmRNA["MSMEG_0001"])

    num_replicate = 2
    num_region = 4
    num_strain = 2
    average_replicate = AvgRep(coverage_ribo, coverage_totalmRNA, num_strain, num_replicate, num_region)
    average_replicate_ribo = average_replicate.get_avg_coverage(average_replicate.ribo_)
    # print(average_replicate_ribo["MSMEG_0001"])
    average_replicate_totalmRNA = average_replicate.get_avg_coverage(average_replicate.totalmRNA_)
    # print(average_replicate_totalmRNA["MSMEG_0001"])

    ### get TPMTotalribo/TPMTotalmRNA for each strain
    ### final go without log2
    num_col = num_strain * num_region
    ribo_total = {}
    ribo_total_avg = {}
    for k_gene, v_coverage in average_replicate_ribo.items():
        ribo_total[k_gene] = []
        ribo_total_avg[k_gene] = []
        for i in range(num_col):
            ribo_temp = average_replicate_ribo[k_gene][i] / average_replicate_totalmRNA[k_gene][i]
            ribo_total[k_gene].append(ribo_temp)
        for j in range(num_region):
            ribo_temp_avg = (ribo_total[k_gene][j] + ribo_total[k_gene][j + 4]) / 2
            ribo_total_avg[k_gene].append(round(ribo_temp_avg, 4))
        print(k_gene + "\t" + "\t".join(str(x) for x in ribo_total_avg[k_gene]))

    # print(ribo_total["MSMEG_0001"])
    # print(ribo_total_avg["MSMEG_0001"])


if __name__ == "__main__":
    main()
