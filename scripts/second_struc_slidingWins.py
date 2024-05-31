
"""get the mean mfe from slicing window mfe"""

import sys
import numpy as np


class LoadMFE:

    def __init__(self, filename_mfe):
        self.filename_mfe_ = filename_mfe
        self.mfe_global_CDS_ = {}
        self.mfe_region_CDS_ = {}
        self.mfe_mean_global_CDS_ = {}
        self.mfe_mean_region_CDS_ = {}

    def get_mfe(self):
        """load mfe for both cds and cds local region (5', mid, 3')"""
        for line in self.filename_mfe_:
            line = line.rstrip().split('\t')
            target_id_raw = line[0].split('_')
            # print(target_id_raw)
            target_mfe = float(line[1])
            if len(target_id_raw) == 3:
                target_id_global_temp = target_id_raw[0]
                target_id_region_temp = target_id_raw[0] + '_' + target_id_raw[2]
            else:
                target_id_global_temp = target_id_raw[0] + '_' + target_id_raw[1]
                target_id_region_temp = target_id_raw[0] + '_' + target_id_raw[1] + '_' + target_id_raw[3]
            # print(target_id_temp)
            if target_id_global_temp not in self.mfe_global_CDS_.keys():
                self.mfe_global_CDS_[target_id_global_temp] = []
            if target_id_region_temp not in self.mfe_region_CDS_.keys():
                self.mfe_region_CDS_[target_id_region_temp] = []
            self.mfe_global_CDS_[target_id_global_temp].append(target_mfe)
            self.mfe_region_CDS_[target_id_region_temp].append(target_mfe)

    def get_mfe_mean(self):
        for k_target, v_mfe in self.mfe_global_CDS_.items():
            self.mfe_mean_global_CDS_[k_target] = round(np.mean(v_mfe), 4)
        for k_target, v_mfe in self.mfe_region_CDS_.items():
            self.mfe_mean_region_CDS_[k_target] = round(np.mean(v_mfe), 4)


def main():
    f_mfe = open(sys.argv[1], 'r')
    mfe = LoadMFE(f_mfe)

    mfe.get_mfe()
    # print(len(mfe.mfe_global_CDS_))
    # print(len(mfe.mfe_region_CDS_))
    # print(mfe.mfe_global_CDS_)
    # print(mfe.mfe_region_CDS_)

    mfe.get_mfe_mean()
    # print(len(mfe.mfe_mean_global_CDS_))
    # print(len(mfe.mfe_mean_region_CDS_))
    # print(mfe.mfe_mean_global_CDS_)
    # print(mfe.mfe_mean_region_CDS_)

    ## keep NA, for region with no subsequence starting from
    # for k_gene, v_mfe in mfe.mfe_mean_global_CDS_.items():
    #     print_all = []
    #     print_all.append(v_mfe)
    #     for region in ['5p', 'mid', '3p']:
    #         k_gene_region = k_gene + '_' + region
    #         v_mfe_region = mfe.mfe_mean_region_CDS_.get(k_gene_region, 'NA')
    #         print_all.append(v_mfe_region)
    #     print(k_gene + '\t' + '\t'.join(str(x) for x in print_all))


    ## fill NA with last available window
    for k_gene, v_mfe in mfe.mfe_mean_global_CDS_.items():
        mfe_global = v_mfe

        k_gene_5p = k_gene + '_5p'
        k_gene_mid = k_gene + '_mid'
        k_gene_3p = k_gene + '_3p'
        mfe_5p = mfe.mfe_mean_region_CDS_[k_gene_5p]

        try:
            mfe_mid = mfe.mfe_mean_region_CDS_[k_gene_mid]
        except:
            mfe_mid = mfe.mfe_region_CDS_[k_gene_5p][-1]
        try:
            mfe_3p = mfe.mfe_mean_region_CDS_[k_gene_3p]
        except:
            try:
                mfe_3p = mfe.mfe_region_CDS_[k_gene_mid][-1]
            except:
                mfe_3p = mfe.mfe_region_CDS_[k_gene_5p][-1]

        print("%s\t%s\t%s\t%s\t%s" % (k_gene, mfe_global, mfe_5p, mfe_mid, mfe_3p))


if __name__ == '__main__':
    main()
