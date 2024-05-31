
import sys
import statistics


class LoadRNAfold:
    """load RNAfold output file from ViennaRNA"""

    def __init__(self, filename_fold):
        self.filename_fold_ = filename_fold
        self.gene_id_ = []
        self.structure_ = []
        self.mfe_ = {}
        self.unpaired_num_ = {}

    def get_structure(self):
        """extract the structure line"""
        lines = [line.rstrip() for line in self.filename_fold_]
        self.gene_id_ = [line[1:] for line in lines[0::3]]
        self.structure_ = lines[2::3]
        # print(self.structure_)

    def get_mfe(self):
        """extract mfe from structure line"""
        for i in range(len(self.gene_id_)):
            id = self.gene_id_[i]
            start_temp = self.structure_[i].rfind('(')
            end_temp = self.structure_[i].rfind(')')
            mfe_num = float(self.structure_[i][start_temp + 1 : end_temp].strip())
            self.mfe_[id] = mfe_num

    def get_unpaired_num(self):
        """extract the unpaired nucleotide number at 5' end"""
        for i in range(len(self.gene_id_)):
            id = self.gene_id_[i]
            start_paired = self.structure_[i].find('(')
            # unpaired_nt = len(self.structure_[i][:start_paired].strip())
            unpaired_nt = len(self.structure_[i][:start_paired])
            self.unpaired_num_[id] = unpaired_nt


class LoadRNAplfold:
    """load unpaired probability of each base.
    Results from RNAplfold output file ViennaRNA. (After extract column in terminal)"""

    def __init__(self, filename_plfold):
        self.filename_plfold_ = filename_plfold
        self.gene_id_ = []
        self.unpaired_prob_raw_ = []
        self.unpaired_prob_mean_global_ = {}
        self.unpaired_prob_mean_local_ = {}
        self.first5_unpaired_ = {}
        self.first3_unpaired_ = {}
        self.first_n_all_unpaired_ = {}
        self.start_codon_unpaired_ = {}

    def get_unpaired_prob(self):
        lines = [line.rstrip() for line in self.filename_plfold_]
        self.gene_id_ = [line for line in lines[0::2]]
        self.unpaired_prob_raw_ = lines[1::2]

    def get_unpaired_prob_mean_global(self):
        for i in range(len(self.gene_id_)):
            id_i = self.gene_id_[i]
            prob_i = self.unpaired_prob_raw_[i].split(' ')
            prob_i_float = [float(i) for i in prob_i]
            self.unpaired_prob_mean_global_[id_i] = statistics.mean(prob_i_float)

    def get_unpaired_prob_mean_local(self, start, end):
        """target region is coded from 3' to 5' e.g. -6 to -16.
        if extract from 5' to 3' needs to modify [start : end]"""
        for i in range(len(self.gene_id_)):
            sum_id_i = 0.
            id_i = self.gene_id_[i]
            prob_i = self.unpaired_prob_raw_[i].split(' ')
            prob_i_float = [float(i) for i in prob_i]
            self.unpaired_prob_mean_local_[id_i] = statistics.mean(prob_i_float[start : end])

    # def get_first5_unpaired(self):
    #     """get the probability that first 5nt are all unpaired.
    #     5th row, 5th column with parameter -u 5, 5nt is 4, 3nt is 2"""
    #     for i in range(len(self.gene_id_)):
    #         id_i = self.gene_id_[i]
    #         prob_i = self.unpaired_prob_raw_[i].split(' ')[4]
    #         self.first5_unpaired_[id_i] = float(prob_i)
    #
    # def get_first3_unpaired(self):
    #     """get the probability that first 3nt are all unpaired.
    #     5th row, 5th column with parameter -u 5, 5nt is 4, 3nt is 2"""
    #     for i in range(len(self.gene_id_)):
    #         id_i = self.gene_id_[i]
    #         prob_i = self.unpaired_prob_raw_[i].split(' ')[2]
    #         self.first3_unpaired_[id_i] = float(prob_i)

    def get_first_n_all_unpaired(self, n):
        """get the probability that first N nt are all unpaired.
        5th row, 5th column with parameter -u 5, 5nt is 4, 3nt is 2"""
        for i in range(len(self.gene_id_)):
            id_i = self.gene_id_[i]
            prob_i = self.unpaired_prob_raw_[i].split(' ')[n - 1]
            # print(id_i)
            self.first_n_all_unpaired_[id_i] = float(prob_i)

    def get_start_codon_all_unpaired(self, n):
        """get the probability that the nt of start codon are all unpaired.
        3th column with parameter -u 3, use n to locate the start codon position in folding sequence"""
        for i in range(len(self.gene_id_)):
            id_i = self.gene_id_[i]
            prob_i = self.unpaired_prob_raw_[i].split(' ')[-n]
            # print(id_i)
            self.start_codon_unpaired_[id_i] = float(prob_i)


def main():
    f_vienna = open(sys.argv[1], 'r')

    if sys.argv[2] == 'RNAfold':
        # extract from RNAfold output file
        rna_fold = LoadRNAfold(f_vienna)
        rna_fold.get_structure()

        if sys.argv[3] == 'mfe':
            rna_fold.get_mfe()
            # print(rna_fold.mfe_)
            for k_gene, v_mfe in rna_fold.mfe_.items():
                print("%s\t%.2f" % (k_gene, v_mfe))
        elif sys.argv[3] == 'unpaired':
            rna_fold.get_unpaired_num()
            for k_gene, v_unpaired in rna_fold.unpaired_num_.items():
                print("%s\t%d" % (k_gene, v_unpaired))
    elif sys.argv[2] == 'RNAplfold':
        rna_plfold = LoadRNAplfold(f_vienna)
        rna_plfold.get_unpaired_prob()
        # print(rna_plfold.prob_)
        # print(len(rna_plfold.prob_))
        # print(rna_plfold.gene_id_)
        if sys.argv[3] == 'mean_global':
            rna_plfold.get_unpaired_prob_mean_global()
            for k_gene, v_prob_mean in rna_plfold.unpaired_prob_mean_global_.items():
                print("%s\t%.4f" % (k_gene, v_prob_mean))
        elif sys.argv[3] == 'mean_local':
            if sys.argv[4] == 'start_codon':
                ### folding region is 5'utr + first 18 nt of cds
                rna_plfold.get_unpaired_prob_mean_local(-18, -15)
            elif sys.argv[4] == 'sd':
                # folding region is 5'utr + first 18 nt of cds
                # -6 to -16
                rna_plfold.get_unpaired_prob_mean_local(-34, -23)
            elif sys.argv[4] == 'sd_50':
                # folding region is last 30 of 5'utr + first 20 nt of cds
                # -6 to -14
                rna_plfold.get_unpaired_prob_mean_local(-34, -25)
            elif sys.argv[4] == 'sd_33':
                # folding region is last 30 of 5'utr + first 3 nt of cds
                # -6 to -14
                rna_plfold.get_unpaired_prob_mean_local(-17, -8)
            elif sys.argv[4] == '5p_5_avg':
                rna_plfold.get_unpaired_prob_mean_local(0, 5)
            elif sys.argv[4] == '5p_3_avg':
                rna_plfold.get_unpaired_prob_mean_local(0, 3)
            for k_gene, v_prob_sum in rna_plfold.unpaired_prob_mean_local_.items():
                print("%s\t%.4f" % (k_gene, v_prob_sum))
        # elif sys.argv[3] == 'first5_unpaired':
        #     rna_plfold.get_first5_unpaired()
        #     for k_gene, v_unpaired_prob in rna_plfold.first5_unpaired_.items():
        #         print("%s\t%.4f" % (k_gene, v_unpaired_prob))
        # elif sys.argv[3] == 'first3_unpaired':
        #     rna_plfold.get_first3_unpaired()
        #     for k_gene, v_unpaired_prob in rna_plfold.first3_unpaired_.items():
        #         print("%s\t%.4f" % (k_gene, v_unpaired_prob))
        elif sys.argv[3] == '5p_5_all':
            rna_plfold.get_first_n_all_unpaired(5)
            for k_gene, v_unpaired_prob in rna_plfold.first_n_all_unpaired_.items():
                print("%s\t%.4f" % (k_gene, v_unpaired_prob))
        elif sys.argv[3] == '5p_3_all':
            rna_plfold.get_first_n_all_unpaired(3)
            for k_gene, v_unpaired_prob in rna_plfold.first_n_all_unpaired_.items():
                print("%s\t%.4f" % (k_gene, v_unpaired_prob))
        elif sys.argv[3] == 'start_codon_all':
            position = int(sys.argv[4])
            rna_plfold.get_start_codon_all_unpaired(position)
            for k_gene, v_unpaired_prob in rna_plfold.start_codon_unpaired_.items():
                print("%s\t%.4f" % (k_gene, v_unpaired_prob))


if __name__ == "__main__":
    main()
