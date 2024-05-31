
import sys
import os


class LoadUTRlenght:
    """load gene 5'UTR length"""

    def __init__(self, filename_fUTR_length):
        self.filename_fUTR_length_ = filename_fUTR_length
        self.gene_fUTR_length_ = {}

    def get_fUTR_length(self):
        next(self.filename_fUTR_length_)
        for line in self.filename_fUTR_length_:
            line = line.rstrip().split("\t")
            gene = line[0]
            self.gene_fUTR_length_[gene] = line[1]


class LoadRNAfold:
    """load deltaG_mRNA RNAfold output file"""

    def __init__(self, filename_fold):
        self.filename_fold_ = filename_fold
        self.gene_id_ = []
        self.structure_ = {}

    def get_structure(self):
        """extract the structure line"""
        lines = [line.rstrip() for line in self.filename_fold_]
        self.gene_id_ = [line[1:] for line in lines[0::3]]
        structures = lines[2::3]
        for i in range(len(self.gene_id_)):
            id = self.gene_id_[i]
            cut_loci_i = structures[i].rfind('(')
            structure_i = structures[i][: cut_loci_i].strip()
            # print(structure_i + 'dot')
            self.structure_[id] = structure_i


def main():
    f_fUTR_length = open(sys.argv[1], 'r')
    f_deltaG_mRNA = open(sys.argv[2], 'r')

    fUTR = LoadUTRlenght(f_fUTR_length)
    fUTR.get_fUTR_length()
    # print(fUTR.gene_fUTR_length_)

    deltaG_mRNA = LoadRNAfold(f_deltaG_mRNA)
    deltaG_mRNA.get_structure()
    # print(deltaG_mRNA.structure_)
    # print(deltaG_mRNA.structure_['MSMEG_6947'][0:5])

    for filename in os.listdir('../feature/translation_init/deltaG_mRNA_ctFiles'):
        if filename.endswith('.ct'):
            with open(os.path.join('../feature/translation_init/deltaG_mRNA_ctFiles', filename)) as f_ct:
                f_ct_index = {}
                for line in f_ct:
                    line = line.rstrip().lstrip().split(" ")
                    index = [x for x in line if x]
                    # print(index)
                    if len(index) == 2:
                        gene = index[1]
                        f_ct_index[gene] = []
                        # print(gene)
                    else:
                        f_ct_index[gene].append(index)
                        # print(index)
                # print(f_ct_index)
                # print(fUTR.gene_fUTR_length_[gene])
                modified_dot = {}
                for i in range(1, 51):
                     modified_dot[i] = []
                     modified_dot[i].append(f_ct_index[gene][i - 1][1])
                     modified_dot[i].append('.')
                # print(modified_dot)

                #### first version (only break base-pairing within regions;
                #### and didn't keep predicted structure outside regions)
                # if fUTR.gene_fUTR_length_[gene] != 'NA':
                #     if int(fUTR.gene_fUTR_length_[gene]) < 12:
                #         for j in range(1, 26):
                #             if f_ct_index[gene][j - 1][4] != '0':
                #                 # print(j, f_ct_index[gene][j - 1][4])
                #                 modified_dot[j][1] = 'x'
                #                 modified_dot[int(f_ct_index[gene][j - 1][4])][1] = 'x'
                #         # print(modified_dot)
                #     else:
                #         for j in range(14, 39):
                #             if f_ct_index[gene][j - 1][4] != '0':
                #                 # print(j, f_ct_index[gene][j - 1][4])
                #                 modified_dot[j][1] = 'x'
                #                 modified_dot[int(f_ct_index[gene][j - 1][4])][1] = 'x'
                #     # print(modifed_dot)
                #     dot_constraints = []
                #     dot_sequence = []
                #     for i in range(1, 51):
                #         dot_constraints.append(modified_dot[i][1])
                #         dot_sequence.append(modified_dot[i][0])
                #     print('>' + gene + '\n' + ''.join(x for x in dot_sequence) + '\n' +
                #     ''.join(x for x in dot_constraints))

                if int(fUTR.gene_fUTR_length_[gene]) == 0 or int(fUTR.gene_fUTR_length_[gene]) >= 12:
                    if int(fUTR.gene_fUTR_length_[gene]) == 0:
                        for j in range(1, 51):
                            if j in range(1, 14):
                                modified_dot[j][1] = 'x'
                                if f_ct_index[gene][j - 1][4] != '0':
                                    modified_dot[int(f_ct_index[gene][j - 1][4])][1] = 'x'
                            else:
                                if f_ct_index[gene][j - 1][4] != '0':
                                    if int(f_ct_index[gene][j - 1][4]) not in range(1, 14):
                                        modified_dot[j][1] = deltaG_mRNA.structure_[gene][j - 1]
                        # print(modified_dot)
                    elif int(fUTR.gene_fUTR_length_[gene]) >= 25:
                        for j in range(1, 51):
                            if j in range(14, 39):
                                modified_dot[j][1] = 'x'
                                if f_ct_index[gene][j - 1][4] != '0':
                                    modified_dot[int(f_ct_index[gene][j - 1][4])][1] = 'x'
                            else:
                                if f_ct_index[gene][j - 1][4] != '0':
                                    if int(f_ct_index[gene][j - 1][4]) not in range(14, 39):
                                        modified_dot[j][1] = deltaG_mRNA.structure_[gene][j - 1]
                    else:
                        fUTR_len = int(fUTR.gene_fUTR_length_[gene])
                        for j in range(1, 51):
                            if j in range(fUTR_len - 11, fUTR_len + 14):
                                modified_dot[j][1] = 'x'
                                if f_ct_index[gene][j - 1][4] != '0':
                                    modified_dot[int(f_ct_index[gene][j - 1][4])][1] = 'x'
                            else:
                                if f_ct_index[gene][j - 1][4] != '0':
                                    if int(f_ct_index[gene][j - 1][4]) not in range(fUTR_len - 11, fUTR_len + 14):
                                        modified_dot[j][1] = deltaG_mRNA.structure_[gene][j - 1]
                    # print(modified_dot)
                    dot_constraints = []
                    dot_sequence = []
                    for i in range(1, 51):
                        dot_constraints.append(modified_dot[i][1])
                        dot_sequence.append(modified_dot[i][0])
                    print('>' + gene + '\n' + ''.join(x for x in dot_sequence) + '\n' +
                    ''.join(x for x in dot_constraints))


if __name__ == "__main__":
    main()
