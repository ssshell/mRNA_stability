
#### modify dot file from RNAfold result for RNAstructure dot2ct to get CT file

import sys

class LoadRNAfold:
    """load RNAfold output file from ViennaRNA"""

    def __init__(self, filename_fold):
        self.filename_fold_ = filename_fold
        self.gene_id_ = []
        self.sequence_ = {}
        self.structure_ = {}

    def get_split(self):
        """separate ID, sequence and structure line"""
        lines = [line.rstrip() for line in self.filename_fold_]
        # print(lines)
        self.gene_id_ = lines[0::3]
        sequences = lines[1::3]
        structures = lines[2::3]
        # print(self.gene_id_)
        # print(sequences)
        # print(structures)

        for i in range(len(self.gene_id_)):
            id = self.gene_id_[i]
            self.sequence_[id] = sequences[i]
            cut_loci_i = structures[i].rfind('(')
            structure_i = structures[i][: cut_loci_i].strip()
            # print(structure_i + 'dot')
            self.structure_[id] = structure_i
        # print(self.sequence_)
        # print(self.structure_)


def main():
    f_vienna = open(sys.argv[1], 'r')

    # extract from RNAfold output file
    rna_fold = LoadRNAfold(f_vienna)
    rna_fold.get_split()

    for gene in rna_fold.gene_id_:
        gene_i = open(gene[1:] + '.dot', 'w')
        gene_i.write("%s\n%s\n%s" % (gene, rna_fold.sequence_[gene], rna_fold.structure_[gene]))
        gene_i.close()


if __name__ == "__main__":
    main()
