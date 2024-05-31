
import sys

class LoadGene:
    """load gene annotation"""

    def __init__(self, filename_gene):
        self.filename_gene_ = filename_gene
        self.gene_ = []
        self.gene_start_ = {}
        self.gene_end_ = {}
        self.gene_strand_ = {}

    def get_annotation(self):
        for line in self.filename_gene_:
            line = line.rstrip().split("\t")
            gene = line[0]
            self.gene_.append(gene)
            self.gene_start_[gene] = int(line[1])
            self.gene_end_[gene] = int(line[2])
            self.gene_strand_[gene] = line[3]


class FillGene:
    """Get feature list with all the genes"""

    def __init__(self, filename_feature, list_gene):
        self.gene_ = list_gene
        self.filename_feature_ = filename_feature
        self.feature_ = {}
        self.feature_length_ = 0
        self.feature_filled_ = {}

    def get_feature(self):
        for line in self.filename_feature_:
            line = line.rstrip().split('\t')
            gene = line[0]
            self.feature_[gene] = []
            self.feature_length_ = len(line) - 1
            for i in range(1, len(line)):
                self.feature_[gene].append(line[i])

    def complete_feature(self):
        """update feature list with all the genes;
        gene with no such feature filled with NA"""
        for gene in self.gene_:
            self.feature_filled_[gene] = self.feature_.get(gene, ['NA'] * self.feature_length_)


def main():
    f_annotation = open(sys.argv[1], 'r')
    f_feature = open(sys.argv[2], 'r')

    target_gene = LoadGene(f_annotation)
    target_gene.get_annotation()

    feature = FillGene(f_feature, target_gene.gene_)
    feature.get_feature()
    feature.complete_feature()
    # print(feature.feature_filled)
    # print(len(feature.feature_filled))
    for k, v in feature.feature_filled_.items():
        print(k + "\t" + "\t".join(str(x) for x in v))


if __name__ == "__main__":
    main()
