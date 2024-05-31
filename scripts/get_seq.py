
import sys


class LoadGenome:
    """load genome sequence"""

    def __init__(self, filename_genome):
        self.filename_genome_ = filename_genome
        self.genome_ = ''

    def get_seq(self):
        sequence = ''
        next(self.filename_genome_)
        for line in self.filename_genome_:
            line = line.rstrip()
            sequence += line
        self.genome_ = sequence


class LoadGene:
    """load gene annotation and tss"""

    def __init__(self, filename_gene, filename_tss):
        self.filename_gene_ = filename_gene
        self.filename_tss_ = filename_tss
        self.gene_ = []
        self.gene_start_ = {}
        self.gene_end_ = {}
        self.gene_strand_ = {}
        self.gene_tss_ = {}

    def get_annotation(self):
        for line in self.filename_gene_:
            line = line.rstrip().split("\t")
            gene = line[0]
            self.gene_.append(gene)
            self.gene_start_[gene] = int(line[1])
            self.gene_end_[gene] = int(line[2])
            self.gene_strand_[gene] = line[3]

    def get_tss(self):
        for line in self.filename_tss_:
            line = line.rstrip().split("\t")
            gene = line[0]
            self.gene_tss_[gene] = line[1]


class ExtractSeq:
    """extract target sequence from genome based on annotation"""

    def __init__(self, genome):
        self.genome_ = genome

    def reverse_complement(self, sequence):
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        return ''.join([complement[base] for base in sequence[::-1]])

    def extract_seq(self, gene, start, end, strand):
        """start and end are 1 based coordinate"""
        if strand == '+':
            result = self.genome_[start - 1 : end]
        else:
            result = self.reverse_complement(self.genome_[start - 1 : end])
        return result


def main():
    f_genome = open(sys.argv[1], 'r')
    f_annotation = open(sys.argv[2], 'r')
    f_tss = open(sys.argv[3], 'r')
    target_genome = LoadGenome(f_genome)
    target_genome.get_seq()
    # print(len(target_genome.genome_))

    # check sequence example
    # print(target_genome.genome_[14110:14130])

    target_gene = LoadGene(f_annotation, f_tss)
    target_gene.get_annotation()
    target_gene.get_tss()
    # print(len(target_gene.gene_tss_))

    extract_seq = ExtractSeq(target_genome.genome_)

    # get 5'UTR sequence
    target_seq_5utr = {}
    target_seq_5utr_length = {}

    # get 3'UTR sequence, N nt after stop codon
    target_seq_3utr_20nt = {}
    target_seq_3utr_50nt = {}
    target_seq_3utr_60nt = {}

    # get sequence between start and end position & length
    target_seq_gene = {}
    target_seq_gene_length = {}

    # get nonstop sequence (between start and end position excluding last 3nts)
    target_seq_gene_nstop = {}

    # get sequence between tss and end position
    target_seq_tss_gene = {}
    # get 5'/first 20 nt of the transcript for genes with known tss
    target_seq_tss_gene_5p20 = {}

    # 5'UTR + 18 cds
    target_seq_5utr_18cds = {}

    # 5' UTR with at least 30 nt + 20 cds
    target_seq_5utr_last30_20cds = {}

    # 5' UTR with at least 30 nt + start codon
    target_seq_5utr_last30_start_codon = {}

    # get 30 bases upstream of start position of gene
    target_seq_30up = {}

    # get 25 bases upstream of start position of gene
    target_seq_25up = {}

    # get N bases upstream of start position of gene
    target_seq_Nup = {}

    # get first 18 bases
    target_seq_gene_5p18 = {}

    # get last 18 bases
    target_seq_gene_3p18 = {}

    # get last 21 bases
    target_seq_gene_3p21 = {}

    # get last 18 nonstop bases
    target_seq_gene_3p18_nstop = {}

    # get TIR (50nt), defined as:
    # 5'UTR >= 25nt:  25nt upstream of start codon and 25nt downstream from start codon
    # 5'UTR < 25nt: 50nt from TSS
    target_seq_TIR = {}

    for gene in target_gene.gene_:
        gene_seq = extract_seq.extract_seq(gene, target_gene.gene_start_[gene],
        target_gene.gene_end_[gene], target_gene.gene_strand_[gene])
        target_seq_gene[gene] = gene_seq
        target_seq_gene_nstop[gene] = gene_seq[: -3]
        target_seq_gene_length[gene] = len(gene_seq)

        target_seq_gene_5p18[gene] = gene_seq[: 18]
        target_seq_gene_3p18[gene] = gene_seq[-18 :]
        target_seq_gene_3p21[gene] = gene_seq[-21 :]
        target_seq_gene_3p18_nstop[gene] = gene_seq[-21 : -3]


        if target_gene.gene_strand_[gene] == '+':
            utr_seq_3p_20 = extract_seq.extract_seq(gene, int(target_gene.gene_end_[gene]) + 1,
            int(target_gene.gene_end_[gene]) + 20, target_gene.gene_strand_[gene])

            utr_seq_3p_50 = extract_seq.extract_seq(gene, int(target_gene.gene_end_[gene]) + 1,
            int(target_gene.gene_end_[gene]) + 50, target_gene.gene_strand_[gene])

            utr_seq_3p_60 = extract_seq.extract_seq(gene, int(target_gene.gene_end_[gene]) + 1,
            int(target_gene.gene_end_[gene]) + 60, target_gene.gene_strand_[gene])
        else:
            utr_seq_3p_20 = extract_seq.extract_seq(gene, int(target_gene.gene_start_[gene]) - 20,
            int(target_gene.gene_start_[gene]) - 1, target_gene.gene_strand_[gene])

            utr_seq_3p_50 = extract_seq.extract_seq(gene, int(target_gene.gene_start_[gene]) - 50,
            int(target_gene.gene_start_[gene]) - 1, target_gene.gene_strand_[gene])

            utr_seq_3p_60 = extract_seq.extract_seq(gene, int(target_gene.gene_start_[gene]) - 60,
            int(target_gene.gene_start_[gene]) - 1, target_gene.gene_strand_[gene])
        target_seq_3utr_20nt[gene] = utr_seq_3p_20
        target_seq_3utr_50nt[gene] = utr_seq_3p_50
        target_seq_3utr_60nt[gene] = utr_seq_3p_60

        # get utr related sequence
        if target_gene.gene_tss_[gene] != 'NA':
            if target_gene.gene_strand_[gene] == '+':
                utr_seq = extract_seq.extract_seq(gene, int(target_gene.gene_tss_[gene]),
                int(target_gene.gene_start_[gene]) - 1, target_gene.gene_strand_[gene])

                if len(utr_seq) >= 30:
                    utr_seq_last30_plus_20 = extract_seq.extract_seq(gene, int(target_gene.gene_start_[gene]) - 30,
                    int(target_gene.gene_start_[gene]) - 1 + 20, target_gene.gene_strand_[gene])

                    utr_seq_last30_plus_start_codon = extract_seq.extract_seq(gene, int(target_gene.gene_start_[gene]) - 30,
                    int(target_gene.gene_start_[gene]) - 1 + 3, target_gene.gene_strand_[gene])
                else:
                    utr_seq_last30_plus_20 = ''
                    utr_seq_last30_plus_start_codon = ''

                utr_plus_18_seq = extract_seq.extract_seq(gene, int(target_gene.gene_tss_[gene]),
                int(target_gene.gene_start_[gene]) - 1 + 18, target_gene.gene_strand_[gene])

                tss_gene_seq = extract_seq.extract_seq(gene, int(target_gene.gene_tss_[gene]),
                int(target_gene.gene_end_[gene]), target_gene.gene_strand_[gene])

                if len(utr_seq) >= 25:
                    TIR_seq = extract_seq.extract_seq(gene, int(target_gene.gene_start_[gene]) - 25,
                    int(target_gene.gene_start_[gene]) - 1 + 25, target_gene.gene_strand_[gene])
                else:
                    TIR_seq = extract_seq.extract_seq(gene, int(target_gene.gene_tss_[gene]),
                    int(target_gene.gene_tss_[gene]) - 1 + 50, target_gene.gene_strand_[gene])

            else:
                utr_seq = extract_seq.extract_seq(gene, int(target_gene.gene_end_[gene]) + 1,
                int(target_gene.gene_tss_[gene]), target_gene.gene_strand_[gene])

                if len(utr_seq) >= 30:
                    utr_seq_last30_plus_20 = extract_seq.extract_seq(gene, int(target_gene.gene_end_[gene]) + 1 - 20,
                    target_gene.gene_end_[gene] + 30, target_gene.gene_strand_[gene])

                    utr_seq_last30_plus_start_codon = extract_seq.extract_seq(gene, int(target_gene.gene_end_[gene]) + 1 - 3,
                    target_gene.gene_end_[gene] + 30, target_gene.gene_strand_[gene])
                else:
                    utr_seq_last30_plus_20 = ''
                    utr_seq_last30_plus_start_codon = ''

                utr_plus_18_seq = extract_seq.extract_seq(gene, int(target_gene.gene_end_[gene]) + 1 - 18,
                int(target_gene.gene_tss_[gene]), target_gene.gene_strand_[gene])

                tss_gene_seq = extract_seq.extract_seq(gene, int(target_gene.gene_start_[gene]),
                int(target_gene.gene_tss_[gene]), target_gene.gene_strand_[gene])

                if len(utr_seq) >= 25:
                    TIR_seq = extract_seq.extract_seq(gene, int(target_gene.gene_end_[gene]) + 1 - 25,
                    target_gene.gene_end_[gene] + 25, target_gene.gene_strand_[gene])
                else:
                    TIR_seq = extract_seq.extract_seq(gene, int(target_gene.gene_tss_[gene]) + 1 - 50,
                    int(target_gene.gene_tss_[gene]), target_gene.gene_strand_[gene])

            target_seq_5utr[gene] = utr_seq
            target_seq_5utr_length[gene] = len(utr_seq)
            target_seq_5utr_18cds[gene] = utr_plus_18_seq
            target_seq_tss_gene[gene] = tss_gene_seq
            target_seq_tss_gene_5p20[gene] = tss_gene_seq[: 20]
            target_seq_5utr_last30_20cds[gene] = utr_seq_last30_plus_20
            target_seq_5utr_last30_start_codon[gene] = utr_seq_last30_plus_start_codon
            target_seq_TIR[gene] = TIR_seq

        else:
            target_seq_5utr[gene] = 'NA'
            target_seq_5utr_length[gene] = 'NA'
            target_seq_5utr_18cds[gene] = 'NA'
            target_seq_tss_gene[gene] = 'NA'
            target_seq_tss_gene_5p20[gene] = 'NA'
            target_seq_5utr_last30_20cds[gene] = 'NA'
            target_seq_5utr_last30_start_codon[gene] = 'NA'


        # get 30 bases upstream of gene start position
        if target_gene.gene_strand_[gene] == '+':
            seq_30 = extract_seq.extract_seq(gene, target_gene.gene_start_[gene] - 30,
            target_gene.gene_start_[gene] - 1, target_gene.gene_strand_[gene])
        else:
            seq_30 = extract_seq.extract_seq(gene, target_gene.gene_end_[gene] + 1,
            target_gene.gene_end_[gene] + 30, target_gene.gene_strand_[gene])
        target_seq_30up[gene] = seq_30

        # get N bases upstream of gene start position
        N = 25
        if target_gene.gene_strand_[gene] == '+':
            seq_N = extract_seq.extract_seq(gene, target_gene.gene_start_[gene] - N,
            target_gene.gene_start_[gene] - 1, target_gene.gene_strand_[gene])
        else:
            seq_N = extract_seq.extract_seq(gene, target_gene.gene_end_[gene] + 1,
            target_gene.gene_end_[gene] + N, target_gene.gene_strand_[gene])
        target_seq_Nup[gene] = seq_N

    if sys.argv[4] == 'gene':
        for gene in target_gene.gene_:
            # fasta format
            print('>' + gene + '\n' + target_seq_gene[gene])
            # to check start codon
            # print('>' + gene + '\n' + target_seq_gene[gene][: 3])

            # bedfile format
            # print("%s\t%s\t%s\t%s\t%s" % (gene, target_gene.gene_start_[gene], target_gene.gene_end_[gene],
            #     target_gene.gene_strand_[gene], target_seq_gene[gene]))
    if sys.argv[4] == 'gene_nstop':
        for gene in target_gene.gene_:
            # fasta format
            print('>' + gene + '\n' + target_seq_gene_nstop[gene])
    elif sys.argv[4] == 'tss_gene':
        # get sequence from tss to the end of gene for leadered genes
        for gene in target_gene.gene_:
            if target_gene.gene_tss_[gene] != 'NA' and target_seq_5utr[gene] != '':
                print('>' + gene + '\n' + target_seq_tss_gene[gene])
    elif sys.argv[4] == 'tss_gene_5p20':
        # get 5'/first 20 nt of the transcript for genes with known tss
        # either 5' of UTR (leadered) or CDS (leaderless) or maybe UTR + CDS if UTR shorter than 20 nt
        for gene in target_gene.gene_:
            if target_gene.gene_tss_[gene] != 'NA':
                print('>' + gene + '\n' + target_seq_tss_gene_5p20[gene])
    elif sys.argv[4] == 'utr':
        for gene in target_gene.gene_:
            if target_gene.gene_tss_[gene] != 'NA' and target_seq_5utr[gene] != '':
                # fasta format
                print('>' + gene + '\n' + target_seq_5utr[gene])
    elif sys.argv[4] == 'utr_exclude_last15':
        for gene in target_gene.gene_:
            if target_gene.gene_tss_[gene] != 'NA' and target_seq_5utr[gene] != '' and target_seq_5utr_length[gene] >= 35:
                # fasta format
                print('>' + gene + '\n' + target_seq_5utr[gene][: -15])
    elif sys.argv[4] == 'utr_18cds':
        for gene in target_gene.gene_:
            if target_gene.gene_tss_[gene] != 'NA' and target_seq_5utr[gene] != '':
                # fasta format
                print('>' + gene + '\n' + target_seq_5utr_18cds[gene])
    elif sys.argv[4] == 'utr_last30_20cds':
        for gene in target_gene.gene_:
            if target_gene.gene_tss_[gene] != 'NA' and target_seq_5utr[gene] != '' and target_seq_5utr_last30_20cds[gene] != '':
                # fasta format
                print('>' + gene + '\n' + target_seq_5utr_last30_20cds[gene])
    elif sys.argv[4] == 'utr_last30_start_codon':
        for gene in target_gene.gene_:
            if target_gene.gene_tss_[gene] != 'NA' and target_seq_5utr[gene] != '' and target_seq_5utr_last30_start_codon[gene] != '':
                # fasta format
                print('>' + gene + '\n' + target_seq_5utr_last30_start_codon[gene])
    elif sys.argv[4] == 'gene_length':
        print(" " + "\t" + 'gene_length')
        for gene in target_gene.gene_:
            print(gene + '\t' + str(target_seq_gene_length[gene]))
    elif sys.argv[4] == 'utr_length':
        print(" " + "\t" + 'fpr_utr_length')
        for gene in target_gene.gene_:
            print(gene + '\t' + str(target_seq_5utr_length[gene]))
    elif sys.argv[4] == 'gene_5p18':
        for gene in target_gene.gene_:
            print('>' + gene + '\n' + target_seq_gene_5p18[gene])
    elif sys.argv[4] == 'gene_3p18':
        for gene in target_gene.gene_:
            print('>' + gene + '\n' + target_seq_gene_3p18[gene])
    elif sys.argv[4] == 'gene_3p21':
        for gene in target_gene.gene_:
            print('>' + gene + '\n' + target_seq_gene_3p21[gene])
    elif sys.argv[4] == 'gene_3p18_nstop':
        for gene in target_gene.gene_:
            print('>' + gene + '\n' + target_seq_gene_3p18_nstop[gene])
    elif sys.argv[4] == 'gene_up30':
        for gene in target_gene.gene_:
            # bedfile format
            print("%s\t%s\t%s\t%s\t%s\t%s" % (gene, target_gene.gene_start_[gene], target_gene.gene_end_[gene],
                target_gene.gene_strand_[gene], target_seq_30up[gene], target_seq_gene[gene]))
    elif sys.argv[4] == 'gene_upN':
        for gene in target_gene.gene_:
            print('>' + gene + '\n' + target_seq_Nup[gene])
            # print("%s\t%s\t%s\t%s\t%s" % (gene, target_gene.gene_start_[gene], target_gene.gene_end_[gene],
            #     target_gene.gene_strand_[gene], target_seq_Nup[gene]))
    elif sys.argv[4] == 'gene_shift':
        n_shift = int(sys.argv[5])
        for gene in target_gene.gene_:
            # fasta format
            print('>' + gene + '\n' + target_seq_gene_nstop[gene][n_shift :])
    elif sys.argv[4] == 'gene_shift_5p18':
        n_shift = int(sys.argv[5])
        for gene in target_gene.gene_:
            # fasta format
            print('>' + gene + '\n' + target_seq_gene[gene][n_shift : n_shift + 18])
    # elif sys.argv[4] == 'gene_excl_18':
    #     n_shift = int(sys.argv[5])
    #     window = n_shift * 18
    #     for gene in target_gene.gene_:
    #         # fasta format
    #         print('>' + gene + '\n' + target_seq_gene[gene][window :])
    elif sys.argv[4] == 'gene_excl_18_5p18':
        n_shift = int(sys.argv[5])
        window = n_shift * 18
        for gene in target_gene.gene_:
            # fasta format
            print('>' + gene + '\n' + target_seq_gene[gene][window : window + 18])

    elif sys.argv[4] == 'utr_3p_20':
        for gene in target_gene.gene_:
            print('>' + gene + '\n' + target_seq_3utr_20nt[gene])

    elif sys.argv[4] == 'utr_3p_50':
        for gene in target_gene.gene_:
            print('>' + gene + '\n' + target_seq_3utr_50nt[gene])

    elif sys.argv[4] == 'utr_3p_60':
        for gene in target_gene.gene_:
            print('>' + gene + '\n' + target_seq_3utr_60nt[gene])

    elif sys.argv[4] == 'TIR':
        for gene in target_gene.gene_:
            if target_gene.gene_tss_[gene] != 'NA':
                # fasta format
                print('>' + gene + '\n' + target_seq_TIR[gene])


if __name__ == "__main__":
    main()
