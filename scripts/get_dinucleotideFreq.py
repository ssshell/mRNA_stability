
import sys
from itertools import product


class LoadSeqFasta:
    """load target sequence in fasta file"""

    def __init__(self, filename_target):
        self.filename_target_ = filename_target
        self.target_ = {}

    def load_fa(self):
        for line in self.filename_target_:
            line = line.rstrip()
            if line.startswith('>'):
                id = line.lstrip('>')
                seq = ''
            else:
                seq += line
                self.target_[id] = seq


class CountSeq:
    """count adjacent nucleotide usage of target sequence"""

    def __init__(self, dict_target):
        self.dict_target_ = dict_target
        self.target_dinucleotide_ = {}

    def get_dinucleotide(self):
        """get all possible dinucleotide combination"""
        nucleotide = ['A', 'C', 'G', 'U']
        self.dinucleotide_ = [''.join(diNt) for diNt in product(nucleotide, repeat = 2)]
        return self.dinucleotide_

    def dinucleotide_count(self, sequence):
        count_seq = {}
        count_diNt = {}
        sequence_U = sequence.replace("T", "U")
        for i in range(0, len(sequence_U) - 1):
            subseq = sequence_U[i : i + 2]

            ####### test session
            # print(subseq)

            count_seq[subseq] = count_seq.get(subseq, 0) + 1

        ####### test session
        # print(count_seq)

        for diNt in self.dinucleotide_:
            if sum(count_seq.values()) == 0:
                count_diNt[diNt] = 'NA'
            else:
                count_diNt[diNt] = round(float(count_seq.get(diNt, 0) / sum(count_seq.values())), 4)
        return count_diNt

    def dinucleotide_count_all(self):
        for k_gene, v_seq in self.dict_target_.items():
            self.target_dinucleotide_[k_gene] = self.dinucleotide_count(v_seq)


def main():
    f_seq = open(sys.argv[1], 'r')
    region = sys.argv[2]

    target_seq = LoadSeqFasta(f_seq)
    target_seq.load_fa()
    # print(target_seq.target_)

    count_dinucleotide = CountSeq(target_seq.target_)
    count_dinucleotide.get_dinucleotide()
    # print(count_dinucleotide.dinucleotide_)

    ####### test session
    # print(count_dinucleotide.dinucleotide_count(target_seq.target_['1']))

    count_dinucleotide.dinucleotide_count_all()
    for k_gene, v_count in count_dinucleotide.target_dinucleotide_.items():
        count_dinucleotide_print = []
        for k_diNt, v_percent in v_count.items():
            col_name = 'adja' + '_' + k_diNt + '_' + region
            count_dinucleotide_print.append(col_name)
    print(" " + "\t" + "\t".join(str(x) for x in count_dinucleotide_print))

    for k_gene, v_count in count_dinucleotide.target_dinucleotide_.items():
        count_percent_print = []
        for k_diNt, v_percent in v_count.items():
            count_percent_print.append(v_percent)
        print(k_gene + "\t" + "\t".join(str(x) for x in count_percent_print))


if __name__ == "__main__":
    main()
