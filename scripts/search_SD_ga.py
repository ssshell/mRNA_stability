
import sys

class LoadSeqFasta:
    """load target sequence in fasta file"""

    def __init__(self, filename_target):
        self.filename_target_ = filename_target
        self.target_ = {}
        self.minus_17to4_ = {}

    def load_fa(self):
        """filter out sequence less than 17nt"""
        for line in self.filename_target_:
            line = line.rstrip()
            if line.startswith('>'):
                id = line.lstrip('>')
                seq = ''
            else:
                seq += line
                if len(seq) >= 17:
                    self.target_[id] = seq

    def get_minus_17to4(self):
        for k_gene, v_seq in self.target_.items():
            self.minus_17to4_[k_gene] = v_seq[-17 : -3]
            # print(v_seq[-17 : -3])


class SearchSDGA:
    """search SD_ga content (G+A di-nts, GG, AA, AG, and GA)
    in SD region (-17 to -4) with 2 sequence frame.
    e.g. AGGATCTCTAAAGG,
    frame1: AG GA TC TC TA AA GG, 4 G+A di-nts.
    frame2: GG AT CT CT AA AG, 3 G+A di-nts.
    Total for both frames = 7 G+A di-nts
    """

    def __init__(self, dict_seq):
        self.seq_ = dict_seq
        self.SD_ga_percent_ = {}
        self.SD_ga_di_count_ = {}

    def get_ga_percent(self):
        for k_gene, v_seq in self.seq_.items():
            ga_num = v_seq.count('G') + v_seq.count('A')
            ga_norm = round(float(ga_num / len(v_seq)), 4)
            self.SD_ga_percent_[k_gene] = ga_norm

    def get_ga_di_nt(self):
        for k_gene, v_seq in self.seq_.items():
            count_seq = {}
            total = 0
            for i in range(2):
                temp_seq = v_seq[i: ]
                for j in range(0, len(temp_seq), 2):
                    subseq = temp_seq[j : j + 2]
                    count_seq[subseq] = count_seq.get(subseq, 0) + 1
            for di_nt in ["GG", "AA", "AG", "GA"]:
                di_nt_count = count_seq.get(di_nt, 0)
                total += int(di_nt_count)
            self.SD_ga_di_count_[k_gene] = total


def main():
    f_seq = open(sys.argv[1], 'r')
    target_seq = LoadSeqFasta(f_seq)
    target_seq.load_fa()
    target_seq.get_minus_17to4()
    # print(target_seq.minus_17to4_['MSMEG_0001'])

    sd = SearchSDGA(target_seq.minus_17to4_)
    sd.get_ga_percent()
    sd.get_ga_di_nt()

    print(" " + "\t" + "SD_GA_percent" + "\t" + "SD_GA_count")
    for k_gene, v_seq in target_seq.target_.items():
        print("%s\t%s\t%s" % (k_gene, str(sd.SD_ga_percent_[k_gene]), str(sd.SD_ga_di_count_[k_gene])))


if __name__ == "__main__":
    main()
