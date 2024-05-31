
import sys


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


class SplitSeqWin:
    """split sequence into N nt window, each time with M nt overlap"""

    def __init__(self, dict_target, window_size, overlap_size):
        self.target_ = dict_target
        self.win_size_ = window_size
        self.overlap_size_ = overlap_size

    def split_seq(self):
        for target_id, target_seq in self.target_.items():
            index_count = 1
            for i in range(0, len(target_seq), self.overlap_size_):
                seq_len = len(target_seq) - i
                headline = '>' + target_id + '_' + str(index_count)
                print(headline)
                index_count += 1
                if seq_len <= self.win_size_:
                    temp_seq = target_seq[i :]
                    print(target_seq[i :])
                    break
                else:
                    temp_seq = target_seq[i : i + self.win_size_]
                    print(target_seq[i : i + self.win_size_])

    def split_seq(self):
        for target_id, target_seq in self.target_.items():
            index_count = 1
            for i in range(0, len(target_seq), self.overlap_size_):
                seq_len = len(target_seq) - i
                five2mid = (len(target_seq) / 3) - 1
                mid2three = ((len(target_seq) / 3) * 2) - 1
                if i < five2mid:
                    region = '5p'
                elif i >= five2mid and i < mid2three:
                    region = 'mid'
                else:
                    region = '3p'
                headline = '>' + target_id + '_' + str(index_count) + '_' + region
                print(headline)
                index_count += 1
                if seq_len <= self.win_size_:
                    temp_seq = target_seq[i :]
                    print(target_seq[i :])
                    break
                else:
                    temp_seq = target_seq[i : i + self.win_size_]
                    print(target_seq[i : i + self.win_size_])



def main():
    f_seq = open(sys.argv[1], 'r')
    target_seq = LoadSeqFasta(f_seq)
    target_seq.load_fa()

    window_size = int(sys.argv[2])
    overlap_size = int(sys.argv[3])
    split_seq = SplitSeqWin(target_seq.target_, window_size, overlap_size)
    split_seq.split_seq()


if __name__ == "__main__":
    main()
