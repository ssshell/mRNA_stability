
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


class SearchSD:
    """search SD with different degree within given sequence"""

    def __init__(self, dict_seq, list_SD):
        self.seq_ = dict_seq
        self.SD_ = list_SD
        self.SD_count_ = {}

    def search_sd(self):
        for k_gene, v_seq in self.seq_.items():
            self.SD_count_[k_gene] = {}
            for sd in self.SD_:
                self.SD_count_[k_gene][sd] = v_seq.count(sd)


def main():
    f_seq = open(sys.argv[1], 'r')
    location = '5p25ntUpStream'
    # location = sys.argv[1].split('.fa')[0].split('_')[-1]
    # print(location)
    target_seq = LoadSeqFasta(f_seq)
    target_seq.load_fa()

    # SDs are substrings of 'AGAAAGGAGGT'
    SD_id_raw = ['AGAAAGGAGGT', 'AGGA', 'AGGAG', 'AGGAGG', 'GAAAGG', 'GAGG', 'GGAG',
        'GGAGG', 'AAGGAG', 'AAGGA', 'AAAGGA', 'AAGG', 'AAAGG', 'GAAAG', 'GAAA',
        'AGGAA', 'GGAA']
    SD_id = []
    for id in SD_id_raw:
        new_id = id + '_' + location
        SD_id.append(new_id)
    # print(SD_id)

    sd = SearchSD(target_seq.target_, SD_id_raw)
    sd.search_sd()
    # print(sd.SD_count_)

    # print columns name and binary degree
    print(' ' + '\t' + '\t'.join(str(x) for x in SD_id))

    for k_gene, v_sd in sd.SD_count_.items():
        sd_freq = []
        # for k_id, v_freq in v_sd.items():
        for sd in SD_id_raw:
            sd_freq.append(v_sd[sd])
        print(k_gene + "\t" + "\t".join(str(x) for x in sd_freq))


if __name__ == "__main__":
    main()
