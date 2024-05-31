
import sys
import math


codon_code_all = """UUU F      CUU L      AUU I      GUU V
UUC F      CUC L      AUC I      GUC V
UUA L      CUA L      AUA I      GUA V
UUG L      CUG L      AUG M      GUG V
UCU S      CCU P      ACU T      GCU A
UCC S      CCC P      ACC T      GCC A
UCA S      CCA P      ACA T      GCA A
UCG S      CCG P      ACG T      GCG A
UAU Y      CAU H      AAU N      GAU D
UAC Y      CAC H      AAC N      GAC D
UAA Stop   CAA Q      AAA K      GAA E
UAG Stop   CAG Q      AAG K      GAG E
UGU C      CGU R      AGU S      GGU G
UGC C      CGC R      AGC S      GGC G
UGA Stop   CGA R      AGA R      GGA G
UGG W      CGG R      AGG R      GGG G"""

codon_code_nonstop = """UUU F      CUU L      AUU I      GUU V
UUC F      CUC L      AUC I      GUC V
UUA L      CUA L      AUA I      GUA V
UUG L      CUG L      AUG M      GUG V
UCU S      CCU P      ACU T      GCU A
UCC S      CCC P      ACC T      GCC A
UCA S      CCA P      ACA T      GCA A
UCG S      CCG P      ACG T      GCG A
UAU Y      CAU H      AAU N      GAU D
UAC Y      CAC H      AAC N      GAC D
CAA Q      AAA K      GAA E
CAG Q      AAG K      GAG E
UGU C      CGU R      AGU S      GGU G
UGC C      CGC R      AGC S      GGC G
CGA R      AGA R      GGA G
UGG W      CGG R      AGG R      GGG G"""


class LoadCDSFasta:
    """load CDS sequence in fasta file"""

    def __init__(self, filename_CDS):
        self.filename_CDS_ = filename_CDS
        self.CDS_ = {}

    def load_CDS(self):
        for line in self.filename_CDS_:
            line = line.rstrip()
            if line.startswith('>'):
                gene = line.lstrip('>')
                seq = ''
            else:
                seq += line
                self.CDS_[gene] = seq


class PairCode:
    """get codon & amino acid pairs code"""

    def __init__(self):
        self.codon_all_ = []
        self.codon_nonstop_ = []
        self.codon_pairs_ = []
        self.amino_acid_code_ = {}
        self.amino_acid_pairs_ = []

    def get_codon_code(self, code_all, code_nonstop):
        codon_all = code_all.split()[0::2]
        codon_nonstop = code_nonstop.split()[0::2]
        for letter in codon_all:
            self.codon_all_.append(letter)
        for letter in codon_nonstop:
            self.codon_nonstop_.append(letter)

    def get_codon_pairs(self):
        for codon_i in self.codon_nonstop_:
            for codon_j in self.codon_all_:
                codon_pair = codon_i + '-' + codon_j
                self.codon_pairs_.append(codon_pair)

    def get_amino_acid_code(self, code_all):
        amino_acid_code = code_all.split()[1::2]
        codon_code = code_all.split()[0::2]
        for i in range(len(codon_code)):
            self.amino_acid_code_[codon_code[i]] = amino_acid_code[i]

    def get_amino_acid_pairs(self):
        for codon_pair in self.codon_pairs_:
            codon_i = codon_pair.split('-')[0]
            codon_j = codon_pair.split('-')[1]
            amino_acid_i = self.amino_acid_code_[codon_i]
            amino_acid_j = self.amino_acid_code_[codon_j]
            amino_acid_pair = amino_acid_i + '-' + amino_acid_j
            if amino_acid_pair not in self.amino_acid_pairs_:
                self.amino_acid_pairs_.append(amino_acid_pair)


class CountPairFreq:
    """count codon & amino acid pairs frequency & relative frequency for CDS sequences"""

    def __init__(self, dict_CDS, codon_pairs, amino_acid_pairs, amino_acid_code):
        self.CDS_ = dict_CDS
        self.codon_pairs_ = codon_pairs
        self.amino_acid_pairs_ = amino_acid_pairs
        self.amino_acid_code_ = amino_acid_code
        self.CDS_codon_pairs_freq_ = {}
        self.CDS_amino_acid_pairs_freq_ = {}
        self.CDS_codon_pairs_relative_freq_ = {}
        self.CDS_amino_acid_pairs_relative_freq_ = {}

    def count_pairs_single_seq(self, sequence):
        count_subseq_codon = {}
        count_subseq_amino_acid = {}
        count_codon_pairs_freq = {}
        count_amino_acid_pairs_freq = {}
        count_codon_pairs_relative_freq = {}
        count_amino_acid_pairs_relative_freq = {}
        subseq = []
        codon_pairs = []
        amino_acid_pairs = []
        sequence_U = sequence.replace("T", "U")
        for i in range(0, len(sequence_U), 3):
            subseq_i = sequence_U[i : i + 3]
            subseq.append(subseq_i)
        ####### test session
        # print(subseq)

        for j in range(len(subseq) - 1):
            codon_pairs_j = subseq[j] + '-' + subseq[j + 1]
            codon_pairs.append(codon_pairs_j)
            amino_acid_pairs_j = self.amino_acid_code_[subseq[j]] + '-' + self.amino_acid_code_[subseq[j + 1]]
            amino_acid_pairs.append(amino_acid_pairs_j)
        ####### test session
        # print(codon_pairs)
        # print(amino_acid_pairs)

        for k_codon_pair in codon_pairs:
            count_subseq_codon[k_codon_pair] = count_subseq_codon.get(k_codon_pair, 0) + 1
        for k_amino_acid_pair in amino_acid_pairs:
            count_subseq_amino_acid[k_amino_acid_pair] = count_subseq_amino_acid.get(k_amino_acid_pair, 0) + 1
        ####### test session
        # print(count_subseq_codon)
        # print(len(count_subseq_codon))
        # print(count_subseq_amino_acid)
        # print(len(count_subseq_amino_acid))

        for code in self.codon_pairs_:
            count_codon_pairs_freq[code] = count_subseq_codon.get(code, 0)
            count_codon_pairs_relative_freq[code] = round(float(count_subseq_codon.get(code, 0) / (len(codon_pairs))), 4)
        for code in self.amino_acid_pairs_:
            count_amino_acid_pairs_freq[code] = count_subseq_amino_acid.get(code, 0)
            count_amino_acid_pairs_relative_freq[code] = round(float(count_subseq_amino_acid.get(code, 0) / (len(amino_acid_pairs))), 4)
        ####### test session
        # print(count_codon_pairs_freq['GCC-AAG'])
        # print(count_amino_acid_pairs_freq['K-A'])
        # print(count_codon_pairs_relative_freq['GCC-AAG'])
        # print(count_amino_acid_pairs_relative_freq['K-A'])

        return count_codon_pairs_freq, count_amino_acid_pairs_freq, count_codon_pairs_relative_freq, count_amino_acid_pairs_relative_freq

    def count_pairs_all_CDS(self):
        for k_gene, v_seq in self.CDS_.items():
            self.CDS_codon_pairs_freq_[k_gene], self.CDS_amino_acid_pairs_freq_[k_gene], \
            self.CDS_codon_pairs_relative_freq_[k_gene], self.CDS_amino_acid_pairs_relative_freq_[k_gene] = self.count_pairs_single_seq(v_seq)


class CountSingleFreq:
    """count single codon & amino acid frequency for CDS sequences"""

    def __init__(self, dict_CDS, codon_code, amino_acid_code):
        self.CDS_ = dict_CDS
        self.codon_code_ = codon_code
        self.amino_acid_code_ = amino_acid_code
        self.CDS_codon_freq_ = {}
        self.CDS_amino_acid_freq_ = {}

    def count_single_seq(self, sequence):
        count_subseq_codon = {}
        count_subseq_amino_acid = {}
        count_codon_freq = {}
        count_amino_acid_freq = {}
        subseq = []
        sequence_U = sequence.replace("T", "U")
        for i in range(0, len(sequence_U), 3):
            subseq_i = sequence_U[i : i + 3]
            count_subseq_codon[subseq_i] = count_subseq_codon.get(subseq_i, 0) + 1
            subseq_amino_acid_i = self.amino_acid_code_[subseq_i]
            count_subseq_amino_acid[subseq_amino_acid_i] = count_subseq_amino_acid.get(subseq_amino_acid_i, 0) + 1
        ####### test session
        # print(count_subseq_codon)
        # print(count_subseq_amino_acid)

        for code_codon in self.codon_code_:
            count_codon_freq[code_codon] = count_subseq_codon.get(code_codon, 0)
        for code_amino_acid in set([val for val in self.amino_acid_code_.values()]):
            count_amino_acid_freq[code_amino_acid] = count_subseq_amino_acid.get(code_amino_acid, 0)
        ####### test session
        # print(count_codon_freq['GCC'])
        # print(count_amino_acid_freq['K'])

        return count_codon_freq, count_amino_acid_freq

    def count_single_all_CDS(self):
        for k_gene, v_seq in self.CDS_.items():
            self.CDS_codon_freq_[k_gene], self.CDS_amino_acid_freq_[k_gene] = self.count_single_seq(v_seq)


class CodonPairScore:
    """get codon-pair scores (CPS) for each codon pair"""

    def __init__(self, codon_code, amino_acid_code, codon_pair_code, amino_acid_pair_code, \
    codon_freq, amino_acid_freq, codon_pair_freq, amino_acid_pair_freq):
        self.codon_code_ = codon_code
        self.amino_acid_code_ = amino_acid_code
        self.codon_pair_code_ = codon_pair_code
        self.amino_acid_pair_code_ = amino_acid_pair_code
        self.CDS_codon_freq_ = codon_freq
        self.CDS_amino_acid_freq_ = amino_acid_freq
        self.CDS_codon_pair_freq_ = codon_pair_freq
        self.CDS_amino_acid_pair_freq_ = amino_acid_pair_freq
        self.Freq_codon_ = {}
        self.Freq_amino_acid_ = {}
        self.Freq_codon_pair_ = {}
        self.Freq_amino_acid_pair_ = {}
        self.CPS_ = {}

    def get_Freq_all_CDS(self):
        for k_gene, v_freq in self.CDS_codon_freq_.items():
            for codon in self.codon_code_:
                self.Freq_codon_[codon] = self.Freq_codon_.get(codon, 0) + int(v_freq[codon])
        for k_gene, v_freq in self.CDS_amino_acid_freq_.items():
            # print(k_gene)
            for amino_acid in set([val for val in self.amino_acid_code_.values()]):
                # print(amino_acid)
                # print(self.Freq_amino_acid_.get(amino_acid, 0) + int(v_freq[amino_acid]))
                self.Freq_amino_acid_[amino_acid] = self.Freq_amino_acid_.get(amino_acid, 0) + int(v_freq[amino_acid])
        for k_gene, v_freq in self.CDS_codon_pair_freq_.items():
            for codon_pair in self.codon_pair_code_:
                self.Freq_codon_pair_[codon_pair] = self.Freq_codon_pair_.get(codon_pair, 0) + int(v_freq[codon_pair])
        for k_gene, v_freq in self.CDS_amino_acid_pair_freq_.items():
            for amino_acid_pair in self.amino_acid_pair_code_:
                self.Freq_amino_acid_pair_[amino_acid_pair] = self.Freq_amino_acid_pair_.get(amino_acid_pair, 0) + int(v_freq[amino_acid_pair])


    def get_CPS(self):
        for codon_pair in self.codon_pair_code_:
            codon_pair_ij = codon_pair.split('-')
            amino_acid_pair = self.amino_acid_code_[codon_pair_ij[0]] + '-' + self.amino_acid_code_[codon_pair_ij[1]]
            codon_pair_freq = self.Freq_codon_pair_[codon_pair]
            codon_i_freq = self.Freq_codon_[codon_pair_ij[0]]
            codon_j_freq = self.Freq_codon_[codon_pair_ij[1]]
            amino_acid_pair_freq = self.Freq_amino_acid_pair_[amino_acid_pair]
            amino_acid_i_freq = self.Freq_amino_acid_[self.amino_acid_code_[codon_pair_ij[0]]]
            amino_acid_j_freq = self.Freq_amino_acid_[self.amino_acid_code_[codon_pair_ij[1]]]
            try:
                codon_pair_expected = (codon_i_freq / amino_acid_i_freq) * (codon_j_freq / amino_acid_j_freq) * amino_acid_pair_freq
                codon_pair_CPS = math.log(codon_pair_freq / codon_pair_expected)
            except:
                codon_pair_expected = 'NA'
                codon_pair_CPS = 'NA'
            ####### test session
            # print(codon_pair, amino_acid_pair, codon_pair_freq, amino_acid_pair_freq, codon_i_freq, codon_j_freq, codon_pair_expected, codon_pair_CPS)

            self.CPS_[codon_pair] = codon_pair_CPS


class CodonPairBias:
    """get codon-pair bias (CPB) for each CDS"""

    def __init__(self, dict_CDS, codon_pair_CPS):
        self.CDS_ = dict_CDS
        self.CPS_ = codon_pair_CPS
        self.CDS_CPB_ = {}

    def get_CPB(self):
        for k_gene, v_seq in self.CDS_.items():
            subseq = []
            codon_pairs = []
            v_seq_U = v_seq.replace("T", "U")
            CPS_v_seq = []
            for i in range(0, len(v_seq_U), 3):
                subseq_i = v_seq_U[i : i + 3]
                subseq.append(subseq_i)

            for j in range(len(subseq) - 1):
                codon_pairs_j = subseq[j] + '-' + subseq[j + 1]
                codon_pairs.append(codon_pairs_j)
            # print(codon_pairs)

            for codon_pair in codon_pairs:
                CPS_v_seq.append(self.CPS_[codon_pair])
            CDS_CPB = sum(CPS_v_seq) / len(CPS_v_seq)
            self.CDS_CPB_[k_gene] = CDS_CPB


def main():
    f_CDS = open(sys.argv[1], 'r')
    feature = sys.argv[2]


    ############## load CDS sequence
    CDS_seq = LoadCDSFasta(f_CDS)
    CDS_seq.load_CDS()
    ####### test session
    # print(CDS_seq.CDS_)
    # for gene, seq in CDS_seq.CDS_.items():
    #     if len(seq) % 3 == 0:
    #         print(gene)


    ############## get codon & amino acid pairs code
    pairs_code = PairCode()
    pairs_code.get_codon_code(codon_code_all, codon_code_nonstop)
    ####### test session
    # print(pairs_code.codon_all_)
    # print(len(pairs_code.codon_all_))
    # print(pairs_code.codon_nonstop_)
    # print(len(pairs_code.codon_nonstop_))
    # diff = [codon for codon in pairs_code.codon_all_ if codon not in pairs_code.codon_nonstop_]
    # print(diff)

    pairs_code.get_codon_pairs()
    ####### test session
    # print(pairs_code.codon_pairs_)
    # print(len(pairs_code.codon_pairs_))
    # codon_stop = [codon_pair for codon_pair in pairs_code.codon_pairs_ if 'UAG' in codon_pair]
    # codon_stop = [codon_pair for codon_pair in pairs_code.codon_pairs_ if 'UAA' in codon_pair]
    # codon_stop = [codon_pair for codon_pair in pairs_code.codon_pairs_ if 'UGA' in codon_pair]
    # print(len(codon_stop))

    pairs_code.get_amino_acid_code(codon_code_all)
    ####### test session
    # print(pairs_code.amino_acid_code_)
    # print(len(pairs_code.amino_acid_code_))
    # print(len(set([val for val in pairs_code.amino_acid_code_.values()])))

    pairs_code.get_amino_acid_pairs()
    ####### test session
    # print(pairs_code.amino_acid_pairs_)
    # print(len(pairs_code.amino_acid_pairs_))
    # print(len([pairs_code for pairs_code in pairs_code.amino_acid_pairs_ if 'Stop' in pairs_code]))


    ############## get frequency & relative frequency of codon & amino acid pairs
    pairs_code_count = CountPairFreq(CDS_seq.CDS_, pairs_code.codon_pairs_, pairs_code.amino_acid_pairs_, pairs_code.amino_acid_code_)
    ####### test session
    # pairs_code_count.count_pairs_single_seq(CDS_seq.CDS_['4'])

    pairs_code_count.count_pairs_all_CDS()
    ####### test session
    # print(len(pairs_code_count.CDS_codon_pairs_freq_))
    # print(len(pairs_code_count.CDS_amino_acid_pairs_freq_))
    # print(len(pairs_code_count.CDS_codon_pairs_relative_freq_))
    # print(len(pairs_code_count.CDS_amino_acid_pairs_relative_freq_))
    # count_dict = [pairs_code_count.CDS_codon_pairs_freq_, pairs_code_count.CDS_amino_acid_pairs_freq_, \
    #     pairs_code_count.CDS_codon_pairs_relative_freq_, pairs_code_count.CDS_amino_acid_pairs_relative_freq_]
    # for count in count_dict:
    #     for k, v in count.items():
    #         print(len(v))
    # print(pairs_code_count.CDS_codon_pairs_freq_['4']['GCC-AAG'])
    # print(pairs_code_count.CDS_amino_acid_pairs_freq_['4']['K-A'])
    # print(pairs_code_count.CDS_codon_pairs_relative_freq_['4']['GCC-AAA'])
    # print(pairs_code_count.CDS_amino_acid_pairs_relative_freq_['4']['K-A'])
    # print(pairs_code_count.CDS_amino_acid_pairs_freq_['3']['R-D'])
    # print(pairs_code_count.CDS_amino_acid_pairs_relative_freq_['3']['R-D'])
    # print(pairs_code_count.CDS_amino_acid_pairs_freq_['1']['P-V'])
    # print(pairs_code_count.CDS_amino_acid_pairs_relative_freq_['1']['P-V'])


    ############## get frequency of single codon & amino acid
    single_code_count = CountSingleFreq(CDS_seq.CDS_, pairs_code.codon_all_, pairs_code.amino_acid_code_)
    ####### test session
    # single_code_count.count_single_seq(CDS_seq.CDS_['1'])

    single_code_count.count_single_all_CDS()
    ####### test session
    # print(single_code_count.CDS_codon_freq_)
    # print(single_code_count.CDS_amino_acid_freq_)
    # print(len(single_code_count.CDS_codon_freq_))
    # print(len(single_code_count.CDS_amino_acid_freq_))
    # count_dict = [single_code_count.CDS_codon_freq_, single_code_count.CDS_amino_acid_freq_]
    # for count in count_dict:
    #     for k, v in count.items():
    #         print(len(v))
    # print(single_code_count.CDS_codon_freq_['4']['GCC'])
    # print(single_code_count.CDS_amino_acid_freq_['4']['K'])


    ############## get codon-pair scores (CPS)
    CPS = CodonPairScore(pairs_code.codon_all_, pairs_code.amino_acid_code_, pairs_code.codon_pairs_, pairs_code.amino_acid_pairs_, \
    single_code_count.CDS_codon_freq_, single_code_count.CDS_amino_acid_freq_, pairs_code_count.CDS_codon_pairs_freq_, \
    pairs_code_count.CDS_amino_acid_pairs_freq_)

    CPS.get_Freq_all_CDS()
    ####### test session
    # print(CPS.Freq_codon_)
    # print(CPS.Freq_amino_acid_)
    # print(CPS.Freq_codon_pair_['CGG-GAU'])
    # print(CPS.Freq_amino_acid_pair_)
    # print(len(CPS.Freq_codon_))
    # print(len(CPS.Freq_amino_acid_))
    # print(len(CPS.Freq_codon_pair_))
    # print(len(CPS.Freq_amino_acid_pair_))

    CPS.get_CPS()
    ####### test session
    # print(CPS.CPS_)
    # print(len(CPS.CPS_))


    ############## get codon-pair bias (CPB)
    CPB = CodonPairBias(CDS_seq.CDS_, CPS.CPS_)
    CPB.get_CPB()
    ####### test session
    # print(CPB.CDS_CPB_)
    # print(len(CPB.CDS_CPB_))


    ############## print out results
    if feature == 'relative_frequency_codon_pair':
        print(" " + "\t" + "\t".join(codon_pair for codon_pair in pairs_code.codon_pairs_))

        for k_gene, v_count in pairs_code_count.CDS_codon_pairs_relative_freq_.items():
            count_gene = []
            for codon_pair in pairs_code.codon_pairs_:
                count_gene.append(v_count[codon_pair])
            print(k_gene + "\t" + "\t".join(str(x) for x in count_gene))

    if feature == 'relative_frequency_amino_acid_pair':
        print(" " + "\t" + "\t".join(amino_acid_pair for amino_acid_pair in pairs_code.amino_acid_pairs_))

        for k_gene, v_count in pairs_code_count.CDS_amino_acid_pairs_relative_freq_.items():
            count_gene = []
            for amino_acid_pair in pairs_code.amino_acid_pairs_:
                count_gene.append(v_count[amino_acid_pair])
            print(k_gene + "\t" + "\t".join(str(x) for x in count_gene))

    if feature == 'CPB':
        for k_gene, v_CPB in CPB.CDS_CPB_.items():
            print(k_gene + "\t" + str(v_CPB))

    if feature == 'CPS':
        for k_codon_pair, v_CPs in CPS.CPS_.items():
            print(k_codon_pair + "\t" + str(v_CPs))


if __name__ == "__main__":
    main()
