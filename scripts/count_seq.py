
import sys

codon_code_string = """UUU F      CUU L      AUU I      GUU V
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

codon_code_string_nstop = """UUU F      CUU L      AUU I      GUU V
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
    """count usage (nucleotide, codon, GC, amino acid) of target sequence"""

    def __init__(self, dict_target):
        self.dict_target_ = dict_target
        self.codon_ = []
        self.amino_ = {}
        self.target_base_ = {}
        self.target_codon_ = {}
        self.target_gc_ = {}
        self.target_amino_ = {}

    def base_count(self, sequence):
        raw = {}
        count_seq = {}
        for base in sequence:
            raw[base] = raw.get(base, 0) + 1
        for base, count in raw.items():
            count_seq[base] = round(float(raw.get(base) / len(sequence)), 4)
        return count_seq

    def base_count_all(self):
        for k_gene, v_seq in self.dict_target_.items():
            self.target_base_[k_gene] = self.base_count(v_seq)

    def get_codon_code(self, code):
        codon_code = code.split()[0::2]
        for letter in codon_code:
            self.codon_.append(letter)
        # print(self.codon_)
        # print(len(self.codon_))

    def codon_count(self, sequence):
        count_seq = {}
        count_codon = {}
        for i in range(0, len(sequence), 3):
            subseq = sequence[i : i + 3]
            count_seq[subseq] = count_seq.get(subseq, 0) + 1
        for code in self.codon_:
            count_codon[code] = round(float(count_seq.get(code, 0) * 3 / (len(sequence))), 4)
        return count_codon

    def codon_count_all(self):
        for k_gene, v_seq in self.dict_target_.items():
            self.target_codon_[k_gene] = self.codon_count(v_seq.replace("T", "U"))

    def gc_count(self):
        for k_gene, v_seq in self.dict_target_.items():
            gc_num = v_seq.count('G') + v_seq.count('C')
            gc_norm = round(float(gc_num / len(v_seq)), 4)
            self.target_gc_[k_gene] = gc_norm

    def get_amino_code(self, code):
        amino_code = code.split()[1::2]
        codon_code = code.split()[0::2]
        for i in range(len(codon_code)):
            self.amino_[codon_code[i]] = amino_code[i]

    def amino_count(self, sequence):
        # entire CDS or 5'/3' region of CDS
        count_seq = {}
        count_amino = {}
        for i in range(0, len(sequence), 3):
            try:
                subseq_amino = self.amino_[sequence[i : i + 3]]
            except:
                subseq_amino = 'NA'
            count_seq[subseq_amino] = count_seq.get(subseq_amino, 0) + 1
        # print(count_seq)
        for amino in self.amino_.values():
            # print(amino)
            count_amino[amino] = round(float(count_seq.get(amino, 0) * 3 / (len(sequence))), 4)
        return count_amino

    def amino_count_all(self):
        for k_gene, v_seq in self.dict_target_.items():
            try:
                self.target_amino_[k_gene] = self.amino_count(v_seq.replace("T", "U"))
            except:
                print(k_gene)


def main():
    f_seq = open(sys.argv[1], 'r')
    feature = sys.argv[2]
    region = sys.argv[3]

    target_seq = LoadSeqFasta(f_seq)
    target_seq.load_fa()
    # print(target_seq.target_)

    count_usage = CountSeq(target_seq.target_)

    if feature == 'base':
        count_usage.base_count_all()
        # print(count_usage.target_base_)
        print(" " + "\t" + "\t".join("Nucl_" + str(x) for x in ["A" + '_' + region, "T" + '_' + region,
        "G" + '_' + region, "C" + '_' + region]))
        for k_gene, v_count in count_usage.target_base_.items():
            print("%s\t%s\t%s\t%s\t%s" % (k_gene, count_usage.target_base_[k_gene].get('A', 0),
            count_usage.target_base_[k_gene].get('T', 0), count_usage.target_base_[k_gene].get('G', 0),
            count_usage.target_base_[k_gene].get('C', 0)))

    elif feature == 'codon':
        # print codon usage
        count_usage.get_codon_code(codon_code_string)
        # print(count_usage.codon_)
        count_usage.codon_count_all()
        # print(count_usage.target_codon_)
        for k_gene, v_count in count_usage.target_codon_.items():
            count_codon_print = []
            for k_codon, v_percent in v_count.items():
                col_name = k_codon + '_' + region
                count_codon_print.append(col_name)
        print(" " + "\t" + "\t".join(str(x) for x in count_codon_print))

        for k_gene, v_count in count_usage.target_codon_.items():
            count_percent_print = []
            for k_codon, v_percent in v_count.items():
                count_percent_print.append(v_percent)
            print(k_gene + "\t" + "\t".join(str(x) for x in count_percent_print))

    elif feature == 'codon_nstop':
        # print codon usage
        count_usage.get_codon_code(codon_code_string_nstop)
        # print(count_usage.codon_)
        count_usage.codon_count_all()
        # print(count_usage.target_codon_)
        for k_gene, v_count in count_usage.target_codon_.items():
            count_codon_print = []
            for k_codon, v_percent in v_count.items():
                col_name = k_codon + '_' + region
                count_codon_print.append(col_name)
        print(" " + "\t" + "\t".join(str(x) for x in count_codon_print))

        for k_gene, v_count in count_usage.target_codon_.items():
            count_percent_print = []
            for k_codon, v_percent in v_count.items():
                count_percent_print.append(v_percent)
            print(k_gene + "\t" + "\t".join(str(x) for x in count_percent_print))

    elif feature == 'gc':
        count_usage.gc_count()
        print(" " + "\t" + 'GC' + '_' + region)
        for k_gene, v_count in count_usage.target_gc_.items():
            print(k_gene + "\t" + str(v_count))

    elif feature == 'amino':
        count_usage.get_amino_code(codon_code_string_nstop)
        # print(count_usage.amino_)
        # print(len(count_usage.amino_))

        count_usage.amino_count_all()

        for k_gene, v_count in count_usage.target_amino_.items():
            count_amino_print = []
            for k_amino, v_percent in v_count.items():
                # print(k_amino)
                col_name = k_amino + '_' + region
                count_amino_print.append(col_name)
        print(" " + "\t" + "\t".join("Amino_" + str(x) for x in count_amino_print))

        for k_gene, v_count in count_usage.target_amino_.items():
            count_percent_print = []
            for k_amino, v_percent in v_count.items():
                count_percent_print.append(v_percent)
            print(k_gene + "\t" + "\t".join(str(x) for x in count_percent_print))


if __name__ == "__main__":
    main()
