
import sys

f_deltaG_mRNA = open(sys.argv[1], 'r')
f_deltaG_init = open(sys.argv[2], 'r')

def load_deltaG(f_deltaG):
    dict_deltaG = {}
    for line in f_deltaG:
        line = line.rstrip().split('\t')
        gene = line[0]
        deltaG = float(line[1])
        dict_deltaG[gene] = deltaG
    return dict_deltaG

deltaG_mRNA = load_deltaG(f_deltaG_mRNA)
deltaG_init = load_deltaG(f_deltaG_init)
# print(deltaG_mRNA)
# print(deltaG_init)

for k_gene, v_deltaG_init in deltaG_init.items():
    deltaG_unfold = v_deltaG_init - deltaG_mRNA[k_gene]
    print("%s\t%.4f" % (k_gene, deltaG_unfold))
