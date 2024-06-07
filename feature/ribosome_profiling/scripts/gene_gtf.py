# create gtf file

import sys

f_annotation = open(sys.argv[1], 'r')
f_leaderless = open(sys.argv[2], 'r')
region = sys.argv[3]
# leadered_ad_others = []
leaderless = []

for line in f_leaderless:
    line = line.rstrip().split("\t")
    leaderless.append(line[0])
# print(leaderless)
# print(len(leaderless))

for line in f_annotation:
    line = line.rstrip().split("\t")
    locus = line[0]
    # feature = line[1]
    feature = 'CDS'
    start = line[1]
    end = line[2]
    strand = line[3]
    seqname = 'NC_008596.1'
    if region == "CDS":
        region_start = start
        region_end = end
    elif region == "CDS_20up":
        if strand == "+":
            region_start = int(start) - 20
            region_end = end
        else:
            region_start = start
            region_end = int(end) + 20
    elif region == '5p_end':
        if locus in leaderless:
            if strand == "+":
                region_start = start
                region_end = int(start) + 18
            else:
                region_start = int(end) - 18
                region_end = end
        else:
            if strand == "+":
                region_start = int(start) - 20
                region_end = int(start) + 18
            else:
                region_start = int(end) - 18
                region_end = int(end) + 20
    elif region == '5p_end_excl':
        if strand == "+":
            region_start = int(start) + 19
            region_end = end
        else:
            region_start = start
            region_end = int(end) - 19

    source = '.'
    score = '.'
    phase = '.'
    attributes = 'gene_id' + ' ' + "\"" + locus + "\"" + ';' + ' ' + 'transcript_id' \
        + ' ' + "\"" + locus + "\"" + ';'
    print("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % \
        (seqname, source, feature, region_start, region_end, score, strand, phase, attributes))
