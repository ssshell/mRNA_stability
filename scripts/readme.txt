
############## get partial sequence coordinates annotation
####### split by either number of portions at the end OR by number of nucleotides from the 5'end starting position
####### CDS
###### by nucleotides
##### first 300 nucleotides; from 5'end starting position
python get_partialAnnotation.py ../index/msmeg_CombinedAnnotation_CDS.bed nucleotide 300 > ../index/msmeg_CombinedAnnotation_CDS_splitBy5pNucleotide_5p300nt.bed

##### first 300 nucleotides; from 3'end ending position
python get_partialAnnotation.py ../index/msmeg_CombinedAnnotation_noTypeField.bed nucleotide 300 3p > ../index/msmeg_CombinedAnnotation_splitBy3pNucleotide_3p300nt_p1.bed
paste ../index/msmeg_CombinedAnnotation.bed ../index/msmeg_CombinedAnnotation_splitBy3pNucleotide_3p300nt_p1.bed | cut -f1-2,7- > ../index/msmeg_CombinedAnnotation_splitBy3pNucleotide_3p300nt.bed


############## get partial sequence coordinates annotation
############# split into three specific regions. 5'end, middle and 3'end of CDS
python get_splitAnnotation.py ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CDS_length.txt 5p > ../index/msmeg_CombinedAnnotation_CDS_oneThird_5p.bed
python get_splitAnnotation.py ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CDS_length.txt mid > ../index/msmeg_CombinedAnnotation_CDS_oneThird_mid.bed
python get_splitAnnotation.py ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CDS_length.txt 3p > ../index/msmeg_CombinedAnnotation_CDS_oneThird_3p.bed


############## extract sequence from genome file
####### 5'UTR sequence
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt utr > ../index/msmeg_CombinedAnnotation_5pUTR.fa

####### 5'end 25nt upstream of leader undefined CDS
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS_leaderUndefined.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt gene_upN > ../index/msmeg_CombinedAnnotation_CDS_leaderUndefined_5pUp25nt.fa

####### 5'end 25nt upstream of leadered and leader_undefined CDS start position
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS_leaderedAndleaderUndefined.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt gene_upN > ../index/msmeg_CombinedAnnotation_CDS_leaderedAndleaderUndefined_5pUp25nt.fa

####### 5'UTR sequence 
###### sequence length >= 35 nt
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS_leaderedAndleaderUndefined.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt utr_exclude_last15 > ../index/msmeg_CombinedAnnotation_5pUTR_longThan35ExclLast15.fa

####### 5'UTR 20nt window with 10nt overlap
###### for 5'UTR longer than 35nt
python seq_windows.py ../index/msmeg_CombinedAnnotation_5pUTR_longThan35ExclLast15.fa 20 10 > ../feature/secondary_structure/UTR_5p_MFE_win/msmeg_CombinedAnnotation_5pUTR_20_10nt.fa

####### 5'UTR sequence plus 18nt of 5'end CDS
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt utr_18cds > ../index/msmeg_CombinedAnnotation_5pUTRplus18CDS.fa

####### 5'UTR sequence last 30nt plus 20nt of 5'end CDS
###### for leadered genes with at least 30nt long 5'UTR
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt utr_last30_20cds > ../index/msmeg_CombinedAnnotation_5pUTRlast30ntPlus20CDS.fa

####### 5'UTR sequence last 30nt plus start codon
###### for leadered genes with at least 30nt long 5'UTR
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt utr_last30_start_codon > ../index/msmeg_CombinedAnnotation_5pUTRlast30ntPlusStartCodon.fa

####### 5'end 20nt of the transcript
###### either 5'UTR (leadered) or CDS (leaderless) or 5'UTR + CDS if 5'UTR is shorter than 20nt (leadered)
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt tss_gene_5p20 > ../index/msmeg_CombinedAnnotation_transcript_5p20nt.fa

####### translation initiation regions sequence (TIR)
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt TIR > ../index/msmeg_CombinedAnnotation_TIR.fa

####### CDS sequence
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt gene > ../index/msmeg_CombinedAnnotation_CDS.fa

####### 5'end 18nt of CDS sequence 
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt gene_5p18 > ../index/msmeg_CombinedAnnotation_CDS_5p18.fa

####### 3'end 18nt of CDS sequence
###### including last 3nt
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt gene_3p18 > ../index/msmeg_CombinedAnnotation_CDS_3p18.fa

####### CDS sequence without stop codon
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt gene_nstop > ../index/msmeg_CombinedAnnotation_CDS_nstop.fa

####### 3'end 18nt of CDS sequence without stop codon
###### 3'end 21nt excluding last 3nt
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt gene_3p18_nstop > ../index/msmeg_CombinedAnnotation_CDS_nstop_3p18.fa

####### CDS sequence 20nt window with 10nt overlap
python seq_windows.py ../index/msmeg_CombinedAnnotation_CDS.fa 20 10 > ../feature/secondary_structure/CDS_MFE_win/msmeg_CombinedAnnotation_CDS_20_10nt.fa

####### CDS sequence 50nt window with 25nt overlap
python seq_windows.py ../index/msmeg_CombinedAnnotation_CDS.fa 50 25 > ../feature/secondary_structure/CDS_MFE_win/msmeg_CombinedAnnotation_CDS_50_25nt.fa

####### CDS sequence 100nt window with 50nt overlap
python seq_windows.py ../index/msmeg_CombinedAnnotation_CDS.fa 100 50 > ../feature/secondary_structure/CDS_MFE_win/msmeg_CombinedAnnotation_CDS_100_50nt.fa

####### 3'UTR sequence
###### use 60nt after stop codon as 3'UTR sequence
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt utr_3p_60 > ../index/msmeg_CombinedAnnotation_3pUTR60nt.fa

####### 3'UTR 20nt window with 10nt overlap
python seq_windows.py ../index/msmeg_CombinedAnnotation_3pUTR60nt.fa 20 10 > ../feature/secondary_structure/UTR_3p_MFE_win/msmeg_CombinedAnnotation_3pUTR_20_10nt.fa


############## get mRNA features
############# 5'UTR
####### 5'UTR length
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt utr_length > ../index/msmeg_CombinedAnnotation_5pUTR_length_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_5pUTR_length_p1.txt > ../index/msmeg_CombinedAnnotation_5pUTR_length_p2.txt
head -n1 ../index/msmeg_CombinedAnnotation_5pUTR_length_p1.txt | sed 's/utr/UTR/' | cat - ../index/msmeg_CombinedAnnotation_5pUTR_length_p2.txt > ../index/msmeg_5pUTR_length.txt

####### 5'UTR nucleotide usage
python count_seq.py ../index/msmeg_CombinedAnnotation_5pUTR.fa base 5pUTR > ../feature/nucleotide/msmeg_nucleotide_5pUTR_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_nucleotide_5pUTR_p1.txt > ../feature/nucleotide/msmeg_nucleotide_5pUTR_p2.txt
head -n1 ../feature/nucleotide/msmeg_nucleotide_5pUTR_p1.txt | cat - ../feature/nucleotide/msmeg_nucleotide_5pUTR_p2.txt > ../feature/nucleotide/msmeg_nucleotide_5pUTR.txt

####### 5'UTR GC percentage
python count_seq.py ../index/msmeg_CombinedAnnotation_5pUTR.fa gc 5pUTR > ../feature/nucleotide/msmeg_nucleotideGC_5pUTR_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_nucleotideGC_5pUTR_p1.txt > ../feature/nucleotide/msmeg_nucleotideGC_5pUTR_p2.txt
head -n1 ../feature/nucleotide/msmeg_nucleotideGC_5pUTR_p1.txt | cat - ../feature/nucleotide/msmeg_nucleotideGC_5pUTR_p2.txt > ../feature/nucleotide/msmeg_nucleotideGC_5pUTR.txt

####### 5'UTR adjacent nucleotide
python get_dinucleotideFreq.py ../index/msmeg_CombinedAnnotation_5pUTR.fa 5pUTR > ../feature/nucleotide/msmeg_dinucleotide_5pUTR_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_dinucleotide_5pUTR_p1.txt > ../feature/nucleotide/msmeg_dinucleotide_5pUTR_p2.txt
head -n1 ../feature/nucleotide/msmeg_dinucleotide_5pUTR_p1.txt | cat - ../feature/nucleotide/msmeg_dinucleotide_5pUTR_p2.txt > ../feature/nucleotide/msmeg_dinucleotide_5pUTR.txt

####### 5'UTR SD strength
###### GA feature in -17 to -4 for leadered and leader_undefined genes
python search_SD_ga.py ../index/msmeg_CombinedAnnotation_5pUTR.fa > ../feature/shine_dalgarno/msmeg_SD_GA_leadered.txt
python search_SD_ga.py ../index/msmeg_CombinedAnnotation_CDS_leaderUndefined_5pUp25nt.fa > ../feature/shine_dalgarno/msmeg_SD_GA_leaderedUndefined.txt
sed 1d ../feature/shine_dalgarno/msmeg_SD_GA_leadered.txt | cat ../feature/shine_dalgarno/msmeg_SD_GA_leaderedUndefined.txt - > ../feature/shine_dalgarno/msmeg_SD_GA_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/shine_dalgarno/msmeg_SD_GA_p1.txt > ../feature/shine_dalgarno/msmeg_SD_GA_p2.txt
head -n1 ../feature/shine_dalgarno/msmeg_SD_GA_p1.txt | cat - ../feature/shine_dalgarno/msmeg_SD_GA_p2.txt > ../feature/shine_dalgarno/msmeg_SD_GA.txt

###### motif frequency within 5'end 25nt upstream of leadered and leader_undefined CDS start position
python search_SD.py ../index/msmeg_CombinedAnnotation_CDS_leaderedAndleaderUndefined_5pUp25nt.fa > ../feature/shine_dalgarno/msmeg_SD_freq_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/shine_dalgarno/msmeg_SD_freq_p1.txt > ../feature/shine_dalgarno/msmeg_SD_freq_p2.txt
head -n1 ../feature/shine_dalgarno/msmeg_SD_freq_p1.txt | cat - ../feature/shine_dalgarno/msmeg_SD_freq_p2.txt > ../feature/shine_dalgarno/msmeg_SD_freq.txt

####### 5'UTR secondary structure
###### MFE of 20nt sliding windows sequences with 10nt overlap 
RNAfold --noPS < ../feature/secondary_structure/UTR_5p_MFE_win/msmeg_CombinedAnnotation_5pUTR_20_10nt.fa > ../feature/secondary_structure/UTR_5p_MFE_win/msmeg_5pUTR_20_10nt_MFE_p1.txt
python second_struc.py ../feature/secondary_structure/UTR_5p_MFE_win/msmeg_5pUTR_20_10nt_MFE_p1.txt RNAfold mfe > ../feature/secondary_structure/UTR_5p_MFE_win/msmeg_5pUTR_20_10nt_MFE_p2.txt
python second_struc_slidingWins.py ../feature/secondary_structure/UTR_5p_MFE_win/msmeg_5pUTR_20_10nt_MFE_p2.txt > ../feature/secondary_structure/UTR_5p_MFE_win/msmeg_5pUTR_20_10nt_MFE_p3.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/secondary_structure/UTR_5p_MFE_win/msmeg_5pUTR_20_10nt_MFE_p3.txt > ../feature/secondary_structure/UTR_5p_MFE_win/msmeg_5pUTR_20_10nt_MFE_p4.txt
echo -e ' \tfprUTR_MFE_20_10nt_5pUTR\tfprUTR_MFE_20_10nt_5p\tfprUTR_MFE_20_10nt_mid\tfprUTR_MFE_20_10nt_3p' | cat - ../feature/secondary_structure/UTR_5p_MFE_win/msmeg_5pUTR_20_10nt_MFE_p4.txt > ../feature/secondary_structure/UTR_5p_MFE_win/msmeg_5pUTR_20_10nt_MFE.txt

###### number of unpaired nt at 5'end
##### folding 5'UTR sequence
RNAfold --noPS < ../index/msmeg_CombinedAnnotation_5pUTR.fa > ../feature/secondary_structure/UTR_5p_unpairedNt/msmeg_5pUTR_MFE_p1.txt
python second_struc.py ../feature/secondary_structure/UTR_5p_unpairedNt/msmeg_5pUTR_MFE_p1.txt RNAfold unpaired > ../feature/secondary_structure/UTR_5p_unpairedNt/msmeg_5pUTR_5pUnpairedNt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/secondary_structure/UTR_5p_unpairedNt/msmeg_5pUTR_5pUnpairedNt_p1.txt > ../feature/secondary_structure/UTR_5p_unpairedNt/msmeg_5pUTR_5pUnpairedNt_p2.txt
echo -e ' \tfprUTR_5pUnpairedNt' | cat - ../feature/secondary_structure/UTR_5p_unpairedNt/msmeg_5pUTR_5pUnpairedNt_p2.txt > ../feature/secondary_structure/UTR_5p_unpairedNt/msmeg_5pUTR_5pUnpairedNt.txt

##### folding 5'UTR sequence plus 18nt of 5'end CDS
RNAfold --noPS < ../index/msmeg_CombinedAnnotation_5pUTRplus18CDS.fa > ../feature/secondary_structure/UTR_5pPlus18ntCDS_unpairedNt/msmeg_5pUTRplus18ntCDS_MFE_p1.txt
python second_struc.py ../feature/secondary_structure/UTR_5pPlus18ntCDS_unpairedNt/msmeg_5pUTRplus18ntCDS_MFE_p1.txt RNAfold unpaired > ../feature/secondary_structure/UTR_5pPlus18ntCDS_unpairedNt/msmeg_5pUTRplus18ntCDS_5pUnpairedNt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/secondary_structure/UTR_5pPlus18ntCDS_unpairedNt/msmeg_5pUTRplus18ntCDS_5pUnpairedNt_p1.txt > ../feature/secondary_structure/UTR_5pPlus18ntCDS_unpairedNt/msmeg_5pUTRplus18ntCDS_5pUnpairedNt_p2.txt
echo -e ' \tfprUTR_plus18ntCDS_5pUnpairedNt' | cat - ../feature/secondary_structure/UTR_5pPlus18ntCDS_unpairedNt/msmeg_5pUTRplus18ntCDS_5pUnpairedNt_p2.txt > ../feature/secondary_structure/UTR_5pPlus18ntCDS_unpairedNt/msmeg_5pUTRplus18ntCDS_5pUnpairedNt.txt

###### unpaired probability
##### folding 5'UTR sequence last 30nt plus 20nt of 5'end CDS
#### for leadered genes with at least 30nt long 5'UTR
### run followings in ~/feature/secondary_structure/UTR_5pLast30Plus20ntCDS_unpairProb/
RNAplfold -W 50 -u 3 < ../../../index/msmeg_CombinedAnnotation_5pUTRlast30ntPlus20CDS.fa
./extract_UnpairProb_2ntUnpaired.sh
./get_CombinedNt.sh > msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_p1.txt

## 2nt averaged unpaired probability of entire sequence
python ../../../scripts/second_struc.py msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_p1.txt RNAplfold mean_global > msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_WholeSeq_p1.txt
python ../../../scripts/feature_table.py ../../../index/msmeg_CombinedAnnotation_CDS.bed msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_WholeSeq_p1.txt > msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_WholeSeq_p2.txt
echo -e ' \tfprUTRlast30nt_plus20ntCDS_2ntAvg_WholeSeq_UnpairedProb' | cat - msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_WholeSeq_p2.txt > msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_WholeSeq.txt

## 2nt averaged unpaired probability of SD region
python ../../../scripts/second_struc.py msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_p1.txt RNAplfold mean_local sd_50 > msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_SDseq_p1.txt
python ../../../scripts/feature_table.py ../../../index/msmeg_CombinedAnnotation_CDS.bed msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_SDseq_p1.txt > msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_SDseq_p2.txt
echo -e ' \tfprUTRlast30nt_plus20ntCDS_2ntAvg_SDseq_UnpairedProb' | cat - msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_SDseq_p2.txt > msmeg_5pUTRlast30ntPlus20CDS_2ntUnpairedProb_SDseq.txt

## unpaired probability of all start codon nt
./extract_UnpairProb_StatCodonUnpaired.sh
./get_CombinedStartCodonNt.sh > msmeg_5pUTRlast30ntPlus20CDS_StartCodonUnpairedProb_p1.txt
python ../../../scripts/second_struc.py msmeg_5pUTRlast30ntPlus20CDS_StartCodonUnpairedProb_p1.txt RNAplfold start_codon_all 18 > msmeg_5pUTRlast30ntPlus20CDS_StartCodonUnpairedProb_p2.txt
python ../../../scripts/feature_table.py ../../../index/msmeg_CombinedAnnotation_CDS.bed msmeg_5pUTRlast30ntPlus20CDS_StartCodonUnpairedProb_p2.txt > msmeg_5pUTRlast30ntPlus20CDS_StartCodonUnpairedProb_p3.txt
echo -e ' \tfprUTRlast30nt_plus20ntCDS_all_StartCodon_UnpairedProb' | cat - msmeg_5pUTRlast30ntPlus20CDS_StartCodonUnpairedProb_p3.txt > msmeg_5pUTRlast30ntPlus20CDS_StartCodonUnpairedProb.txt

##### folding 5'UTR sequence last 30nt plus start codon
#### for leadered genes with at least 30nt long 5'UTR
### run followings in ~/feature/secondary_structure/UTR_5pLast30PlusStartCodon_unpairProb/
RNAplfold -W 33 -u 3 < ../../../index/msmeg_CombinedAnnotation_5pUTRlast30ntPlusStartCodon.fa
./extract_UnpairProb_2ntUnpaired.sh
./get_CombinedNt.sh > msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_p1.txt

## 2nt averaged unpaired probability of entire sequence
python ../../../scripts/second_struc.py msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_p1.txt RNAplfold mean_global > msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_WholeSeq_p1.txt
python ../../../scripts/feature_table.py ../../../index/msmeg_CombinedAnnotation_CDS.bed msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_WholeSeq_p1.txt > msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_WholeSeq_p2.txt
echo -e ' \tfprUTRlast30nt_plusStartCodon_2ntAvg_WholeSeq_UnpairedProb' | cat - msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_WholeSeq_p2.txt > msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_WholeSeq.txt

## 2nt averaged unpaired probability of SD region
python ../../../scripts/second_struc.py msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_p1.txt RNAplfold mean_local sd_33 > msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_SDseq_p1.txt
python ../../../scripts/feature_table.py ../../../index/msmeg_CombinedAnnotation_CDS.bed msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_SDseq_p1.txt > msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_SDseq_p2.txt
echo -e ' \tfprUTRlast30nt_plusStartCodon_2ntAvg_SDseq_UnpairedProb' | cat - msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_SDseq_p2.txt > msmeg_5pUTRlast30ntPlusStartCodon_2ntUnpairedProb_SDseq.txt

## unpaired probability of all start codon nt
./extract_UnpairProb_StatCodonUnpaired.sh
./get_CombinedStartCodonNt.sh > msmeg_5pUTRlast30ntPlusStartCodon_StartCodonUnpairedProb_p1.txt
python ../../../scripts/second_struc.py msmeg_5pUTRlast30ntPlusStartCodon_StartCodonUnpairedProb_p1.txt RNAplfold start_codon_all 1 > msmeg_5pUTRlast30ntPlusStartCodon_StartCodonUnpairedProb_p2.txt
python ../../../scripts/feature_table.py ../../../index/msmeg_CombinedAnnotation_CDS.bed msmeg_5pUTRlast30ntPlusStartCodon_StartCodonUnpairedProb_p2.txt > msmeg_5pUTRlast30ntPlusStartCodon_StartCodonUnpairedProb_p3.txt
echo -e ' \tfprUTRlast30nt_plusStartCodon_all_StartCodon_UnpairedProb' | cat - msmeg_5pUTRlast30ntPlusStartCodon_StartCodonUnpairedProb_p3.txt > msmeg_5pUTRlast30ntPlusStartCodon_StartCodonUnpairedProb.txt

############# 5'transcript
####### 5'transcript secondary structure
###### MFE of 5'end 20nt of the transcript
##### 5’ UTR for Leadered, CDS for Leaderless or 5’ UTR plus CDS if UTR is shorter than 20 nt
RNAfold --noPS < ../index/msmeg_CombinedAnnotation_transcript_5p20nt.fa > ../feature/secondary_structure/transcript_5p20nt_MFE/msmeg_transcript_5p20nt_MFE_p1.txt
python second_struc.py ../feature/secondary_structure/transcript_5p20nt_MFE/msmeg_transcript_5p20nt_MFE_p1.txt RNAfold mfe > ../feature/secondary_structure/transcript_5p20nt_MFE/msmeg_transcript_5p20nt_MFE_p2.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/secondary_structure/transcript_5p20nt_MFE/msmeg_transcript_5p20nt_MFE_p2.txt > ../feature/secondary_structure/transcript_5p20nt_MFE/msmeg_transcript_5p20nt_MFE_p3.txt
echo -e ' \ttranscript_5p20nt_MFE' | cat - ../feature/secondary_structure/transcript_5p20nt_MFE/msmeg_transcript_5p20nt_MFE_p3.txt > ../feature/secondary_structure/transcript_5p20nt_MFE/msmeg_transcript_5p20nt_MFE.txt

###### number of unpaired nt at 5'end
python second_struc.py ../feature/secondary_structure/transcript_5p20nt_MFE/msmeg_transcript_5p20nt_MFE_p1.txt RNAfold unpaired > ../feature/secondary_structure/transcript_5p20nt_unpairedNt/msmeg_transcript_5p20nt_5pUnpairedNt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/secondary_structure/transcript_5p20nt_unpairedNt/msmeg_transcript_5p20nt_5pUnpairedNt_p1.txt > ../feature/secondary_structure/transcript_5p20nt_unpairedNt/msmeg_transcript_5p20nt_5pUnpairedNt_p2.txt
echo -e ' \ttranscript_5p20nt_5pUnpairedNt' | cat - ../feature/secondary_structure/transcript_5p20nt_unpairedNt/msmeg_transcript_5p20nt_5pUnpairedNt_p2.txt > ../feature/secondary_structure/transcript_5p20nt_unpairedNt/msmeg_transcript_5p20nt_5pUnpairedNt.txt

###### unpaired probability
##### run followings in ~/feature/secondary_structure/transcript_5p20nt_unpairProb
RNAplfold -W 20 -u 5 < ../../../index/msmeg_CombinedAnnotation_transcript_5p20nt.fa
./extract_UnpairProb_AllNt.sh
./get_combined3nt.sh > msmeg_transcript_5p20nt_First3ntUnpairedProb_p1.txt
./get_combined5nt.sh > msmeg_transcript_5p20nt_First5ntUnpairedProb_p1.txt
./extract_UnpairProb_EachNt.sh
./get_combinedNt.sh > msmeg_transcript_5p20nt_EachNtUnpairedProb_p1.txt

#### first 3nt of 5'end
### probability of all 3nt are unpaired
python ../../../scripts/second_struc.py msmeg_transcript_5p20nt_First3ntUnpairedProb_p1.txt RNAplfold 5p_3_all > msmeg_transcript_5p20nt_First3ntUnpairedProb_AllNt_p1.txt
python ../../../scripts/feature_table.py ../../../index/msmeg_CombinedAnnotation_CDS.bed msmeg_transcript_5p20nt_First3ntUnpairedProb_AllNt_p1.txt > msmeg_transcript_5p20nt_First3ntUnpairedProb_AllNt_p2.txt
echo -e ' \ttranscript_5p20nt_all_First3nt_UnpairedProb' | cat - msmeg_transcript_5p20nt_First3ntUnpairedProb_AllNt_p2.txt > msmeg_transcript_5p20nt_First3ntUnpairedProb_AllNt.txt

### average probability of each of the 3nt is unpaired
python ../../../scripts/second_struc.py msmeg_transcript_5p20nt_EachNtUnpairedProb_p1.txt RNAplfold mean_local 5p_3_avg > msmeg_transcript_5p20nt_First3ntUnpairedProb_AvgNt_p1.txt
python ../../../scripts/feature_table.py ../../../index/msmeg_CombinedAnnotation_CDS.bed msmeg_transcript_5p20nt_First3ntUnpairedProb_AvgNt_p1.txt > msmeg_transcript_5p20nt_First3ntUnpairedProb_AvgNt_p2.txt
echo -e ' \ttranscript_5p20nt_Avg_First3nt_UnpairedProb' | cat - msmeg_transcript_5p20nt_First3ntUnpairedProb_AvgNt_p2.txt > msmeg_transcript_5p20nt_First3ntUnpairedProb_AvgNt.txt

#### first 5nt of 5'end
### probability of all 5nt are unpaired
python ../../../scripts/second_struc.py msmeg_transcript_5p20nt_First5ntUnpairedProb_p1.txt RNAplfold 5p_5_all > msmeg_transcript_5p20nt_First5ntUnpairedProb_AllNt_p1.txt
python ../../../scripts/feature_table.py ../../../index/msmeg_CombinedAnnotation_CDS.bed msmeg_transcript_5p20nt_First5ntUnpairedProb_AllNt_p1.txt > msmeg_transcript_5p20nt_First5ntUnpairedProb_AllNt_p2.txt
echo -e ' \ttranscript_5p20nt_all_First5nt_UnpairedProb' | cat - msmeg_transcript_5p20nt_First5ntUnpairedProb_AllNt_p2.txt > msmeg_transcript_5p20nt_First5ntUnpairedProb_AllNt.txt

### average probability of each of the 5nt is unpaired
python ../../../scripts/second_struc.py msmeg_transcript_5p20nt_EachNtUnpairedProb_p1.txt RNAplfold mean_local 5p_5_avg > msmeg_transcript_5p20nt_First5ntUnpairedProb_AvgNt_p1.txt
python ../../../scripts/feature_table.py ../../../index/msmeg_CombinedAnnotation_CDS.bed msmeg_transcript_5p20nt_First5ntUnpairedProb_AvgNt_p1.txt > msmeg_transcript_5p20nt_First5ntUnpairedProb_AvgNt_p2.txt
echo -e ' \ttranscript_5p20nt_Avg_First5nt_UnpairedProb' | cat - msmeg_transcript_5p20nt_First5ntUnpairedProb_AvgNt_p2.txt > msmeg_transcript_5p20nt_First5ntUnpairedProb_AvgNt.txt

####### 5'transcript unfold of translation initiation region MFE
###### native TIR RNA structure (deltaG_mRNA)
RNAfold --noPS < ../index/msmeg_CombinedAnnotation_TIR.fa > ../feature/translation_init/msmeg_deltaG_mRNA_p1.txt
python second_struc.py ../feature/translation_init/msmeg_deltaG_mRNA_p1.txt RNAfold mfe > ../feature/translation_init/msmeg_deltaG_mRNA.txt

###### TIR bound by an initiating ribosome  (deltaG_init)
##### run followings in ~/feature/translation_init
python ../../../scripts/get_dot4TIR.py ../msmeg_deltaG_mRNA_p1.txt

##### run followings back in ./scripts
./get_dot2ct.sh > dot2ct.log 

##### get sequence file with base-pairing constraints
python get_dotByCt.py ../index/msmeg_5pUTR_length.txt ../feature/translation_init/msmeg_deltaG_mRNA_p1.txt > ../index/msmeg_deltaG_init.fa
RNAfold -C --noPS < ../index/msmeg_deltaG_init.fa > ../feature/translation_init/msmeg_deltaG_init_p1.txt 
python second_struc.py ../feature/translation_init/msmeg_deltaG_init_p1.txt RNAfold mfe > ../feature/translation_init/msmeg_deltaG_init.txt

###### energy required by the ribosome to unfold the mRNA at the TIR (deltaG_unfold)
python get_deltaG_unfold.py ../feature/translation_init/msmeg_deltaG_mRNA.txt ../feature/translation_init/msmeg_deltaG_init.txt > ../feature/translation_init/msmeg_deltaG_unfold_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/translation_init/msmeg_deltaG_unfold_p1.txt > ../feature/translation_init/msmeg_deltaG_unfold_p2.txt
echo -e ' \tMFE_unfold' | cat - ../feature/translation_init/msmeg_deltaG_unfold_p2.txt > ../feature/translation_init/msmeg_deltaG_unfold.txt

####### 5'transcript ribosome occupancy
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/ribosome_profiling/ribo_normalizedBy_totalmRNA.txt > ../feature/ribosome_profiling/msmeg_riboProfiling_p1.txt
echo -e ' \tribo_CDS\tribo_CDSplus20up\tribo_5p\tribo_5p_excl' | cat - ../feature/ribosome_profiling/msmeg_riboProfiling_p1.txt > ../feature/ribosome_profiling/msmeg_riboProfiling.txt
cut -f1,4 ../feature/ribosome_profiling/msmeg_riboProfiling.txt > ../feature/ribosome_profiling/msmeg_riboProfiling_5pTranscript.txt

####### first and second nt of 5'end of the transcript
###### either 5'UTR (leadered) or CDS (leaderless)
sed -E '/^[^>]/ s/.{18}$//' ../index/msmeg_CombinedAnnotation_transcript_5p20nt.fa | paste - - | awk '{$2 = substr($2, 1, 1) "\t" substr($2, 2)}1' | sed 's/^.//' | tr ' ' \\t > ../index/msmeg_transcript_5pFirstAndSecond_nt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_transcript_5pFirstAndSecond_nt_p1.txt > ../index/msmeg_transcript_5pFirstAndSecond_nt_p2.txt
echo -e ' \t5p_1stNT\t5p2ndNT' | cat - ../index/msmeg_transcript_5pFirstAndSecond_nt_p2.txt > ../index/msmeg_transcript_5pFirstAndSecond_nt.txt

############# CDS
####### CDS length
python get_seq.py ../index/Mycobacterium_smegmatis_MC2-155_genome_v4.fasta ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_TSS.txt gene_length > ../index/msmeg_CombinedAnnotation_CDS_length_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../index/msmeg_CombinedAnnotation_CDS_length_p1.txt > ../index/msmeg_CombinedAnnotation_CDS_length_p2.txt
head -n1 ../index/msmeg_CombinedAnnotation_CDS_length_p1.txt | sed 's/gene/CDS/' | cat - ../index/msmeg_CombinedAnnotation_CDS_length_p2.txt > ../index/msmeg_CDS_length.txt

####### CDS nucleotide usage
###### usage of entire CDS (including last 3nt)
python count_seq.py ../index/msmeg_CombinedAnnotation_CDS.fa base CDS > ../feature/nucleotide/msmeg_nucleotide_CDS_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_nucleotide_CDS_p1.txt > ../feature/nucleotide/msmeg_nucleotide_CDS_p2.txt
head -n1 ../feature/nucleotide/msmeg_nucleotide_CDS_p1.txt | cat - ../feature/nucleotide/msmeg_nucleotide_CDS_p2.txt > ../feature/nucleotide/msmeg_nucleotide_CDS.txt

###### usage of 5'end 18nt of CDS
python count_seq.py ../index/msmeg_CombinedAnnotation_CDS_5p18.fa base CDS_5p18nt > ../feature/nucleotide/msmeg_nucleotide_CDS_5p18nt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_nucleotide_CDS_5p18nt_p1.txt > ../feature/nucleotide/msmeg_nucleotide_CDS_5p18nt_p2.txt
head -n1 ../feature/nucleotide/msmeg_nucleotide_CDS_5p18nt_p1.txt | cat - ../feature/nucleotide/msmeg_nucleotide_CDS_5p18nt_p2.txt > ../feature/nucleotide/msmeg_nucleotide_CDS_5p18nt.txt

###### usage of 3'end 18nt of CDS
##### including last 3nt
python count_seq.py ../index/msmeg_CombinedAnnotation_CDS_3p18.fa base CDS_3p18nt > ../feature/nucleotide/msmeg_nucleotide_CDS_3p18nt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_nucleotide_CDS_3p18nt_p1.txt > ../feature/nucleotide/msmeg_nucleotide_CDS_3p18nt_p2.txt
head -n1 ../feature/nucleotide/msmeg_nucleotide_CDS_3p18nt_p1.txt | cat - ../feature/nucleotide/msmeg_nucleotide_CDS_3p18nt_p2.txt > ../feature/nucleotide/msmeg_nucleotide_CDS_3p18nt.txt

####### CDS GC percentage
###### usage of entire CDS (including last 3nt)
python count_seq.py ../index/msmeg_CombinedAnnotation_CDS.fa gc CDS > ../feature/nucleotide/msmeg_nucleotideGC_CDS_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_nucleotideGC_CDS_p1.txt > ../feature/nucleotide/msmeg_nucleotideGC_CDS_p2.txt
head -n1 ../feature/nucleotide/msmeg_nucleotideGC_CDS_p1.txt | cat - ../feature/nucleotide/msmeg_nucleotideGC_CDS_p2.txt > ../feature/nucleotide/msmeg_nucleotideGC_CDS.txt

###### usage of 5'end 18nt of CDS
python count_seq.py ../index/msmeg_CombinedAnnotation_CDS_5p18.fa gc CDS_5p18nt > ../feature/nucleotide/msmeg_nucleotideGC_CDS_5p18nt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_nucleotideGC_CDS_5p18nt_p1.txt > ../feature/nucleotide/msmeg_nucleotideGC_CDS_5p18nt_p2.txt
head -n1 ../feature/nucleotide/msmeg_nucleotideGC_CDS_5p18nt_p1.txt | cat - ../feature/nucleotide/msmeg_nucleotideGC_CDS_5p18nt_p2.txt > ../feature/nucleotide/msmeg_nucleotideGC_CDS_5p18nt.txt

###### usage of 3'end 18nt of CDS
##### including last 3nt
python count_seq.py ../index/msmeg_CombinedAnnotation_CDS_3p18.fa gc CDS_3p18nt > ../feature/nucleotide/msmeg_nucleotideGC_CDS_3p18nt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_nucleotideGC_CDS_3p18nt_p1.txt > ../feature/nucleotide/msmeg_nucleotideGC_CDS_3p18nt_p2.txt
head -n1 ../feature/nucleotide/msmeg_nucleotideGC_CDS_3p18nt_p1.txt | cat - ../feature/nucleotide/msmeg_nucleotideGC_CDS_3p18nt_p2.txt > ../feature/nucleotide/msmeg_nucleotideGC_CDS_3p18nt.txt

####### CDS adjacent nucleotide
###### usage of entire CDS (including last 3nt)
python get_dinucleotideFreq.py ../index/msmeg_CombinedAnnotation_CDS.fa CDS > ../feature/nucleotide/msmeg_dinucleotide_CDS_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_dinucleotide_CDS_p1.txt > ../feature/nucleotide/msmeg_dinucleotide_CDS_p2.txt
head -n1 ../feature/nucleotide/msmeg_dinucleotide_CDS_p1.txt | cat - ../feature/nucleotide/msmeg_dinucleotide_CDS_p2.txt > ../feature/nucleotide/msmeg_dinucleotide_CDS.txt

###### usage of 5'end 18nt of CDS
python get_dinucleotideFreq.py ../index/msmeg_CombinedAnnotation_CDS_5p18.fa CDS_5p18nt > ../feature/nucleotide/msmeg_dinucleotide_CDS_5p18nt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_dinucleotide_CDS_5p18nt_p1.txt > ../feature/nucleotide/msmeg_dinucleotide_CDS_5p18nt_p2.txt
head -n1 ../feature/nucleotide/msmeg_dinucleotide_CDS_5p18nt_p1.txt | cat - ../feature/nucleotide/msmeg_dinucleotide_CDS_5p18nt_p2.txt > ../feature/nucleotide/msmeg_dinucleotide_CDS_5p18nt.txt

###### usage of 3'end 18nt of CDS
##### including last 3nt
python get_dinucleotideFreq.py ../index/msmeg_CombinedAnnotation_CDS_3p18.fa CDS_3p18nt > ../feature/nucleotide/msmeg_dinucleotide_CDS_3p18nt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/nucleotide/msmeg_dinucleotide_CDS_3p18nt_p1.txt > ../feature/nucleotide/msmeg_dinucleotide_CDS_3p18nt_p2.txt
head -n1 ../feature/nucleotide/msmeg_dinucleotide_CDS_3p18nt_p1.txt | cat - ../feature/nucleotide/msmeg_dinucleotide_CDS_3p18nt_p2.txt > ../feature/nucleotide/msmeg_dinucleotide_CDS_3p18nt.txt

####### CDS start & stop codon usage
###### start codon
paste - - < ../index/msmeg_CombinedAnnotation_CDS.fa | sed 's/>//' | awk '{print $1"\t"substr($2,1,3)}' > ../feature/codon/msmeg_CDS_StartCodon_p1.txt 
python convert2binary.py ../feature/codon/msmeg_CDS_StartCodon_p1.txt start_codon > ../feature/codon/msmeg_CDS_StartCodon_p2.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/codon/msmeg_CDS_StartCodon_p2.txt > ../feature/codon/msmeg_CDS_StartCodon_p3.txt
head -n1 ../feature/codon/msmeg_CDS_StartCodon_p2.txt | sed 's/T/U/g' | cat - ../feature/codon/msmeg_CDS_StartCodon_p3.txt > ../feature/codon/msmeg_CDS_StartCodon.txt

###### stop codon
paste - - < ../index/msmeg_CombinedAnnotation_CDS.fa | sed 's/>//'| awk '{print $1"\t"substr($2,length($2) - 2, length($2))}' > ../feature/codon/msmeg_CDS_StopCodon_p1.txt
python convert2binary.py ../feature/codon/msmeg_CDS_StopCodon_p1.txt stop_codon > ../feature/codon/msmeg_CDS_StopCodon_p2.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/codon/msmeg_CDS_StopCodon_p2.txt > ../feature/codon/msmeg_CDS_StopCodon_p3.txt
head -n1 ../feature/codon/msmeg_CDS_StopCodon_p2.txt | sed 's/T/U/g' | cat - ../feature/codon/msmeg_CDS_StopCodon_p3.txt > ../feature/codon/msmeg_CDS_StopCodon.txt

####### CDS codon usage
###### usage of entire CDS (excluding last 3nt)
python count_seq.py ../index/msmeg_CombinedAnnotation_CDS_nstop.fa codon_nstop CDSnstop > ../feature/codon/msmeg_codon_CDSnstop_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/codon/msmeg_codon_CDSnstop_p1.txt > ../feature/codon/msmeg_codon_CDSnstop_p2.txt
head -n1 ../feature/codon/msmeg_codon_CDSnstop_p1.txt | cat - ../feature/codon/msmeg_codon_CDSnstop_p2.txt > ../feature/codon/msmeg_codon_CDSnstop.txt

###### usage of 5'end 18nt of CDS 
python count_seq.py ../index/msmeg_CombinedAnnotation_CDS_5p18.fa codon_nstop CDS_5p18nt > ../feature/codon/msmeg_codon_CDS_5p18nt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/codon/msmeg_codon_CDS_5p18nt_p1.txt > ../feature/codon/msmeg_codon_CDS_5p18nt_p2.txt
head -n1 ../feature/codon/msmeg_codon_CDS_5p18nt_p1.txt | cat - ../feature/codon/msmeg_codon_CDS_5p18nt_p2.txt > ../feature/codon/msmeg_codon_CDS_5p18nt.txt

###### usage of 3'end 18nt of CDS excluding stop codon
##### last 21nt excluding last 3nt
python count_seq.py ../index/msmeg_CombinedAnnotation_CDS_nstop_3p18.fa codon_nstop CDSnstop_3p18nt > ../feature/codon/msmeg_codon_CDSnstop_3p18nt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/codon/msmeg_codon_CDSnstop_3p18nt_p1.txt > ../feature/codon/msmeg_codon_CDSnstop_3p18nt_p2.txt
head -n1 ../feature/codon/msmeg_codon_CDSnstop_3p18nt_p1.txt | cat - ../feature/codon/msmeg_codon_CDSnstop_3p18nt_p2.txt > ../feature/codon/msmeg_codon_CDSnstop_3p18nt.txt

####### CDS ajacent codon pair bias
python get_codonPairScoreBias.py ../index/msmeg_CombinedAnnotation_CDS.fa CPB > ../feature/codon/msmeg_CPB_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/codon/msmeg_CPB_p1.txt > ../feature/codon/msmeg_CPB_p2.txt
echo -e ' \tCodonPairBias' | cat - ../feature/codon/msmeg_CPB_p2.txt > ../feature/codon/msmeg_CPB.txt

####### CDS secondary structure
###### MFE of sliding windows sequences
##### 20nt window with 10nt overlap
RNAfold --noPS < ../feature/secondary_structure/CDS_MFE_win/msmeg_CombinedAnnotation_CDS_20_10nt.fa > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_20_10nt_MFE_p1.txt
python second_struc.py ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_20_10nt_MFE_p1.txt RNAfold mfe > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_20_10nt_MFE_p2.txt
python second_struc_slidingWins.py ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_20_10nt_MFE_p2.txt > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_20_10nt_MFE_p3.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_20_10nt_MFE_p3.txt > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_20_10nt_MFE_p4.txt
echo -e ' \tCDS_MFE_20_10nt_CDS\tCDS_MFE_20_10nt_5p\tCDS_MFE_20_10nt_mid\tCDS_MFE_20_10nt_3p' | cat - ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_20_10nt_MFE_p4.txt > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_20_10nt_MFE.txt

##### 50nt window with 25nt overlap
RNAfold --noPS < ../feature/secondary_structure/CDS_MFE_win/msmeg_CombinedAnnotation_CDS_50_25nt.fa > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_50_25nt_MFE_p1.txt
python second_struc.py ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_50_25nt_MFE_p1.txt RNAfold mfe > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_50_25nt_MFE_p2.txt
python second_struc_slidingWins.py ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_50_25nt_MFE_p2.txt > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_50_25nt_MFE_p3.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_50_25nt_MFE_p3.txt > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_50_25nt_MFE_p4.txt
echo -e ' \tCDS_MFE_50_25nt_CDS\tCDS_MFE_50_25nt_5p\tCDS_MFE_50_25nt_mid\tCDS_MFE_50_25nt_3p' | cat - ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_50_25nt_MFE_p4.txt > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_50_25nt_MFE.txt

##### 100nt window with 50nt overlap
RNAfold --noPS < ../feature/secondary_structure/CDS_MFE_win/msmeg_CombinedAnnotation_CDS_100_50nt.fa > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_100_50nt_MFE_p1.txt
python second_struc.py ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_100_50nt_MFE_p1.txt RNAfold mfe > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_100_50nt_MFE_p2.txt
python second_struc_slidingWins.py ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_100_50nt_MFE_p2.txt > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_100_50nt_MFE_p3.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_100_50nt_MFE_p3.txt > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_100_50nt_MFE_p4.txt
echo -e ' \tCDS_MFE_100_50nt_CDS\tCDS_MFE_100_50nt_5p\tCDS_MFE_100_50nt_mid\tCDS_MFE_100_50nt_3p' | cat - ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_100_50nt_MFE_p4.txt > ../feature/secondary_structure/CDS_MFE_win/msmeg_CDS_100_50nt_MFE.txt

###### number of unpaired nt at 5'end
##### folding entire CDS sequence, including last 3nt
RNAfold --noPS < ../index/msmeg_CombinedAnnotation_CDS.fa > ../feature/secondary_structure/CDS_5p_unpairedNt/msmeg_CDS_MFE_p1.txt
python second_struc.py ../feature/secondary_structure/CDS_5p_unpairedNt/msmeg_CDS_MFE_p1.txt RNAfold unpaired > ../feature/secondary_structure/CDS_5p_unpairedNt/msmeg_CDS_5pUnpairedNt_p1.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/secondary_structure/CDS_5p_unpairedNt/msmeg_CDS_5pUnpairedNt_p1.txt > ../feature/secondary_structure/CDS_5p_unpairedNt/msmeg_CDS_5pUnpairedNt_p2.txt
echo -e ' \tCDS_5pUnpairedNt' | cat - ../feature/secondary_structure/CDS_5p_unpairedNt/msmeg_CDS_5pUnpairedNt_p2.txt > ../feature/secondary_structure/CDS_5p_unpairedNt/msmeg_CDS_5pUnpairedNt.txt

####### CDS ribosome occupancy
cut -f1-3,5 ../feature/ribosome_profiling/msmeg_riboProfiling.txt > ../feature/ribosome_profiling/msmeg_riboProfiling_CDS.txt

############# 3'UTR
####### 3'UTR secondary structure
###### use 60nt after stop codon as 3'UTR sequence
##### MFE of 20nt sliding windows sequences with 10nt overlap
RNAfold --noPS < ../feature/secondary_structure/UTR_3p_MFE_win/msmeg_CombinedAnnotation_3pUTR_20_10nt.fa > ../feature/secondary_structure/UTR_3p_MFE_win/msmeg_3pUTR_20_10nt_MFE_p1.txt
python second_struc.py ../feature/secondary_structure/UTR_3p_MFE_win/msmeg_3pUTR_20_10nt_MFE_p1.txt RNAfold mfe > ../feature/secondary_structure/UTR_3p_MFE_win/msmeg_3pUTR_20_10nt_MFE_p2.txt
python second_struc_slidingWins.py ../feature/secondary_structure/UTR_3p_MFE_win/msmeg_3pUTR_20_10nt_MFE_p2.txt > ../feature/secondary_structure/UTR_3p_MFE_win/msmeg_3pUTR_20_10nt_MFE_p3.txt
python feature_table.py ../index/msmeg_CombinedAnnotation_CDS.bed ../feature/secondary_structure/UTR_3p_MFE_win/msmeg_3pUTR_20_10nt_MFE_p3.txt > ../feature/secondary_structure/UTR_3p_MFE_win/msmeg_3pUTR_20_10nt_MFE_p4.txt
echo -e ' \tthrprUTR_MFE_20_10nt_3pUTR\tthrprUTR_MFE_20_10nt_5p\tthrprUTR_MFE_20_10nt_mid\tthrprUTR_MFE_20_10nt_3p' | cat - ../feature/secondary_structure/UTR_3p_MFE_win/msmeg_3pUTR_20_10nt_MFE_p4.txt > ../feature/secondary_structure/UTR_3p_MFE_win/msmeg_3pUTR_20_10nt_MFE.txt

