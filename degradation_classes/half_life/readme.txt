
####### get half_life values in separate files
###### log phase
##### cd /half_life_logPhase/con_no_atc\ half-life\ calls (export to .CSV files from .xlsx)
#### check number of genes & get half_life values
cut -d',' -f1,5,10 Half-life\ calls\ 9-19-20-Table\ 1.csv | awk -F',' '$2==$3{print$1}'|wc -l

cut -d',' -f1,5 Half-life\ calls\ 9-19-20-Table\ 1.csv | sed '1d'| sed "s/\,/\\t/g" > ../msmeg_halfLife_logPhase.txt

cut -d',' -f1,6 Half-life\ calls\ 9-19-20-Table\ 1.csv | sed '1d'| sed "s/\,/\\t/g" > ../msmeg_halfLife_MSE_logPhase.txt

###### hypoxia
##### cd /half_life_hypoxia
cut -d',' -f1,3 Hypoxia\ half-life\ calls.csv | sed '1d' | sed "s/\,/\\t/g" > msmeg_halfLife_hypoxia.txt

cut -d',' -f1,6 Hypoxia\ half-life\ calls.csv | sed '1d' | sed "s/\,/\\t/g" > msmeg_halfLife_MSE_hypoxia.txt

###### fold change of hypoxia to log phase
##### cd /half_life_hyp2logPhaseFC
#### same results using con_no_atc_metrics_half-lives.xlsx or Hypoxia half-life calls.xlsx
### get half_life_fold_increase_in_hyp.csv from con_no_atc_metrics_half-lives.xlsx
tr ',' '\t' < half_life_fold_increase_in_hyp.csv | awk '$2!=0' > msmeg_halfLife_hyp2logPhaseFC.txt

### get hypoxia MSE from Hypoxia half-life calls.xlsx; to check value with MSE from hypoxia folder
## checked. they are the same.
cut -d',' -f1,7 1_5\ MSE\ \<0.5\ fil\ 0\ and\ pos\ slop-Table\ 1.csv | sed '1d' | sed "s/\,/\\t/g" > ../half_life_summary_MSE_hypoxia.txt

####### get half_life values and fold change in hypoxia from one file
###### get half_life_summary.csv from file ./half_life_hyp2logPhaseFC/con_no_atc_metrics_half-lives.xlsx -> Combined
##### checked. the same values as from seperate files above
#### this file only have genes with fold change available in hypoxia, not complete list of hypoxia half-life available genes

###### get log phase MSE from file ./half_life_hyp2logPhaseFC/con_no_atc_metrics_half-lives.xlsx -> Combined; to check value with MSE from log phase folder
##### checked. they are the same

############## at the end summary
####### genes in the above single file only have hypoxia half-life for those have fold change in hypoxia
####### for global half-life distribution in individual condition use separate files with the most complete available genes  
####### half-life & fold change are the same for 3988 genes from these two ways of loading 
####### half-life log phase is the same; hypoxia round to 2 digits is the same; fold change is the same
####### CDS_halfLife_completeClass.txt contains all the half-life values and classes
####### use this file for downstream analysis. remove NAs when extract for each condition/class
####### all the value checking process can be found in get_halfLifeValueClass.R

############## adding RNaseE half-life
####### get the half-life file from file "RNase E KD half-lives.xlsx"