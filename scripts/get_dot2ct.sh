#!/bin/bash

cd /Users/huamingsun/Desktop/PhD/shelllab/mRNA_stability/msmeg_ThisIsTheFinal/feature/translation_init/RNAstructure/exe

for file in ../../deltaG_mRNA_dotFiles/*dot;do
	base1="${file##*/}"
	base2="${base1%.*}"
	./dot2ct $file ../../deltaG_mRNA_ctFiles/$base2.ct
done
