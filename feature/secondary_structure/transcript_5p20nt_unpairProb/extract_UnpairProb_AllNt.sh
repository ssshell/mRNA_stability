#!/bin/bash
for file in *lunp
do
	id="${file%_*}"
	cut -f6 $file|sed '1,2d' > $id.5ntUnpaired.txt
	cut -f4 $file|sed '1,2d' > $id.3ntUnpaired.txt
done
