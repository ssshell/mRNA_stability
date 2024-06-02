#!/bin/bash
for file in *lunp
do
	id="${file%_*}"
	cut -f3 $file|sed '1,3d' > $id.2ntUnpaired.txt
done
