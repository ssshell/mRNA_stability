#!/bin/bash

for file in *lunp
do
	id="${file%_*}"
	cut -f2 $file|sed '1,2d' > $id.EachNtUnpaired.txt
done
