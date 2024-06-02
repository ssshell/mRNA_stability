#!/bin/bash

for file in *.EachNtUnpaired.txt;do
	id="${file%%.*}"
	echo $id
	paste -sd" " $file
done
