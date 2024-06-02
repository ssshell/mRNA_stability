#!/bin/bash

for file in *.2ntUnpaired.txt;do
	id="${file%%.*}"
	echo $id
	paste -sd" " $file
done
