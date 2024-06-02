#!/bin/bash

for file in *3ntUnpaired.txt;do
	id="${file%%.*}"
	echo $id
	paste -sd" " $file
done
