#!/bin/bash

for file in *5ntUnpaired.txt;do
	id="${file%%.*}"
	echo $id
	paste -sd" " $file
done
