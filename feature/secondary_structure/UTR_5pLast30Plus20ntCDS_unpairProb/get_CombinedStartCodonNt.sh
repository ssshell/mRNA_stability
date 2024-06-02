#!/bin/bash

for file in *.StartCodonUnpaired.txt;do
	id="${file%%.*}"
	echo $id
	paste -sd" " $file
done
