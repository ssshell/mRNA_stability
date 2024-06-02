#!/bin/bash
for file in *lunp
do
	id="${file%_*}"
	cut -f4 $file|sed '1,4d' > $id.StartCodonUnpaired.txt
done
