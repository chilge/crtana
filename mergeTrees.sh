#!/bin/bash

if [ -d "./Trees" ]
then
	echo "Tree directory found."
	nfile=`ls -1 ./Trees/orm*.root | wc -l`
	echo $nfile " orm test trees found. Merging..."
	hadd ./Trees/merged"$nfile".root ./Trees/orm*.root
else 
	echo "Tree directory could not be found! Exiting..."
	exit 9999

fi
