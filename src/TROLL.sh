#!/bin/bash

clear

echo "BUILD"
g++  main.cpp  -o TROLL.out
echo "...done"
echo ""

echo "PARAMETERS"
echo "...done"
echo ""

echo "RUN"
read -p "Confirm run (yes/no)  : " confirm
if [ $confirm == yes ]
	then ./TROLL.out -i'input.txt' -o./OUTPUT/'test2'
	else echo "You cancelled the run."
fi
echo "...done"
echo ""
