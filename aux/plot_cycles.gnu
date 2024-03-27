#!/usr/bin/gnuplot -c


#This script takes in three input paramters
#It can be invoked with the command
#gnuplot -c script_name filein fileout lineno
#The input arguments have the following meaning
#filein: name of the file with input data
#fileout: name of the file to which the plot is meant to be written
#titlename: The title of the plot to be output
#lineno: Number of lines in the infile


filein=ARG1
fileout=ARG2
title_name=ARG3
lineno=ARG4

set terminal pdf
set output fileout
set xrange [lineno:0]
set grid
set title title_name
set xlabel 'Number of Clusters'
set ylabel "Min Distance Between Clusters [{A}]"
set y2label 'Inter-Cluster distance increase from previous %'
#unset key
set xtics lineno/10
set ytics nomirror
set y2tics
plot filein using 1:5 axes x1y1 with steps lw 2 lt 1 title 'distance',\
filein using 1:6 axes x1y2 with points pt 6 ps 0.4 lt 3 title 'increase%'
