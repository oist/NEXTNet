#!/usr/bin/env gnuplot

set terminal cairolatex pdf #GOOD ONE
#set term tikz standalone size 10cm, 5cm
set output "figure.tex"
#set encoding iso_8859_1

set format x "$%g$"
set format y "10^{%L}"

#set title "Mixture of Poisson  $(0.75\\times \\mathbf{3} + 0.25 \\times \\mathbf{50})$"
set xlabel "$t$"
set ylabel "Number of infected"
#rotate by 0

set yrange [1:100000]


set logscale y
set key

set key bottom right

set datafile separator ","
plot 0.7*exp(0.2934*x) dashtype '.-_' lt -1 title " EL",\
0.7*exp(0.351101*x) dashtype 2 lt -1 title "extended EL",\
"datahetero.csv" lt -1 notitle
