#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output("decay.pdf")

set xlabel "Time [s]"
set ylabel "Decay Power [MW]"

set format y "10^{%L}"

set log y

plot "decay_power_new.dat" using ($1/86400):($3+$4) with lines notitle lt 1 lw 3.0 lc rgb "blue"
