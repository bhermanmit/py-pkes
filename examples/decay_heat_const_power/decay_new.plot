#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output("decay_new.pdf")

set xlabel "Time [s]"
set ylabel "Decay Power [MW]"

plot "decay_power_new.dat" using 1:($3+$4) with lines lt 1 lw 3.0 lc rgb "blue"
