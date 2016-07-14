#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output("decay.pdf")

set xlabel "Time [s]"
set ylabel "Decay Power [W]"

plot "decay_power.dat" using 1:($3+$4) with lines lt 1 lw 3.0 lc rgb "blue"
