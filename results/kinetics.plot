#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output("kinetics.pdf")

set ytics nomirror
set y2tics autofreq

set xlabel "Time [s]"
set ylabel "Power [MW]"
set y2label "Reactivity [$]"

set log y

set yrange [1.e-6:5000e6]

plot \
     "power.dat" using ($2/86400):4 with lines lt 2 lw 3.0 lc rgb "blue" axes x1y1 title "Power Trace", \
     "power.dat" using ($2/86400):($3/0.006108) with lines lt 1 lw 3.0 lc rgb "red" axes x1y2 title "Reactivity"
