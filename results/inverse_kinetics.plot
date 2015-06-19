#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output("inverse_kinetics.pdf")

set ytics nomirror
set y2tics autofreq

set xlabel "Time [s]"
set ylabel "Power [W]"
set y2label "Reactivity [$]"

set log y

plot "reactivity.dat" using 1:3 with lines lt 1 lw 3.0 lc rgb "blue" axes x1y1 title "Input Power Trace", \
     "reactivity.dat" using 1:($2/0.006108) with lines lt 1 lw 3.0 lc rgb "red" axes x1y2 title "Reactivity using Inverse Kinetics"
