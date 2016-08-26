#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output("inverse_kinetics.pdf")

set ytics nomirror
set y2tics autofreq

set xlabel "Time [hr]"
set ylabel "Power [W]"
set y2label "Reactivity [$]"

set log y

# set y2range [-5.e-6:5.5e-5]

plot "reactivity.dat" using 2:4 with lines lt 1 lw 3.0 lc rgb "blue" axes x1y1 title "Input Power Trace", \
     "reactivity.dat" using 2:($3/0.006108) with lines lt 2 lw 3.0 lc rgb "red" axes x1y2 title "Reactivity using Inverse Kinetics"
