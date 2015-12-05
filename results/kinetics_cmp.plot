#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output("kinetics_cmp.pdf")

set ytics nomirror
set y2tics autofreq

set xlabel "Time [s]"
set ylabel "Power [W]"
set y2label "Reactivity [$]"

set log y

set yrange [1:5000e6]
set xrange [0:72]

plot "power.dat" using 2:4 with lines lt 1 lw 3.0 lc rgb "blue" axes x1y1 title "Calculated Power Trace", \
     "reactivity.dat" using 2:4 with lines lt 1 lw 3.0 lc rgb "orange" axes x1y1 title "Target Power Trace", \
     "power.dat" using 2:($3/0.006108) with lines lt 1 lw 3.0 lc rgb "red" axes x1y2 title "Reactivity using Inverse Kinetics"
