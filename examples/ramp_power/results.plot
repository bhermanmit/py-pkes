#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output("power_comparision.pdf")

set ytics nomirror
set yrange [1:1300]
set y2tics autofreq
set y2range [-5:3]

set xlabel "Time [s]"
set ylabel "Power [W]"
set y2label "Reactivity [$]"

plot "reactivity.dat" using 1:3 with lines lt 1 lw 3.0 lc rgb "blue" axes x1y1 title "Input Power Trace", \
     "reactivity.dat" using 1:($2/0.006108) with lines lt 1 lw 3.0 lc rgb "red" axes x1y2 title "Reactivity using Inverse Kinetics", \
     "power.dat"      using 1:3 with lines lt 2 lw 3.0 lc rgb "orange" axes x1y1 title "Power using Point Kinetics"
