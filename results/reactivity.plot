#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output("reactivity.pdf")

set ytics nomirror
set y2tics autofreq

set xlabel "Time [hr]"
set ylabel "Reactivity [%dk/k]"
set y2label "Reactivity [$]"

set format y "%.2tx10^{%L}"
set format y2 "%.2tx10^{%L}"

set yrange [0:1e-4]
set y2range [0:1e-6/0.006108]

set key font ",10"

plot "reactivity.dat" using 2:($3*100) with lines lt 2 lw 3.0 lc rgb "red" axes x1y1 title "Reactivity using Inverse Kinetics", \
     '' using 2:($3/0.006108) with lines lt 2 lw 3.0 lc rgb "red" axes x1y2 notitle
