#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output("all.pdf")

set xlabel "Time [s]"
set ylabel "Power [MW]"

set log y
set format y "10^{%L}"

set yrange [1.e-6:5000e6]

plot \
     "power_interp.dat" using ($1/86400):2 with lines lt 1 lw 3.0 lc rgb "blue" axes x1y1 title "Fission Power", \
     "decay_power_new.dat" using ($1/86400):3 with lines lt 1 lw 3.0 lc rgb "orange" title "U-235 Decay Power", \
     "decay_power_new.dat" using ($1/86400):4 with lines lt 1 lw 3.0 lc rgb "red" title "U-238 Decay Power", \
     "<paste  power_interp.dat decay_power_new.dat" using ($1/86400):($2+$4+$5) with lines lt 1 lw 3.0 lc rgb "#006400" title "Total Power"
