#!/usr/bin/env gnuplot

set terminal pdf enhanced dashed
set output("all.pdf")

set xlabel "Time [hr]"
set ylabel "Power [MW]"

set log y
set format y "10^{%L}"

set yrange [1.e-6:5000e6]
set xrange [0.0:100.0]

plot \
     "<paste  power_interp.dat decay_power_new.dat" using ($1/86400):($2+$4+$5) with lines lt 1 lw 3.0 lc rgb "#006400" title "Total Power", \
     "power_interp.dat" using ($1/86400):($2*(1-0.006108)) with lines lt 1 lw 3.0 lc rgb "red" axes x1y1 title "Prompt Fission Power", \
     "power_interp.dat" using ($1/86400):($2*0.006108) with lines lt 1 lw 3.0 lc rgb "blue" axes x1y1 title "Delayed Fission Power", \
     "decay_power_new.dat" using ($1/86400):($3+$4) with lines lt 1 lw 3.0 lc rgb "orange" title "U-235 and U-238 Decay Power", \
     "actinide_decay.dat" using ($1/86400+72.0):2 with lines lt 1 lw 3.0 lc rgb "black" title "Actinide Decay Power"
