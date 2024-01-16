reset
set log y
set grid
unset key
set yrange[0.9:1.1]
set format x "%.1te%02T"
set format y "%.1te%02T"
set ylabel '{/Symbol d}N(a)'
set xlabel 'a [m]'
set title " Dust mixture 001; "
plot '-' with points lc rgb "#0000F0"
2e-07	1
e
