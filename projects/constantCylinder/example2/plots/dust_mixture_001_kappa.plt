reset
set log y
set grid
set yrange[916360:8.82809e+06]
set format x "%.1te%02T"
set format y "%.1te%02T"
set ylabel '{/Symbol k}_{avg}  [m^2/kg]'
set xlabel '{/Symbol l} [m]'
set title " Dust mixture 001; "
plot '-' with points title '{/Symbol k}_{ext,x}' lc rgb "#0000F0",'-' with points title '{/Symbol k}_{ext,y}' lc rgb "#000090",'-' with points title '{/Symbol k}_{abs,x}' lc rgb "#FF0000",'-' with points title '{/Symbol k}_{abs,y}' lc rgb "#900000",'-' with points title '{/Symbol k}_{sca,x}' lc rgb "#FFFF00",'-' with points title '{/Symbol k}_{sca,y}' lc rgb "#909000"
6.6e-07	8.02554e+06
e
6.6e-07	8.02554e+06
e
6.6e-07	1.01818e+06
e
6.6e-07	1.01818e+06
e
6.6e-07	7.00736e+06
e
6.6e-07	7.00736e+06
e
