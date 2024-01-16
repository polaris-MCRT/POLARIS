reset
set log y
set grid
set yrange[1.37625e+06:1.33231e+07]
set format x "%.1te%02T"
set format y "%.1te%02T"
set ylabel '{/Symbol k}_{avg}  [m^2/kg]'
set xlabel '{/Symbol l} [m]'
set title " Dust mixture 001; "
plot '-' with points title '{/Symbol k}_{ext,x}' lc rgb "#0000F0",'-' with points title '{/Symbol k}_{ext,y}' lc rgb "#000090",'-' with points title '{/Symbol k}_{abs,x}' lc rgb "#FF0000",'-' with points title '{/Symbol k}_{abs,y}' lc rgb "#900000",'-' with points title '{/Symbol k}_{sca,x}' lc rgb "#FFFF00",'-' with points title '{/Symbol k}_{sca,y}' lc rgb "#909000"
6.6e-07	1.21119e+07
e
6.6e-07	1.21119e+07
e
6.6e-07	1.52917e+06
e
6.6e-07	1.52917e+06
e
6.6e-07	1.05828e+07
e
6.6e-07	1.05828e+07
e
