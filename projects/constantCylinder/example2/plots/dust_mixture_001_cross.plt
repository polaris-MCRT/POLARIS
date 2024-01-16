reset
set log y
set grid
set yrange[1.11275e-13:1.07201e-12]
set format x "%.1te%02T"
set format y "%.1te%02T"
set ylabel 'C_{avg} [m^{2}]'
set xlabel '{/Symbol l} [m]'
set title " Dust mixture 001; "
plot '-' with points title 'C_{ext,x}' lc rgb "#0000F0",'-' with points title 'C_{ext,y}' lc rgb "#000090",'-' with points title 'C_{abs,x}' lc rgb "#FF0000",'-' with points title 'C_{abs,y}' lc rgb "#900000",'-' with points title 'C_{sca,x}' lc rgb "#FFFF00",'-' with points title 'C_{sca,y}' lc rgb "#909000"
6.6e-07	9.74552e-13
e
6.6e-07	9.74552e-13
e
6.6e-07	1.23639e-13
e
6.6e-07	1.23639e-13
e
6.6e-07	8.50913e-13
e
6.6e-07	8.50913e-13
e
