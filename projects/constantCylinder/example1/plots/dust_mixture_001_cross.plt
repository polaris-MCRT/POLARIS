reset
set log y
set grid
set yrange[4.61186e-14:4.46462e-13]
set format x "%.1te%02T"
set format y "%.1te%02T"
set ylabel 'C_{avg} [m^{2}]'
set xlabel '{/Symbol l} [m]'
set title " Dust mixture 001; "
plot '-' with points title 'C_{ext,x}' lc rgb "#0000F0",'-' with points title 'C_{ext,y}' lc rgb "#000090",'-' with points title 'C_{abs,x}' lc rgb "#FF0000",'-' with points title 'C_{abs,y}' lc rgb "#900000",'-' with points title 'C_{sca,x}' lc rgb "#FFFF00",'-' with points title 'C_{sca,y}' lc rgb "#909000"
6.6e-07	4.05875e-13
e
6.6e-07	4.05875e-13
e
6.6e-07	5.12429e-14
e
6.6e-07	5.12429e-14
e
6.6e-07	3.54632e-13
e
6.6e-07	3.54632e-13
e
