reset
set log y
set grid
set yrange[0.367001:3.55284]
set format x "%.1te%02T"
set format y "%.1te%02T"
set ylabel 'Q_{avg}'
set xlabel '{/Symbol l} [m]'
set title " Dust mixture 001; "
plot '-' with points title 'Q_{ext,x}' lc rgb "#0000F0",'-' with points title 'Q_{ext,y}' lc rgb "#000090",'-' with points title 'Q_{abs,x}' lc rgb "#FF0000",'-' with points title 'Q_{abs,y}' lc rgb "#900000",'-' with points title 'Q_{sca,x}' lc rgb "#FFFF00",'-' with points title 'Q_{sca,y}' lc rgb "#909000"
6.6e-07	3.22985
e
6.6e-07	3.22985
e
6.6e-07	0.407778
e
6.6e-07	0.407778
e
6.6e-07	2.82207
e
6.6e-07	2.82207
e
