reset
set log y
set grid
set yrange[9e+99:0]
set format x "%.1te%02T"
set format y "%.1te%02T"
set ylabel 'C_{avg} [m^{2}]'
set xlabel '{/Symbol l} [m]'
set title " Dust mixture 001; "
plot '-' with points title '|dC_{ext}| (Cpol)' lc rgb "#0000FF",'-' with points title '|dC_{abd}|'  lc rgb "#FF0000",'-' with points title '|dC_{phas}| (C_{circ})' lc rgb "#800080"
e
e
e
