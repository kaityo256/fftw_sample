set term png
set out "acf.png"
set log y
set yra [0.01:1]
set xra [0:46.05]
set pointsize 1.0
set xla "t"
set yla "<p(0)p(t)>"
p "acf.dat" t "Data"  pt 7\
, exp(-0.1*x) lt 1 lc 0

