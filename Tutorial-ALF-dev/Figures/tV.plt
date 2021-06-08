set terminal postscript eps enhanced solid color 'Times-Roman' 18
set size 0.8,0.8
set title "{L=28 tV, {/Symbol b}t=20}"
set out 'tV.eps'
set fit errorvariables
set xlabel "r"
set ylabel "<n(r)n(0)>"
data_file = "tV.dat"

plot "tV.dat" i 0  u 1:2:3 w e lc rgb "black"      lt 1 pt 7 lw 2 t "V/t=1" ,  \
      ''      i 0  u 1:2   w l lc rgb "black"      lt 1      lw 2 t "" ,       \
      ''      i 1  u 1:2:3 w e lc rgb "red"        lt 1 pt 7 lw 2 t "V/t=2" ,  \
      ''      i 1  u 1:2   w l lc rgb "red"        lt 1      lw 2 t "" ,       \
      ''      i 2  u 1:2:3 w e lc rgb "royalblue"  lt 1 pt 7 lw 2 t "V/t=2.5" ,\
      ''      i 2  u 1:2   w l lc rgb "royalblue"  lt 1      lw 2 t ""

!epstopdf tV.eps
