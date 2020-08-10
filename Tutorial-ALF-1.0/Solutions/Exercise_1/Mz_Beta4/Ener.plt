set terminal postscript eps enhanced solid color 'Times-Roman' 18
set size 0.8,0.8
set out 'Ener_B4.eps'
set fit errorvariables
tmp='Ener_fit.dat'
set print tmp
set xlabel "{/Symbol  D t}^2"
set ylabel "Energy"
f(x) = a   + b * x

fit f(x) "Ener.dat" u ($1*$1):2:3 via a,b
print "0.0  ", a,  a_err, b, b_err,  "\n\n"

plot  "Ener.dat" u ($1*$1):2:3  w e lt 11  t "", \
      "Ener_exact.dat"  u 1:2 w p lt  1 t "Exact", \
      "Ener_fit.dat"  u 1:2:3 w p lt  11 t "Fit", \
      f(x)     w l lt 1 t ""

!epstopdf Ener_B4.eps
!open Ener_B4.pdf
