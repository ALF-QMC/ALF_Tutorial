set terminal postscript eps enhanced solid color 'Times-Roman' 18
set size 0.8,0.8
set title "{L=28 tV, {/Symbol b}t=20}"
set out 'tV.eps'
set fit errorvariables
#tmp='Ener_fit.dat'
#set print tmp
set xlabel "r"
set ylabel "<n(r)n(0)>"
f(x) = a   + b * x

#print "0.0  ", a,  a_err, b, b_err,  "\n\n"

plot  "Den.dat" i 0  u 1:2:3 w e lt 2 t "V/t=1" ,\
      "Den.dat" i 1  u 1:2:3 w e lt 3 t "V/t=2" ,\
      "Den.dat" i 2  u 1:2:3 w e lt 4 t "V/t=2.5" ,\
      "Den.dat" i 0  u 1:2 w l lt 2 t "" ,\
      "Den.dat" i 1  u 1:2 w l lt 3 t "" ,\
      "Den.dat" i 2  u 1:2 w l lt 4 t ""

!epstopdf tV.eps
!open tV.pdf
