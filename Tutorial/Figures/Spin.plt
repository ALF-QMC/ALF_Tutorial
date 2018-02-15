set terminal postscript eps enhanced solid color 'Times-Roman' 18
set size 0.8,0.8
set title "{L=14 Hubbard Ladder, {/Symbol b}t=10, U/t=4 }"
set out 'Ladder.eps'
set fit errorvariables
#tmp='Ener_fit.dat'
#set print tmp
set xlabel "r"
set ylabel "S(r,0)"
f(x) = a   + b * x

#print "0.0  ", a,  a_err, b, b_err,  "\n\n"

plot  "Spin.dat" i 0  u 1:2:3 w e lt 2 t "t_y=0" ,\
      "Spin.dat" i 1  u 1:2:3 w e lt 3 t "t_y=1" ,\
      "Spin.dat" i 2  u 1:2:3 w e lt 4 t "t_y=2" ,\
      "Spin.dat" i 0  u 1:2 w l lt 2 t "" ,\
      "Spin.dat" i 1  u 1:2 w l lt 3 t "" ,\
      "Spin.dat" i 2  u 1:2 w l lt 4 t ""

!epstopdf Ladder.eps
!open Ladder.pdf
