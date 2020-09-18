set terminal postscript eps enhanced solid color 'Times-Roman' 18
set size 0.8,0.8
set title "{L=14 Hubbard Ladder, {/Symbol b}t=10, U/t=4 }"
set out 'ladder.eps'
set fit errorvariables
#tmp='Ener_fit.dat'
#set print tmp
set xlabel "r"
set ylabel "S(r,0)"
f(x) = a   + b * x

#print "0.0  ", a,  a_err, b, b_err,  "\n\n"

plot  "ladder.dat" i 0  u 1:2:3 w e lt 2 t "t_y=0" ,\
      "ladder.dat" i 1  u 1:2:3 w e lt 3 t "t_y=1" ,\
      "ladder.dat" i 2  u 1:2:3 w e lt 4 t "t_y=2" ,\
      "ladder.dat" i 0  u 1:2 w l lt 2 t "" ,\
      "ladder.dat" i 1  u 1:2 w l lt 3 t "" ,\
      "ladder.dat" i 2  u 1:2 w l lt 4 t ""

!epstopdf ladder.eps
!open ladder.pdf
