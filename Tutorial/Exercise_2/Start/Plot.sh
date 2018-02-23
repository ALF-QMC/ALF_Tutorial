rm -f $1"_plot"
while read x y
do
  if (( $(echo "$x < 0.00" | bc -l) )) ; then
    read line
  else
    read o p best err dump 
    if (( $(echo "$y == 0.00" | bc -l) )) ; then
      echo $x $best $err >> $1"_plot"
    fi
  fi
done < $1
