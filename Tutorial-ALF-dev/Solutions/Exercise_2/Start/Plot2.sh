echo "Deleting $1_plot if exists."
rm -f $1"_plot"
echo "Reading $1 and printing some of its lines to $1_plot."
while read x y
do
    echo $x $y
    if (( $(echo "$x < 0.00" | sed 's/E./*10^/g' | bc -l) )) ; then
	read line
    else
	read o p best err dump
	if (( $(echo "$y == 0.00" | sed 's/E./*10^/g' | bc -l) )) ; then
            echo $x $best $err >> $1"_plot"
	fi
    fi
done < $1
