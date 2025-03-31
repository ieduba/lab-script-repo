for loop in scr-combo-mapped-filt_500.mcool-10000-loops dH1-combo-mapped-filt_500.mcool-10000-loops; do
	tail -n +2 $loop\.csv > $loop\.tmp
	python loop-convert.py -l $loop\.tmp
	bedtools sort -i $loop\.tmp.longrange > $loop\.longrange
	bgzip $loop\.longrange
	tabix -0 -s 1 -b 2 -e 3 $loop\.longrange.gz
done
rm *.tmp*

