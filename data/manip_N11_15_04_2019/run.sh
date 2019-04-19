export SWAP=1
for f in *.txt;
do
	echo $f;
	../../bin/cbridge4 $f R0=82.6 
done
