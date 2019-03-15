export SWAP=1
for f in *.txt;
do
#R0=76.5
	echo $f;
	../bin/cbridge4 $f R0=76.5
done
