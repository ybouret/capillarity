PROG=./bin/cbridge4
for f in `find verre_teflon -name "*.txt"`;
do
	echo $f;
	$PROG $f
done
