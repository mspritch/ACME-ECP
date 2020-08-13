#!/bin/csh
cp control.txt test.txt
echo "BEFORE:"
cat test.txt
/bin/sed -i 's/\!\$acc parallel/\!\$omp parallel/g' ./test.txt
/bin/sed -i 's/async(asyncid)//g' ./test.txt
echo "--------"
echo "AFTER:"
cat test.txt
