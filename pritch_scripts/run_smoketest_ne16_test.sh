#!/bin/bash

for COUNTER in 2080 2112 2176 2304 2560
#
do
sed -e "s/NPPP/${COUNTER}/g;" smoketest_ne16_test.csh > smoketest_ne16_test_NP${COUNTER}.csh

chmod 700 smoketest_ne16_test_NP${COUNTER}.csh
./smoketest_ne16_test_NP${COUNTER}.csh

done
mv smoketest_ne16_test_NP*.csh /home1/07088/tg863871/temp/

