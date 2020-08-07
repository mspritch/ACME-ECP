#!/bin/bash

for COUNTER in 4 3 2 1 0.5 0.2 0.1
#
do
sed -e "s/NTIME/${COUNTER}/g;" test_NT_ne4.csh > Process_ne4_nx8_32_NTIME_${COUNTER}.csh

chmod 700 Process_ne4_nx8_32_NTIME_${COUNTER}.csh
./Process_ne4_nx8_32_NTIME_${COUNTER}.csh

mv Process_ne4_nx8_32_NTIME_* /home1/07088/tg863871/temp/

done
