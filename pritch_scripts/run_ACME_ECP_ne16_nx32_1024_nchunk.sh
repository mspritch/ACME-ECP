#!/bin/bash

TEMP_PATH="/home1/07088/tg863871/temp/"
cd ${TEMP_PATH}

for COUNTER in 2065 2066 2067 2068 2069 2070 2071 2072 2073 2074 2075 2076 
#
do
sed -e "s/NCHUNK/${COUNTER}/g;" Process_NCHUNK_ne16_nx32_1024.csh > Process_ne16_nx32_1024_NCHUNK_${COUNTER}.csh

chmod 700 Process_ne16_nx32_1024_NCHUNK_${COUNTER}.csh
./Process_ne16_nx32_1024_NCHUNK_${COUNTER}.csh

done
