#!/bin/bash

if [ -e "output.sh" ]
then
rm "output.sh"
fi

COUNTER=1
while [  $COUNTER -lt 2 ]; do

cat covHRMf_test2.slurm | sed "s/nodelist.txt/nodelist_$COUNTER.txt/g;s/*/$COUNTER/g" > run_$COUNTER.slurm

echo "sbatch run_$COUNTER.slurm" >> output.sh
let COUNTER=COUNTER+1
done
