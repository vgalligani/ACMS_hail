#!/bin/bash

for i in $(seq 1 10); do # snow arts loop
      	python main_run_efficient.py $i 0 3 999 &   # this is for WSM6
done 
wait


for i in $(seq 1 10); do # snow arts loop
        python main_run_efficient.py $i 0 4 999 &   # this is for WSM6
done 
wait



echo "All done completed"


