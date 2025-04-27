#!/bin/bash

for j in $(seq 0 10); do # grau arts loop 
      	python main_run_efficient.py 0 $j &
done 
wait


echo "All done completed"


