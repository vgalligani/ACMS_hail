#!/bin/bash

for i in $(seq 1 10); do # snow arts loop
	for j in $(seq 0 10); do # grau arts loop 
      		python main_run_efficient.py $i $j &
      		#echo $i'and '$j 
     
	done
 	wait 
done 

echo "All done completed"


