#!/bin/bash

for i in $(seq 0 10); do # snow arts loop
	for j in $(seq 0 10); do # grau arts loop 
      		python main_run_efficient.py $i $j 0 999 &
      		#echo $i'and '$j 
     
	done
 	wait 
done 

for i in $(seq 0 10); do # snow arts loop
        for j in $(seq 0 10); do # grau arts loop 
                python main_run_efficient.py $i $j 1 999 &
                #echo $i'and '$j 

        done
        wait
done

echo "All done completed"


