#!/bin/bash
# use ./run_bash.sh 

MAX_JOBS=2
# =============================================================================
# For graupelsp and liu for snow (and aggregates)
#count=0
#for i in $(seq 0 10); do # snow arts loop
#for i in 2 5; do   
#    python main_process_efficient_grau.py '' $i & 
#    echo $i
#    ((count++))
#    if  (( count % MAX_JOBS == 0 )); then
#        echo 'waiting at '$i
#        wait
#    fi    
#done 
#wait
#echo "All done completed for liu and grau soft sphere"

# =============================================================================
#count=0
#for i in $(seq 0 10); do # snow arts loop
#    python main_process_efficient_grau.py 'halfgrau' $i & 
#    echo $i      
#    ((count++))
#    if  (( count % MAX_JOBS == 0 )); then
#        wait
#    fi  
#done
#wait 
#echo "All done completed for liu and grau soft sphere but halfgrau"

# =============================================================================
for i in $(seq 4 10); do 
	for j in 8 9 10; do   #	for j in 2 5 8 9 10; do     
    	python main_process_efficient.py '' $i $j & 
    	wait
	done
	echo 'all graupel liu combinations for '$i
#	wait 
done 
echo "All done completed for liu liu "

# =============================================================================
#for i in $(seq 0 10); do # snow arts loop
#	for j in 2 5 8 9 10; do     
#    	python main_process_efficient.py 'halfgrau' $i $j & 
#	done
#	echo 'all graupel liu combinations for '$i
# 	wait 
#done 
#echo "All done completed for liu liu but half grau" 



