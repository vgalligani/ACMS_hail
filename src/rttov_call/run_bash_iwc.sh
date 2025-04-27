#!/bin/bash

#for i in $(seq 1 10); do # snow arts loop
#      	python main_run_efficient.py $i 0 3 &   # this is for WSM6
#done 
#wait
#python main_run_efficient.py 9 9 0 'grau_iwc'  & #  # this is for eqMass WSM6
#python main_run_efficient.py 9 9 0 'snow_iwc'  & #  # this is for eqMass WSM6
#python main_run_efficient.py 3 3 0 'grau_iwc'  & #  # this is for eqMass WSM6
#python main_run_efficient.py 3 3 0 'snow_iwc'  & #  # this is for eqMass WSM6
python main_run_efficient.py 3 3 1 'grausnow_iwc'  & #  # this is for eqMass WSM6
python main_run_efficient.py 9 9 1 'grausnow_iwc'  & #  # this is for eqMass WSM6


#python main_run_efficient.py 9 9 1 'grau_iwc' &  #  # this is for WSM6
#python main_run_efficient.py 9 9 1 'snow_iwc' &  #  # this is for WSM6
#python main_run_efficient.py 3 3 1 'grau_iwc' &  #  # this is for WSM6
#python main_run_efficient.py 3 3 1 'snow_iwc' &  #  # this is for WSM6

wait

echo "All done completed"


