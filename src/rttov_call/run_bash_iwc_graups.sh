#!/bin/bash

#for i in $(seq 1 10); do # snow arts loop
#      	python main_run_efficient.py $i 0 3 &   # this is for WSM6
#done 
#wait
python main_run_efficient.py 9 99 0 'grau_iwc'  & #  # this is for eqMass WSM6
python main_run_efficient.py 9 99 0 'snow_iwc'  & #  # this is for eqMass WSM6
python main_run_efficient.py 3 99 0 'grau_iwc'  & #  # this is for eqMass WSM6
python main_run_efficient.py 3 99 0 'snow_iwc'  & #  # this is for eqMass WSM6
python main_run_efficient.py 9 99 0 'grausnow_iwc'  & #  # this is for eqMass WSM6
python main_run_efficient.py 3 99 0 'grausnow_iwc'  & #  # this is for eqMass WSM6


python main_run_efficient.py 9 99 1 'grau_iwc' &  #  # this is for WSM6
python main_run_efficient.py 9 99 1 'snow_iwc' &  #  # this is for WSM6
python main_run_efficient.py 3 99 1 'grau_iwc' &  #  # this is for WSM6
python main_run_efficient.py 3 99 1 'snow_iwc' &  #  # this is for WSM6
python main_run_efficient.py 9 99 1 'grausnow_iwc'  & #  # this is for eqMass WSM6
python main_run_efficient.py 3 99 1 'grausnow_iwc'  & #  # this is for eqMass WSM6

wait

echo "All done completed"


