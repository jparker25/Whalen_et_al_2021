import numpy as np
from matplotlib import pyplot as plt
import os, sys

def run_cmd(str):
    print(f"{str}\n")
    os.system(str)

t1 = 500; t2 = 2000; T = 50;
runs = [0,6,7,8,9,10];


T2 = [1000,2000]; T1 = [100,500];
for i in [0,1]:
    t1 = T1[i]; t2 = T2[i];
    for run in runs:
        for freq in [2,15,25]:
            data_direc = f"/home/jep/Whalen2021/jitter_profile_tests/varying_shift/1_percent/t1_{t1}_t2_{t2}_freq_{freq}/T{T}/random_run_{run}"
            run_cmd(f"cp {data_direc}/frequency_tracker.txt /home/jep/Whalen2021/jitter_profile_tests/results_summary/NEW/frequency_profile_t1_{t1}_t2_{t2}_freq_{freq}_run_{run}.txt")
