import os, sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.signal import find_peaks
import time
import datetime

def run_cmd(cmd):
    print(cmd)
    os.system(cmd)
    print()

def plot_stn(matrix,save_dir,tt):
    fig, ax = plt.subplots(2,1,figsize=(8,6),dpi=300)
    for i in range(matrix.shape[0]):
        for k in range(matrix.shape[1]):
            if matrix[i,k] == 1:
                ax[0].scatter(tt[k],i+1,color='k',s=2)
    delta = delta_multiplier(matrix)
    ax[1].plot(tt,delta)
    fig.savefig(f'{save_dir}/spikes.pdf')
    plt.close()

def delta_multiplier(matrix):
    return np.sum(matrix,axis=0)

def spike_matrix(tt,source,rand,Ts):
    n = len([sd for sd in os.listdir(source) if 'spikes_' in sd])
    matrix = np.zeros((n,len(tt)))
    randT = np.random.rand()*(Ts*1000-tt[-1]);
    for i in range(1,n+1):
        spikes = np.loadtxt(f'{source}/spikes_{i}.txt')*1000
        if rand:
            spikes = spikes[np.where((spikes >= randT) & (spikes < randT+tt[-1]))] - randT
        for t in range(len(tt)-1):
            if len(np.where((spikes>=tt[t]) & (spikes < tt[t+1]))[0]) > 0:
                matrix[i-1,t]=len(np.where((spikes>=tt[t]) & (spikes < tt[t+1]))[0])
    return matrix, randT


dt = .025; # millisecons
T = 50; # simtime seconds
delay = 10 # seconds
npop = 100;

compile_sim = True;
run_sim = True;
gen_gstn_figs = True;
remove_stn = False
no_connections = False;

sim_time_series=np.arange(0,(T+delay)*1000+dt,dt)

animals = [2,3,4,5]
batches = [5,4,7,9]
record_time = [246,242,213,239];

wgs=[6]
ranges = np.array([[0.13,0.23]])
animal_samples = 4

wgs = [6]


rands = [True]
sim = 1
tot_time = 0;
tot_sims = len(animals)*len(ranges)*len(wgs)*animal_samples*len(rands)
animal_sims = len(ranges)*len(wgs)*animal_samples*len(rands)
dir = "large_run";
run_cmd(f"mkdir {dir}")
for animal in range(len(animals)):
    animal_dir = f"{dir}/animal_{animals[animal]}"
    run_cmd(f"mkdir {animal_dir}")
    stn_spike_path = f"DD_STN_animals/animal_{animals[animal]}/batch_{batches[animal]}"
    print("Animal: ",animals[animal])
    animal_sim = 1

    for rand in rands:
        save_dir = f"{animal_dir}/rand_{rand}"
        run_cmd(f"mkdir {save_dir}")
        for samples in range(animal_samples):
            sample_dir = f"{save_dir}/sample_{samples+1}"
            run_cmd(f"mkdir {sample_dir}")
            if gen_gstn_figs:
                matrix, randT = spike_matrix(sim_time_series,stn_spike_path,rand,record_time[animal])
                delta = delta_multiplier(matrix)
                np.savetxt(f'delta_multiplier.txt',delta,fmt='%d',delimiter=" ",newline="\n")
                np.savetxt(f'{sample_dir}/delta_multiplier.txt',delta,fmt='%d',delimiter=" ",newline="\n",header=f"randT: {randT:.4f}")
                plot_stn(matrix,sample_dir,sim_time_series)
            else:
                run_cmd(f"cp {sample_dir}/delta_multiplier.txt ./delta_multiplier.txt")
            for rng in ranges:
                for wg in wgs:
                    start_time = time.time()
                    sim_dir = f"{sample_dir}/sim_gstn_rng_{rng[0]}_{rng[1]}_wstn_{wg}"
                    run_cmd(f"mkdir {sim_dir}")
                    if compile_sim and sim == 1:
                        run_cmd("g++ -std=c++17 -O2 snr_stn_input_full.cc -o snr_stn_input_full")

                    if run_sim and not no_connections:
                        run_cmd(f"./snr_stn_input_full -arch 2 -dt {dt} -T {T} -sd 749141 -dep_off -fac_off -gstnton {rng[0] if not remove_stn else 0} {rng[1] if not remove_stn else 0} -wgstn_dyn {wg} -stn_delay {delay} -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .05 .25 -wg .05 .25 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop {npop} -bip_snr 0.5 -write_syn 0 {sim_dir}/competitive -o >{sim_dir}/competitive.hst")

                    elif no_connections:
                        run_cmd(f"./snr_stn_input_full -arch 2 -dt {dt} -T {T} -sd 749141 -dep_off -fac_off -gstnton {rng[0] if not remove_stn else 0} {rng[1] if not remove_stn else 0} -wgstn_dyn {wg} -stn_delay {delay} -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w 0 0 -wg 0 0 -wstn 0 0 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop {npop} -bip_snr 0.5 -write_syn 0 {sim_dir}/competitive -o >{sim_dir}/competitive.hst") # only using for special cases

                    run_cmd(f"cp stn_data_ON.txt {sim_dir}/")
                    end_time = time.time()
                    elapsed = end_time - start_time
                    tot_time += elapsed
                    est_rem = ((tot_time)/sim)*(tot_sims-sim)

                    print(f"Completed {sim}/{tot_sims} total; {animal_sim}/{animal_sims} animal")
                    print(f"Sim time: {str(datetime.timedelta(seconds=elapsed))}, Total: {str(datetime.timedelta(seconds=tot_time))}")
                    print(f"Est Rem: {str(datetime.timedelta(seconds=est_rem))}, Time Done: {str(datetime.datetime.now()+datetime.timedelta(seconds=est_rem))}")
                    sim+= 1
                    animal_sim += 1
