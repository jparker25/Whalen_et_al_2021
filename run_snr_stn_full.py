# John Parker, last edited Jan 2023
# Simulates STN related results (Figure 5 data)
# Refer to git repo for more information
import numpy as np
from matplotlib import pyplot as plt
import time, os, datetime

# Edits below may lead to instability
def run_cmd(cmd):
    # simple function to run commands on the terminal through python
    print(cmd)
    os.system(cmd)

def plot_stn(matrix,save_dir,tt):
    # Saves spikes for STN data
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
    # Generates matrix for delta spikes
    return np.sum(matrix,axis=0)

def spike_matrix(tt,source,rand,Ts):
    # Give sproperly formated spike matrix
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


# Default simulation parameters 
dt = .025; 
T = 50; 
delay = 10; 
npop = 100;

compile_sim = True;
run_sim = True;
gen_gstn_figs = True;
remove_stn = False; # used in test cases
no_connections = False; # used in test cases

sim_time_series=np.arange(0,(T+delay)*1000+dt,dt)

animals = [2,3,4,5] # animals to analyze
batches = [5,4,7,9] # batches from corresponding animals to analyze
record_time = [246,242,213,239]; # Recording times to analyze for corresponding animals

wgs=[6] # default dynamic stn weight
ranges = np.array([[0.13,0.23]]) # params for gstnton for c++ sim
animal_samples = 4 # number of samples for each animal

# params to set up simulations and print progress
rands = [True]
sim = 1
tot_time = 0;
tot_sims = len(animals)*len(ranges)*len(wgs)*animal_samples*len(rands)
animal_sims = len(ranges)*len(wgs)*animal_samples*len(rands)
# Create directory for storage
dir = "stn_sims";
run_cmd(f"mkdir {dir}")
# Run sims
for animal in range(len(animals)): # iterate thorugh all animals
    animal_dir = f"{dir}/animal_{animals[animal]}" # create each animal directory
    run_cmd(f"mkdir {animal_dir}")
    stn_spike_path = f"DD_STN_animals/animal_{animals[animal]}/batch_{batches[animal]}" # Find STN spike path for each experimental results
    print("Animal: ",animals[animal])
    animal_sim = 1

    for rand in rands: # iterate through random or not
        save_dir = f"{animal_dir}/rand_{rand}" # create storage directory
        run_cmd(f"mkdir {save_dir}")
        for samples in range(animal_samples): # itreate through all samples
            sample_dir = f"{save_dir}/sample_{samples+1}"
            run_cmd(f"mkdir {sample_dir}") # create storage directory
            if gen_gstn_figs: # Generate figures if specified
                matrix, randT = spike_matrix(sim_time_series,stn_spike_path,rand,record_time[animal])
                delta = delta_multiplier(matrix)
                np.savetxt(f'delta_multiplier.txt',delta,fmt='%d',delimiter=" ",newline="\n")
                np.savetxt(f'{sample_dir}/delta_multiplier.txt',delta,fmt='%d',delimiter=" ",newline="\n",header=f"randT: {randT:.4f}")
                plot_stn(matrix,sample_dir,sim_time_series)
            else:
                run_cmd(f"cp {sample_dir}/delta_multiplier.txt ./delta_multiplier.txt")
            for rng in ranges: # iterate through all gstnton ranges
                for wg in wgs: # iterate through all dynamic stn values
                    start_time = time.time() # record start time
                    sim_dir = f"{sample_dir}/sim_gstn_rng_{rng[0]}_{rng[1]}_wstn_{wg}"
                    run_cmd(f"mkdir {sim_dir}") # create sim directory

                    # Run simulations based on params
                    if compile_sim and sim == 1:
                        run_cmd("g++ -std=c++17 -O2 snr_stn_input_full.cc -o snr_stn_input_full")

                    if run_sim and not no_connections:
                        run_cmd(f"./snr_stn_input_full -arch 2 -dt {dt} -T {T} -sd 749141 -dep_off -fac_off -gstnton {rng[0] if not remove_stn else 0} {rng[1] if not remove_stn else 0} -wgstn_dyn {wg} -stn_delay {delay} -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w .05 .25 -wg .05 .25 -wstn .025 .025 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop {npop} -bip_snr 0.5 -write_syn 0 {sim_dir}/competitive -o >{sim_dir}/competitive.hst")

                    elif no_connections:
                        run_cmd(f"./snr_stn_input_full -arch 2 -dt {dt} -T {T} -sd 749141 -dep_off -fac_off -gstnton {rng[0] if not remove_stn else 0} {rng[1] if not remove_stn else 0} -wgstn_dyn {wg} -stn_delay {delay} -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w 0 0 -wg 0 0 -wstn 0 0 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop {npop} -bip_snr 0.5 -write_syn 0 {sim_dir}/competitive -o >{sim_dir}/competitive.hst") # only using for special cases


                    run_cmd(f"cp stn_data_ON.txt {sim_dir}/") # Copy data to correpsonding directory

                    # Record end time and print progress of all simulations
                    end_time = time.time()
                    elapsed = end_time - start_time
                    tot_time += elapsed
                    est_rem = ((tot_time)/sim)*(tot_sims-sim)

                    print(f"Completed {sim}/{tot_sims} total; {animal_sim}/{animal_sims} animal")
                    print(f"Sim time: {str(datetime.timedelta(seconds=elapsed))}, Total: {str(datetime.timedelta(seconds=tot_time))}")
                    print(f"Est Rem: {str(datetime.timedelta(seconds=est_rem))}, Time Done: {str(datetime.datetime.now()+datetime.timedelta(seconds=est_rem))}")
                    sim+= 1
                    animal_sim += 1
