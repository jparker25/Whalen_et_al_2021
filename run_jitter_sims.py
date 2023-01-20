# John Parker, last edited Jan 2023
# Runs all jitter related simulations (Figure 4 data)
# Refer to git repo for more information
import numpy as np
from matplotlib import pyplot as plt
import os

def run_cmd(str):
    # simple function to run commands on the terminal through python
    print(f"{str}\n")
    os.system(str)

data_direc = "/home/jep/Whalen_et_al_2021/jitter_sims" # change to corresponding directory on local machine

# Edits below may lead to instability
dt = 0.025; T = 50; shift =0.01; # default parameters for sims

run_cmd("g++ -std=c++17 -O2 SNr_jitter.cc -o SNr_jitter") # compile c++ code
for T in [50]:
    for tj in [[100,1000],[500,2000]]: # jitter levels
        Tj1 = tj[0]; Tj2 = tj[1];
        jitt1 = int(Tj1/dt); jitt2 = int(Tj2/dt);
        for implement_jitter in [0,1]: # run cases with and without jitter
            for f in [2,15,25]: # run cases with different GPe frequencies
                for r in [0,6,7,8,9,10]: # run cases with different random seeds
                    store_direc = f"{data_direc}/t1_{Tj1}_t2_{Tj2}_freq_{f}/T{T}/random_run_{r}" if implement_jitter == 1 else f"{data_direc}/no_jitter_freq_{f}/T{T}/random_run_{r}" # specify storage directory
                    run_cmd(f"mkdir -p {store_direc}") # make the storage directory

                    seed = 749141 if r == 0 else r # specify the seed (0 gives original default seed)
                    run_cmd(f"./SNr_jitter -arch 2 -dt {dt} -T {T} -sd {seed} -dep_off -fac_off -gstnton 0.15 0.25 -iapp -0.000 -0.000 -glk 0.04 0.04 "+
                    f"-c3 0 0 -w 0.05 0.25 -wg .05 .25 -wstn .025 .025 -fixedGABA -70 -osc_freq {f} -fracup_gpe 0.55 -osc_shape_gpe rect " +
                    f"-osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 "
                    f"-write_syn 0 data/competitive -o >data/competitive.hst -jitt1 {jitt1} -jitt2 {jitt2} -shift {shift} -implement_jitter {implement_jitter}") # run the c++ code with jitter params

                    freqs = np.loadtxt("frequency_tracker.txt") # read in the jitter frequency

                    # Plot and store the jitter for the sim
                    fig, axe = plt.subplots(2,1,figsize=(8,6),dpi=300)
                    axe[0].plot(freqs[:,0],freqs[:,1],label="Jitter Freq")
                    ylims = axe[0].get_ylim()
                    axe[0].vlines(np.arange(0,T*1000+jitt1*dt,jitt1*dt),ylims[0],ylims[1],color="k",linestyle="dashed",linewidth=0.5,label="Freq Change")
                    axe[0].vlines(np.arange(0,T*1000+jitt2*dt,jitt2*dt),ylims[0],ylims[1],color="r",linestyle="dashed",linewidth=0.5,label="New Shift")
                    axe[0].hlines(np.mean(freqs[:,1]),0,T*1000,color="g",linewidth=2,linestyle="dotted",label="Avg Freq")
                    axe[0].set_ylabel("$f$"); axe[0].set_xlabel("$t$")
                    axe[0].set_xlim(0,T*1000);
                    axe[0].legend()

                    axe[1].plot(freqs[1:,0],np.diff(freqs[:,1]))
                    axe[1].set_ylabel("$\Delta f$"); axe[1].set_xlabel("$t$")
                    axe[1].set_xlim(0,T*1000);


                    plt.suptitle(f"Frequency ({f}) - Jitter over time ")
                    figure_title = f"t1_{Tj1}_t2_{Tj2}_osc_freq_{f}_run_{r}.pdf" if implement_jitter == 1 else f"no_jitter_osc_freq_{f}_run_{r}.pdf"
                    plt.savefig(f"{store_direc}/{figure_title}")
                    plt.close()

                    # Plot and store the jitter for the sim
                    fig, axe = plt.subplots(1,1,figsize=(8,6),dpi=300)
                    axe.plot(freqs[:,0],freqs[:,1],label="Jitter Freq")
                    axe.hlines(np.mean(freqs[:,1]),0,T*1000,color="g",linewidth=2,linestyle="dotted",label="Avg Freq")
                    axe.set_ylabel("Frequency"); axe.set_xlabel("Time")
                    axe.set_xlim(0,T*1000);
                    axe.spines['right'].set_visible(False)
                    axe.spines['top'].set_visible(False)
                    plt.savefig(f"{store_direc}/freq_profile.pdf")
                    plt.close()

                    # Copy output data from c++ code to designated directory
                    run_cmd(f"cp frequency_tracker.txt {store_direc}/frequency_tracker.txt")
                    run_cmd(f"cp data/competitive.* {store_direc}/")