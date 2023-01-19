import numpy as np
from matplotlib import pyplot as plt
import os, sys

def run_cmd(str):
    print(f"{str}\n")
    os.system(str)

dt = 0.025; T = 50;
shift =0.01;
Tj1 = 100; Tj2 = 1000;
implement_jitter = 0;
jitt1 = int(Tj1/dt); jitt2 = int(Tj2/dt);

data_direc = "/home/jep/Whalen2021/jitter_profile_tests/varying_shift/1_percent"

run_cmd("g++ -std=c++17 -O2 SNr_jitter.cc -o SNr_jitter")
for T in [50]:
    for implement_jitter in [0,1]:
        for f in [2,15,25]:
            for r in [0,6,7,8,9,10]:
                store_direc = f"{data_direc}/t1_{Tj1}_t2_{Tj2}_freq_{f}/T{T}/random_run_{r}" if implement_jitter == 1 else f"{data_direc}/no_jitter_freq_{f}/T{T}/random_run_{r}"
                run_cmd(f"mkdir -p {store_direc}")

                seed = 749141 if r == 0 else r
                run_cmd(f"./SNr_jitter -arch 2 -dt {dt} -T {T} -sd {seed} -dep_off -fac_off -gstnton 0.15 0.25 -iapp -0.000 -0.000 -glk 0.04 0.04 "+
                f"-c3 0 0 -w 0.05 0.25 -wg .05 .25 -wstn .025 .025 -fixedGABA -70 -osc_freq {f} -fracup_gpe 0.55 -osc_shape_gpe rect " +
                f"-osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop 100 -bip_snr 0.5 "
                f"-write_syn 0 data/competitive -o >data/competitive.hst -jitt1 {jitt1} -jitt2 {jitt2} -shift {shift} -implement_jitter {implement_jitter}")

                freqs = np.loadtxt("frequency_tracker.txt")

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

                fig, axe = plt.subplots(1,1,figsize=(8,6),dpi=300)
                axe.plot(freqs[:,0],freqs[:,1],label="Jitter Freq")
                axe.hlines(np.mean(freqs[:,1]),0,T*1000,color="g",linewidth=2,linestyle="dotted",label="Avg Freq")
                axe.set_ylabel("Frequency"); axe.set_xlabel("Time")
                axe.set_xlim(0,T*1000);
                axe.spines['right'].set_visible(False)
                axe.spines['top'].set_visible(False)
                plt.savefig(f"{store_direc}/freq_profile.pdf")
                plt.close()

                run_cmd(f"cp frequency_tracker.txt {store_direc}/frequency_tracker.txt")
                run_cmd(f"cp data/competitive.* {store_direc}/")