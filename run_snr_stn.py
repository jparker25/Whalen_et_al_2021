import os, sys
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.signal import find_peaks

def run_cmd(cmd):
    print(cmd)
    os.system(cmd)

def plot_stn(matrix,save_dir):
    fig, ax = plt.subplots(2,1,figsize=(8,6),dpi=300)
    for i in range(matrix.shape[0]):
        for k in range(matrix.shape[1]):
            if matrix[i,k] == 1:
                ax[0].scatter(time[k],i+1,color='k',s=2)
    delta = delta_multiplier(matrix)
    ax[1].plot(time,delta)
    fig.savefig(f'{save_dir}/spikes.pdf')
    plt.close()

def delta_multiplier(matrix):
    return np.sum(matrix,axis=0)

def spike_matrix(time,source,rand,Ts):
    n = len([sd for sd in os.listdir(source) if 'spikes_' in sd])
    matrix = np.zeros((n,len(time)))
    randT = np.random.rand()*(Ts*1000-time[-1]);
    for i in range(1,n+1):
        spikes = np.loadtxt(f'{source}/spikes_{i}.txt')*1000
        if rand:
            spikes = spikes[np.where((spikes >= randT) & (spikes < randT+time[-1]))] - randT
        for t in range(len(time)-1):
            if len(np.where((spikes>=time[t]) & (spikes < time[t+1]))[0]) > 0:
                matrix[i-1,t]=len(np.where((spikes>=time[t]) & (spikes < time[t+1]))[0])
    return matrix


dt = .025; # millisecons
T = 20; # simtime seconds
delay = 10 # seconds
npop = 6;
rand2 = False;

compile_sim = False; run_sim = True;
gen_gstn_figs = True;
remove_stn = False

gstn_static = 0.025;
wgstn = 2;

gstatics = [0.02,0.025,0.03,0.035]
wgs = [1.5,2,2.5,3]
gstatics=[0.2]
wgs = [2]
time=np.arange(0,(T+delay)*1000+dt,dt)

animal = 5;
batch = 9;
choice = np.random.randint(0,high=3)
animals = [2,3,4,5]
batches = [5,4,7,9]
animals=[3]
batches=[4]
record_time = [246,242,213,239];
record_time=[242];

main_save_dir = "single_neurons_test"
run_cmd(f"mkdir {main_save_dir}")
for animal in range(len(animals)):
    stn_spike_path = f"DD_STN_animals/animal_{animals[animal]}/batch_{batches[animal]}"


    tot_sims = 0
    animal_save_dir = f"{main_save_dir}/animal_{animal+2}"
    run_cmd(f"mkdir {animal_save_dir}")

    for rand in [True]:
        save_dir = f"{animal_save_dir}/rand_{rand}"
        run_cmd(f"mkdir {save_dir}")
        sims = 0
        for gstn in gstatics:
            for wg in wgs:
                sims += 1
                tot_sims += 1
                if gen_gstn_figs and sims == 1:
                    matrix = spike_matrix(time,stn_spike_path,rand,record_time[animal])
                    delta = delta_multiplier(matrix)
                    np.savetxt(f'delta_multiplier.txt',delta,fmt='%d',delimiter=" ",newline="\n")
                    np.savetxt(f'{save_dir}/delta_multiplier.txt',delta,fmt='%d',delimiter=" ",newline="\n")
                    plot_stn(matrix,save_dir)

                if compile_sim and sims == 1:
                    run_cmd("g++ -std=c++17 -O2 snr_stn_input.cc -o snr_stn_input")

                if run_sim:
                    run_cmd(f"./snr_stn_input -arch 2 -dt {dt} -T {T} -sd 749141 -dep_off -fac_off -gstnton {0.15 if not remove_stn else 0} {0.25 if not remove_stn else 0} -gton_stn_static {gstn} -wgstn_dyn {wg} -stn_delay {delay} -iapp -0.000 -0.000 -glk 0.04 0.04 -c3 0 0 -w 0 0 -wg 0 0 -wstn 0 0 -fixedGABA -70 -osc_freq 2 -fracup_gpe 0.55 -osc_shape_gpe rect -osc_mod_gpe 24 24 -osc_mod_gpe2 24 24 -osc_cent_gpe 25 -n_delays_gpe 1 -delay_std_gpe 34.6164 -n_pop {npop} -bip_snr 0.5 -write_syn 0 data/competitive -o >data/competitive.hst")

                    run_cmd(f'mv stn_data_{"ON" if not remove_stn else "OFF"}.txt {save_dir}/gstn_{gstn}_wg_{wg}.txt')
                    print(f"Sim {tot_sims}/{2*len(animals)*len(gstatics)*len(wgs)} completed")


        freqs = np.zeros((len(gstatics),len(wgs)))
        for i in range(freqs.shape[0]):
            for j in range(freqs.shape[1]):
                data = np.loadtxt(f"{save_dir}/gstn_{gstatics[i]}_wg_{wgs[j]}.txt")
                data = data[np.where(data[:,0] >= 0)[0],:]
                vd = data[:,2]
                pks, _ = find_peaks(vd,height=(-25))
                freqs[i,j] = len(pks)/T

        fig,ax = plt.subplots(1,1,figsize=(8,6),dpi=300)
        hm=sns.heatmap(freqs,xticklabels=wgs,vmin=4,vmax=16,yticklabels=gstatics,annot=True,fmt=".2f",ax=ax,cmap="jet",cbar_kws={'label': r'$f$ (Hz)','orientation':'horizontal'})
        ax.set_yticklabels(gstatics,rotation=0)
        ax.set_ylabel(r"$g_{STN}^{static}$",rotation=0)
        ax.set_xlabel(r"$W_{STN}$")

        fig.savefig(f"{save_dir}/gstn_wg_freq_heatmap.pdf")
        plt.close()


sys.exit()


stn_data = np.loadtxt(f'stn_data_{"ON" if not remove_stn else "OFF"}.txt')
t = stn_data[:,0]
#t = np.arange(-delay*1000,T*1000,dt)
print(stn_data.shape)
for neuron in range(npop):
    v = stn_data[:,1+5*neuron];
    vd = stn_data[:,2+5*neuron];
    gstn = stn_data[:,3+5*neuron];
    gstn_static = stn_data[:,4+5*neuron];
    wgstn = stn_data[:,5+5*neuron];

    pks, _ = find_peaks(v,height=(-25))
    print(len(pks)/T)
    fig, ax = plt.subplots(3,1,figsize=(8,4),dpi=300)
    ax[0].plot(t,v,color="blue")
    ax[1].plot(t,vd,color="blue")
    ax[2].plot(t,(gstn+gstn_static)*10,color="blue")

    print(np.mean(10*(gstn+gstn_static)))

    ax[0].set_ylabel(r"$V^S_m$ (mV)"); ax[0].set_xlabel(r"$t$ (ms)")
    ax[1].set_ylabel(r"$V^D_m$ (mV)"); ax[1].set_xlabel(r"$t$ (ms)")
    ax[0].set_ylim([-70,20]); ax[1].set_ylim([-70,20])
    ax[2].set_ylabel(r"$g_{STN}$ (nS/pF)"); ax[2].set_xlabel(r"$t$ (ms)")
    #ax[2].set_ylim([0,3);
    ax[2].hlines(0.5,0,t[-1],color="gray",linestyle="dashed")
    ax[2].hlines(2.5,0,t[-1],color="gray",linestyle="dashed")
    sns.despine()


    for i in ["left","right","top","bottom"]:
        for k in range(3):
            ax[k].spines[i].set_linewidth(3);
            ax[k].tick_params("both",width=3)
            #ax[k].set_xlim([0,time[-1]])
    plt.tight_layout()
    fig.savefig(f"v_gpe_stn_{'ON' if not remove_stn else 'OFF'}_neuron_{neuron+1}.pdf")
    plt.close()

run_cmd(f"xdg-open v_gpe_stn_{'ON' if not remove_stn else 'OFF'}_neuron_6.pdf")

'''
gstn =np.zeros(len(time));
gstn[0] = np.random.rand()*0.25;
wstn = 2; taustn = 3;
for i in range(1,len(time)):
    gstn[i] = gstn[i-1]+dt*(-gstn[i-1]/taustn + wstn*delta[i-1])

fig, ax = plt.subplots(1,1,figsize=(8,6),dpi=300)
ax.plot(time,10*gstn)
ax.set_ylabel(r"$g_{STN}$ (nS/pF)")
ax.set_xlabel(r"$t$ (ms)");
ax.set_ylim([0,3]);
ax.hlines(0.5,0,time[-1],color="gray",linestyle="dashed")
ax.hlines(2.5,0,time[-1],color="gray",linestyle="dashed")
fig.savefig("gstn_de_sim.pdf")
plt.close()
'''
