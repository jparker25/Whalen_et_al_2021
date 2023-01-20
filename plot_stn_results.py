import numpy as np
from scipy.signal import find_peaks
from matplotlib import pyplot as plt
import seaborn as sns
import os, sys
import matplotlib.gridspec as gridspec
import scipy.integrate as integrate
import pandas as pd

def run_cmd(cmd):
    print(cmd)
    os.system(cmd)

def get_freq(data):
    freq = np.zeros(data.shape[1]-1)
    for n in range(len(freq)):
        pks,_ = find_peaks(data[:,n+1],height=(-20))
        freq[n] = len(pks)/(data[-1,0]);
    return np.mean(freq)

def paper_figure(sim_csv,real_csv,w):
    sims=pd.read_csv(sim_csv)
    real=pd.read_csv(real_csv)

    stats = ['mean_rate','frac_ap','frac_ip','apip_power_ratio']
    fig_labels = ["FR (Hz)","Frac. AP","Frac. IP","Power Ratio"]
    fig_alpha = ["B","C","D","E"]

    fig = plt.figure(figsize=(8,6),dpi=300)
    gs = fig.add_gridspec(4,8)
    hist_axes = [fig.add_subplot(gs[2:,0:2]),fig.add_subplot(gs[2:,2:4]),fig.add_subplot(gs[2:,4:6]),fig.add_subplot(gs[2:,6:])]
    summary_axes = [fig.add_subplot(gs[0:2,i]) for i in range(8)]
    for stat in range(len(stats)):
        clow = real[real["type"]=="conf_low"][stats[stat]]
        chigh = real[real["type"]=="conf_high"][stats[stat]]
        rmean = real[real["type"]=="mean"][stats[stat]]
        sns.histplot(data=sims[sims["wstn_static"]==w],x=stats[stat],ax=hist_axes[stat],color="gray")
        hist_axes[stat].set_xlabel(fig_labels[stat])
        ylims = hist_axes[stat].get_ylim()
        hist_axes[stat].vlines(rmean,ylims[1]-0.25,ylims[1]+0.25,color="b")
        hist_axes[stat].vlines(clow,ylims[1]-0.25,ylims[1]+0.25,color="k")
        hist_axes[stat].vlines(chigh,ylims[1]-0.25,ylims[1]+0.25,color="k")
        hist_axes[stat].hlines(ylims[1],clow,chigh,color="k")

        for k in ['left','right','top','bottom']:
            hist_axes[stat].spines[k].set_visible(False if ((k== 'right') or (k=='top')) else True)
            hist_axes[stat].spines[k].set_linewidth(3)
            hist_axes[stat].tick_params('both',width=3)
            hist_axes[stat].text(0.05,1,fig_alpha[stat],transform=hist_axes[stat].transAxes,fontweight="bold",fontsize=16)

    labels = list(real.columns)
    fig_labels = ["FR (Hz)","CV2 ISI","Frac. AP","Frac. IP", "Frac. No Osc.","Power Ratio","AP Power CV","IP Power CV"]

    lims=[[0,25],[0,1],[0,0.7], [0,0.7],  [0,0.7],  [0,5], [0,1],[0,1]]

    animal_color={2:"r",3:"g",4:"orange",5:"purple"}
    for i in range(len(summary_axes)):

        summary_axes[i].hlines(real[real["type"]=="mean"][labels[i]],0.8,1.2,color="b",linewidth=2.5)
        summary_axes[i].hlines(real[real["type"]=="conf_low"][labels[i]],0.85,1.15,color="k",linewidth=1.5)
        summary_axes[i].hlines(real[real["type"]=="conf_high"][labels[i]],0.85,1.15,color="k",linewidth=1.5)
        summary_axes[i].vlines(1,real[real["type"]=="conf_low"][labels[i]],real[real["type"]=="conf_high"][labels[i]],color="k",linewidth=1.5)
        summary_axes[i].set_ylim(lims[i])
        summary_axes[i].set_xlim([0.7,1.3])
        summary_axes[i].set_xlabel(fig_labels[i],rotation=15)
        summary_axes[i].set_xticklabels([])

        for animal in range(2,6):
            data = np.mean(sims[(sims["animal"]==animal) & (sims["wstn_static"]==w)][labels[i]])
            summary_axes[i].scatter(1.1,data,color=animal_color[animal],marker="<",s=80)

        for k in ['left','right','top','bottom']:
            summary_axes[i].spines[k].set_visible(False if k!= 'left' else True)
            summary_axes[i].spines[k].set_linewidth(3)
            summary_axes[i].tick_params('y',width=3)
            summary_axes[i].tick_params('x',width=0)

        if i == 0:
            summary_axes[i].text(0.05,1,"A",transform=summary_axes[i].transAxes,fontweight="bold",fontsize=16)

    fig.tight_layout()
    fig.savefig(f"paper_fig_{w}.eps")
    plt.close()

    run_cmd(f"xdg-open paper_fig_{w}.eps")

def supplemental_figure(sim_dir,w):
    spks = np.loadtxt(f"{sim_dir}/sim_gstn_rng_0.13_0.23_wstn_{w}/competitive.hst")
    gstns = np.loadtxt(f"{sim_dir}/sim_gstn_rng_0.13_0.23_wstn_{w}/stn_data_ON.txt")
    rates = np.array([spks[(spks[:,1]==n) & (spks[:,0] > 3)].shape[0]/(50-3) for n in range(100)])
    gstatics = gstns[0,4::5]
    fig,ax = plt.subplots(1,1,figsize=(8,6),dpi=300)

    ax.scatter(gstatics,rates,color="b")
    ax.set_xlabel(r"$g_{STN}^{static}$ (nS/pF)")
    ax.set_ylabel(r"$f$ (Hz)")

    sns.despine()
    fig_labels = ["A","B"]
    for k in ['left','right','top','bottom']:
            ax.spines[k].set_linewidth(3)
            ax.tick_params('both',width=3)
    fig.tight_layout()
    fig.savefig("supplemental_figure1.eps")
    plt.close()

    run_cmd("xdg-open supplemental_figure1.eps")


paper_figure("stn_sims_results.csv","real_data_results.csv",6)
supplemental_figure("stn_sims/animal_2/rand_True/sample_3",6)

sys.exit()
