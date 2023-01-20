# John Parker, last edited Jan 2023
# Plots jitter related simulations for 15Hz cases (Figure 4)
# No additional analysis is performed
# Refer to git repo for more information
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import os

def run_cmd(str):
    # simple function to run commands on the terminal through python
    print(f"{str}\n")
    os.system(str)


jitter_sims_direc = "/home/jep/Whalen_et_al_2021/jitter_sims" # change to corresponding directory on local machine

# Edits below may lead to instability
results_file = f"{jitter_sims_direc}/jitter_results.csv"  # where CSV file is stored from MATLAB script
T1 = [100,500]; T2 = [1000,2000]; # default jitter params
df = pd.read_csv(results_file) # read in CSV

cols = df.columns # read in pandas data

# EXP data bounds
mfr = [16.6597, 18.2274, 19.9521];
mcv2 = [0.661594, 0.714409, 0.770275];
fap = [0.310595, 0.410526, 0.516192];
fip = [0.24306, 0.336842, 0.441104];
nonOsc = [0.169055, 0.252632, 0.352189];
ratio = [2.01356, 2.76119, 3.81588];
apcv = [0.378699, 0.559011, 0.746587];
ipcv = [0.512115, 0.666985, 0.825455];

# organizing data
exps = [mfr, fap, fip, nonOsc, ratio]
exps_lims = [25, 0.7, 0.7,0.7, 5]
labels = ["FR (Hz)", "Frac. AP", "Frac. IP", "Frac. No Osc.", "Power Ratio"]

# Create overall figure and adjust for visibility
fig = plt.figure(figsize=(16,8),dpi=300)
sub1 = 4; sub2 = 10; sub3 = 8; sub2offset = 2; 
gs = fig.add_gridspec(12,sub1+sub2+sub3+2) 
ax00 = fig.add_subplot(gs[0:6,0:sub1])
ax10 = fig.add_subplot(gs[6:,0:sub1])

# Create panels A and D
axes = [ax00,ax10]
fig_label = ["A","D"]
scale = 1.0;
for i in [0,1]:
    t1 = T1[i]; t2 = T2[i];
    jitter15 = df[(df['Jitter?'] == 1) & (df['Frequency'] == 15) & (df['T1'] == t1) & (df['T2'] == t2)]
    jitter2 = df[(df['Jitter?'] == 1) & (df['Frequency'] == 2) & (df['T1'] == t1) & (df['T2'] == t2)]
    noJ25 = df[(df['Jitter?'] == 0) & (df['Frequency'] == 15)]
    noJ2 = df[(df['Jitter?'] == 0) & (df['Frequency'] == 2)]

    markers, caps, bars = axes[i].errorbar(0,np.mean(jitter15['AP'])*scale,yerr=np.std(jitter15['AP'])*scale,color="tab:orange",marker=">",capsize=2,alpha=1,markersize=12)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]
    markers, caps, bars = axes[i].errorbar(0.15,np.mean(jitter15['IP'])*scale,yerr=np.std(jitter15['IP'])*scale,color="tab:orange",marker=">",capsize=2,alpha=1,markersize=12)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]


    markers, caps, bars = axes[i].errorbar(0,np.mean(jitter2['AP'])*scale,yerr=np.std(jitter2['AP'])*scale,color="g",marker=">",capsize=2,alpha=1,markersize=12)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]
    markers, caps, bars = axes[i].errorbar(0.15,np.mean(jitter2['IP'])*scale,yerr=np.std(jitter2['IP'])*scale,color="g",marker=">",capsize=2,alpha=1,markersize=12)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]

    markers, caps, bars = axes[i].errorbar(0,np.mean(noJ25['AP'])*scale,yerr=np.std(noJ25['AP'])*scale,color="gray",marker="<",capsize=2,alpha=1,markersize=12)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]
    markers, caps, bars = axes[i].errorbar(0.15,np.mean(noJ25['IP'])*scale,yerr=np.std(noJ25['IP'])*scale,color="gray",marker="<",capsize=2,alpha=1,markersize=12)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]


    markers, caps, bars = axes[i].errorbar(0,np.mean(noJ2['AP'])*scale,yerr=np.std(noJ2['AP'])*scale,color="tab:purple",marker="<",capsize=2,alpha=1,markersize=12)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]
    markers, caps, bars = axes[i].errorbar(0.15,np.mean(noJ2['IP'])*scale,yerr=np.std(noJ2['IP'])*scale,color="tab:purple",marker="<",capsize=2,alpha=1,markersize=12)
    [bar.set_alpha(0.5) for bar in bars]
    [cap.set_alpha(0.5) for cap in caps]

    axes[i].scatter(-5,-5,marker=">",color="tab:orange",label="15Hz Jitter",s=120)
    axes[i].scatter(-5,-5,marker=">",color="g",label="2Hz Jitter",s=120)
    axes[i].scatter(-5,-5,marker="<",color="gray",label="15Hz No Jitter",s=120)
    axes[i].scatter(-5,-5,marker="<",color="tab:purple",label="2Hz No Jitter",s=120)


    axes[i].set_xlim(0,50*1000);
    axes[i].spines['right'].set_visible(False)
    axes[i].spines['top'].set_visible(False)
    axes[i].spines['bottom'].set_visible(False)
    axes[i].spines['left'].set_linewidth(3)
    axes[i].tick_params(axis="x",length=0)
    axes[i].tick_params(axis="y",width=3,labelsize=12)
    axes[i].set_xlim([-0.1,0.25]);
    axes[i].set_ylim([0,1.5e-5]);

    axes[i].set_xticks([0,0.15],["SNr AP","SNr IP"],fontsize=12)
    axes[i].set_ylabel("Power",fontsize=12)

    axes[i].text(-0.18,0.95,fig_label[i],fontsize=20,fontweight="bold",transform=axes[i].transAxes)
    axes[0].legend(prop={'size':8})

# Create panels B and E
for i in [0,1]:
    t1 = T1[i]; t2 = T2[i];
    jitter15 = df[(df['Jitter?'] == 1) & (df['Frequency'] == 15) & (df['T1'] == t1) & (df['T2'] == t2)]
    jitter2 = df[(df['Jitter?'] == 1) & (df['Frequency'] == 2) & (df['T1'] == t1) & (df['T2'] == t2)]
    noJ25 = df[(df['Jitter?'] == 0) & (df['Frequency'] == 15)]
    noJ2 = df[(df['Jitter?'] == 0) & (df['Frequency'] == 2)]

    axe = []
    if i == 0:
        axe = [fig.add_subplot(gs[0:6,k:k+sub2offset]) for k in range(sub1+1,sub1+sub2,sub2offset)]
    else:
        axe = [fig.add_subplot(gs[6:,k:k+sub2offset]) for k in range(sub1+1,sub1+sub2,sub2offset)]

    scale = [1,1,1,1,1,1,1,1,1]
    exp_cols = ["MFR","FAP","FIP","NonOsc","Ratio"]
    for j in range(len(exp_cols)):

        axe[j].hlines(exps[j][1],0.85,1.15,color="b",linewidth=2.5)
        axe[j].hlines(exps[j][0],0.90,1.1,color="k",linewidth=1.5)
        axe[j].hlines(exps[j][2],0.90,1.1,color="k",linewidth=1.5)
        axe[j].vlines(1,exps[j][0],exps[j][2],color="k",linewidth=1.5)

        axe[j].scatter(0.95,np.mean(jitter15[exp_cols[j]])*scale[j],color="tab:orange",marker=">",s=160)
        axe[j].scatter(1.05,np.mean(noJ25[exp_cols[j]])*scale[j],color="gray",marker="<",s=160)
        axe[j].scatter(0.95,np.mean(jitter2[exp_cols[j]])*scale[j],color="g",marker=">",s=160)
        axe[j].scatter(1.05,np.mean(noJ2[exp_cols[j]])*scale[j],color="tab:purple",marker="<",s=160)


        axe[j].set_xlim([0.8,1.4]);
        axe[j].set_xticklabels([]);
        axe[j].set_ylim([0,exps_lims[j]]);
        axe[j].spines['right'].set_visible(False)
        axe[j].spines['top'].set_visible(False)
        axe[j].spines['bottom'].set_visible(False)
        axe[j].tick_params(axis="x",length=0)
        axe[j].set_xlabel(labels[j],rotation=10,fontsize=12)
        axe[j].spines['left'].set_linewidth(3)
        axe[j].tick_params(axis="y",width=3,labelsize=10)

        axe[j].xaxis.set_label_coords(0.05, 0)

        axe[0].text(-0.7,0.95,"B" if i == 0 else "E",fontsize=20,fontweight="bold",transform=axe[0].transAxes)


# Create panels C and F
fig_label = ["C","F"]
colors = ["dimgray","rosybrown","darkturquoise","seagreen","darkmagenta","royalblue"]
for i in [0,1]:
    t1 = T1[i]; t2 = T2[i];
    freq = 15;

    if i == 0:
        axe = fig.add_subplot(gs[0:6,sub1+sub2+2:])
    else:
        axe = fig.add_subplot(gs[6:,sub1+sub2+2:])

    count = 0;
    for run in [0,6,7,8,9,10]:
        freqs = np.loadtxt(f"{jitter_sims_direc}/t1_{t1}_t2_{t2}_freq_{freq}/T50/random_run_{run}/frequency_tracker.txt") 

        axe.plot(freqs[:,0],freqs[:,1],label="Jitter Freq",color=colors[count])
        ylims = axe.get_ylim()
        axe.set_ylabel("Frequency",fontsize=12); axe.set_xlabel("Time (s)",fontsize=12)
        axe.set_xlim(0,50*1000);
        axe.spines['right'].set_visible(False)
        axe.spines['top'].set_visible(False)
        axe.tick_params(axis="both",width=3,labelsize=12)
        axe.spines['left'].set_linewidth(3)
        axe.spines['bottom'].set_linewidth(3)
        axe.set_xticks(list(range(0,60000,10000)))
        axe.set_xticklabels([int(x) for x in list(range(0,60,10))])
        axe.text(-0.1,0.95,fig_label[i],fontsize=20,fontweight="bold",transform=axe.transAxes)
        count += 1

# Adjust figure and save
plt.subplots_adjust(hspace=10)
fig.savefig("jitter_figure.eps",bbox_inches='tight')
plt.close()

# Open figure
run_cmd("xdg-open jitter_figure.eps")
