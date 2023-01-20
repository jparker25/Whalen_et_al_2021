% Tim C. Whalen, last edited May 2021
% Minor changes by John Parker, Jan 2023
% Plots results of SNr.cc simulations and compares to in vivo data
% Refer to git repo for more information
% Ignores SNr neurons which fire slower than 5 Hz (should be none with
% proper tuning) to match experimental methods from Whalen et al. 2020 JNP

%% TO RUN
% To compare against in vivo data, requires the following data from the
% jparker25/Whalen_et_al_2021 github repository
%   Whalen2021_SNr_pow_by_ecog.mat
%   Whalen2021_SNr_conf_depleted.mat
%   Whalen2021_SNr_conf_control.mat
% Also requires the output obtained by running SNr.cc, see file for details
%
% Specify which figure you want to recreate by uncommenting the section
% below
% Or, set the values yourself:
% infile: the filename with no extension) outputted by SNr.cc
% control: 1 if control (pacemaking, not osc GPe) simulation (-control in .cc)
% T: time length of recording (-T in .cc)
% force_freq: the frequency of oscillating GPe neurons (-osc_freq in .cc)
% search_freqs: [lo, hi] frequencies to search over for significant
%   oscillations (should include force_freq!)
% n_snr: number of SNr neurons (also number of GPe neurons, which are
%   typically split into two equal-sized populations. -n_pop in .cc)
% scatter_colors: 1 to plot two SNr populations with different colors (only
%   makes sense for full adn partial segregated models)

%% Choose a figure to re-create by uncommenting a block

% FIGURE 2C-D: Fully segregated pathways
% infile = 'full_seg';
% control = 0;
% scatter_colors = 1;

% % %% FIGURE 2F-G: Partially-segregated pathways
% infile = 'partial_seg';
% control = 0;
% scatter_colors = 1;

% %% FIGURE 3B-C: Competitive model 
% Also run this before running Whalen2021_plot_balance)
infile = 'competitive';
control = 0;
scatter_colors = 0;

% %% FIGURE 3D: Competitive model in healthy conditions
% % Note: the oscillation detection method has ~5% false positive rate. As
% % such, ~5 oscillating neurons may be detected even if they exhibit no
% % underlying oscillatory process. Checking if each neuron is oscillating at
% % the same frequency in the detection band (here, 0.5-4Hz) is a good sanity
% % check as to whether a consistent oscillatory process is occurring across
% % the entire population
% infile = 'competitive_control';
% control = 1;
% scatter_colors = 0;

% %% FIGURE 4A-B: Competitive model with QIF neurons
% % If running Whalen2021_plot_qif_noise for Figure 2C, keep ALL input
% % variable sets (including this one) commented.
% infile = 'competitive_qif';
% control = 0;
% scatter_colors = 0;

%% OTHER PARAMETERS TO CHANGE if you've run SNr.cc with custom inputs
T = 50;
force_freq = 2; % Comment out for certain MATLAB files
search_freqs = [0.5, 4]; % Comment out for certain MATLAB files
n_snr = 100; % assumes there are two GPe populations (oscillaitng and poisson) each of size n_snr/2


%% Main code - edits below may lead to instability

[ts] = import_cc_output([infile '.hst'],3, n_snr);
[ts_gpe] = import_cc_output([infile '.gpe'],3, n_snr);

sdf_std = .05; % sec, smoothing factor for SDF
dt = 0.001; % downsample to this timestep for SDF
osc_only = 1; % still calculates both w/ (pttype_all) or w/o (pttype_osc) non-oscillating neurons, but sets the main pttype variable to desired version
npop = n_snr/2; % size of each SNr population
gpe_inds = 1:npop; % index of GPe neuron to determine AP/IP from - average over all

% create data struct for inut to oscillation detector
data = struct();
data.files = {infile};
data.ts{1}=ts; data.T=T; data.rates{1}=cellfun(@(x) size(x,1)/(T-3), ts);
data.nfiles=1; data.animalcodes=1;
data.files_trunc = ['SNr_' infile(1:end)];

data.osc = struct();
data.osc.min_rate = 5;
data.osc.to_plot = 0;
data.osc.srch_lo = search_freqs(1);
data.osc.srch_hi = search_freqs(2);
[ data ] = renewalPSD_phaseShift_batch( data);

[~, freq_ind] = min(abs(data.osc.freqs-force_freq)); % index closest to forcing frequency
ts_min = data.ts{1}(data.rates{1}>data.osc.min_rate);
rates_min = data.rates{1}(data.rates{1}>data.osc.min_rate);

pows = data.osc.psd_unc{1}(:,freq_ind)./rates_min;
pop_inds = [zeros(npop,1); ones(npop,1)]; % which population each SNr neuron in
pop_inds_min = pop_inds(data.rates{1}>data.osc.min_rate);

% find AP vs IP
[~, gpe_sdf] = gaussSpikeDensity( cell2mat(ts_gpe(gpe_inds)), sdf_std, 0, T);

lo = 0; hi = T;
nu = size(ts_min,1);
sdftimes = (lo+sdf_std*3):dt:(hi-sdf_std*3); % from gaussSpikeDensity
sdfs = zeros(nu,length(sdftimes));

maxlag = floor((1/force_freq)/1.2*1000);
xcs = cell(nu,1);
sig1 = gpe_sdf-mean(gpe_sdf); % GPe SDF for comparison against all SNr neurons

% pttype stands for peak-trough-type - 1 for peak relationship with GPe, -1
% for trough, 0 for no significant relationship
pttype_osc = zeros(nu,1);
pttype_all = zeros(nu,1); % like above, but excludes 0's - finds peak or trough even if non-significant
for u=1:nu
    [~, sdfs(u,:)] = gaussSpikeDensity( ts_min{u}, sdf_std, 0, T);
    sig2 = sdfs(u,:)-mean(sdfs(u,:));
    xcs{u} = xcorr(sig1,sig2,maxlag,'normalized');
    [ma, maind] = max(xcs{u});
    [mi, miind] = min(xcs{u});
    domax = abs(maind-(maxlag+1)) < abs(miind-(maxlag+1)); % 1 if peak xcorr w/ GPe, - if trough
    pttype_osc(u) = (~domax) - domax;
    pttype_all(u) = (~domax) - domax;
    if ~(data.osc.has_osc{1}(u))
        pttype_osc(u)=0;
    end
end

if osc_only
    pttype = pttype_osc;
else
    pttype = pttype_all;
end

pows_ap = pows(pttype==1);
pows_ip = pows(pttype==-1);
mean_ap = mean(pows_ap);
mean_ip = mean(pows_ip);

ap_osc = pows(pop_inds_min==0 & pttype==1);
ip_osc = pows(pop_inds_min==0 & pttype==-1);
ap_pois = pows(pop_inds_min==1 & pttype==1);
ip_pois = pows(pop_inds_min==1 & pttype==-1);

%% Plot powers by AP/IP
if scatter_colors
    dg = [0 .7 0];
else
    dg = 'k';
end
[mean_ap mean_ip]
figure
hold on
if ~isnan(mean_ap)
    plot([-.05 .05], [mean_ap mean_ap],'k','LineWidth',2)
end
if ~isnan(mean_ip)
    plot([.95 1.05], [mean_ip mean_ip],'Color',dg,'LineWidth',2)
end
if ~isempty(ap_osc)
    scatter(rand(size(ap_osc))*.1-.05,ap_osc,200,'k','.')
end
if ~isempty(ap_pois)
    scatter(rand(size(ap_pois))*.1-.05,ap_pois,200,dg,'.')
end
if ~isempty(ip_osc)
    scatter(rand(size(ip_osc))*.1-.05+1,ip_osc,200,'k','.')
end
if ~isempty(ip_pois)
    scatter(rand(size(ip_pois))*.1-.05+1,ip_pois,200,dg,'.')
end

ylabel('Power')
set(gca,'xtick',[0 1],'xticklabels',{'SNr AP';'SNr IP'})
makeNice(gca)
set(gcf,'Position',[10 400 430 300])

% some stats on power distirbutions
meanpow_ratio_sim = mean_ap/mean_ip;
CV_AP_sim = std(pows_ap)/mean_ap;
CV_IP_sim = std(pows_ip)/mean_ip;

mean_rate_sim = mean(data.rates{1});
cv_rate_sim = std(data.rates{1})./mean_rate_sim;

cv2isi_sim = cellfun(@CV2ISI,getit(data.ts))';
mean_cv2isi_sim = mean(cv2isi_sim);
cv_cv2isi_sim = std(cv2isi_sim)/mean_cv2isi_sim;

%% Bootstrap stats and plot confidence intervals
seed = 3458890;
rng(seed)
nMC = 1000;
alpha = .025; % two sides, alpha*2 gives 1-conf
i1 = floor(nMC*alpha)+1; i2 = ceil(nMC*(1-alpha))+1;
tickangle = 20; % xlabel plotting angle

if control
    load Whalen2021_SNr_conf_control
    rate_boot_real = bootsample(rates_control, nMC);
    cv2isi_boot_real = bootsample(cv2isi_control,nMC);
    
    means_rate_real = mean(rate_boot_real,2);
    cvs_rate_real = sort(std(rate_boot_real,0,2)./means_rate_real);
    means_rate_real = sort(means_rate_real);
    
    means_cv2isi_real = mean(cv2isi_boot_real,2);
    cvs_cv2isi_real = sort(std(cv2isi_boot_real,0,2)./means_cv2isi_real);
    means_cv2isi_real = sort(means_cv2isi_real);
    
    mean_rate_conf_real = means_rate_real([i1 i2]);
    cvs_rate_conf_real = cvs_rate_real([i1 i2]);
    mean_cv2isi_conf_real = means_cv2isi_real([i1 i2]);
    cvs_cv2isi_conf_real = cvs_cv2isi_real([i1 i2]);
    
    figure
    subplot(1,2,1)
    plot_conf(mean(rates_control), mean_rate_sim, mean_rate_conf_real, [0 35], 'Mean Firing Rate (Hz)', tickangle)
    
    subplot(1,2,2)
    plot_conf(mean(cv2isi_control), mean_cv2isi_sim, mean_cv2isi_conf_real, [0 1], 'Mean CV2 ISI', tickangle)

else
    load Whalen2021_SNr_pow_by_ecog
    load Whalen2021_SNr_conf_depleted
    ap_boot_real = bootsample(peakpows_r, nMC);
    ip_boot_real = bootsample(troughpows_r, nMC);
    rate_boot_real = bootsample(rates_real, nMC);
    cv2isi_boot_real = bootsample(cv2isi_real,nMC);
    
    means_rate_real = mean(rate_boot_real,2);
    cvs_rate_real = sort(std(rate_boot_real,0,2)./means_rate_real);
    means_rate_real = sort(means_rate_real);
    
    means_cv2isi_real = mean(cv2isi_boot_real,2);
    cvs_cv2isi_real = sort(std(cv2isi_boot_real,0,2)./means_cv2isi_real);
    means_cv2isi_real = sort(means_cv2isi_real);
    
    means_ap_real = mean(ap_boot_real,2);
    means_ip_real = mean(ip_boot_real,2);
    cvs_ap_real = sort(std(ap_boot_real,[],2)./means_ap_real);
    cvs_ip_real = sort(std(ip_boot_real,[],2)./means_ip_real);
    meanrat_real = sort(means_ap_real./means_ip_real);
    
    meanrat_conf_real = meanrat_real([i1 i2]);
    cvs_ap_conf_real = cvs_ap_real([i1 i2]);
    cvs_ip_conf_real = cvs_ip_real([i1 i2]);
    mean_rate_conf_real = means_rate_real([i1 i2]);
    cvs_rate_conf_real = cvs_rate_real([i1 i2]);
    mean_cv2isi_conf_real = means_cv2isi_real([i1 i2]);
    cvs_cv2isi_conf_real = cvs_cv2isi_real([i1 i2]);
    
    zz=[mean_rate_sim mean_cv2isi_sim   sum(pttype==1)/size(pttype,1) sum(pttype==-1)/size(pttype,1) sum(pttype==0)/size(pttype,1) meanpow_ratio_sim CV_AP_sim CV_IP_sim]
    figure
    subplot(1,8,1)
    plot_conf(mean(rates_real), mean_rate_sim, mean_rate_conf_real, [0 25], 'Mean Firing Rate (Hz)', tickangle)
    
    subplot(1,8,2)
    plot_conf(mean(cv2isi_real), mean_cv2isi_sim, mean_cv2isi_conf_real, [0 1], 'Mean CV2 ISI', tickangle)
    
    subplot(1,8,3)
    plot_conf(pn_AP, sum(pttype==1)/size(pttype,1), pn_AP_conf, [0 0.7], 'Fraction of AP Units', tickangle)
    
    subplot(1,8,4)
    plot_conf(pn_IP, sum(pttype==-1)/size(pttype,1), pn_IP_conf, [0 0.7], 'Fraction of IP Units', tickangle)
    
    subplot(1,8,5)
    plot_conf(pn_noosc, sum(pttype==0)/size(pttype,1), pn_noosc_conf, [0 0.7], 'Fraction of Non-Osc Units', tickangle)
    
    subplot(1,8,6)
    plot_conf(meanpow_ratio, meanpow_ratio_sim, meanrat_conf_real, [0 5], 'AP/IP Mean Power Ratio', tickangle)
    
    subplot(1,8,7)
    plot_conf(CV_AP, CV_AP_sim, cvs_ap_conf_real, [0 1], 'AP Power Pop. CV', tickangle)
    
    subplot(1,8,8)
    plot_conf(CV_IP, CV_IP_sim, cvs_ip_conf_real, [0 1], 'IP Power Pop. CV', tickangle)

    set(gcf,'Position',[500 400 600 300])
end

% get an nMC x sample matrix of bootstrapped samples
function [mat] = bootsample(input, nMC)
    len = length(input);
    mat = reshape(datasample(input,nMC*len,'replace',true),nMC,len);
end

% plotting each panel of confidence interval comparison figure
function [] = plot_conf(real, sim, conf, yrange, label, tickangle)
    hold on
    plot([1 1], [conf(1) conf(2)],'k','LineWidth',1.5)
    plot([.85 1.15], [conf(1) conf(1)], 'k','LineWidth',1.5)
    plot([.85 1.15], [conf(2) conf(2)], 'k','LineWidth',1.5)
    plot([0.8 1.2],ones(1,2)*real,'b','LineWidth',2.5)
    plot(1.1,sim,'r<','MarkerSize',10,'MarkerFaceColor','r')
    set(gca,'xtick',[1 2],'xticklabels',{label})
    xtickangle(gca,tickangle)
    xlim([0.7 1.3])
    ylim(yrange)
    makeNice(gca)
end