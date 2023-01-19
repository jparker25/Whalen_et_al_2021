% Tim C. Whalen, last edited May 2021
% Plots Whalen 2021 Figure 3C, fitness of QIF model with varying noise
% levels

%% TO RUN
% First, run the 7 QIF models in SNr.cc (noiseless and the 6 levels of
% non-zero noise), and make sure the outputs are in your MATLAB path.
%
% This script runs the Whalen2021_plot_fits script several times, storing
% its results each iteration
% To run this script, do not set any input variables (that is, comment out
% each set of inputs at the top of the file).


%% Main code - edits below may lead to instability

files = {'competitive_qif';
    'competitive_qif_noise015';
    'competitive_qif_noise030';
    'competitive_qif_noise045';
    'competitive_qif_noise060';
    'competitive_qif_noise075';
    'competitive_qif_noise090'};
control = 0;
scatter_colors = 0;
nfiles = size(files,1);
noises = 0:.015:.09;

cv = zeros(nfiles,1);
apfrac = zeros(nfiles,1);
ipfrac = zeros(nfiles,1);
nofrac = zeros(nfiles,1);
powrat = zeros(nfiles,1);
for f = 1:nfiles
    infile = files{f};
    Whalen2021_plot_fits
    cv(f) = mean_cv2isi_sim;
    apfrac(f) = sum(pttype==1)/size(pttype,1);
    ipfrac(f) = sum(pttype==-1)/size(pttype,1);
    nofrac(f) = sum(pttype==0)/size(pttype,1);
    powrat(f) = meanpow_ratio_sim;
end

ys = [cv apfrac ipfrac nofrac powrat];
confs = [mean_cv2isi_conf_real, pn_AP_conf', pn_IP_conf', pn_noosc_conf', meanrat_conf_real];
confmeans = [mean(cv2isi_real), pn_AP, pn_IP, pn_noosc, meanpow_ratio];
figure
for i = 1:size(ys,2)
    subplot(5,1,i)
    hold on
    rectangle('Position',[noises(1) confs(1,i) noises(end)-noises(1) confs(2,i)-confs(1,i)],'FaceColor',[.95 .95 .95],'LineStyle','none');
    plot([noises(1), noises(end)],ones(1,2)*confmeans(i),'Color',[.7 .7 1],'LineWidth',2)
    plot(noises, ys(:,i),'r','LineWidth',2)
    xlim([noises(1),noises(end)])
end