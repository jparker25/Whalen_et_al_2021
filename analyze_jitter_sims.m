% analyze_jitter_sims.m
% Comment out force_freq and search_freq initial variables in
% Whalen2021_plot_fits before running.

data_direc = 'jitter_sims'; % directory where all jitter results are stored

tj1 = [100 500];
tj2 = [1000 2000];
ffreq = [2 15 25];
sfreqs = [[0.5 4]; [8 30]; [12 35]];
T = 50;
runs = [0 6 7 8 9 10];

output_data = zeros(length(runs)*size(sfreqs,1)*length(tj1)+length(runs)*length(ffreq),20);

count = 0;
for i = 1:length(tj1)
    for f = 1:length(ffreq)
        for r=1:length(runs)
            dir = sprintf('%s/t1_%g_t2_%g_freq_%g/T%g/random_run_%g/competitive.*',data_direc,tj1(i),tj2(i),ffreq(f),T,runs(r))
            copyfile(dir,'data/')
            force_freq = ffreq(f);
            search_freqs = sfreqs(f,:);
            Whalen2021_plot_fits;
            close all;
            count = count + 1;
            seed = runs(r);
            if r == 1
                seed = 749141;
            end
            output_data(count,:) = [ffreq(f),1,tj1(i),tj2(i),T,seed,1e-5,mean_ap,mean_ip,mean_rate_sim,mean_cv2isi_sim, sum(pttype==1)/size(pttype,1), sum(pttype==-1)/size(pttype,1) ,sum(pttype==0)/size(pttype,1) ,meanpow_ratio_sim,CV_AP_sim,CV_IP_sim,search_freqs(1),search_freqs(2),force_freq];
        end
    end
end

for f = 1:length(ffreq)
    for r = 1:length(runs)
        dir = sprintf('%s/no_jitter_freq_%g/T%g/random_run_%g/competitive.*',data_direc,ffreq(f),T,runs(r))
        copyfile(dir,'data/')
        force_freq = ffreq(f);
        search_freqs = sfreqs(f,:);
        Whalen2021_plot_fits;
        close all;
        count = count + 1;
        seed = runs(r);
        if r == 1
            seed = 749141;
        end
        output_data(count,:) = [ffreq(f),0,-1,-1,T,seed,1e-5,mean_ap,mean_ip,mean_rate_sim,mean_cv2isi_sim, sum(pttype==1)/size(pttype,1), sum(pttype==-1)/size(pttype,1) ,sum(pttype==0)/size(pttype,1) ,meanpow_ratio_sim,CV_AP_sim,CV_IP_sim,search_freqs(1),search_freqs(2),force_freq];
    end
end

output_table = array2table(output_data);
output_table.Properties.VariableNames(1:20) = {'Frequency','Jitter?','T1','T2','total Time','Random Seed run','APIP Scale','AP','IP','MFR','MCV2','FAP','FIP','NonOsc','Ratio','APCV','IPCV','Low Freq','Hi Freq','Force Freq'};
writetable(output_table,sprintf('%s/jitter_results.csv',data_direc));

