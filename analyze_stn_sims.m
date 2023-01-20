% analyze_stn_sims.m
% John Parker, last edited Jan 2023
% script to analyze several STN simulations with Whalen2021_plot_fits.m
% search_freq = [0.5, 4] and force_freq = 2;

% Changes below may lead to instability

% Default STN sim parameters from python script
animals = [2,3,4,5];
wstn = ["6"];
save_dir = "stn_sims";
samples = 4;
rngs = [0.13 0.23];
output_data = zeros(length(animals)*length(wstn)*samples,13); % data storage variable

% Iterate through all parameters and store outputs
count = 0;
for animal=1:length(animals)
    for sample=1:samples
        for ws=1:length(wstn)
            dir = sprintf("%s/animal_%g/rand_True/sample_%g/sim_gstn_rng_%g_%g_wstn_%s/competitive.*",save_dir,animals(animal),sample,rngs(1),rngs(2),wstn(ws));
            copyfile(dir,"data/")
            Whalen2021_plot_fits;
            close all;
            count = count + 1;
            output_data(count,:) = [animals(animal),sample,rngs(1),rngs(2),eval(wstn(ws)),zz];
        end
    end
end

% Save data to CSV for further analysis
T = array2table(output_data);
T.Properties.VariableNames(1:13) = {'animal','sample','gstn_low','gstn_high','wstn_static','mean_rate','mean_cv2','frac_ap','frac_ip','frac_non_osc','apip_power_ratio','ap_power_pop_cv','ip_power_pop_cv'};
writetable(T,sprintf('%s_results.csv',save_dir));

% Read in real data and record in CSV for comparisons
real_data = [mean(rates_real),mean(cv2isi_real),pn_AP,pn_IP,pn_noosc,meanpow_ratio,CV_AP,CV_IP];
confs = [mean_rate_conf_real;mean_cv2isi_conf_real;pn_AP_conf'; pn_IP_conf'; pn_noosc_conf';meanrat_conf_real;cvs_ap_conf_real;cvs_ip_conf_real];
conf_res = zeros(3,length(confs)/2);

for i=1:length(confs)/2
    conf_res(1,i) = confs(2*i-1);
    conf_res(2,i) = confs(2*i);
    conf_res(3,i) = real_data(i);
end

headers = {'mean_rate','mean_cv2','frac_ap','frac_ip','frac_non_osc','apip_power_ratio','ap_power_pop_cv','ip_power_pop_cv'};
rowNames={'conf_low' 'conf_high' 'mean'};
 
confTable = array2table(conf_res);
confTable.Properties.VariableNames(1:length(confs)/2) = headers;
confTable.("type")=rowNames';
writetable(confTable,'real_data_results.csv');