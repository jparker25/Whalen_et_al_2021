% supplementary_figure2.m
% John Parker, last edited Jan 2023
% Comment out force_freq and search_freq initial variables in
% Whalen2021_plot_fits before running.
clear all;

% Changes below may lead to instability

% Analyze 2Hz simulation
twoHz_dir = "jitter_sims/no_jitter_freq_2/T50/random_run_0"; 
dir = sprintf('%s/competitive.*',twoHz_dir);
copyfile(dir,'data/')
force_freq = 2;
search_freqs = [0.5, 4];
Whalen2021_plot_fits;
Whalen2021_plot_balance;
close all;

% Store 2Hz simulation data
twoHz.cents = cents;
twoHz.non_osc_cent = non_osc_cent;
twoHz.xnoG = xnoG;
twoHz.ynoS = ynoS;
twoHz.pttype_osc = pttype_osc;
twoHz.strG_total = strG_total;
twoHz.strS_total = strS_total;

% Analyze 15Hz simulation
fifteenHz_dir = "jitter_sims/no_jitter_freq_15/T50/random_run_8";
dir = sprintf('%s/competitive.*',fifteenHz_dir);
copyfile(dir,'data/')
force_freq = 15;
search_freqs = [8, 22];
Whalen2021_plot_fits;
Whalen2021_plot_balance;
close all;

% Store 15Hz simulation data
fifteenHz.cents = cents;
fifteenHz.non_osc_cent = non_osc_cent;
fifteenHz.xnoG = xnoG;
fifteenHz.ynoS = ynoS;
fifteenHz.pttype_osc = pttype_osc;
fifteenHz.strG_total = strG_total;
fifteenHz.strS_total = strS_total;


% Plot figure
figure
subplot(1,2,1)
hold on
scatter(twoHz.cents(1,1), twoHz.cents(1,2), 2000, [1 .75 .75], '.') 
scatter(twoHz.cents(2,1), twoHz.cents(2,2), 2000, [.75 .75 1], '.')
scatter(twoHz.non_osc_cent(1), twoHz.non_osc_cent(2), 2000, 'k', '.')

scatter(twoHz.xnoG,twoHz.ynoS,200,[0.8 0.8 0.8],'.')
scatter(twoHz.strG_total(twoHz.pttype_osc==-1),twoHz.strS_total(twoHz.pttype_osc==-1),200,'b','.')
scatter(twoHz.strG_total(twoHz.pttype_osc==1),twoHz.strS_total(twoHz.pttype_osc==1),200,'r','.')

plot([0 800],[400 100],'k','LineWidth',3)
line([400 400],[0 250],'Color','k','LineWidth',3)
hold off

[hLg, icons]=legend({'AP','IP',"Non-Osc"})
icons = findobj(icons,'Type','patch');
icons = findobj(icons,'Marker','none','-xor');
set(icons(1:3),'MarkerSize',20);
xlabel('Total Power From GPe')
ylabel('Total Power From SNr')
xlim([0 800])
ylim([0 700])
makeNice(gca)
set(gcf,'Position',[10 400 450 400])

subplot(1,2,2)
hold on
scatter(fifteenHz.cents(1,1), fifteenHz.cents(1,2), 2000, [1 .75 .75], '.') 
scatter(fifteenHz.cents(2,1), fifteenHz.cents(2,2), 2000, [.75 .75 1], '.')
scatter(fifteenHz.non_osc_cent(1), fifteenHz.non_osc_cent(2), 2000, 'k', '.')
scatter(fifteenHz.xnoG,fifteenHz.ynoS,200,[0.8 0.8 0.8],'.')
plot([0 800],[400 100],'k','LineWidth',3)
line([400 400],[0 250],'Color','k','LineWidth',3)
for type = [1 -1]
    col = cols(1+(type<0),:);
    scatter(fifteenHz.strG_total(fifteenHz.pttype_osc==type),fifteenHz.strS_total(fifteenHz.pttype_osc==type),200,col,'.')
end

hold off
[hLg, icons]=legend({'AP','IP',"Non-Osc"})
icons = findobj(icons,'Type','patch');
icons = findobj(icons,'Marker','none','-xor');
set(icons(1:3),'MarkerSize',20);
xlabel('Total Power From GPe')
ylabel('Total Power From SNr')
xlim([0 800])
ylim([0 700])
makeNice(gca)
set(gcf,'Position',[10 400 450 400])


figure;
hold on
fifteen_all = []
for i = 1:length(fifteenHz.xnoG)
    fifteen_all = [fifteen_all; [fifteenHz.xnoG(i) fifteenHz.ynoS(i);]];
end

xx15 = fifteenHz.strG_total(fifteenHz.pttype_osc==-1);
yy15 = fifteenHz.strS_total(fifteenHz.pttype_osc==-1);
for i = 1:length(xx15)
    fifteen_all = [fifteen_all; [xx15(i) yy15(i);]];
end

xx15 = fifteenHz.strG_total(fifteenHz.pttype_osc==1);
yy15 = fifteenHz.strS_total(fifteenHz.pttype_osc==1);
for i = 1:length(xx15)
    fifteen_all = [fifteen_all; [xx15(i) yy15(i);]];
end

scatter(twoHz.cents(1,1), twoHz.cents(1,2), 2000, [1 .75 .75], '.') 
scatter(twoHz.cents(2,1), twoHz.cents(2,2), 2000, [.75 .75 1], '.')
scatter(mean(fifteen_all(:,1)),mean(fifteen_all(:,2)),2000,[0.9 0.9 0.9],'.')
scatter(fifteen_all(:,1),fifteen_all(:,2),200,[0.7 0.7 0.7],'.')

scatter(twoHz.strG_total(twoHz.pttype_osc==-1),twoHz.strS_total(twoHz.pttype_osc==-1),200,'b','.')
scatter(twoHz.strG_total(twoHz.pttype_osc==1),twoHz.strS_total(twoHz.pttype_osc==1),200,'r','.')


hold off

legend({'AP','IP',"15Hz"})
xlabel('Total Power From GPe')
ylabel('Total Power From SNr')
xlim([0 800])
ylim([0 700])
makeNice(gca)
set(gcf,'Position',[10 400 450 400])

