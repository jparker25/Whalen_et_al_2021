%%% plot powers around force_freq and GPe average freq
%%% run Whalen2021_plot_fits first with competitive model

gpe_mean = 24;
[~,freq_ind_gpe] = min(abs(data.osc.freqs-gpe_mean));
pows_gpe = data.osc.psd_unc{1}(:,freq_ind_gpe)./rates_min;

fig=figure;
hold on;
scatter(1:size(pows,1),pows,'b',"filled");
scatter(1:size(pows_gpe,1),pows_gpe,'r',"filled");
plot([1 size(pows,1)],[1 1].*mean(pows),":k","LineWidth",2);
plot([1 size(pows,1)],[1 1].*mean(pows_gpe),"--k","LineWidth",2);
xlabel("SNr Neuron");
ylabel("Power");
legend("2 Hz Power","24 Hz Power","Mean 2Hz Power","Mean 24Hz Power");

hold off
makeNice(gca)

saveas(fig,"snr_power_at_gpe_mean.eps","epsc");
saveas(fig,"snr_power_at_gpe_mean.fig");