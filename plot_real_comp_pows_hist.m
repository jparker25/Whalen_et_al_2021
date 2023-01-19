%%% plot powers of competitive model and compare with real data as
%%% histogram
%%% run Whalen2021_plot_fits first with competitive model

vivo = load('Whalen2021_SNr_depleted.mat');

rates_min_vivo = vivo.data.rates{1}(vivo.data.rates{1}>vivo.data.osc.min_rate);
pows_vivo = vivo.data.osc.psd_unc{1,:}(:,freq_ind)./rates_min_vivo;

figure;
hold on;
h=histogram(pows,"BinWidth",0.25e-5);
histogram(pows_vivo,"BinWidth",0.25e-5);
hold off;
legend("Synaptic Competition","In Vivo");
xlabel("Delta Power");
ylabel("Count");
makeNice(gca);

