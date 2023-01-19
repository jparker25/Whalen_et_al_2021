%% hist plot of ap/ip in vivo with competitive 
%% run Whalen2021_plot_fits with competitive model first

in_vivo = load('Whalen2021_SNr_pow_by_ecog.mat');

bins = 20;
redge = 5e-5;
edges = 0:redge/(bins+1):redge;

figure;
subplot(2,1,1)
hold on
histogram(in_vivo.peakpows_r,'BinEdges',edges,'FaceColor','k')
histogram(pows_ap,'BinEdges',edges,'FaceColor','r')
hold off
legend({'in vivo','syn. comp.'})
xlabel('AP Delta Power')
ylabel('Count')
makeNice(gca)

subplot(2,1,2)
hold on
histogram(in_vivo.troughpows_r,'BinEdges',edges,'FaceColor','k')
histogram(pows_ip,'BinEdges',edges,'FaceColor','b')
hold off
legend({'in vivo','syn. comp.'})
xlabel('IP Delta Power')
ylabel('Count')

%ylabel('Power')
%set(gca,'xtick',[0 1],'xticklabels',{'SNr AP';'SNr IP'})
makeNice(gca)
%set(gcf,'Position',[10 400 430 300])