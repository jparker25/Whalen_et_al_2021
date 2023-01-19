% Tim C. Whalen, last edited May 2021
% Plots Whalen 2021 Figure 4 for one simulation run
% Note that Fig 4C pools 5 simulation runs; this code only strS_total one.

%% TO RUN
% First, run Whalen2021_plot_fits with the desired output from SNr.cc
% Then, set scale_by_power = 1 for power*strength (Fig 4B) or 0 for
% strength only (Fig A)

scale_by_power = 1; % to multiply strengths by osc power

% JEP
force_freq = 2;
[~, freq_ind] = min(abs(data.osc.freqs-force_freq)); % index closest to forcing frequency
%

%% Main code - edits below may lead to instability

%% First fig: find each neuron position in strength plane
[recG, recS, strG, strS] = import_prp_file([infile '.prp'], npop*2);

if scale_by_power
    data_gpe = struct();
    data_gpe.ts{1}=ts_gpe; data_gpe.T=T; data_gpe.rates{1}=cellfun(@(x) size(x,1)/T, ts_gpe);
    data_gpe.nfiles=1; data_gpe.animalcodes=1;
    data_gpe.files_trunc = ['GPe_' infile(1:end)];
    [ data_gpe ] = renewalPSD_phaseShift_batch( data_gpe);
    
    for i = 1:npop*2
        strG{i} = strG{i}.*data_gpe.osc.psd{1}(recG{i},freq_ind);
        strS{i} = strS{i}.*data.osc.psd{1}(recS{i},freq_ind);
    end
end

strG_total = cellfun(@sum,strG);
strS_total = cellfun(@sum,strS);

xG = strG_total(pttype_osc~=0);
yS = strS_total(pttype_osc~=0);
lab = pttype_osc(pttype_osc~=0);

points_all = [xG yS];

% find centroids and line between them
cents = [mean(points_all(lab==1,:)); mean(points_all(lab==-1,:))];
% JEP
xnoG =  strG_total(pttype_osc==0);
ynoS = strS_total(pttype_osc==0);
nopoints_all = [xnoG ynoS];
non_osc_cent = mean(nopoints_all);

%
cent_slope = (cents(1,2)-cents(2,2))/(cents(1,1)-cents(2,1));
cent_int = cents(1,2)-cent_slope*cents(1,1);

% plotting plane
cols = [1 .2 .2; .2 .2 1];
figure
hold on
plot([0,-cent_int/cent_slope],[cent_int,0],'--','Color',[.9 .9 .9],'LineWidth',5)
scatter(cents(1,1), cents(1,2), 2000, [1 .75 .75], '.') 
scatter(cents(2,1), cents(2,2), 2000, [.75 .75 1], '.')
scatter(non_osc_cent(1), non_osc_cent(2), 2000, 'k', '.')
scatter(xnoG,ynoS,200,[0.8 0.8 0.8],'.')

for type = [1 -1]
    col = cols(1+(type<0),:);
    scatter(strG_total(pttype_osc==type),strS_total(pttype_osc==type),200,col,'.')
end

legend({'Centroid Axis','AP','IP',"Non-Osc"})
xlabel('Total Power From GPe')
ylabel('Total Power From SNr')
makeNice(gca)
%xlim([0 200])
set(gcf,'Position',[10 400 450 400])

%% Second fig: distance to centroids along axis defined by centroids

% rotate space s.t. line connecting centroids is horizontal
slope_cents = (cents(1,2)-cents(2,2))/(cents(2,1)-cents(1,1));
theta = atan(-slope_cents);
rotmat = [cos(theta) -sin(theta); sin(theta) cos(theta)];

points_rot = points_all*rotmat;
cents_rot = cents*rotmat;

dists2ap = -(points_rot(:,1)-cents_rot(1,1));
dists2ip = (points_rot(:,1)-cents_rot(2,1));

scale = mean([mean(abs(dists2ap(lab==1))) mean(abs(dists2ip(lab==-1)))]);
dists2ap_scaled = dists2ap/scale;
dists2ip_scaled = dists2ip/scale;

% plotting histogram
if scale_by_power
    title_append = ' (Weight * Power)';
else
    title_append = ' (Weight)';
end
edges = -2:.5:6;
figure
subplot(2,1,1)
histogram(dists2ap_scaled(lab==-1),edges,'FaceColor',[.5 .5 .5])
ylabel('Neuron Count')
title(['IP Scaled Displacement from AP Cluster', title_append])
makeNice(gca)
subplot(2,1,2)
histogram(dists2ip_scaled(lab==1),edges,'FaceColor',[.5 .5 .5])
title(['AP Scaled Displacement from IP Cluster', title_append])
xlabel('Scaled Displacement')
makeNice(gca)