% ap_ip_examples_phase_offset.m
% John Parker, last edited Jan 2023
% script to generate example of AP/IP in comparison with GPe
% Plots phase offsets from GPe with Whalen2021_plot_fits.m
% search_freq = [0.5, 4] and force_freq = 2;
% Details in manuscript

% Changes below may lead to instability
gpe1 = gpe_sdf-mean(gpe_sdf); % Read in GPe data

% Arrays to be filled for figures
phaseShiftsAP = []; phaseShiftsIP = [];
phaseShiftsAPabs = []; phaseShiftsIPabs = [];
actPhaseAP = []; actPhaseIP = [];
redHisto = []; blueHisto = [];

% Iterate through SNr neurons and record necessary info for offsets
count = 0;
for i = 1:100
    if pttype(i) == -1 || pttype(i) == 1
        snr1 = sdfs(i,:)-mean(sdfs(i,:));
            
        [r,lags] = xcorr(gpe1,snr1,maxlag,'normalized');
        
        % find tp/tt by max 
        [ma, maind] = max(r);
        [mi, miind] = min(r);

        tp = sdftimes(abs(maind-(maxlag+1)));
        tt = sdftimes(abs(miind-(maxlag+1)));

        if pttype(i) == 1
            if maind-(maxlag+1) < 0
                phaseShiftsAPabs = [phaseShiftsAPabs; -tp];
            else
                phaseShiftsAPabs = [phaseShiftsAPabs; tp];
            end
        end

        if pttype(i) == -1
            if miind-(maxlag+1) < 0
                phaseShiftsIPabs = [phaseShiftsIPabs; -tt];
            else
                phaseShiftsIPabs = [phaseShiftsIPabs; tt];
            end
        end

        % DIFF TP - TT
        diffTpTt = abs(tp) - abs(tt);
        if pttype(i) == -1
            [max1,I1] = max(data.osc.psd{1}(i,:));
            actPhaseIP = [actPhaseIP; 2*pi*diffTpTt*data.osc.freqs(I1)];
            phaseShiftsIP = [phaseShiftsIP; diffTpTt];
            if tt < tp
                blueHisto = [blueHisto; tt];
            else
                redHisto = [redHisto; tp];
                count = count + 1;
            end
        elseif pttype(i) == 1
            [max1,I1] = max(data.osc.psd{1}(i,:));
            actPhaseAP = [actPhaseAP; 2*pi*(abs(tt)-abs(tp))*data.osc.freqs(I1)];
            phaseShiftsAP = [phaseShiftsAP; abs(tt)-abs(tp)];
            if tt < tp
                blueHisto = [blueHisto; tt];
                
            else
                redHisto = [redHisto; tp];
                
            end
        end
    end
end

% Specify AP/IP unit for ideal example
ap_unit = 73;
ip_unit = 69;

% Shift appropriately for best visualization
gpe_fig_sdf = ((gpe_sdf-0*mean(gpe_sdf)) / (max(gpe_sdf-0*mean(gpe_sdf)))) - 2;
ap_fig_sdf = ((sdfs(ap_unit,:)-0*mean(sdfs(ap_unit,:))) / (max(sdfs(ap_unit,:)-0*mean(sdfs(ap_unit,:))))) - 1;
ip_fig_sdf = (sdfs(ip_unit,:)-0*mean(sdfs(ip_unit,:))) / (max(sdfs(ip_unit,:)-0*mean(sdfs(ip_unit,:))));

ap_spikes = ts_min{ap_unit}; ip_spikes = ts_min{ip_unit};

% Generate example figure
figure; 
hold on; 
h1 = plot(sdftimes,gpe_fig_sdf,"Color",[0.5 0.5 0.5 0.5],"LineWidth",5); 
h2 = plot(sdftimes, ap_fig_sdf, "Color",[0.8 0 0 0.5],"LineWidth",5); 
h3 = plot(sdftimes, ip_fig_sdf, "Color", [0 0 0.8 0.5],"LineWidth",5); 
scatter(ip_spikes,ones(length(ip_spikes),1)*1.4,500,[0 0 0.8],'|',"MarkerFaceAlpha",0.5,"LineWidth",2);
scatter(ap_spikes,ones(length(ap_spikes),1),500,[0.8 0 0],'|',"MarkerFaceAlpha",0.5,"LineWidth",2);

xlim([12.75 15.75])
ylim([min(gpe_fig_sdf) 1.6])
xlabel('Time (s)')
set(gca,'YTick',[mean(gpe_fig_sdf) mean(ap_fig_sdf) mean(ip_fig_sdf) 1 1.4],'YTickLabels',{'GPe' 'AP Unit' 'IP Unit' 'AP Unit' 'IP Unit'})
set(gca,'XTick',[12:1:15])
set(gca,'TickLength',[0, 0.01])
set(gca,'LineWidth',80)
makeNice(gca)
hold off;

% Generate phase offset figure
figure
edges = -pi:pi/16:pi;
peakbar = histcounts(actPhaseAP,edges);
troughbar = histcounts(actPhaseIP,edges);
hold on
h = bar(edges(1:end-1),peakbar,'histc');
h.FaceColor = [.8 0 0];
h.LineStyle = 'none';
h.FaceAlpha = 0.5;
hold on
h = bar(edges(1:end-1),troughbar,'histc');
h.FaceColor = [0 0 0.8];
h.LineStyle = 'none';
h.FaceAlpha = 0.5;
set(gca,'XTick',[-pi,-pi/2,0,pi/2],'XTickLabels',{'-\pi' '-\pi/2' '0' '\pi/2'})
xlim([-pi pi/2])
xlabel('Offset from GPe (Phase)')
ylabel('Unit Count')
set(gca,'LineWidth',80)
makeNice(gca)
