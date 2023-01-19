% Tim C. Whalen, last edited May 2021
% Plots Whalen 2021 Figure 2E-H, pairwise phase histograms

%% TO RUN
% For simulated data, set dataset to 0 and run Whalen2021_plot_fits using
% your desired output from SNr.cc, then run. The default parameters for the
% "competitive model" should replicate the bottom panels of Fig 2F-H

% For real data, two different datasets are used depending on whether we
% need AP/IP labeled data (which is only possible in a subset of data which
% contains a M1 ECoG reference)
% To replicate the top panels of Fig 2G-H (not requiring ECoG), set dataset
% to 1 and ensure Whalen2021_SNr_depleted.mat is in your path, then run.
% To replicate Fig 2E (requiring ECoG), set dataset to 2 and ensure
% Whalen2021_SNr_depleted_ecog.mat is in your path, then run (this file may
% not be in the repository due to its size). This will also produce figures
% similar to G and H but needlessly using only a subset of the data, so 
% they are not in the paper.

dataset = 1; % 0 for simulation, 1 for all in vivo SNr data, 2 for SNr with ECoG

%% Main code - edits below may lead to instability

%% Load data into consistent format
if dataset==0 % Simulation data
    pttypes = {pttype};
elseif dataset==1 % all SNr
    force_freq = 0;
    load('Whalen2021_SNr_depleted','data');
    nfiles = data.nfiles;
    pttypes = cell(nfiles,1);
    for f = 1:data.nfiles % ignore pttypes when no ecog
        rates = data.rates{f};
        pttypes{f} = zeros(size(data.rates{f}));
    end
elseif dataset==2 % ecog SNr
    force_freq = 0;
    load('Whalen2021_SNr_depleted_ecog','data');
    pttypes = data.ecog_reg.pttype;
else
    error('Invalid value for "dataset"; must be 0, 1 or 2')
end

%% Main analysis
dt = 0.001; % sec
sdf_std = .05; % sec

nfiles = data.nfiles;
sdfs = cell(nfiles,1);
xcs = cell(nfiles,1);
pairlags = cell(nfiles,1);
ptmat = cell(nfiles,1);

apap_delays = cell(nfiles,1);
apip_delays = cell(nfiles,1);
ipip_delays = cell(nfiles,1);

for f = 1:nfiles
    ts = data.ts{f}(data.rates{f}>data.osc.min_rate);
    ts = ts(data.osc.has_osc{f}==1);
    pttype_defined = pttypes{f}(data.osc.has_osc{f}==1);
    T = data.T(f);
    lo = 0; hi = T;
    nu = size(ts,1);
    sdftimes = (lo+sdf_std*3):dt:(hi-sdf_std*3); % from gaussSpikeDensity
    sdfs{f} = zeros(nu,length(sdftimes));
    xcs{f} = cell(nu,nu);
    
    for u = 1:nu
        [~, sdfs{f}(u,:)] = gaussSpikeDensity( ts{u}, sdf_std, 0, T);
    end
    pairlags{f} = zeros(nu,nu)+NaN;
    ptmat{f} = zeros(nu,nu)+NaN;
    
    nap = sum(pttype_defined==1);
    nip = sum(pttype_defined==-1);
    apap_delays{f} = zeros((nap^2-nap)/2,1)+NaN;
    apip_delays{f} = zeros(nap*nip,1)+NaN;
    ipip_delays{f} = zeros((nip^2-nip)/2,1)+NaN;
    
    % counts for filling abvoe arrays
    napap = 0; 
    napip = 0;
    nipip = 0;
    
    for i = 1:nu
        for j = i+1:nu
            sig1 = sdfs{f}(i,:)';
            sig2 = sdfs{f}(j,:)';
            
            FS = 1/dt; slen = 16000; blen = 20000; step = 8000; cut = 4000; sdf_to_plot = 0;
            maxlag = (blen-slen)/2;
            lags = -maxlag:maxlag;
            xco_avg = windowedXCorr(sig1, sig2, slen, blen, 'step_size', step, 'to_plot', sdf_to_plot, 'FS', FS);
            
            % find whether peak or trough closer to zero
            if force_freq % for simulations with strong oscillations which do not attenuate correlation at long lags
                newlag = 600/force_freq; % .6 of period
                
                [maxvals, maxinds] = findpeaks(xco_avg(maxlag+1-newlag:maxlag+newlag));
                [~, whichmax] = min(abs(maxinds-newlag));
                maxval = maxvals(whichmax); 
                maxind = maxinds(whichmax)+maxlag-newlag;
                
                [minvals, mininds] = findpeaks(-xco_avg(maxlag+1-newlag:maxlag+newlag));
                [~, whichmin] = min(abs(mininds-newlag));
                minval = -minvals(whichmin); 
                minind = mininds(whichmin)+maxlag-newlag;
            else
                [maxval, maxind] = max(xco_avg);
                [minval, minind] = min(xco_avg);
            end
            domax = abs(maxind-maxlag)<abs(minind-maxlag); % 1 if max ccloser to zero, 0 othw.
            ptmat{f}(i,j) = domax+(~domax*-1); % becomes 1 for peak, -1 for trough
            xcs{f}{i,j} = xco_avg'; % not used but preserved for visualization
            
            if domax
                exind = maxind;
            else
                exind = minind;
            end
            pairlags{f}(i,j) = lags(exind);
            
            if pttype_defined(i)==1
                if pttype_defined(j)==1
                    apap_delays{f}(napap+1) = lags(maxind);
                    napap = napap+1;
                elseif pttype_defined(j)==-1
                    apip_delays{f}(napip+1) = -lags(minind); % flip so ap/ip always wrt ap
                    napip = napip+1;
                end
            elseif pttype_defined(i)==-1
                if pttype_defined(j)==1
                    apip_delays{f}(napip+1) = lags(minind);
                    napip = napip+1;
                elseif pttype_defined(j)==-1
                    ipip_delays{f}(nipip+1) = lags(maxind);
                    nipip = nipip+1;
                end
            end
        end
    end
end

% All lags
pairlags_all = cell2mat(cellfun(@singlecolon,pairlags,'UniformOutput',0));
pairlags_all = (pairlags_all(~isnan(pairlags_all)));
edges = (0:15:300)-0.5;

% Within and between population lags
ptmat_all = cell2mat(cellfun(@singlecolon,ptmat,'UniformOutput',0));
ptmat_all = ptmat_all(~isnan(ptmat_all));

%% Plot pairwise phase histograms
figure
subplot(2,1,1)
hold on
histogram(abs(pairlags_all(ptmat_all==1)),edges)
title('Pair lags (same pop)')
ylabel('Pair Count')
xlim([0 300])
makeNice(gca)
subplot(2,1,2)
hold on
histogram(abs(pairlags_all(ptmat_all==-1)),edges)
title('Pair lags (opp pop)')
xlabel('Lag (ms)')
ylabel('Pair Count')
xlim([0 300])
makeNice(gca)

%% Plot AP-IP directed lags (when possible)
apip_lags = cell2mat(apip_delays);
if ~isempty(apip_lags)
    figure
    histogram(apip_lags,[-125:25:125],'FaceColor',[.5 .5 .5])
    hold on
    plot([0 0], [0 12],'--','Color',[0 0 0],'LineWidth',2)
    xlabel('<- AP Leads                   Delay (ms)                   IP Leads ->')
    ylabel('Count')
    title('AP vs IP Delays')
    p = signrank(cell2mat(apip_delays));
end

%% helpers

function [Acol] = singlecolon(A)
% Functional version of (:)
Acol=A(:);
end