function [recG, recS, strG, strS] = import_prp_file(file, n_snr)
% Opens a .prp file from a run of SNr.cc to get presynaptic GPe and SNr
% neuron indices and strengths for each postsynaptic SNr
% Returns:
% recG: {:,1} cell where ith element is (:,1) double of ONE-BASED indices 
%   of each presynaptic GPe neuron for the ith SNr neuron
% strG: same as recG, but the synaptic strengths of each synapse
% recS and strS are the same for presynaptic SNr neurons instead of GPe

fid = fopen(file);
line = fgetl(fid);
try
    if sum(line~='SNr neurons receiving GPe input from:')>0
        disp('Warning: unexepected GPe-SNr target header in .prp')
    end
catch
	disp('Warning: unexepected GPe-SNr target header in .prp')
end

recG = cell(n_snr,1);
for i = 1:n_snr
    recG{i} = str2num(fgetl(fid))'+1;
end

line = fgetl(fid);
try
    if sum(line~='SNr neurons receiving SNr input from:')>0
        disp('Warning: unexepected SNr-SNr target header in .prp')
    end
catch
	disp('Warning: unexepected SNr-SNr target header in .prp')
end

recS = cell(n_snr,1);
for i = 1:n_snr
    recS{i} = str2num(fgetl(fid))'+1;
end

line = fgetl(fid); % stnton title
line = fgetl(fid); % stnton values
line = fgetl(fid); % blank
line = fgetl(fid); % 
try
    if sum(line~='SNr-SNr Synapse Strengths:')>0
        disp('Warning: unexepected SNr-SNr strength header in .prp')
    end
catch
	disp('Warning: unexepected SNr-SNr strength header in .prp')
end

strS = cell(n_snr,1);
for i = 1:n_snr
    strS{i} = str2num(fgetl(fid))';
end

line = fgetl(fid); % blank
line = fgetl(fid); % 
try
    if sum(line~='SNr-GPe Synapse Strengths:')>0
        disp('Warning: unexepected SNr-GPe strength header in .prp')
    end
catch
	disp('Warning: unexepected SNr-GPe strength header in .prp')
end

strG = cell(n_snr,1);
for i = 1:n_snr
    strG{i} = str2num(fgetl(fid))';
end

fclose(fid);
end

