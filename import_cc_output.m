function [ts] = import_cc_output(file, sub, N)
    % Import simulated SNr population spikes and returns cell of spike times
    % Inputs:
    % file: a .txt (e.g. .hst from SNr.cc output) with two \t separated columns
    %   where 1st column is time (in sec), 2nd is spiked neuron index (0-based).
    % sub: seconds to discard from beginning
    % N: number of neurons in simulation
    % Output:
    % {N,1} cell of (:,1) spike times for each neuron
    
    k = importdata(file);
    time = k(:,1);
    nind = round(k(:,2));
    
    nind = nind(time>=sub);
    time = time(time>=sub)-sub;
    
    ts = cell(N,1);
    for i = 0:max(nind)
        ts{i+1} = time(nind==i);
    end
end

