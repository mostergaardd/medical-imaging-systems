function [single_line] = simulate_single_line(...
    vessel_angle, vessel_diameter, f_prf, fs, vz, c, n_emissions, pulse)
%SIMULATE_SINGLE_LINE Summary of this function goes here
%   Detailed explanation goes here
T_prf = 1 / f_prf;

% Estimate time shift
ts = (2*vz*cos(vessel_angle) / c) * T_prf;
n_ts = ceil(ts * fs);

% Calculate number of samples in vessel
vessel_length = vessel_diameter / sin(vessel_angle); %m
depth2time = (2*vessel_length) / c; % s
N_scatter = floor(depth2time * fs);

% Initialize some random scatters that we will time shift
scatter = randn(N_scatter,1);
received_signal = conv(pulse, scatter);

% housekeeping
N_received = length(received_signal);
single_line = zeros(n_emissions,N_received);
single_line(1, :) = received_signal;
cnt = 1;
for i=(2:n_emissions)     % Shift the received signal
    index = cnt * n_ts;
    
    signal = received_signal(1:end-index)';
    pad = N_received - length(signal);
    single_line(i, :) = [zeros(1, pad), signal];
    
    cnt = cnt + 1;
end

end

