function [data, N_scatter] = simulate_single_line(...
    vessel_angle, vessel_diameter, f_prf, fs, vz, c, ...
    n_emissions, pulse, seed, err_std, add_noise, add_st_sig)
%SIMULATE_SINGLE_LINE Summary of this function goes here
%   Detailed explanation goes here
rng(seed);
T_prf = 1 / f_prf;

% Amplitude of vessel boundary
if add_st_sig
    vessel_amp = 100;
else
    vessel_amp = 0;
end

% Estimate time shift
ts = ((2*vz*cos(vessel_angle)) / c) * T_prf;
n_ts = round(ts * fs);

% Calculate number of samples in vessel
vessel_length = vessel_diameter / sin(vessel_angle); %m
depth2time = (2*vessel_length) / c; % s
N_scatter = round(depth2time * fs);

% Initialize some random scatters that we will time shift
scatter = randn(N_scatter,1);
stationary_part = repmat(randn(1, N_scatter), n_emissions, 1);

% Maybe just use zeros outside vessel
if ~add_st_sig
    stationary_part = stationary_part .* 0;
end

% Time-shift the scatters in the vessel
data = zeros(n_emissions, 3*N_scatter);
data(1,N_scatter+1:2*N_scatter) = scatter(:);
for i = 2:n_emissions
   data(i,N_scatter+1:2*N_scatter) =  circshift(scatter(:), (i-1) * n_ts);
end

% Surround the vessel with stationary signal and add interface scatter
data(:, [1:N_scatter, 2*N_scatter+1:end]) = [stationary_part, stationary_part];
data(:, [N_scatter, 2*N_scatter]) = ones(n_emissions, 2) .* vessel_amp;

% Convolute phantom with US pulse
for i=1:n_emissions
    data(i,:) = conv(data(i,:), pulse, 'same');
end

% Maybe add random noise to the measurement
if add_noise
    data = data + err_std * randn(size(data));
end

end

