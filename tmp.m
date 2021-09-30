clear; close all; clc;

load('pulse.mat');

%>> whos
%>> Name       Size               Bytes  Class 
%>> fs           1x1                 8  double              
%>> pulse      351x1              2808  double    

vessel_angle = deg2rad(45); % angle between beam and vessel [rad]
vessel_diameter = 10e-3;    % diameter of vessel [m]
start_depth = 20e-3;        % vessel start depth [m]
dx = 2e-4;                  % axial resolution [m]
lines = 64;                 % number of RF lines 

f_prf = 5e3;                %Hz
vz = 0.15;                  %m/s
T_prf = 1 / f_prf;          %s
c = 1500;                   %m/s
f0 = 3.2e6;                 % Hz                  
err_std = 0.5;

ts = (2*vz*cos(vessel_angle) / c) * T_prf;
n_ts = round(ts * fs);
n_signals = 500;

vessel_length = vessel_diameter / sin(vessel_angle); %m
delta_depth = dx / tan(pi - vessel_angle);
delta_depth = floor(((2*delta_depth) / 1500)*fs);

depth2time = (2*vessel_length) / 1500; % s
N_scatter = floor(depth2time * fs);

depth = floor(( 2 * start_depth / 1500 ) * fs);
total_depth = depth + N_scatter + floor(depth / 2);

scatter = randn(N_scatter,1);
%scatter(1) = 50; 
%scatter(end) = 50;
phantom = vertcat(zeros(depth, 1), scatter, ...
            zeros(total_depth - depth - N_scatter, 1));
received_signal = conv(pulse, scatter);
N_received = length(received_signal);
received_signal = received_signal;% + err_std * randn(size(received_signal));

single_line = zeros(n_signals,N_received);
single_line(1, :) = received_signal;
cnt = 1;
for i=(2:n_signals)
    
    index = cnt * n_ts;
    
    signal = received_signal(1:end-index)';
    
    pad = N_received - length(signal);
    
    single_line(i, :) = [zeros(1, pad), signal];
    
    cnt = cnt + 1;
end

x = (0:length(single_line)-1) ./ f_prf;
y = (0:size(single_line,1)-1);
[X,Y] = meshgrid(x,y);
figure;
imagesc(x, y, single_line);
set(gca,'view',[-90, -90]);
colormap(hot);

%%

factor = - (c * f_prf) / (4 * pi * f0);

complex_sig = hilbert(single_line);
[r, lags] = xcorr(complex_sig(50,:));
idx = find(lags==1);

v_ = factor * atan(imag(r(idx)) / real(r(idx)));

%%
received_matrix = [];
scatter(1) = 50; 
scatter(end) = 50;
for j=(1:lines)
    phantom = vertcat(zeros(depth, 1), scatter, ...
        zeros(total_depth - depth - N_scatter, 1));
    received_signal = conv(pulse, phantom);
    received_signal = received_signal + err_std * randn(size(received_signal));
    received_matrix = [received_matrix, received_signal];

    depth = depth + delta_depth;
end

received_matrix_h = abs(hilbert(received_matrix));
received_matrix_drc = drc(unit_normalize(received_matrix_h), 60);

[rf_samples, n_lines] = size(received_matrix_drc);

width_x = (0:n_lines) .* dx;
width_x = (0:dx:n_lines);
width_x = width_x - width_x(end)/2;

t = (0:rf_samples-1) ./ fs;
depth = (t .* 1500) ./ 2;
figure;
imagesc(width_x, depth*1e2, int32(received_matrix_drc * 127));
colormap(gray(128));
ylabel('Depth [cm]');
xlabel('Axial displacement [m]')
%set(gca, 'view', [0, -90]);

%%



%%

f_prf = 5e3; %Hz
vz = 0.15; %m/s
T_prf = 1 / f_prf; %s
c = 1500; %m/s

ts = (2*vz*cos(vessel_angle) / c) * T_prf;
n_ts = round(ts * fs);

n_signals = 100;

received_over_time = zeros(n_signals,rf_samples);

received_over_time(1, :) = received_signal;
cnt = 1;

for i=2:n_signals
    index = cnt * n_ts;
    
    signal = received_signal(1:end-index)';
    
    pad = rf_samples - length(signal);
    
    received_over_time(i, :) = [zeros(1, pad), signal];
    
    cnt = cnt + 1;
end

%%

function v = acorr_estimate(rf_0, rf_1, c, f_prf, f0)
    rf_0_c = hilbert(rf_0);
    rf_1_c = hilbert(rf_1);
    
    

end

function y = unit_normalize(x)
    y = (x - min(min(x))) ./ (max(max(x)) - min(min(x)));
end

function y = drc(x, dynamic_range)
%DRC Apply dynamic range compression to a signal
%   The input should be unit-normalized

    x_max = 10^(-0/20);
    x_min = 10^(-dynamic_range/20);

    x(x < x_min) = x_min;
    x(x > x_max) = x_max;

    C1 = 1.0 / log10( x_max / x_min);
    C3 = 1.0 / x_min;

    y = C1 .* log10(C3 .* x);
end

