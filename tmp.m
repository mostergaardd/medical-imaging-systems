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
n_signals = 100;

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

winsize = round(1/f0 * fs);
v_est = acorr_estimate(single_line, 1, c, f_prf, f0, vessel_angle);

mean(v_est)
std(v_est)

figure;
plot(v_est(200:end-200));

%%

function v = acorr_estimate(rf, winsize, c, f_prf, f0, rad_angle)
    rf_iq = hilbert(rf);
    N_rf = size(rf_iq, 2);
    
    factor = - (c * f_prf) / (4 * pi * f0);
    
    v = zeros(1, N_rf - winsize);
    
    for n = 1:N_rf - winsize
        r = rf_iq(:, n:n+winsize-1);
        
        vs = 0;
        for m = 1:winsize
            rm = r(:, m);
            
            y = imag(rm);
            x = real(rm);
            
            Ri = 0;
            Rr = 0;
            for i = (1:length(rm)-1)
                Ri = Ri + (y(i+1)*x(i) - x(i+1)*y(i));
                Rr = Rr + (x(i+1)*x(i) + y(i+1)*y(i));
            end
            Ri = Ri / length(rm);
            Rr = Rr / length(rm);
            
            vs_m = factor * atan2(Ri, Rr);
            vs_m = vs_m / cos(rad_angle);
            
            vs = vs + vs_m; 
        end
        
        v(n) = vs / winsize; 
    end
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

