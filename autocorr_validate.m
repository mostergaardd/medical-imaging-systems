clear; close all; clc;

% Set the random seed
rng(0);

% Pulse parameters
f0 = 3.2e6;         % Probe center frequency [hz]
M = 6;              % Number of cycles in pulse [n]
fs = 1e8;           % System sample rate [hz]
n_fft = 1024;       % Number of samples to use in fft

[pulse, t, pulse_F, f, f0_est] = generate_pulse(f0, M, fs, n_fft);
f0_est = abs(f0_est);

% Plotting
figure;
subplot(211);
plot((0:length(pulse)-1)./fs, pulse);
title(['Pulse for CFM. ', num2str(M), ' cycles.']);
ylabel('Amplitude [v]');
xlabel('Time [s]');
axis tight; grid on;

subplot(212);
plot(f,abs(pulse_F));
title('Amplitude spectrum of pulse');
xlabel('Frequency [hz]');
ylabel('Amplitude');
axis tight; grid on;

%% Initialize some simulated data

vessel_angle = deg2rad(45); % angle between beam and vessel [rad]
vessel_diameter = 10e-3;    % diameter of simulated vessel [m]
f_prf = 5e3;                % pulse repetition frequency [hz]
T_prf = 1 / f_prf;          % pulse repetition time [s]
vz = 0.15;                  % simulated velocity [m/s]
c = 1500;                   % speed of sound [m/s]
err_std = 0.5;              % standard dev of noise
n_emissions = 50;           % number of emissions to generate [n]

single_line = simulate_single_line(vessel_angle, vessel_diameter, ...
    f_prf, fs, vz, c, n_emissions, pulse);

% Plotting
x = (0:length(single_line)-1) ./ fs;
x = (x.*c)./2;
y = (0:size(single_line,1)-1);
[X,Y] = meshgrid(x,y);
figure;
imagesc(x, y, single_line);
ylabel('Number of emissions (i)');
xlabel('Depth in vessel [m]');
title('Simulated emissions')
set(gca,'view',[-90, -90]);
colormap(hot);
ax = get(gca);
ax.XAxis.Exponent = -3;

%% Validate on the simulated data

winsize = abs(round(1/f0_est * fs));
[v, depth] = autocorr_estimator(single_line, winsize, c, f_prf, f0_est, ...
    vessel_angle, fs);

% Plotting
figure;
plot(depth(200:end-200), v(200:end-200));
xlabel('Depth [m]');
ylabel('Velocity [m/s]');
title({['Estimated velocities. Target v=', num2str(vz), ' m/s'], ...
    ['Winsize: ', num2str(winsize), ...
    ', f0: ', num2str(f0_est, '%.3G'), 'Hz']});
axis tight; grid on;

