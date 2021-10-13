clear; close all; clc;
%%
all_biases = {};
all_stdevs = {};
%%

% System parameters
f0 = 3.2e6;                 % Probe center frequency [hz]
M = 4;                      % Number of cycles in pulse [n]
fs = 1e8;                   % System sample rate [hz]
n_fft = 1024;               % Number of samples to use in fft
vessel_angle = deg2rad(45); % angle between beam and vessel [rad]
vessel_diameter = 10e-3;    % diameter of simulated vessel [m]
f_prf = 5e3;                % pulse repetition frequency [hz]
T_prf = 1 / f_prf;          % pulse repetition time [s]
vz = 0.15;                  % simulated velocity [m/s]
c = 1500;                   % speed of sound [m/s]
err_std = 0.5;              % standard dev of noise
seed = 0;                   % random seed

biases = [];
stdevs = [];

% Hyperparams
apply_mf                = false;
apply_ec                = false;
add_noise               = false;
add_stationary_signal   = false;
modify_pulse            = true;
doplot                  = true;
overwrite_winsize       = false;
winsize_                = 1;

for n_emissions=4:20
    %%
    load('pulse.mat');

    if modify_pulse
        figure; plot(pulse); hold on; 
        mod_sig = sin(2*pi*f0.*(0:1/fs:M/f0)) .* 0.1*max(pulse);
        %mod_sig = mod_sig(:) .* hanning(length(mod_sig));
        %mod_sig = (mod_sig - min(mod_sig)) ./ (max(mod_sig) - min(mod_sig));
        pulse = conv(pulse, mod_sig);
        plot(pulse);
    end

    t = (0:length(pulse))./fs;   % time
    f = (-n_fft/2:n_fft/2-1)./n_fft.*fs;
    [f0_est, pulse_F] = estimate_f0(pulse, n_fft, f, doplot);
    f0_est = abs(f0_est);

    matched_h = matched_filter(pulse, doplot);

    if doplot
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
    end

    % Initialize some simulated data

    [single_line, vessel_depth] = simulate_single_line(... 
        vessel_angle, vessel_diameter, f_prf, fs, vz, c,...
        n_emissions, pulse, seed, err_std, add_noise, add_stationary_signal);

    if doplot
        % Plotting
        x = (0:length(single_line)-1) ./ fs;
        x = (x.*c)./2;
        y = (0:size(single_line,1)-1);
        [X,Y] = meshgrid(x,y);
        figure();
        imagesc(x, y, single_line);
        ylabel('Number of emissions (i)');
        xlabel('Depth in vessel [m]');
        title('Simulated emissions')
        set(gca,'view',[-90, -90]);
        colormap(hot);
        ax = get(gca);
        ax.XAxis.Exponent = -3;
    end

    % Validate on the simulated data

    winsize = abs(round(1/f0_est * fs));
    if overwrite_winsize
        winsize = winsize_;
    end

    [v, depth] = autocorr_estimator(single_line, winsize, c, f_prf, f0_est, ...
        vessel_angle, fs, matched_h, apply_ec, apply_mf);

    if doplot
        % Plotting
        figure();
        plot(depth(), v());
        plot(depth(vessel_depth:end-vessel_depth), v(vessel_depth:end-vessel_depth));
        xlabel('Depth [m]');
        ylabel('Velocity [m/s]');
        title({['Estimated velocities. Target v=', num2str(vz), ' m/s'], ...
            ['Winsize: ', num2str(winsize), ...
            ', f0: ', num2str(f0_est, '%.3G'), 'Hz']});
        axis tight; grid on;
    end

    [bias, stdev] = evaluate_estimate(vz, v(vessel_depth:end-vessel_depth));
    biases = [biases, bias];
    stdevs = [stdevs, stdev];
    
end

all_stdevs{end+1} = stdevs(:);
all_biases{end+1} = biases(:);
%%

figure;
ax1=subplot(4,1,(1:2));
for i = 1:size(all_biases,2)
    plot(4:20, all_biases{i});
    ylabel('MSE [m/s]');
    grid on; hold on;
end
set(gca, 'XTickLabel', [])
title('Estimate precision as a function of emissions');
legend({'Winsize=1', 'Winsize=T0', 'Modify pulse', 'EC', 'Noise', 'MF'});

ax2=subplot(4,1,(3:4));
for i = 1:size(all_stdevs,2)
    plot(4:20, all_stdevs{i});
    grid on; hold on;
    ylabel('Std [m/s]')
end

linkaxes([ax1, ax2], 'x');
xlim([3, 21]);
xlabel('Number of lines in estimate');

%%

function [bias, stdev] = evaluate_estimate(target_vz, vz)
    bias = mean((target_vz - vz).^2);
    stdev = std((target_vz - vz).^2);
end

