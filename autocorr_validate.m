clear; close all; clc;

global f0 M fs n_fft vessel_angle vessel_diameter f_prf
global T_prf vz c seed winsize_

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
seed = 0;                   % random seed

%% Evaluate windowsize on-off

all_biases = {};
all_stdevs = {};

apply_mf                = false;
apply_ec                = false;
add_noise               = false;
add_stationary_signal   = false;
modify_pulse            = false;
doplot                  = false;
do_estimate_f0          = true;
winsize_                = 1;
n_emissions = 20;

biases = []; stdevs = [];

for overwrite_winsize=[1, 0]
    
    [biases, stdevs] = run_experiment(modify_pulse, do_estimate_f0, err_std, add_noise, add_stationary_signal, overwrite_winsize, apply_mf, apply_ec, doplot);
    
    all_stdevs{end+1} = stdevs(:);
    all_biases{end+1} = biases(:);
end

display_bias_plot(all_biases, all_stdevs, {'Impact of window size'}, {'Winsize=1', 'Winsize=T_0'}, 'estimate_winsize.png');

%% Evaluate importance of modifying the pulse

all_biases = {};
all_stdevs = {};

apply_mf                = false;
apply_ec                = false;
add_noise               = false;
add_stationary_signal   = false;
doplot                  = false;
overwrite_winsize       = false;
do_estimate_f0          = true;
n_emissions = 20;

for modify_pulse=[false, true]
    
    [biases, stdevs] = run_experiment(modify_pulse, do_estimate_f0, err_std, add_noise, add_stationary_signal, overwrite_winsize, apply_mf, apply_ec, doplot);
    
    all_stdevs{end+1} = stdevs(:);
    all_biases{end+1} = biases(:);
end

display_bias_plot(all_biases, all_stdevs, {'Impact of modifying pulse', 'Winsize: T_0'}, {'Original', 'Modified'}, 'estimate_pulse_modify.png');


%% Evaluate importance of estimating f0

all_biases = {};
all_stdevs = {};

apply_mf                = false;
apply_ec                = false;
add_noise               = false;
add_stationary_signal   = false;
modify_pulse            = false;
doplot                  = false;
overwrite_winsize       = false;
n_emissions = 20;

for do_estimate_f0=[false, true]
    
    [biases, stdevs] = run_experiment(modify_pulse, do_estimate_f0, err_std, add_noise, add_stationary_signal, overwrite_winsize, apply_mf, apply_ec, doplot);
    
    all_stdevs{end+1} = stdevs(:);
    all_biases{end+1} = biases(:);
end

display_bias_plot(all_biases, all_stdevs, {'Impact of using estimate of f_0', 'Winsize: T_0'}, {'f_0', 'f_0-est'}, 'estimate_do_estimate_f0.png');

%% Evaluate importance stationary signal

all_biases = {};
all_stdevs = {};

apply_mf                = false;
apply_ec                = false;
add_noise               = false;
add_stationary_signal   = true;
modify_pulse            = false;
doplot                  = false;
overwrite_winsize       = false;
do_estimate_f0          = false;
n_emissions = 20;

[biases, stdevs] = run_experiment(modify_pulse, do_estimate_f0, err_std, add_noise, add_stationary_signal, overwrite_winsize, apply_mf, apply_ec, doplot);
    
all_stdevs{end+1} = stdevs(:);
all_biases{end+1} = biases(:);

% Add echo cancelling
apply_ec                = true;

[biases, stdevs] = run_experiment(modify_pulse, do_estimate_f0, err_std, add_noise, add_stationary_signal, overwrite_winsize, apply_mf, apply_ec, doplot);
    
all_stdevs{end+1} = stdevs(:);
all_biases{end+1} = biases(:);

display_bias_plot(all_biases, all_stdevs, {'Impact of echo cancelling', 'Winsize: T_0'}, {'No EC', 'With EC'}, 'estimate_EC.png');

%% SNR estimate

add_noise             = false;
add_stationary_signal = false;
apply_ec              = false;
overwrite_winsize     = false;
n_emissions           = 8;

load('pulse.mat');
matched_h = matched_filter(pulse, 0);

[single_line, vessel_depth] = simulate_single_line(... 
            vessel_angle, vessel_diameter, f_prf, fs, vz, c,...
            n_emissions, pulse, seed, err_std, add_noise, add_stationary_signal);
        
winsize = abs(round(1/f0 * fs));
signal_power = max(max(single_line));

all_biases = {};
all_stdevs = {};

for apply_mf = [0 1]
    
    biases = []; stdevs = [];
    
    for snr = [-10 -5 0 5 10 15 20]
        p_noise = 10^(signal_power/20) / 10^(snr/20);

        noise = p_noise * randn(size(single_line));
        data = single_line + noise;

        [v, depth] = autocorr_estimator(data, winsize, c, f_prf, f0, ...
                vessel_angle, fs, matched_h, apply_ec, apply_mf);

        [bias, stdev] = evaluate_estimate(vz, v(vessel_depth:end-(vessel_depth)));
        
        biases = [biases, bias];
        stdevs = [stdevs, stdev];
    end
        
    all_stdevs{end+1} = stdevs(:);
    all_biases{end+1} = biases(:);
    
end

SNRS = [-10 -5 0 5 10 15 20];

figure;
ax1=subplot(2,1,1);
for i = 1:size(all_biases,2)
    plot(SNRS, all_biases{i});
    ylabel('MSE [m/s]');
    grid on; hold on;
end
set(gca, 'XTickLabel', [])
title({'Importance of matched filter in presence of noise', 'N_c: 8'});
legend({'No matched filter', 'With matched filter'});

ax2=subplot(2,1,2);
for i = 1:size(all_stdevs,2)
    plot(SNRS, all_stdevs{i});
    grid on; hold on;
    ylabel('Std [m/s]')
end
legend({'No matched filter', 'With matched filter'});
linkaxes([ax1, ax2], 'x');
xlim([-11, 21]);
xlabel('SNR [dB]');
saveas(gcf, 'estimate_snr.png');

%% Functions

function display_bias_plot(biases, stdevs, title_str, legend_str, save_name)
    figure;
    ax1=subplot(2,1,1);
    for i = 1:size(biases,2)
        plot(4:20, biases{i});
        ylabel('MSE [m/s]');
        grid on; hold on;
    end
    set(gca, 'XTickLabel', [])
    title(title_str);
    legend(legend_str);

    ax2=subplot(2,1,2);
    for i = 1:size(stdevs,2)
        plot(4:20, stdevs{i});
        grid on; hold on;
        ylabel('Std [m/s]')
    end
    legend(legend_str);
    linkaxes([ax1, ax2], 'x');
    xlim([3, 21]);
    xlabel('Number of emissions in estimate');
    saveas(gcf, save_name);
end

function [biases, stdevs] = run_experiment(modify_pulse, do_estimate_f0, err_std, add_noise, add_stationary_signal, overwrite_winsize, apply_mf, apply_ec, doplot)
    global f0 M fs n_fft vessel_angle vessel_diameter f_prf
    global vz c seed winsize_

    biases = [];
    stdevs = [];

    for n_emissions=4:20
        
        load('pulse.mat');

        if modify_pulse
            
            mod_sig = sin(2*pi*f0.*(0:1/fs:M/f0)) .* 0.1*max(pulse);
            pulse = conv(pulse, mod_sig,'same');
            
            if doplot && n_emissions == 20
                t = (0:length(pulse)-1) ./ fs;
                figure; plot(t, pulse); hold on; 
                plot((0:length(pulse)-1) ./ fs,pulse);
                legend('Pulse','Modified pulse');
                xlabel('Time [s]'); ylabel('Amplitude [V]');
                grid on; axis tight;
                title('Modifying a measured ultrasound pulse');
                saveas(gcf,'pulse_modify.png');
            end
        end

        t = (0:length(pulse))./fs;   % time
        f = (-n_fft/2:n_fft/2-1)./n_fft.*fs;
        [f0_est, pulse_F] = estimate_f0(pulse, n_fft, f, ...
            doplot && n_emissions == 20);
        f0_est = abs(f0_est);
        
        if ~do_estimate_f0
            f0_est = f0;
        end 

        matched_h = matched_filter(pulse, doplot && n_emissions == 20);

        if doplot && n_emissions==20
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

        if doplot && n_emissions == 20
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

        if doplot && n_emissions == 20
            % Plotting
            figure();
            %plot(depth(vessel_depth+150:end-(vessel_depth+150)), v(vessel_depth+150:end-(vessel_depth+150)));
            plot(depth(:), v(:));
            xlabel('Depth [m]');
            ylabel('Velocity [m/s]');
            title({['Estimated velocities. Target v=', num2str(vz), ' m/s'], ...
                ['Winsize: ', num2str(winsize), ...
                ', f0: ', num2str(f0_est, '%.3G'), 'Hz']});
            axis tight; grid on;
        end

        [bias, stdev] = evaluate_estimate(vz, v(vessel_depth:end-(vessel_depth)));
        %[bias, stdev] = evaluate_estimate(vz, v);
        biases = [biases, bias];
        stdevs = [stdevs, stdev];

    end
end

function [bias, stdev] = evaluate_estimate(target_vz, vz)
    bias = mean((target_vz - vz).^2);
    stdev = std((target_vz - vz).^2);
end

