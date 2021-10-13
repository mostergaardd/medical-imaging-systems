close all; clc; clear;

load('cfm_carotis.mat');

% Constants
f0      = 5e6;          % Probe center frequency [Hz]
fs      = 40e6;         % System sample rate [Hz]
c       = 1540;         % Speed of sound [m/s]
fprf    = 6e3;          % Pulse repetition frequency [Hz]
M       = 8;            % Number of cycles in pulse
start_d = 1e-3;         % Data start depth [m]
end_d   = 30e-3;        % Data end depth [m]
x_min   = -9.75e-3;     % Data start left [m]
x_max   = 9.75e-3;      % Data start right [m]

% Variables
drc_db   = 60;          % Dynamic range [dB]
scales   = 128;         % Number of gray levels in b-mode 
flow_scales = 2048;     % Number of color levels in cfm
apply_mf = false;       % Don't apply matched filter
apply_ec = true;        % Do apply echo cancelling

% Use window size in estimate of one pulse length
winsize = abs(round(1/f0 * fs));
vessel_angle = deg2rad(45);               % assume 45 deg angle

% Process b-mode image
bmode_image = process_rf(bmode_data(:,1:end-1), drc_db, scales);

% House keeping of dimensions
time = (0:size(bmode_image,1)-1) ./ fs;
depth = (c .* time) ./ 2;                  % convert from t -> depth
depth = depth(:) + start_d;                % adjust for start depth
lateral_dist = linspace(x_min, x_max, size(bmode_image,2));

% Estimate velocity
all_vs = zeros(size(rf_cfm_data, [1,3]));
for i=1:size(rf_cfm_data, 3)
    data = double(rf_cfm_data(:,:,i))';
    
    [v, depth] = autocorr_estimator(data, winsize, c, fprf, f0, ...
            vessel_angle, fs, 0, apply_ec, apply_mf);
    v = [v; zeros(winsize,1)];
    
    all_vs(:,i) = v;
end

% Mask image
v_masked = all_vs .* vessel;

% Interpolate to same size as bmode image using linear interpolation
all_vs_interp = zeros(size(bmode_image));
X = linspace(x_min, x_max, size(all_vs,2));
for i=1:size(all_vs_interp,1)
    all_vs_interp(i,:) = interp1(X, v_masked(i,:), lateral_dist, 'nearest');
end

% Display 2x the largest detected velocity
v_act_max = 2 * min(min(all_vs_interp));

% Initialize an index-image for flow
x = zeros(size(bmode_image));

% Handle pos flow
pos_flow = all_vs_interp > 0;
x(pos_flow) = 1 + 2*flow_scales * ...
    ((all_vs_interp(pos_flow) - min(all_vs_interp(pos_flow))) /...
    (abs(v_act_max) - min(all_vs_interp(pos_flow))));

% Handle neg flow
neg_flow = all_vs_interp < 0;
x(neg_flow) = 1 + flow_scales * ...
    ((all_vs_interp(neg_flow) - v_act_max) ./ ...
    (max(all_vs_interp(neg_flow)) - v_act_max));

% Combine bmode and flow index images
X = bmode_image + uint32(x);

% Create a colormap and convert to an RGB image
range = linspace(1, 0.1, flow_scales);
blues = [zeros(flow_scales,2), range'];
reds = [flipud(range'), zeros(flow_scales,2)];
cmap = colormap([gray(scales); blues; reds]);

% Map indices to RGB for display
RGB = ind2rgb(X,cmap);

% Display the image
imagesc(lateral_dist*1e3, depth*1e3, RGB);
title('Axial flow');
xlabel('Lateral distance [mm]')
ylabel('Axial distance [mm]')
ylim([start_d*1e3, 27]);

% Setup colorbar
c = colorbar;
c.Label.String = 'Velocity [m/s]';
c.Limits = [scales/(scales + 2*flow_scales) 1.0];
c.Ticks = linspace(scales/(scales + 2*flow_scales), 1, 9);
tick_labels = {};
delta = (-v_act_max - v_act_max) / 8;
for i=1:9
    tick_labels{end+1} = num2str(v_act_max + (i-1)*delta,'%.2f');
end
c.TickLabels = tick_labels;