function [pulse, t, PULSE, f, f0_est] = generate_pulse(f0,M,fs,n_fft)
%GENERATE_PULSE Summary of this function goes here
%   Detailed explanation goes here

t = (0:1/fs:M/f0);  % time
f = (-n_fft/2:n_fft/2-1)./n_fft.*fs;

% Generate an ultrasound pulse
pulse = sin(2*pi*f0.*t);

% Apply a window to the pulse
pulse_w = pulse' .* hanning(length(pulse));
pulse = pulse_w(:);

[f0_est, PULSE] = estimate_f0(pulse, n_fft, f, 1);

end

