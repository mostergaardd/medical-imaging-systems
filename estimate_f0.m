function [f0_est, PULSE] = estimate_f0(pulse, n_fft, f, do_plot)
%ESTIMATE_F0 Summary of this function goes here
%   Detailed explanation goes here

pulse = pulse(:);
f = f(:);

% Calculate spectrum
PULSE = fftshift(fft(pulse, n_fft));

% Estimate f0 from pulse
[~, i] = max(abs(PULSE));
PULSE_dx = diff(abs(PULSE),1);

% Search for zero-crossings around max peak
found_fw = false;
ii = i;
while ~found_fw
    found_fw = sign(PULSE_dx(ii)) ~= sign(PULSE_dx(ii+1));
    fw_idx = ii;
    ii = ii + 1;
end

found_bw = false;
ii = i-1;
while ~found_bw
    found_bw = sign(PULSE_dx(ii)) ~= sign(PULSE_dx(ii-1));
    bw_idx = ii;
    ii = ii - 1;
end

% Extract segment of pulse
segment = abs(PULSE(bw_idx:fw_idx+1));
f_segment = f(bw_idx:fw_idx+1);

% Method Jørgen
f0_est = sum(segment .* f_segment) / sum(segment);

% Method MØ
%segment_sum = cumsum(segment);
%[~, I] = min(abs(segment_sum - prctile(segment_sum, 50)));
%f0_est = f_segment(I);

if do_plot
   figure; 
   
   subplot(211);
   plot(f, abs(PULSE));
   title('Full two-sided amplitude spectrum');
   xlabel('Frequency [hz]');
   ylabel('Amplitude [a.u.]');
   grid on; axis tight;
   
   subplot(212);
   plot(f_segment, segment);
   xline(f0_est, 'r', [num2str(f0_est, '%.3G'), ' Hz']);
   title('Segment of amplitude spectrum');
   xlabel('Frequency [hz]')
   ylabel('Amplitude [a.u.]');
   legend('Segment of spectrum', 'Estimated center frequency',...
       'Location', 'northwest');
   grid on; axis tight;
end
end

