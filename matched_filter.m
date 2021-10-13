function [h] = matched_filter(pulse, doplot)
%MATCHED_FILTER Summary of this function goes here
%   Detailed explanation goes here
pulse = pulse(:);
h = flip(pulse);

if doplot
    figure;
    
    subplot(211);
    plot(pulse); title('Ultrasound pulse');
    axis tight; grid on;
    
    subplot(212);
    plot(h); title('Matched filter');
    axis tight; grid on;
end

end

