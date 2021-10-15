function [h] = matched_filter(pulse, doplot)
%MATCHED_FILTER Return the time-reversed matched filter of a signal
%   Arguments
%   pulse: signal template
%   doplot: flag for plotting

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

