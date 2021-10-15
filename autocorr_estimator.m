function [v,depth] = autocorr_estimator...
    (rf, winsize, c, f_prf, f0, rad_angle, fs, m_h, apply_echo_c, apply_mf)
%AUTOCORR_ESTIMATOR Applies autocorrelation flow estimation
%   rf: RF signal to calculate CFM on. Real signal. Should be [emissions,
%   rf_signals]. 
%   winsize: Optional averaging over samples in depth. A good choice is one
%       pulse width                      [samples]
%   c: Speed of sound                    [m/s]
%   f_prf: Pulse repetition frequency    [hz]
%   f0: center frequency of us pulse     [hz]
%   rad_angle: vessel angle in radians   [rads]
%   fs: sample frequency                 [hz]
%   m_h: matched filter                 
%   appply_echo_c: flag for applying ec  [bool]
%   apply_mf: 
%       flag for applying matched filter [bool]

% Apply matched filter
if apply_mf
    for i=1:size(rf,1)
        rf(i,:) = conv(rf(i,:), m_h(:), 'same');
    end
end

% Apply echo cancelling
rf_ec = rf;
if apply_echo_c
    rf_ec = rf - mean(rf, 1);
end

% Obtain the IQ signal
rf_iq = hilbert(rf_ec);
N_rf = size(rf_iq, 2);

factor = - (c * f_prf) / (4 * pi * f0);

v = zeros(1, N_rf - winsize);

% Loop over depth
for n = 1:N_rf - winsize                
    rf_n = rf_iq(:, n:n+winsize-1);

    vs = 0;
    % Average over window in depth
    for m = 1:winsize                   
        rf_nm = rf_n(:, m);
        
        % Unpack IQ signal
        y = imag(rf_nm);
        x = real(rf_nm);

        Ri = 0;
        Rr = 0;
        % Average over emissions
        for i = (1:length(rf_nm)-1)        
            % Calculate acorr at lag 1
            Ri = Ri + (y(i+1)*x(i) - x(i+1)*y(i));
            Rr = Rr + (x(i+1)*x(i) + y(i+1)*y(i));
        end

        % Collect the results and estimate velocity
        vs_m = factor * atan2(Ri, Rr);
        vs_m = vs_m / cos(rad_angle);

        vs = vs + vs_m; 
    end

    % Apply averaging factor
    v(n) = vs / winsize; 
end

% Also return a vector with the depth information
time = ((0:length(v)-1) + round((winsize/2))) ./ fs;
depth = (time * c) / 2;

v = v(:);
depth = depth(:);

end

