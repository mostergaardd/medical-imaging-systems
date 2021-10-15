function [output_lines] = process_rf(lines, dynamic_range, scale)
%PROCESS_RF Process RF data for B-mode imaging
%   Apply envelope detection, dynamic range compression and grayscale
%   mapping on the given RF data
%   Arguments
%   lines: RF data
%   dynamic_range: dynamic rance in decibels
%   scale: number of grayscales to map to

line_env = abs(hilbert(lines));

x_min = min(min(line_env));
x_max = max(max(line_env));
env_scaled = (line_env - x_min) / (x_max - x_min);

rf_db = 20 .* log10(env_scaled);
db_tresh = max(max(rf_db)) - dynamic_range;
rf_db( rf_db <= db_tresh ) = db_tresh;

x_min = min(min(rf_db));
x_max = max(max(rf_db));

output_lines = scale * (rf_db - x_min) / (x_max - x_min);
output_lines = uint32(output_lines);

end

