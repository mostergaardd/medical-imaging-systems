function [output_lines] = process_rf(lines, dynamic_range, scale)
%PROCESS_RF Summary of this function goes here
%   Detailed explanation goes here

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

