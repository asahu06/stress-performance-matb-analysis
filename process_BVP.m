function [BVP_aligned] = process_BVP(patient_folder)
%process_BVP applies a moving mean filter and aligns BVP data
%   Loads BVP and tag data from a patient folder. Uses movmean() with a 10-
%   second window to apply a moving average filter to the raw BVP amplitude
%   envelope. Shifts the time vector based on the tag timestamp and 
%   interpolates the smoothed data to a 10-second master timeline to match 
%   MATB-II logs.

% Load data using dynamic file paths
BVP = load(fullfile(patient_folder, 'BVP.csv'));
tags = load(fullfile(patient_folder, 'tags.csv'));
tag_timestamp = tags(1);

% Extract raw BVP data, rate, and start time
BVP_raw = BVP(3:end); % 64 Hz
BVP_rate = BVP(2);
BVP_start_timestamp = BVP(1);

% Calculate offset and create shifted time vector
offset_seconds = tag_timestamp - BVP_start_timestamp;
t_BVP = (0:length(BVP_raw) - 1) / BVP_rate;
t_BVP_shifted = t_BVP - offset_seconds; % shift t=0 to task start

% Smoothen Data
BVP_amplitude = envelope(BVP_raw);
window_seconds = 10; % 10 second window
window_samples = window_seconds * BVP_rate;
BVP_smoothed = movmean(BVP_amplitude, window_samples);

% Interpolate to the MATB-II Timeline - 30 minutes of task = 1800 s
t_master = 0:10:1800;
BVP_aligned = interp1(t_BVP_shifted, BVP_smoothed, t_master, 'linear');

end