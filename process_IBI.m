function [LF_HF_ratios] = process_IBI(patient_folder)
%%process_IBI Calculates Heart Rate Variability (HRV) metrics
%   Takes a patient data folder as input and loads IBI and tag data.
%   Aligns the uneven inter-beat intervals to the task start time and 
%   interpolates them onto a continuous 4 Hz timeline using spline 
%   interpolation. Subdivides the 30-minute task into 5, 6-minute blocks. 
%   Applies spectral analysis (bandpower) to calculate the Low Frequency 
%   and High Frequency power, outputting the LF/HF ratio for each block.
%   
%   LF_HF_ratios is a vector of 5 values, representing the sympathovagal
%   balance (cognitive stress) for each event period.

% Load data using dynamic file paths
IBI_data = readmatrix(fullfile(patient_folder, 'IBI.csv'));
tags = load(fullfile(patient_folder, 'tags.csv'));
HR_info = load(fullfile(patient_folder, 'HR.csv')); 

% Extract timestamps
tag_timestamp = tags(1);
session_start_timestamp = HR_info(1);

% Extract IBI columns and remove initial timestamp artifact
valid_beats = IBI_data(:, 2) > 0;
beat_times_raw = IBI_data(valid_beats, 1);
ibi_durations = IBI_data(valid_beats, 2);

% Align time to the MATB-II task start
offset_seconds = tag_timestamp - session_start_timestamp;
beat_times_shifted = beat_times_raw - offset_seconds;

% Interpolate uneven IBI data onto a 4 Hz timeline
fs_interp = 4; % 4 Hz resampling rate
t_uniform = 0:(1/fs_interp):1800; % 30 minutes of task time
ibi_uniform = interp1(beat_times_shifted, ibi_durations, t_uniform, 'spline');

% Calculate LF/HF Ratio for all 5 6-min blocks
samples_per_block = fs_interp * 6 * 60; 
LF_HF_ratios = zeros(1, 5);

for i = 1:5
    % Find the start and end of the current block
    current_start = 1 + (i - 1) * samples_per_block;
    current_end = current_start + samples_per_block - 1;
    
    % Prevent exceeding array bounds
    if current_end > length(ibi_uniform)
        current_end = length(ibi_uniform);
    end
    
    % Isolate and detrend the current 6-minute chunk
    ibi_chunk = ibi_uniform(current_start:current_end);
    ibi_chunk_detrended = detrend(ibi_chunk);
    
    % Calculate Power Spectral Density
    LF_power = bandpower(ibi_chunk_detrended, fs_interp, [0.04, 0.15]);
    HF_power = bandpower(ibi_chunk_detrended, fs_interp, [0.15, 0.40]);
    
    % Calculate the ratio
    LF_HF_ratios(i) = LF_power / HF_power;
end
end