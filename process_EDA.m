function [EDA_integrated, EDA_diff_aligned] = process_EDA(patient_folder)
%%process_EDA differentiates and integrates EDA data
%   Takes a patient data folder (e.g. PPG\p01) as input, and extracts 
%   relevant EDA and timestamp data from the folder. Calculates derivative
%   data using gradient() and shifts the alignment to match the MATB-II 
%   data. Then isolates the 30 minute window of the actual experiment, and 
%   subdivides it into 5, 6-minute blocks. Uses trapz() to calculate the 
%   area under the curve for each block.
%   
%   EDA_integrated is a vector of 5 values, representing cumulative
%   sympathetic arousal for each event period.
%   
%   EDA_differentiated flattens baseline EDA data to isolate rapid spikes
%   in electrodermal activity (sweating).

% Load data using dynamic file paths
EDA = load(fullfile(patient_folder, 'EDA.csv'));

% The first tag is the task start.
tags = load(fullfile(patient_folder, 'tags.csv'));
tag_timestamp = tags(1);

% Extract actual EDA physiological data
EDA_raw = EDA(3:end); % 4 Hz

% Extract sampling rates and wristband start time
EDA_rate = EDA(2);
EDA_start_timestamp = EDA(1);

% Create time vector for data based on sampling rate
t_EDA = (0:length(EDA_raw) - 1) / EDA_rate;

% Differentiate EDA data
dt_EDA = 1 / EDA_rate;
EDA_differentiated = gradient(EDA_raw, dt_EDA);

% Calculate the start index for the 30-minute task
offset_seconds = tag_timestamp - EDA_start_timestamp;
start_index = round(offset_seconds * EDA_rate) + 1;

% Find sympathetic arousal - integrate all 5 6-min blocks
samples_per_block = EDA_rate * 6 * 60; % 6 minutes
EDA_integrated = zeros(1, 5);

for i = 1:5
    % Find the start and end of the current block
    current_start = start_index + (i - 1) * samples_per_block;
    current_end = current_start + samples_per_block - 1;
    t_block = current_start:current_end;
    
    % Integrate
    EDA_integrated(i) = trapz(t_EDA(t_block), EDA_raw(t_block));
end

% Interpolate differentiated EDA to the MATB-II Timeline
t_master = 0:10:1800;
t_EDA_shifted = t_EDA - offset_seconds;
EDA_diff_aligned = interp1(t_EDA_shifted, EDA_differentiated, t_master, 'linear');

end