clc; clear; close all;

cd(fileparts(which(mfilename)));

% ============================= Load Data =================================

% Preallocate matrices to hold patient data
all_EDA_integrated = zeros(37, 5); 
all_EDA_diff = zeros(37, 181);
all_smoothed_BVP = zeros(37, 181);
all_LF_HF_ratios = zeros(37, 5);

% Define base directory - Make sure your MATLAB Current Folder is set to 
% the folder where you saved the patient data (PPG folder)
base_dir = 'PPG';

% Loop through every patient ID
for patient_id = 1:37
    
    % Skip the excluded patients
    if patient_id == 14 || patient_id == 23
        continue; 
    end
    
    % Dynamically generate the folder name
    p_folder = sprintf('p%02d', patient_id);
    folder_name = fullfile(base_dir, p_folder);
    
    % Call functions to process EDA, BVP, IBI, and HR data
    [all_EDA_integrated(patient_id, :), all_EDA_diff(patient_id, :)] = process_EDA(folder_name);
    all_smoothed_BVP(patient_id, :) = process_BVP(folder_name);
    all_LF_HF_ratios(patient_id, :) = process_IBI(folder_name);
    
end

% ======================= Regression Analysis =============================

matb_dir = 'MATB-II';
all_task_error_phases = zeros(37, 5);

for patient_id = 1:37
    if patient_id == 14 || patient_id == 23
        continue
    end
    file_name = fullfile(matb_dir, sprintf('p%02dresman.csv', patient_id));
    data = readtable(file_name);
    
    % Convert MM:SS time strings to seconds
    time_strings = data.ELAPSED_TIME;
    t_matb = zeros(length(time_strings), 1);
    for i = 1:length(time_strings)
        parts = split(time_strings{i}, ':');
        t_matb(i) = str2double(parts{1}) * 60 + str2double(parts{2});
    end
    
    % Task error = average absolute deviation of tanks A and B from 2500
    task_error = (abs(data.DIFF_A) + abs(data.DIFF_B)) / 2;
    
    % Average task error within each 6-minute phase
    phase_duration = 360;
    for i = 1:5
        phase_start = (i - 1) * phase_duration;
        phase_end   =  i      * phase_duration;
        in_phase = t_matb >= phase_start & t_matb < phase_end;
        all_task_error_phases(patient_id, i) = mean(task_error(in_phase), 'omitnan');
    end
end
disp('MATB-II phase data loaded.')

t_master = 0:10:1800;
all_BVP_phase_means = zeros(37, 5);
for patient_id = 1:37
    if patient_id == 14 || patient_id == 23
        continue
    end
    for i = 1:5
        phase_start = (i - 1) * 360;
        phase_end   =  i      * 360;
        in_phase = t_master >= phase_start & t_master < phase_end;
        all_BVP_phase_means(patient_id, i) = mean(all_smoothed_BVP(patient_id, in_phase), 'omitnan');
    end
end

valid_patients = true(37, 1);
valid_patients(14) = false;
valid_patients(23) = false;

task_error_vec = reshape(all_task_error_phases(valid_patients, :), [], 1);
EDA_vec        = reshape(all_EDA_integrated(valid_patients, :), [], 1);
BVP_vec        = reshape(all_BVP_phase_means(valid_patients, :), [], 1);
IBI_vec        = reshape(all_LF_HF_ratios(valid_patients, :), [], 1);

% Remove NaN rows
valid_rows = ~isnan(task_error_vec) & ~isnan(EDA_vec) & ~isnan(BVP_vec) & ~isnan(IBI_vec);
task_error_vec = task_error_vec(valid_rows);
EDA_vec        = EDA_vec(valid_rows);
BVP_vec        = BVP_vec(valid_rows);
IBI_vec        = IBI_vec(valid_rows); 

SS_tot = sum((task_error_vec - mean(task_error_vec)).^2);

% EDA
p_EDA  = polyfit(EDA_vec, task_error_vec, 1);
R2_EDA = 1 - sum((task_error_vec - polyval(p_EDA, EDA_vec)).^2) / SS_tot;

% BVP
p_BVP  = polyfit(BVP_vec, task_error_vec, 1);
R2_BVP = 1 - sum((task_error_vec - polyval(p_BVP, BVP_vec)).^2) / SS_tot;

% IBI (HRV)
p_IBI   = polyfit(IBI_vec, task_error_vec, 1);
R2_IBI  = 1 - sum((task_error_vec - polyval(p_IBI, IBI_vec)).^2) / SS_tot;

fprintf('EDA R² = %.3f\n', R2_EDA)
fprintf('BVP R² = %.3f\n', R2_BVP)
fprintf('HRV R² = %.3f\n', R2_IBI)

figure;

% EDA
subplot(1, 3, 1)
scatter(EDA_vec, task_error_vec, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.5)
hold on
x_range = linspace(min(EDA_vec), max(EDA_vec), 100);
plot(x_range, polyval(p_EDA, x_range), 'r', 'LineWidth', 2)
xlabel('EDA Sympathetic Arousal (\muS)')
ylabel('Mean Tank Deviation (units)')
title(sprintf('EDA vs Task Error\nR² = %.3f', R2_EDA))
legend('Data', 'Regression Line', 'Location', 'best')
grid on

% BVP
subplot(1, 3, 2)
scatter(BVP_vec, task_error_vec, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.5)
hold on
x_range = linspace(min(BVP_vec), max(BVP_vec), 100);
plot(x_range, polyval(p_BVP, x_range), 'r', 'LineWidth', 2)
xlabel('BVP Amplitude (normalized)')
ylabel('Mean Tank Deviation (units)')
title(sprintf('BVP vs Task Error\nR² = %.3f', R2_BVP))
legend('Data', 'Regression Line', 'Location', 'best')
grid on

% IBI (HRV)
subplot(1, 3, 3)
scatter(IBI_vec, task_error_vec, 40, 'b', 'filled', 'MarkerFaceAlpha', 0.5)
hold on
x_range = linspace(min(IBI_vec), max(IBI_vec), 100);
plot(x_range, polyval(p_IBI, x_range), 'r', 'LineWidth', 2)
xlabel('LF/HF Ratio (Stress)')
ylabel('Mean Tank Deviation (units)')
title(sprintf('HRV vs Task Error\nR² = %.3f', R2_IBI))
legend('Data', 'Regression Line', 'Location', 'best')
grid on

sgtitle('Phase-Level Physiological Signals vs Task Performance')

% ======================= Statistical Testing =============================

valid_patients = [1:13, 15:22, 24:37];
N = length(valid_patients);

[h_EDA, p_EDA] = ttest(all_EDA_integrated(valid_patients,1), all_EDA_integrated(valid_patients,2));
[h_HRV, p_HRV] = ttest(all_LF_HF_ratios(valid_patients,1), all_LF_HF_ratios(valid_patients,2));
[h_Task, p_Task] = ttest(all_task_error_phases(valid_patients, 1), all_task_error_phases(valid_patients, 2));

figure('Name', 'Statistical Validation of Stress Protocol');

% Subplot 1: EDA Response
subplot(1, 3, 1)
data_eda = [mean(all_EDA_integrated(valid_patients, 1), 'omitnan'), mean(all_EDA_integrated(valid_patients, 2), 'omitnan')];
sd_eda   = [std(all_EDA_integrated(valid_patients, 1), 'omitnan'), std(all_EDA_integrated(valid_patients, 2), 'omitnan')];
sem_eda  = sd_eda ./ sqrt(N);
bar(data_eda, 'FaceColor', [0.4 0.6 0.8]); hold on;
errorbar(1:2, data_eda, sem_eda, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
ylabel('Integrated EDA');
set(gca, 'XTickLabel', {'Baseline', 'Stress'});
title(sprintf('EDA \np = %.4f', p_EDA));

% Subplot 2: HRV Response 
subplot(1, 3, 2)
data_hrv = [mean(all_LF_HF_ratios(valid_patients, 1), 'omitnan'), mean(all_LF_HF_ratios(valid_patients, 2), 'omitnan')];
sd_hrv   = [std(all_LF_HF_ratios(valid_patients, 1), 'omitnan'), std(all_LF_HF_ratios(valid_patients, 2), 'omitnan')];
sem_hrv  = sd_hrv ./ sqrt(N);
bar(data_hrv, 'FaceColor', [0.8 0.4 0.4]); hold on;
errorbar(1:2, data_hrv, sem_hrv, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
ylabel('LF/HF Ratio');
set(gca, 'XTickLabel', {'Baseline', 'Stress'});
title(sprintf('Heart Stress (HRV)\np = %.4f', p_HRV));

% Subplot 3: Task Performance
subplot(1, 3, 3)
data_task = [mean(all_task_error_phases(valid_patients, 1), 'omitnan'), mean(all_task_error_phases(valid_patients, 2), 'omitnan')];
sd_task   = [std(all_task_error_phases(valid_patients, 1), 'omitnan'), std(all_task_error_phases(valid_patients, 2), 'omitnan')];
sem_task  = sd_task ./ sqrt(N);
bar(data_task, 'FaceColor', [0.5 0.8 0.5]); hold on;
errorbar(1:2, data_task, sem_task, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
ylabel('Mean Tank Deviation');
set(gca, 'XTickLabel', {'Baseline', 'Stress'});
title(sprintf('Task Performance Error\np = %.4f', p_Task));

sgtitle('Validation of Stress Induction: Baseline vs. Stress');


