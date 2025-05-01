% Main script for the SANDI analysis using random forest regression
% Reference: Ramos-Llorden, Lee, ..., Huang, Nature BME, 2025
% Based on Palombo M. et al. Neuroimage 2020: https://doi.org/10.1016/j.neuroimage.2020.116835
%
% Requirements:
% - SANDI-Matlab-Toolbox v1.0 (https://github.com/palombom/SANDI-Matlab-Toolbox-v1.0)
%   Please make sure the functions are installed and the path is correctly set.
%
% Author: Hansol Lee (0000-0003-2112-1197)

clear all; close all; clc;

%% Add paths to SANDI functions
restoredefaultpath % Restore default MATLAB path
addpath(genpath('/autofs/cluster/connectome2/Bay8_C2/bids/code/SANDI/functions')); % Add SANDI functions

%% User-defined parameters
sub = 'sub-001';                                 % List of subjects
rf = '/autofs/cluster/connectome2/Bay8_C2/bids'; % Root folder
delta = 13;                                      % Diffusion gradient separation (ms)
smalldelta = 6;                                  % Diffusion gradient duration (ms)
Dsoma = 3;                                       % Soma diffusivity (µm^2/ms)
Din_UB = 3;                                      % Intra-neurite diffusivity UB (µm^2/ms)
Rsoma_UB = 12;                                   % Soma radius UB (µm)
De_UB = 3;                                       % Extra-cellular diffusivity UB (µm^2/ms)
seed_rng = 1;                                    % Random seed
MLmodel = 'RF';                                  % Machine learning model: 'RF', 'MLP', 'GRNN'
Nset = 1e4;                                      % Number of noise samples for SNR estimation

%% Start timer
overallTimer = tic;

%% Define data paths
data_f = fullfile(rf, sub, 'dwi_process/s13_SANDI');        % Path to DWI processed data
noise_f = fullfile(rf, 'derivatives/processed_dwi', sub);   % Path to estimated noise
out_f = fullfile(rf, 'derivatives/SANDI2', sub);            % Output folder for SANDI results
if ~exist(out_f, 'dir')
    mkdir(out_f);                                           % Create output folder if it does not exist
end

%% Load b-values and acquisition parameters
bvalues_filename = fullfile(data_f, 'bvals_D13ms.txt');         % b-values file (FSL format)
data_filename = fullfile(data_f, 'D13ms.nii.gz');               % Preprocessed DWI data file
mask_filename = fullfile(noise_f, [sub '_brain_mask.nii.gz']);  % Brain mask file
noisemap_filename = fullfile(noise_f, [sub '_sigma.nii.gz']);   % Noise map file from MP-PCA denoising

%% Calculate SNR
noisemap = normalize_noisemap(noisemap_filename, data_filename, mask_filename, bvalues_filename); % Normalize noise map by b=0 image

sigma = noisemap(:);
sigma = sigma(sigma > 0);
sigma = min(max(sigma, 0), 1);

sigma_sampled = sigma(randi(numel(sigma), [Nset, 1]));
sigma_median = median(sigma_sampled);
SNR = 1 / sigma_median;

%% Train Machine Learning model
disp(['Training SANDI ML model for ' sub '...']);

load_acquisition_parameters(bvalues_filename, delta, smalldelta, out_f); % Save acquisition parameters

schemefile = fullfile(out_f, 'diravg.scheme');
[Mdl, train_perf] = setup_and_run_model_training_rician( ...
    schemefile, SNR, out_f, Dsoma, Din_UB, Rsoma_UB, De_UB, seed_rng, MLmodel); % Train random forest regression model

%% Preprocess: Compute direction-averaged signal
make_direction_average(out_f, data_filename, bvalues_filename, delta, smalldelta); % Compute direction-averaged DWI

%% Fit SANDI model
disp(['Fitting SANDI model for ' sub '...']);

img_data = fullfile(out_f, 'diravg_signal.nii.gz');
run_model_fitting(img_data, mask_filename, schemefile, out_f, Mdl, MLmodel, sub); % Perform model fitting using trained ML model

%% Display total execution time
elapsedTime = toc(overallTimer);
disp(['Finished SANDI analysis for ' sub ' in ' num2str(elapsedTime, '%.2f') ' seconds.']);