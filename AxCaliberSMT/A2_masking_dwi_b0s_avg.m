close all
clear all
clc

%%

dpSub = '/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-01';

% fpdwi = fullfile(dpSub, 'sub-01_preprocessed_dwi.nii.gz');
fpb0s = fullfile(dpSub, 'sub-01_preprocessed_b0s.nii.gz');
fpb0s_avg = fullfile(dpSub, 'sub-01_preprocessed_b0s_avg.nii.gz');
fpb0s_avg_brain = fullfile(dpSub, 'sub-01_preprocessed_b0s_avg_brain.nii.gz');

% cmd = ['fslroi ' fpdwi ' ' fpb0s ' 0 10'];
% [status, result] = system(cmd, '-echo');

cmd = ['fslmaths ' fpb0s ' -Tmean ' fpb0s_avg];
[status, result] = system(cmd, '-echo');

cmd = ['bet ' fpb0s_avg ' ' fpb0s_avg_brain ' -f 0.1 -g 0 -m'];
[status, result] = system(cmd, '-echo');
