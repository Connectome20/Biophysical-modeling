close all
clear all
clc

%%

addpath(genpath('/autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/fsl.64bit/5.0.7/etc/matlab'));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%PREPARE THE MODEL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dpSub = '/autofs/cluster/connectome2/Bay8_C2/bids/derivatives/processed_dwi/sub-01';

fpBvec = fullfile(dpSub, 'sub-01_preprocessed_dwi.bvecs');
fpBval= fullfile(dpSub, 'sub-01_preprocessed_dwi.bvals');

bval_i = dlmread(fpBval);
bvec_i = dlmread(fpBvec);

bval_is = bval_i(:,[12:116 126:470 471:746 748:887]);
bvec_is = bvec_i(:,[12:116 126:470 471:746 748:887]);

bval_13 = bval_is(:,451:end);
bvec_13 = bvec_is(:,451:end);

bval_30 = bval_is(:,1:450);
bvec_30 = bvec_is(:,1:450);

Bd_13=13*10^(-3).*ones(size(bval_13));
Bd_30=30*10^(-3).*ones(size(bval_30));
Sd_13=6*10^(-3).*ones(size(bval_13));
Sd_30=6*10^(-3).*ones(size(bval_30));

bval = [bval_13 bval_30];
bvec = [bvec_13 bvec_30];
Bd = [Bd_13 Bd_30];
Sd = [Sd_13 Sd_30];

gyro = 2*pi*42.58*1e6; % gyromagnetic ratio (rad/s/T)
Gd = sqrt(bval./((Bd-(Sd/3)).*(gyro^2).*(Sd.^2)))*10^3;

scheme_all = [bvec', Gd', Bd', Sd'];
scheme_13 = scheme_all(1:size(bval_13,2),:);
scheme_30 = scheme_all(size(bval_13,2)+1:end,:);

save scheme.mat scheme_all scheme_13 scheme_30
