% close all
clear all
clc

%%

addpath(genpath('/autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/fsl.64bit/5.0.7/etc/matlab'));
addpath(genpath(pwd));

%%

sag = 44;
axi = 41;
cor = 62;

rf = '/autofs/cluster/Connectome2/Bay8_C2/bids/derivatives/';
of = fullfile(rf,'AxCaliber_SMT/sub-01/');

mask_n = fullfile(rf, '/preprocessed_dwi/sub-01_preprocessed_b0s_avg_brain.nii.gz');
mask_i = read_avw(mask_n);
mask = mask_i(:,:,:);

% wm_n = fullfile(rf, 'data/20230906/wm.nii.gz');
% wm_i = read_avw(wm_n);
% wm = wm_i(sag,:,:);

[r,c,sl] = ind2sub(size(mask),find(mask>0));

files = struct2cell(dir([of, '/mat/whole/SMT_fit_*.mat']));

n_in = 1000;
nJobs = size(files,2);

a = zeros(size(mask));
D_h = zeros(size(mask));
f_r = zeros(size(mask));
f_csf = zeros(size(mask));
f_h = zeros(size(mask));

for iJob=1:nJobs

    load(fullfile(files{2,iJob},files{1,iJob}))

    for ii=1:size(axc_mcmc_est_thin,2)

        jj=(iJob-1)*n_in+ii;

        if size(axc_mcmc_est_thin{ii},1)~=0

            a(r(jj),c(jj),sl(jj)) = mean(axc_mcmc_est_thin{1,ii}.a_samples);
            D_h(r(jj),c(jj),sl(jj)) = mean(axc_mcmc_est_thin{1,ii}.D_h_samples);
            f_r(r(jj),c(jj),sl(jj)) = mean(axc_mcmc_est_thin{1,ii}.f_r_samples);
            f_csf(r(jj),c(jj),sl(jj)) = mean(axc_mcmc_est_thin{1,ii}.f_csf_samples);
            f_h(r(jj),c(jj),sl(jj)) = 1-f_r(r(jj),c(jj),sl(jj))-f_csf(r(jj),c(jj),sl(jj));

        end
    end
end

a(a<0)=0;
D_h(D_h<0)=0;

f_r(f_r<0)=0;
f_r(f_r>1)=1;

f_csf(f_csf<0)=0;
f_csf(f_csf>1)=1;

f_h(f_h<0)=0;
f_h(f_h>1)=1;

% a = a.*wm;
% D_h = D_h.*wm;
% f_r = f_r.*wm;
% f_csf = f_csf.*wm;
% f_h = f_h.*wm;

a_f = zeros(size(mask_i));
D_h_f = zeros(size(mask_i));
f_r_f = zeros(size(mask_i));
f_csf_f = zeros(size(mask_i));
f_h_f = zeros(size(mask_i));

a_f(:,:,:) = a;
D_h_f(:,:,:) = D_h;
f_r_f(:,:,:) = f_r;
f_csf_f(:,:,:) = f_csf;
f_h_f(:,:,:) = f_h;

hdr_ref=load_nifti(fullfile(rf, 'preprocessed_dwi/sub-01_preprocessed_b0s.nii.gz'));

hdr_ref.vol = a_f.*10^6;
hdr_ref.dim(5)=size(hdr_ref.vol, 4);

save_nifti(hdr_ref, [of, '/sub-01_a.nii.gz']);

hdr_ref.vol = D_h_f.*10^9;
hdr_ref.dim(5)=size(hdr_ref.vol, 4);

save_nifti(hdr_ref, [of, '/sub-01_Dh.nii.gz']);

hdr_ref.vol = f_r_f;
hdr_ref.dim(5)=size(hdr_ref.vol, 4);

save_nifti(hdr_ref, [of, '/sub-01_fr.nii.gz']);

hdr_ref.vol = f_csf_f;
hdr_ref.dim(5)=size(hdr_ref.vol, 4);

save_nifti(hdr_ref, [of, '/sub-01_fcsf.nii.gz']);

hdr_ref.vol = f_h_f;
hdr_ref.dim(5)=size(hdr_ref.vol, 4);

save_nifti(hdr_ref, [of, '/sub-01_fh.nii.gz']);

figure
subplot(2,3,1),imagesc(fliplr(flipud(squeeze(a(sag,:,:))'))*10^6,[2 6]),axis square off, colormap jet
subplot(2,3,2),imagesc(fliplr(flipud(squeeze(D_h(sag,:,:))'))*10^9,[0 1]),axis square off, colormap hot
subplot(2,3,4),imagesc(fliplr(flipud(squeeze(f_r(sag,:,:))')),[0 0.8]),axis square off, colormap hot
subplot(2,3,5),imagesc(fliplr(flipud(squeeze(f_h(sag,:,:))')),[0 0.8]),axis square off, colormap hot
subplot(2,3,6),imagesc(fliplr(flipud(squeeze(f_csf(sag,:,:))')),[0 0.8]),axis square off, colormap jet
