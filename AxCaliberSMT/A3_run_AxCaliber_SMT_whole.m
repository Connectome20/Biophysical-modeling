close all
clear all
clc

%%

addpath(genpath('/autofs/cluster/pubsw/2/pubsw/Linux2-2.3-x86_64/packages/fsl.64bit/5.0.7/etc/matlab'));

%%

rf = '/autofs/cluster/Connectome2/Bay8_C2/bids/derivatives/';
of = fullfile(rf,'AxCaliber_SMT/sub-01/');
my_pool = parpool(32);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%PREPARE THE MODEL%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load scheme.mat

scheme = scheme_all;

nVolsPerShell = [69 69 69 69 35 35 35 35 35 35 35 69 69 69 69 69];
nShells = length(nVolsPerShell);
end_vol_num = zeros(nShells,1);
for i = 1:nShells
    end_vol_num(i) = sum(nVolsPerShell(1:i));
end
start_vol_num = end_vol_num-nVolsPerShell.'+1;
idx_n_D = [start_vol_num, end_vol_num];

%%%%SET DIFFUSION TIMES
axc_mod_param.big_delta{1}=13*10^(-3);
axc_mod_param.big_delta{2}=30*10^(-3);

axc_mod_param.small_delta{1}=6*10^(-3);
axc_mod_param.small_delta{2}=6*10^(-3);

gyro=2*pi*42.58*1e6; %rad/s/T
axc_mod_param.gmr=gyro; %rad/s/T

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%CALCULATE ZEROS OF BESSEL FUNCTION FOR REF%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%CALCULATE ZERO CROSSINGS OF BESSEL FUNCTION OF THE FIRST KIND

axc_mod_param.k_beta=10;

beta_mtx=zeros(axc_mod_param.k_beta);

zeros_bess_tmp=zerobess('DJ',1,axc_mod_param.k_beta);
axc_mod_param.beta_mtx=zeros_bess_tmp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%ASS. BRUTE-FORCE DISCRETIZATION FOR "a" and "D_r"and "D_h"%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%SET AXON DIAMETERS TO BE INCLUDED IN THE MODEL

a_max=2e-5; % i.e, maximum axon diameter: 40 um

axc_mod_param.a_vec= linspace(1e-7,a_max,20); % i.e., min. axon diam: 0.1um

%%%%%SET DIFFUSION COEFS TO BE INCLUDED IN THE MODEL
axc_mod_param.D_r=1.7e-9;
axc_mod_param.D_h_vec=linspace(1e-11,3e-9,20);
axc_mod_param.sig_vec=linspace(0.01,0.1,10);
axc_mod_param.D_csf=3e-9; % Set D_csf to free diffusivity of water

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%MCMC ESTIMATE "SINGLE FIBER" MODEL%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ind_start_tmp=4;
axc_mcmc_param.a_start=5e-6;
axc_mcmc_param.D_h_start=1e-9;
axc_mcmc_param.f_r_start=0.5;
axc_mcmc_param.f_csf_start=0.1;
axc_mcmc_param.sig_start=0.05; %0.05;

axc_mod_param.a_old=axc_mcmc_param.a_start;
axc_mod_param.D_h_old=axc_mcmc_param.D_h_start;
axc_mod_param.f_r_old=axc_mcmc_param.f_r_start;
axc_mod_param.f_csf_old=axc_mcmc_param.f_csf_start;
axc_mod_param.sig_old=axc_mcmc_param.sig_start;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%MCMC PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

axc_mcmc_param.n_samples=2e5; 
axc_mcmc_param.mu0=0;

n_burnin=2e4;
n_thin=100;

%%

sag = 44;
axi = 41;
cor = 62;

diff_n = fullfile(rf, 'data/noise/20231025/dwi.nii.gz');
mask_n = fullfile(rf, 'data/noise/20231025/dwi_b0s_avg_brain_mask.nii.gz');
noise_n = fullfile(rf, 'data/noise/20231025/dwi_10b0s_std_median.nii.gz');

diff_i = read_avw(diff_n);
diff_i = diff_i(:,:,:,[12:116 126:470 471:746 748:887]);
diff_i = diff_i(:,:,:,[451:end 1:450]);

mask_i = read_avw(mask_n);
mask = mask_i(:,:,:);

noise_i = read_avw(noise_n);
noise = repmat(noise_i,[1 1 1 size(diff_i,4)]);
noise = noise(:,:,:,:);

diff = sqrt(max((diff_i.^2 - noise.^2), 0));

diff_r = reshape(diff,[size(diff,1)*size(diff,2)*size(diff,3),size(diff,4)]);

[r,c,sl] = ind2sub(size(mask),find(mask>0));

diff_f = diff_r(find(mask>0),:)';
nVoxs = size(diff_f,2);

n_in = 1000;
nJobs = floor(nVoxs./n_in)+1;

for iJob=1:round(nJobs./2)

    if iJob == nJobs
        diff_ft = squeeze(diff_f(:,(iJob-1)*n_in+1:end));
    else
        diff_ft = squeeze(diff_f(:,(iJob-1)*n_in+1:iJob*n_in));
    end

    for ii=1:size(diff_ft,2)

        tmp = []; DELTA_v = zeros(1,nShells); Gvec_v = zeros(1,nShells);
        tmp_bval_v = zeros(1,nShells);

        for ish = 1:nShells
    
            idx = idx_n_D(ish,1:2);
            bvecs = scheme(idx(1):idx(2),1:3);
            gvals = scheme(idx(1):idx(2),4);
            Gvec_v(ish)  = scheme(idx(1)+1,4);
            DELTA_v(ish) = scheme(idx(1)+1,5);
            littledelta = scheme(idx(1)+1,6);

            tmp_bval_v(ish) = (Gvec_v(ish)*gyro*littledelta)^2 ...
                *(DELTA_v(ish)-littledelta/3);

            b0s_idx = find(gvals == 0);
            dwi_idx = find(gvals ~= 0);
        
            s_v = squeeze(diff_ft(idx(1):idx(2),ii));

            s0 = mean(s_v(b0s_idx));
            s_v = s_v(dwi_idx)/s0;
            tmp = cat(1,tmp,mean(s_v));
    
        end

        Gvec_13 = Gvec_v(DELTA_v == 13*10^(-3));
        Gvec_30 = Gvec_v(DELTA_v == 30*10^(-3));
        
        %%%%q_vec from experimental data
        axc_mod_param.q_vec{1}=axc_mod_param.gmr*axc_mod_param.small_delta{1}*Gvec_13;
        axc_mod_param.q_vec{2}=axc_mod_param.gmr*axc_mod_param.small_delta{2}*Gvec_30;
        
        axc_mod_param.g_vec{1}=Gvec_13;
        axc_mod_param.g_vec{2}=Gvec_30;

        axc_mod_param.b_vec{1}=tmp_bval_v(DELTA_v == 13*10^(-3));
        axc_mod_param.b_vec{2}=tmp_bval_v(DELTA_v == 30*10^(-3));
    
        avgsig_over_sphere_m = tmp(DELTA_v == 13*10^(-3))';
        avgsig_over_sphere_m = cat(1,avgsig_over_sphere_m, tmp(DELTA_v == 30*10^(-3))');
    
        data_to_use(:,:,ii) = avgsig_over_sphere_m;

    end
    
    axc_exp_data=[];

    for ii=1:size(diff_ft,2)
            
        axc_exp_data(ii).Y{1}=data_to_use(1,:,ii)';
        axc_exp_data(ii).Y{2}=data_to_use(2,:,ii)';
       
    end

    axc_mcmc_est=cell(size(diff_ft,2),1);

    switch size(diff_ft,2)
        case 1
            %%%SET WALLCLOCK TIME
            disp(['    fitting starts @ ', datestr(now, 'HH:MM:SS mm/dd/yy')]);
            for ii=1:size(diff_ft,2)
                
                try

                    axc_mcmc_est{ii}=mcmc_calc_axcaliber_single_fiber_CSF_SMT_hp(axc_mod_param,axc_exp_data(ii),axc_mcmc_param);
                    axc_mcmc_est_thin{ii}.a_samples=axc_mcmc_est{ii}.a_samples(n_burnin:n_thin:end);
                    axc_mcmc_est_thin{ii}.D_h_samples=axc_mcmc_est{ii}.D_h_samples(n_burnin:n_thin:end);
                    axc_mcmc_est_thin{ii}.f_r_samples=axc_mcmc_est{ii}.f_r_samples(n_burnin:n_thin:end);
                    axc_mcmc_est_thin{ii}.f_csf_samples=axc_mcmc_est{ii}.f_csf_samples(n_burnin:n_thin:end);
                    axc_mcmc_est_thin{ii}.sig_samples=axc_mcmc_est{ii}.sig_samples(n_burnin:n_thin:end);
                
                catch
                end
                
            end
        otherwise    % PARFOR
            %%%SET WALLCLOCK TIME
            disp(['    fitting starts @ ', datestr(now, 'HH:MM:SS mm/dd/yy')]);
            parfor ii=1:size(diff_ft,2)
                
                try

                    axc_mcmc_est{ii}=mcmc_calc_axcaliber_single_fiber_CSF_SMT_hp(axc_mod_param,axc_exp_data(ii),axc_mcmc_param);
                    axc_mcmc_est_thin{ii}.a_samples=axc_mcmc_est{ii}.a_samples(n_burnin:n_thin:end);
                    axc_mcmc_est_thin{ii}.D_h_samples=axc_mcmc_est{ii}.D_h_samples(n_burnin:n_thin:end);
                    axc_mcmc_est_thin{ii}.f_r_samples=axc_mcmc_est{ii}.f_r_samples(n_burnin:n_thin:end);
                    axc_mcmc_est_thin{ii}.f_csf_samples=axc_mcmc_est{ii}.f_csf_samples(n_burnin:n_thin:end);
                    axc_mcmc_est_thin{ii}.sig_samples=axc_mcmc_est{ii}.sig_samples(n_burnin:n_thin:end);
                                
                catch
                end

            end
    end

    save([of, '/mat/whole/SMT_fit_', num2str(iJob,'%03d'), '.mat'], 'axc_mcmc_est_thin')

    clear axc_mcmc_est axc_mcmc_est_thin

end

disp(['    fitting done @ ', datestr(now, 'HH:MM:SS mm/dd/yy')]);
