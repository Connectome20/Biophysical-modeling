function out=mcmc_calc_axcaliber_single_fiber_CSF_SMT_hp(mod_param,obs_data,mcmc_param)

%%%% mcmc_calc_axcaliber_single_fiber_master.m
%%%% Performs Markov chain Monte Carlo sampling for parameters in the
%%%% van Gelderen model given the data
%%%% Uses Gaussian noise model

%%%REFORMAT OBSERVED DATA


Y_obs=cat(1,obs_data.Y{1},obs_data.Y{2});

%%%ALLOCATE SAMPLES
[ns,tmp]=size(Y_obs);

a_samples=zeros(mcmc_param.n_samples,1);
D_h_samples=zeros(mcmc_param.n_samples,1);
f_r_samples=zeros(mcmc_param.n_samples,1);
f_csf_samples=zeros(mcmc_param.n_samples,1);
sig_samples=zeros(mcmc_param.n_samples,1);

logp_samples=zeros(mcmc_param.n_samples,1);

Y_r_pred_tmp=zeros(mcmc_param.n_samples,ns);
Y_h_pred_tmp=zeros(mcmc_param.n_samples,ns);
Y_csf_pred_tmp=zeros(mcmc_param.n_samples,ns);
Y_pred_tmp=zeros(mcmc_param.n_samples,ns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INITIALIZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_param.a_old=mcmc_param.a_start;
test_param.D_h_old=mcmc_param.D_h_start;
test_param.f_r_old=mcmc_param.f_r_start;
test_param.f_csf_old=mcmc_param.f_csf_start;
test_param.sig_old=mcmc_param.sig_start;

test_param.a=test_param.a_old;
test_param.D_h=test_param.D_h_old;
test_param.sig=test_param.sig_old;

axc_fwd_mod=calc_axcaliber_model_CSF_SMT_hp(mod_param,test_param);

Y_r_pred=reshape(axc_fwd_mod.E_r_mtx',ns,1);
Y_h_pred=reshape(axc_fwd_mod.E_h_mtx',ns,1);
Y_csf_pred=reshape(axc_fwd_mod.E_csf_mtx',ns,1);
Y_pred=test_param.f_r_old*Y_r_pred+(1-test_param.f_r_old-test_param.f_csf_old)*Y_h_pred+test_param.f_csf_old*Y_csf_pred;

logp_old=sum(sum((Y_obs-Y_pred-mcmc_param.mu0*ones(size(Y_obs))).^2)/(2*test_param.sig^2))-(length(Y_obs)/2)*log(1/test_param.sig^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%MAIN MCMC LOOP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%PROPOSAL DISTRIBUTION VALUES
sig_a=min(mod_param.a_vec(2)-mod_param.a_vec(1))/2;
sig_D_h=min(mod_param.D_h_vec(2)-mod_param.D_h_vec(1))/4;
sig_f_r=0.1/2;
sig_f_csf=0.1/2;
sig_sig=0.01/2;

for ss=1:mcmc_param.n_samples;
    
%%%%%SAMPLE NEW CANDIDATE STATE S

test_param.a_new=test_param.a_old+sig_a*randn(1,1);
test_param.D_h_new=test_param.D_h_old+sig_D_h*randn(1,1);
while (test_param.D_h_new>=mod_param.D_r)
    test_param.D_h_new=test_param.D_h_old+sig_D_h*randn(1,1);
end
test_param.f_r_new=test_param.f_r_old+sig_f_r*randn(1,1);
test_param.f_csf_new=test_param.f_csf_old+sig_f_csf*randn(1,1);
test_param.sig_new=test_param.sig_old+sig_sig*randn(1,1);

test_param.a=test_param.a_new;
test_param.D_h=test_param.D_h_new;
test_param.sig=test_param.sig_new;

axc_fwd_mod=calc_axcaliber_model_CSF_SMT_hp(mod_param,test_param);

%%%%EVALUATE NEW LIKELIHOOD
Y_r_pred=reshape(axc_fwd_mod.E_r_mtx',ns,1);
Y_h_pred=reshape(axc_fwd_mod.E_h_mtx',ns,1);
Y_csf_pred=reshape(axc_fwd_mod.E_csf_mtx',ns,1);
Y_pred=test_param.f_r_new*Y_r_pred+(1-test_param.f_r_new-test_param.f_csf_new)*Y_h_pred+test_param.f_csf_new*Y_csf_pred;

in_limits_a=min(mod_param.a_vec)<test_param.a_new && test_param.a_new<max(mod_param.a_vec); 
in_limits_D_h=min(mod_param.D_h_vec)<test_param.D_h_new && test_param.D_h_new<max(mod_param.D_h_vec);
in_limits_f_r=0<test_param.f_r_new && test_param.f_r_new <1;
in_limits_f_csf=0<test_param.f_csf_new && test_param.f_csf_new <1;
in_limits_sig=min(mod_param.sig_vec)<test_param.sig_new && test_param.sig_new<max(mod_param.sig_vec);

if in_limits_a && in_limits_D_h && in_limits_f_r && in_limits_sig && in_limits_f_csf;   
    
    logp_new=sum(sum((Y_obs-Y_pred-mcmc_param.mu0*ones(size(Y_obs))).^2)/(2*test_param.sig^2))-(length(Y_obs)/2)*log(1/test_param.sig^2);

else
    
    logp_new=10^20; %%%%LARGE NUMBER

end;

%%%%TEST ACCEPTANCE

% accept_ratio=min(1,exp(-logp_new)/exp(-logp_old)); 
% 06-23-2014:
% exp(-logp_old) may be NaN if logp_old is a large number ==> 
% Matlab ignores this result and sets accept_ratio = 1. 
% Correct problem by taking difference of logp_new and logp_old
accept_ratio=min(1,exp(-(logp_new-logp_old)));

%alpha=unifrnd(0,1);
alpha=rand(1,1);

if  alpha < accept_ratio;
    
    %%UPDATE   
    test_param.a_old=test_param.a_new;
    test_param.D_h_old=test_param.D_h_new;
    test_param.f_r_old=test_param.f_r_new;
    test_param.f_csf_old=test_param.f_csf_new;
    test_param.sig_old=test_param.sig_new;
    logp_old=logp_new;
    
end
    
a_samples(ss)=test_param.a_old;
D_h_samples(ss)=test_param.D_h_old;
f_r_samples(ss)=test_param.f_r_old;
f_csf_samples(ss)=test_param.f_csf_old;
sig_samples(ss)=test_param.sig_old;
logp_samples(ss)=logp_old;


Y_r_pred_tmp(ss,:,:)=Y_r_pred;
Y_h_pred_tmp(ss,:,:)=Y_h_pred;
Y_csf_pred_tmp(ss,:,:)=Y_csf_pred;
Y_pred_tmp(ss,:,:)=Y_pred;

end

out.a_samples=a_samples;
out.D_h_samples=D_h_samples;
out.f_r_samples=f_r_samples;
out.f_csf_samples=f_csf_samples;
out.logp_samples=logp_samples;
out.sig_samples=sig_samples;
% out.seed=s;