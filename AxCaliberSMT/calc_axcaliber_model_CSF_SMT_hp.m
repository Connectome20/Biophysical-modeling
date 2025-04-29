function [ out ] = calc_axcaliber_model_CSF_SMT_hp(mod_param,test_param)

%%% calc_axcaliber_model_master.m
%%% Calculates Gaussian phase distribution model based on van Gelderen JMRB 1994
%%% Incorporates CSF model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%FIRST DEAL WITH RESTRICTED MODEL%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%DO THE ACTUAL CALCULATION

% len = length(mod_param.q_vec);
q_l_1 = length(mod_param.q_vec{1});
q_l_2 = length(mod_param.q_vec{2});

% for i = 1:len
%     cmd = ['q_l_', num2str(i),'=length(mod_param.q_vec{',num2str(i),'});'];
%     eval(cmd);
% end


small_delta=[mod_param.small_delta{1}*ones(1,q_l_1),...
    mod_param.small_delta{2}*ones(1,q_l_2)];
big_delta=[mod_param.big_delta{1}*ones(1,q_l_1),...
    mod_param.big_delta{2}*ones(1,q_l_2)];

q_vec=[mod_param.q_vec{1},mod_param.q_vec{2}];
g_vec=[mod_param.g_vec{1},mod_param.g_vec{2}];


gmr=mod_param.gmr;

E_r_mtx=zeros(size(q_vec));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%HERE THE ACTUAL COMPUTATION HAPPENS%%%%%%%%%
                    
%%%%%%%%%%%%%%%%FIRST CASE N=0%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
ds=size(big_delta);

E0_tmp=zeros(mod_param.k_beta,ds(1),ds(2));
                    
alphm=mod_param.beta_mtx./((test_param.a)/2);

for kk_ind=1:mod_param.k_beta;
    
    alpha_2 = alphm(kk_ind)^2;
    factor_1=2*mod_param.D_r*alpha_2*small_delta;
    factor_2=2*exp(-mod_param.D_r*alpha_2*small_delta);
    factor_3=2*exp(-mod_param.D_r*alpha_2*big_delta);
    factor_4=exp(-mod_param.D_r*alpha_2*(big_delta-small_delta));
    factor_5=exp(-mod_param.D_r*alpha_2*(big_delta+small_delta));
    
    factor_6=mod_param.D_r^2*alphm(kk_ind)^6*(((test_param.a)/2)^2*alpha_2-1);
    
    E0_tmp(kk_ind,:,:)=(factor_1-2+factor_2+factor_3-factor_4-factor_5)./factor_6;
    
end

E0_tmp_sum=reshape(sum(E0_tmp,1),1,ds(2));
L_perp = -2*gmr^2*E0_tmp_sum;   % L+

L_para = -(big_delta-small_delta/3).*(gmr*small_delta).^2*mod_param.D_r; % L//

if all((L_perp-L_para)>0)
    zr = g_vec.*sqrt(L_perp-L_para);
    out.E_r_mtx= sqrt(pi)./(2*zr).*exp(g_vec.^2.*L_perp).*erf(zr);
else
    disp('L_perp-L_para:');
    disp(num2str(L_perp-L_para));
end
% E_r_mtx= (-2*gmr^2*g_vec.^2.*E0_tmp_sum);

%%%%STORE OUTPUT
% out.E_r_mtx=E_r_mtx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%THEN DO THE HINDERED MODEL%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%INITIALIZE THE FORWARD MODEL TO ZEROS
L_perp_h = -(big_delta-small_delta/3).*(gmr*small_delta).^2*test_param.D_h;
if all((L_perp_h-L_para)>0)
    zh = g_vec.*sqrt(L_perp_h-L_para);
    % E_h_mtx=exp(-q_vec.^2*test_param.D_h.*(big_delta-small_delta/3));
    out.E_h_mtx= sqrt(pi)./(2*zh).*exp(g_vec.^2.*L_perp_h).*erf(zh);
else
    disp('test_param.D_h:');
    disp(num2str(test_param.D_h));
    disp('mod_param.D_r:');
    disp(num2str(mod_param.D_r));
end

% out.E_h_mtx=E_h_mtx;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%FINALLY THE CSF MODEL%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out.E_csf_mtx=exp(-q_vec.^2*mod_param.D_csf.*(big_delta-small_delta/3));

%%%%STORE OUTPUT
% out.E_csf_mtx=E_csf_mtx;