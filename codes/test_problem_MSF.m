%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Paul Diaz 
% Date: 5/1/2018 
% Email: pdiaz303@gmail.com or paul.diaz@colorado.edu
% Please cite our paper:
% Diaz P., Doostan A., & Hampton J. "Sparse polynomial chaos expansions via
% compressed sensing and D-optimal design." Computer Methods in Applied
% Mechanics and Engineering. Volume 336, 1 July 2018, Pages 640-666
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Testing  D-Optimal SP (DSP) vs. (SP)
% To run this script in the background use the linux command:
% matlab -nodesktop -nodisplay <test_problem_MSF.m &> file.out & 
% This script uses subspace pursuit (SP) and modified subspace pursuit 
% to investigate the reconstruction of manufactured sparse functions

clear;clc;
addpath(genpath('./'))
% Dimension and order of the problem 
% try either (d,p) = (20,2) or (d,p) = (2,20)
dim = 20    % stochastic dimension d
order = 2 % total order of basis p
index_pc = nD_polynomial_array(dim,order);
P = size(index_pc,1) % size of the PC basis set
M = 20*P % Set as 20*P, candidate pools of different sizes are not provided  
polytype = 'L' % 'h' for hermite, 'L' for legendre
sampling_type = 'C' % 'C' for coherence optimal, 'S' for standard sampling

load_name = ['MSF_data_',num2str(dim),'d_',num2str(order),'p_',num2str(M),'M_',polytype,'_',sampling_type]
load(load_name)

delete(gcp('nocreate'))
pool_data = parpool;
rng('shuffle'); % Shuffles on local worker only
%Shuffle on each parallel worker
seed_offset = randi(floor(intmax/10),1);
parfor kk = 1:pool_data.NumWorkers
  rng(kk + seed_offset);
end

s = 60  %Number of non-zero coefficients in PCE
sample_sizes =  round(linspace(2*s,220,10)) % Values of N
noise = 0.03  % Additive Noise Level
weight = false; %Normalization for SP sensing 
R = 4 %Number of repeated experiments paper uses R = 240

Num_samps = length(sample_sizes); %Specifies an outer loop size

%Initialize statistics to be saved
SP_stats = zeros(Num_samps,3);
GDSP_stats = zeros(Num_samps,3);
DSP_stats = zeros(Num_samps,3);
ORACLE_stats = zeros(Num_samps,3);

OMP_stats = zeros(Num_samps,3);
DOMP_stats = zeros(Num_samps,3);
GDOMP_stats = zeros(Num_samps,3);

%Initialize arrays to store error at each of the R repetitions
sp_err = zeros(R,Num_samps);
gdsp_err = zeros(R,Num_samps);
dsp_err = zeros(R,Num_samps);
oracle_err = zeros(R,Num_samps);

sp_support = zeros(R,Num_samps);
gdsp_support = zeros(R,Num_samps);
dsp_support = zeros(R,Num_samps);


fprintf('Running...\n')
tic

for n = 1:Num_samps %Iterate to next sample size
    N = sample_sizes(n); %Specific sample size
    parfor r = 1:R   % Repeat test R times
        data_inds = datasample(1:M,round(0.5*M),'replace',false);
        psi = sample.lhs(data_inds,:);
        %%% Manufacture exact PC coefficients, c_ref
        indices = datasample(1:P,s,'replace',false);
        c_ref = zeros(P,1);
        c_ref(indices) = randn(s,1);
        u = psi*c_ref;
        %Introduce additive noise
        the_noise = noise*abs(u).*randn(length(u),1);
        u = u + the_noise;
        %Apply weights which correspond to a sampling scheme
        psi = diag(sample.w(data_inds))*psi;
        u = diag(sample.w(data_inds))*u;

        %Construct oracle solution
        c_hat = zeros(P,1);
        oracle_design = rrqr_dopt_adapt(psi(:,indices),N,[]);
        %oracle_design = Designs(n,1:N);
        c_hat(indices) =  pinv(psi(oracle_design,indices))*u(oracle_design);
        e12 = norm(c_hat - c_ref,2)/norm(c_ref,2);

        %Construct a random design matrix of candidate rows and RHS's
        inds = datasample(1:round(0.5*M),N,'replace',false);
        psi_rand = psi(inds,:); u_rand = u(inds);
        %Solve with subspace pursuit MC or Coh-Opt
        [c_hat ] = SP(psi_rand,0,u_rand,weight);
        e1 = norm(c_hat - c_ref,2)/norm(c_ref,2)
        supp = find(c_hat ~= 0);
        sp_support(r,n) = length( intersect(indices, supp))/s;
        %Solve with modified subspace pursuit Seq-D-Coh-Opt
        [c_hat] = DSP(psi,0,N,u,round(0.8*N),weight);
        e2 = norm(c_hat - c_ref,2)/norm(c_ref,2)
        supp = find(c_hat ~= 0);
        dsp_support(r,n) = length( intersect(indices, supp))/s;
        %Solve with supspace pursuit and N-point, D-optimal design D-MC and
        %D-Coh-Opt 
        design = rrqr_dopt_adapt(psi,N,[]);
        psi_dopt = psi(design,:); u_dopt = u(design);
        [c_hat ] = SP(psi_dopt,0,u_dopt,weight);
        e4 = norm( c_hat - c_ref,2)/norm(c_ref,2)
        supp = find(c_hat ~= 0);
        gdsp_support(r,n) = length( intersect(indices, supp))/s;

        %save relative errors
        sp_err(r,n) = e1;
        dsp_err(r,n) =  e2;
        gdsp_err(r,n) = e4;
        oracle_err(r,n) = e12;
    end

    SP_stats(n,1) =  mean(sp_err(:,n)); %Computation of mean relative error
    SP_stats(n,2) =  std(sp_err(:,n)); %Computation of Stdev of relative error
    SP_stats(n,3) =   mean(sp_support(:,n));
    DSP_stats(n,1) =  mean(dsp_err(:,n));
    DSP_stats(n,2) =  std(dsp_err(:,n));
    DSP_stats(n,3) =   mean(dsp_support(:,n));
    GDSP_stats(n,1) =  mean(gdsp_err(:,n)); 
    GDSP_stats(n,2) =  std(gdsp_err(:,n));
    GDSP_stats(n,3) =   mean(gdsp_support(:,n));
    ORACLE_stats(n,1) = mean(oracle_err(:,n));
    ORACLE_stats(n,2) = std(oracle_err(:,n));

    % Display current iteration statistics
    fprintf('Avg rel err | Stdev rel err | percent (supp) successfully identified \n')
    sp_iter =  SP_stats(n,:)
    dsp_iter = DSP_stats(n,:)
    gdsp_iter = GDSP_stats(n,:)
    oracle_iter = ORACLE_stats(n,:)

    save_name = ['noise_03_',num2str(s),'s_',num2str(dim),'d_',num2str(order),'p_',num2str(R),'R_',num2str(M),'M_',polytype,'_',sampling_type];
    %save(save_name) %Save the data to an output file
    fprintf('Complete sample size %i\n', N)
    time = toc
end
%Output data summary
standard_SP = table(sample_sizes',SP_stats(:,1),SP_stats(:,2),SP_stats(:,3),'VariableNames',{'N' 'Rel_err' 'Std_err' 'Perc_supp'})
D_opt_SP = table(sample_sizes',DSP_stats(:,1),DSP_stats(:,2),DSP_stats(:,3), 'VariableNames',{'N' 'Rel_err' 'Std_err' 'Perc_supp'})
global_D_SP =  table(sample_sizes',GDSP_stats(:,1),GDSP_stats(:,2),GDSP_stats(:,3),'VariableNames',{'N' 'Rel_err' 'Std_err' 'Perc_supp'})
ORACLE = table(sample_sizes',ORACLE_stats(:,1),ORACLE_stats(:,2), 'VariableNames',{'N' 'Rel_err' 'Std_err' })
save_name = ['noise_03_',num2str(s),'s_',num2str(dim),'d_',num2str(order),'p_',num2str(R),'R_',num2str(M),'M_',polytype,'_',sampling_type];
%save(save_name)

%% This part of the script is for generating the .mat data files
% The .mat files have already been generated,
% this could block need not be ran.

% clear;clc;
% addpath(genpath('./'))
% delete(gcp('nocreate'))
% pool_data = parpool;
% rng('shuffle'); % Shuffles on local worker only
% %Shuffle on each parallel worker
% seed_offset = randi(floor(intmax/10),1);
% parfor kk = 1:pool_data.NumWorkers
%   rng(kk + seed_offset);
% end
% %%
% format short g;
% % Dimension and order of the problem 
% dim = 2    % stochastic dimension
% order = 20 % total order of basis 
% index_pc = nD_polynomial_array(dim,order);
% P = size(index_pc,1) % size of the PC basis set
% M = 20*P
% polytype = 'L' % 'h' for hermite, 'L' for legendre
% sampling_type = 'C' % 'C' for coherence optimal, 'S' for standard sampling
% 
% % %This code is used for Coherence Optimal Sampling
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basis_opt.dim = dim; % Number of initial dimensionst
% eval_opt.max_dim = dim; % Maximum dimension of problem
% sample_opt.r_dim = dim;
% 
% sample_opt.initial_size = M; % Size of candidate pool
% eval_opt.p_type = repmat(polytype,1,eval_opt.max_dim);  % for hermite polynomials 
% if dim < order && polytype == 'L'
%     sample_opt.prop_handle = @asym_prop; 
% end 
% if dim < order && polytype =='h'
%     sample_opt.prop_handle = @herm_sphere_prop;
% end
% if dim >= order
%    sample_opt.prop_handle = @orth_prop; % input proposals are from orthogonality measure 
% end
% 
% eval_opt.qoi_handle = @zero_eval; % Rhs is zero here, easy enough to change
% 
% % Evaluation Related Options
% eval_opt.grad = false; % Requires compatible RHS evaluation
% eval_opt.grad_dim = 0; % Number of dimensions with gradients. (must correspond to first dimensions)
% sample_opt.w_correction_handle = @l2_w; % Minimize l2-coherence
% 
% % Most of the rest are incidental/ancillary options
% % Sample Related Options
% sample_opt.sr = 0; % Minimal Sampling Rate
% sample_opt.min_sample_percent = 0.1; % Percent new samples to generate
% sample_opt.max_sample_percent = 1;
% sample_opt.n_workers = pool_data.NumWorkers; % Number of parallel workers from pool to use.
% sample_opt.burn_in = 1000; % Number of burn-in samples
% sample_opt.log_col_rate = -4; % Logarithm of collision rate (MCMC samples to be discarded)
% sample_opt.w_handle = @l2_w; % Minimize l2-coherence
% basis_opt.type_handle = @basis_anisotropic_total_order; % Initial basis is 'to' total order
% 
% basis_opt.ord = order*ones(1,basis_opt.dim); % Initial total order of approximation
% basis_opt.pc_flag = false; % Whether or not to identify preconditioner
% basis_opt.validation_strikes = 6; % Number of strikes for basis validation
% basis_opt.expand_coeff = 1.3; % Order expansion is done multiplicatively (with ceiling function) Lower numbers insure+1 effectiveness.
% basis_opt.dim_add = 2; % Number of dimensions to add order 1 approximation to.
% basis_opt.validation_iters = 3; % Number of iterations for basis adaptation
% basis_opt.order_bounds = false;
% basis_opt.basis_identify_handle = @basis_validation_anisotropic_total_order;
% sample_opt.order = basis_opt.ord; 
% basis = basis_init(basis_opt, eval_opt); % Initialization of basis
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% psi = zeros(M,P); 
% sample = struct;
% %Generate Samples
% if sampling_type == 'C' % Coherence optimal sampling 
% sample = sample_init(basis,sample_opt, eval_opt); % Initialization of Sample (coherence optimal)
% W = zeros(M,1);
% parfor m = 1:M
%      B = norm(sample.lhs(m,:),2);
%      W(m) = 1/B;
% end
% sample.w = W;
% clear W 
% end
% if sampling_type == 'S' % Standard sampling 
%   if polytype == 'L'
%     sample.rv = 2*rand(M,dim) - 1;  
%     parfor isim = 1:M
%         psi(isim,:) = piset(sample.rv(isim,:),index_pc);
%     end
%     sample.lhs = psi;
%   end 
%   if polytype =='h'
%     sample.rv = randn(M,dim);  
%     parfor isim = 1:M
%         psi(isim,:) = piset_hermite(sample.rv(isim,:),index_pc); 
%     end
%   end
%   sample.w = ones(M,1);
%   sample.lhs = psi;
% end
% 
% save_name = ['MSF_data_',num2str(dim),'d_',num2str(order),'p_',num2str(M),'M_',polytype,'_',sampling_type]
% %save(save_name)
