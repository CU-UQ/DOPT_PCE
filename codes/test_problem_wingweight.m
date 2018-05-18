%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Paul Diaz 
% Date: 5/1/2018 
% Email: pdiaz303@gmail.com or paul.diaz@colorado.edu
% Please cite our paper:
% Diaz P., Doostan A., & Hampton J. "Sparse polynomial chaos expansions via
% compressed sensing and D-optimal design." Computer Methods in Applied
% Mechanics and Engineering. Volume 336, 1 July 2018, Pages 640-666
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Testing D-Optimal SP (DSP) vs SP
% To run this script in the background use the linux command:
% matlab -nodesktop -nodisplay <test_problem_wingweight.m &> file.out & 
% This script uses subspace pursuit (SP) and modified subspace pursuit 
% to investigate the reconstruction of the ishigami function
clear;clc;
addpath(genpath('./'))
delete(gcp('nocreate'))
pool_data = parpool;
rng('shuffle'); % Shuffles on local worker only
%Shuffle on each parallel worker
seed_offset = randi(floor(intmax/10),1);
parfor kk = 1:pool_data.NumWorkers
  rng(kk + seed_offset);
end

load('wingweight_recon_data')
load('wingweight_validation_data')
sample_sizes =  round(linspace(40,280,10))
R = 1000

Num_samps = length(sample_sizes);
SP_stats = zeros(Num_samps,2);
GSP_stats = zeros(Num_samps,2);
DSP_stats = zeros(Num_samps,2);

sp_err = zeros(R,Num_samps);
gsp_err = zeros(R,Num_samps);
dsp_err = zeros(R,Num_samps);
fprintf('Running...\n')
tic

for n = 1:Num_samps %Iterate to next sample size
    N = sample_sizes(n); %Specific sample size
    parfor r = 1:R   % Repeat test R times
        data_inds = datasample(1:M,round(0.5*M),'replace',false);
        psi_cands = psi(data_inds,:); u_cands = U(data_inds);
        %Construct a random matrix of candidate rows for MC or Coh-Opt
        inds = datasample(1:size(psi_cands,1),N,'replace',false);
        psi_rand = psi_cands(inds,:); u_rand = u_cands(inds);
        %Solve with subspace pursuit MC or Coh-Opt
        [c_hat] = SP(psi_rand,0,u_rand,false);
        e1 = norm(psi_validation*c_hat - u_validation)/norm(u_validation)
        %Solve with supspace pursuit and N-point, D-optimal design D-MC and
        %D-Coh-Opt
        design = rrqr_dopt_adapt(psi_cands,N,[]);
        psi_dopt = psi_cands(design,:); u_dopt = u_cands(design);
        [c_hat] = SP(psi_dopt,0,u_dopt,false);
        e3 = norm(psi_validation*c_hat - u_validation)/norm(u_validation)
        %Solve with modified subspace pursuit Seq-D-Coh-Opt
        [c_hat] = DSP(psi_cands,0,N,u_cands,round(0.8*N),false);
        e2 = norm(psi_validation*c_hat - u_validation)/norm(u_validation)

        sp_err(r,n) = e1;
        dsp_err(r,n) = e2;
        gsp_err(r,n) = e3;
        if mod(r,100)==0
            fprintf('r = %i  \n', r)
        end
    end
    SP_stats(n,1) =  mean(sp_err(:,n)); 
    SP_stats(n,2) =  std(sp_err(:,n));

    DSP_stats(n,1) =  mean(dsp_err(:,n));
    DSP_stats(n,2) =  std(dsp_err(:,n));

    GSP_stats(n,1) =  mean(gsp_err(:,n)); 
    GSP_stats(n,2) =  std(gsp_err(:,n));
    
    clear Psi 
    save_name = ['wingweight',num2str(dim),'d_',num2str(order),'p_',num2str(R),'R'];
    %save(save_name)
    fprintf('Complete sample size %i\n', N)
    disp(table(sample_sizes(1:n)',SP_stats(1:n,1),SP_stats(1:n,2),'VariableNames',{'N' 'Rel_err' 'Std_err'}))
    disp(table(sample_sizes(1:n)',GSP_stats(1:n,1),GSP_stats(1:n,2),'VariableNames',{'N' 'Rel_err' 'Std_err'}))
    disp(table(sample_sizes(1:n)',DSP_stats(1:n,1),DSP_stats(1:n,2),'VariableNames',{'N' 'Rel_err' 'Std_err'}))
    time = toc
end

fprintf('Columns: 1 =  Avg relative error, 2 = Std Dev Rel Err, 3 = Avg Cond #, 4 = Std Dev Cond #\n')
disp(table(sample_sizes',SP_stats(:,1),SP_stats(:,2),'VariableNames',{'N' 'Rel_err' 'Std_err'}))
disp(table(sample_sizes',GSP_stats(:,1),GSP_stats(:,2),'VariableNames',{'N' 'Rel_err' 'Std_err'}))
disp(table(sample_sizes',DSP_stats(:,1),DSP_stats(:,2),'VariableNames',{'N' 'Rel_err' 'Std_err'}))
save_name = ['wingweight',num2str(dim),'d_',num2str(order),'p_',num2str(R),'R'];
%save(save_name)

%% This part of the script is for generating the .mat data files
% The .mat files have already been generated,
% this could block need not be ran.
% Generate Validation Data
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
% 
% tic
% dim = 10
% order = 3
% % Order of expansion 
% index_pc = nD_polynomial_array(dim,order);
% 
% % size of the PC basis set
% P = size(index_pc,1); 
% M = 20000;
% sample.rv = 2*rand(M,dim) - 1; 
% U = zeros(M,1);
% psi = zeros(M,P);
% parfor m = 1:M
%     psi(m,:) = piset(sample.rv(m,:),index_pc);
%     [qoi,~] = wingweight(sample.rv(m,:));
%     U(m) = qoi;
% end
% psi_validation = psi;  
% u_validation = U;
% time = toc
% %save('wingweight_validation_data','psi_validation','u_validation')
% 
% % Generate Reconstruction Data
% tic 
% dim = 10
% order = 3
% % Order of expansion 
% index_pc = nD_polynomial_array(dim,order);
% 
% % size of the PC basis set
% P = size(index_pc,1); 
% 
% M = 20*P
% 
% %%%% Sampling Setup
% basis_opt.dim = dim; % Number of initial dimensions
% eval_opt.max_dim = dim; % Maximum dimension of problem
% sample_opt.r_dim = dim;
% 
% sample_opt.initial_size = M; % Size of candidate pool
% polytype = 'L' % 'h' for hermite, 'L' for legendre
% sampling_type = 'C' % 'C' for coherence optimal, 'S' for standard sampling, 'A' for asymptotic sampling
% eval_opt.p_type = repmat(polytype,1,eval_opt.max_dim);  % for hermite polynomials 
% 
% sample_opt.prop_handle = @orth_prop; % input proposals are from orthogonality measure 
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
% 
% % Create samples and data
% sample = sample_init(basis,sample_opt, eval_opt); % Initialization of Sample (coherence optimal)
% W = zeros(M,1);
% U = zeros(M,1);
% psi = zeros(M,P);
% parfor m = 1:M
%     psi(m,:) = piset(sample.rv(m,:),index_pc);
%     B = norm(psi(m,:),2);
%     W(m) = 1/B;
%     [qoi,~] = wingweight(sample.rv(m,:));
%     U(m) = qoi;
% end
% sample.lhs = psi;
% sample.w = W;
% clear W
% psi = diag(sample.w)*sample.lhs;
% U = diag(sample.w)*U;
% time = toc
% %save('wingweight_recon_data')