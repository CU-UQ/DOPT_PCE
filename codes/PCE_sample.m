function [ sample ] = PCE_sample( sample_data )
% This function generates inputs and constructs a polynomial basis.
% Inputs:
% sample_data.M = number_of_candidate_points;
% sample_data.dim = dimension_of_points;
% sample_data.polytype = polytype;   % 'h' for hermite, 'L' for [-1,1]^dim
% legendre polynomials
% sample_data.sampling_type = sampling_type; % 'C' for coherence optimal, 
% 'S' for standard sampling, 'A' for asymptotic sampling
% sample_data.order = order;
% Outputs:
% sample.w
% sample.lhs
% sample.rvs

M = sample_data.M;
dim = sample_data.dim;
polytype = sample_data.polytype;
sampling_type = sample_data.sampling_type;
order = sample_data.order;
NumWorkers = sample_data.NumWorkers;

% Order of expansion 
index_pc = nD_polynomial_array(dim,order);
% size of the PC basis set

%%%% Sampling Setup
sample_opt.r_dim = dim;
sample_opt.order = order*ones(1,dim); % Initial total order of approximation
sample_opt.initial_size = M; % Size of candidate pool
sample_opt.w_correction_handle = @l2_w; % Minimize l2-coherence
sample_opt.sr = 0; % Minimal Sampling Rate
sample_opt.min_sample_percent = 0.1; % Percent new samples to generate
sample_opt.max_sample_percent = 1;
sample_opt.n_workers = NumWorkers; % Number of parallel workers from pool to use.
sample_opt.burn_in = 1000; % Number of burn-in samples
sample_opt.log_col_rate = -4; % Logarithm of collision rate (MCMC samples to be discarded)
sample_opt.w_handle = @l2_w; % Minimize l2-coherence

% Evaluation Related Options
eval_opt.grad = false; % Requires compatible RHS evaluation
eval_opt.grad_dim = 0; % Number of dimensions with gradients. (must correspond to first dimensions)
eval_opt.max_dim = dim; % Maximum dimension of problem
eval_opt.p_type = repmat(polytype,1,eval_opt.max_dim);  % for hermite polynomials 
eval_opt.qoi_handle = @zero_eval; % Rhs is zero here, easy enough to change



% Basis Related Options
basis_opt.dim = dim; % Number of initial dimensionst
basis_opt.ord = order*ones(1,basis_opt.dim); % Initial total order of approximation
basis_opt.pc_flag = false; % Whether or not to identify preconditioner
basis_opt.validation_strikes = 6; % Number of strikes for basis validation
basis_opt.expand_coeff = 1.3; % Order expansion is done multiplicatively (with ceiling function) Lower numbers insure+1 effectiveness.
basis_opt.dim_add = 2; % Number of dimensions to add order 1 approximation to.
basis_opt.validation_iters = 3; % Number of iterations for basis adaptation
basis_opt.order_bounds = false;
basis_opt.basis_identify_handle = @basis_validation_anisotropic_total_order;
basis_opt.type_handle = @basis_anisotropic_total_order; % Initial basis is 'to' total order

if dim < order && polytype == 'L'
    sample_opt.prop_handle = @asym_prop; 
end 
if dim < order && polytype =='h'
    sample_opt.prop_handle = @herm_sphere_prop;
end
if dim >= order
   sample_opt.prop_handle = @orth_prop; % input proposals are from orthogonality measure 
end

basis = basis_init(basis_opt, eval_opt); % Initialization of basis


%Generate Samples
if sampling_type == 'C' % Coherence optimal sampling 
sample = sample_init(basis,sample_opt, eval_opt); % Initialization of Sample (coherence optimal)
W = zeros(M,1);
for m = 1:M
     B = norm(sample.lhs(m,:),2);
     W(m) = 1/B;
end
sample.w = W;
clear W 
end
if sampling_type == 'S' % Standard sampling 
  if polytype == 'L'
    sample.rv = 2*rand(M,dim) - 1;  
    for isim = 1:M
        psi(isim,:) = piset(sample.rv(isim,:),index_pc);
    end
    sample.lhs = psi;
  end 
  if polytype =='h'
    sample.rv = randn(M,dim);  
    for isim = 1:M
        psi(isim,:) = piset_hermite(sample.rv(isim,:),index_pc); 
    end
  end
  sample.w = ones(M,1);
  sample.lhs = psi;
end
if sampling_type == 'A' % Asymptotic
  if polytype == 'L'
    sample.w = zeros(1,M);
    sample.rv = 2*betarnd(0.5,0.5,M,dim)-1;
    for isim = 1:M
        weight = 1;
        for l = 1:dim
           weight = weight*(1-sample.rv(isim,l)^2)^(1/4);
        end
        psi(isim,:) = piset(sample.rv(isim,:),index_pc);
        sample.w(isim) = weight;
    end
        sample.lhs = psi;
    %psi = diag(sample.w)*psi;
  end 
  if polytype =='h'
        [rv,w]  = unif_ball_sample(M,sqrt(2)*sqrt(2*order+1),dim);
        sample.rv = rv;
        sample.w = w;
        for isim = 1:M
            psi(isim,:) = piset_hermite(sample.rv(isim,:),index_pc); 
        end
        sample.lhs = psi;
  end
end







end

