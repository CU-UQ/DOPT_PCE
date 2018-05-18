function [ k, TOL ] = cross_val_SP( Psi, u, weight )
%Cross-validation function via Subspace Pursuit
% Inputs:
% Psi , NxP measurement matrix 
% u , Nx1 QoI vector
% Tol , tolernace for SP function

R = 4; %Number of 'folds'
N = size(Psi,1);
P = size(Psi,2);
N_r = floor(0.8*N); %Numer of reconstruction samples
N_v = N - N_r; % Number of validation sample
N_k = 10; %Number of sparsities to attempt

K = round(linspace(1,floor(N/2),N_k)); %Sparsities to attempt

K = union(K,K);
N_k = length(K);
if nargin < 3 
    weight = false;
end
for i = 1:N_k
    k = K(i);
    for r = 1:R
        recon_inds = datasample(1:N,N_r,'Replace',false);
        val_inds = setdiff(1:N,recon_inds);
        Psi_r = Psi(recon_inds,:);
        Psi_v = Psi(val_inds,:);
        u_r = u(recon_inds);
        u_v = u(val_inds);
        c_hat = SP(Psi_r,k,u_r,weight);
        val_errs(r) = norm(Psi_v*c_hat - u_v)/norm(u_v);
    end
    SP_errs(i) = mean(val_errs);
end

[~,min_ind] = min(SP_errs);
k = K(min_ind);
TOL =  SP_errs(min_ind);






end

