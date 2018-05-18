%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Paul Diaz 
% Date: 5/1/2018 
% Email: pdiaz303@gmail.com or paul.diaz@colorado.edu
% Please cite our paper:
% Diaz P., Doostan A., & Hampton J. "Sparse polynomial chaos expansions via
% compressed sensing and D-optimal design." Computer Methods in Applied
% Mechanics and Engineering. Volume 336, 1 July 2018, Pages 640-666
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c_hat ] = DSP( Psi_c, K, N, u_c, N_0,weight)
% Modified Subspace Pursuits
% Inputs: Psi_c  N_cxP candidate measurement matrix
%       : K = approximate bound on signal sparsity such that K >= s.
%         set K = 0 if bound is unknown to cross validate K.
%       : N = max # of measurements to be used in reconstruction
%       : u_c is an N_c x 1 vector of measurements (QOI samples) 
%       : weight is either true or false to normalize the support set est.
% Output: c Px1 vector of PCE coefficients
%INITIALIZATION:  
%Construct inital (2*K)xP design matrix and corresponding rows of data
P = size(Psi_c,2); max_iter = size(Psi_c,2);
design = rrqr_dopt_adapt(Psi_c, max(N_0,2*K),[]);

Psi = Psi_c(design,:); u = u_c(design);
if nargin < 6 || weight == false
   W = eye(size(Psi,2));
   weight = false;
else
   W = get_matrix_weights(Psi); 
   Psi = Psi*W;
end

c_hat = zeros(P,1);
stop = 0; num_iter=0; u_stop = Inf;

if K == 0
    [K,~] = cross_val_SP(Psi,u,weight);
    cross_val = true;
end
%Initial support estimation
x = zeros(P,1);
corr = abs(Psi'*u); [vals,~] = sort(corr,'descend');
I = find(corr >= vals(K));
%Initial residual calculation
x(I) = pinv(Psi(:,I))*u;  
ur0 = u - Psi*x;

while stop == 0 %ITERATION:   
   corr = abs(Psi'*ur0); [vals,~] = sort(corr,'descend');
   II = find(corr >= vals(K)); II = union(I,II);
   %LSP iteration.
   x = zeros(P,1); x(II) = pinv(Psi(:,II))*u;
   %Update support estimation   
   I0 = I; [vals,~] = sort(abs(x),'descend'); 
   I =  find( abs(x) >= vals(K));
   %D-optimal design update
   if length(u) < N
     design = rrqr_dopt_adapt(Psi_c(:,I),1,design);
     %design = alpha_optimal(Psi_c,design,I,1);
     if length(design) < length(u)+1
         fprintf('design falure in DSP\n')
         c_hat = 'failure';
         return;
     end
     Psi = Psi_c(design,:); u = u_c(design);
     if weight == true
        W = get_matrix_weights(Psi);
        Psi = Psi*W;
     end
     if cross_val == true 
        [K,~] = cross_val_SP(Psi,u,weight);
     end
   end  
   %Update residual 
   x = zeros(P,1); x(I) = pinv(Psi(:,I))*u; 
   ur = u - Psi*x;
   num_iter=num_iter+1;
   %Check stopping criteria
   if (norm(ur,2) >= norm(u_stop) && length(u) == N)
       I = I0;
       stop = 1;
   end
   if  num_iter == max_iter
       stop = 1;
   end
   ur0 = ur;
   if length(u) == N
       u_stop = ur0;
   end
   
end
c_hat(I) = pinv(Psi(:,I))*u;   
c_hat = W*c_hat;
end