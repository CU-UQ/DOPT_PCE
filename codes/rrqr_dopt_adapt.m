%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Paul Diaz 
% Date: 5/1/2018 
% Email: pdiaz303@gmail.com or paul.diaz@colorado.edu
% Please cite our paper:
% Diaz P., Doostan A., & Hampton J. "Sparse polynomial chaos expansions via
% compressed sensing and D-optimal design." Computer Methods in Applied
% Mechanics and Engineering. Volume 336, 1 July 2018, Pages 640-666
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ design ] = rrqr_dopt_adapt( psi, n_rows, design)
%rrqr_dopt computes the best N-point D-optimal design for a wide matrix psi
% Inputs:
%       : psi, is an NxP matrix of candidate rows
%       : n_rows, is the number of rows to be added to the design
%       : design,  indices of current rows of psi which are in the design
%         If no design exists, set design = []
% Output: design, is an Nx1 vector whose indices correspond to the rows of  
%        psi which produce the optimal design.

if size(psi,1) < size(psi,2)
    design = datasample(1:size(psi,1),n_rows,'replace',false);
else
    if isempty(design) %Subset Selection
        [~,~,V] = svd(psi');
        [~,~,piv] = qr(V(:,1:size(psi,2))','vector');
        [~,~,pivot_vector] = rank_k_exchange(V(:,1:size(psi,2))',piv,size(psi,2),1);
        design = pivot_vector(1:n_rows);
    else %Design Adaptation via RRQR
        [Q,R] = qr(psi(design,:)');
        PP = R\(Q'*psi');
        [~,R,piv] = qr(psi'-psi(design,:)'*PP,'vector');
        rk = length(find(abs(diag(R)) >= 1e-12));
        [~,~,pivot_vector] = rank_k_exchange(psi'-psi(design,:)'*PP, piv,rk, 1);
        n_rows_added = 0;
        for n = 1:size(psi,1)
           if isempty(intersect(pivot_vector(n),design))
              design = [design, pivot_vector(n)];
              n_rows_added = n_rows_added+1;
           end
           if n_rows_added == n_rows
               return
           end
        end
    end
end



end

