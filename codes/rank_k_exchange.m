function [Q,R,piv] = rank_k_exchange(AA, piv, rk, thresh)
%
% This method exchanges pivots in QR alogirhtm to improve determinant of
% the resulting Ain R= [A B; 0 C]);
%
% inputs: (Via any QR factorization) Recommended: 2x modified GS
% AA - matrix to be approximated
% piv - vector keeping track of the pivots AA(:,piv) \approx Q*R
% rk - Number of Rows/Columns for which R is triangular (A is rk x rk)
% thresh >= 1 , exchange until exchanges with disc >= thresh is impossible

% outputs:
% Q - orthogonal matrix (updated)
% R - Upper Triangular matrix, A = QR (updated)
% piv - vector keeping track of the pivots (updated)

% Initialize decomposition/values we need
n_cols = size(AA,2); % Number of potential columns
n_rows = size(AA,1);
[Q,R] = qr(AA(:,piv));
% Now we swap iteratively

max_thresh = thresh+1; % Just to enter loop

while max_thresh > thresh

    % Recompute (or initialize) our A,B,C matrices
    A = R(1:rk,1:rk);                  % R = |A  B|
    B = R(1:rk,(rk+1):n_cols);         %     |0  C|
    C = R((rk+1):n_rows,(rk+1):n_cols);

    invA = inv(A); % Need this inverse
    invABsq = invA*B; %#ok<MINV> % Need this values
    invABsq = invABsq.^2; 

    A_inv_norm_sq = zeros(rk,1);
    for i = 1:rk
        A_inv_norm_sq(i) = norm(invA(i,:))^2;
    end
    C_norm_sq = zeros(rk,1);
    for i = 1:(n_cols-rk)
        C_norm_sq(i) = norm(C(:,i))^2;
    end
    
    % identify column i to swap with column rk+j
    max_thresh = 0;
    best_i = 0;
    best_j = 0;
    for i = 1:rk
        for j =1:(n_cols-rk)
            disc = sqrt(invABsq(i,j) + (A_inv_norm_sq(i)*C_norm_sq(j))); % Make exchange that maximizes this value
            if disc > max_thresh
                max_thresh = disc;
                best_i = i;
                best_j = j;
            end
        end
    end
    
    if max_thresh >= thresh
        i = best_i; % Swap column best_i
        j = best_j; % Swap column best_j
        
        % Swap column i for column rk+j in pivot
        temp = piv(rk+j);
        piv(rk+j) = piv(i);
        piv(i) = temp;
        
        % Now we recompute R, and update Q accordingly to new pivots
        % This can likely be done efficiently with givens rotations
        % Important: qr method should not pivot
        [Q,R] = qr(AA(:,piv));
        return
    end
 
end

end