function [W] = get_matrix_weights(Phi)
%   Compute a matrix W for the linear system Phi*W*W^(-1)*c = v;
%   such that the columns of Phi*W have 2-norm equal to 1.

P = size(Phi,2);
W = zeros(size(Phi,2));

for p=1:P
    W(p,p) = 1/norm(Phi(:,p)); 
end

