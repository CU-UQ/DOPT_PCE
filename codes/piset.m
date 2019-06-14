% Evaluates a multi_D PC basis at a given point (xi_1,...,xi_d)

function pc_xi = piset(xi,index_pc)

pc_xi = ones(size(index_pc,1),1);

p = sum(index_pc(size(index_pc,1),:));
Legendre = my_legendre_1d_reg(p,xi);

for id=1:size(index_pc,2)
    nnz_index = find(index_pc(:,id)>0);
    if find(index_pc(:,id)>0)
        pc_xi(nnz_index) = pc_xi(nnz_index).*Legendre(index_pc(nnz_index,id)+1,id);
    end
end
