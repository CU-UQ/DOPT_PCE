function [ design, score ] = design_select_algorithm(number_of_design_points, sample_data )


sample = PCE_sample(sample_data);
Psi_c = diag(sample.w)*sample.lhs;
P = size(sample.lhs,2);


design = rrqr_dopt_adapt(Psi_c,number_of_design_points,[]);
if number_of_design_points < P
    Psi_temp = Psi_c(design,1:length(design));
    
else
    Psi_temp = Psi_c(design,:);
end

D = Psi_temp'*Psi_temp;
D = D/norm(D,'fro');
[~, E, ~] = svd(D,'econ');

log_det = sum(log(diag(E)));
score = exp(log_det*(1/P));

dims = size(sample.rv(design,:));
design = sample.rv(design,:);
design = reshape(design',1,dims(1)*dims(2));

end
