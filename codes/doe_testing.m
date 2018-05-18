%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Paul Diaz 
% Date: 5/1/2018 
% Email: pdiaz303@gmail.com or paul.diaz@colorado.edu
% Please cite our paper:
% Diaz P., Doostan A., & Hampton J. "Sparse polynomial chaos expansions via
% compressed sensing and D-optimal design." Computer Methods in Applied
% Mechanics and Engineering. Volume 336, 1 July 2018, Pages 640-666
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code will not run as it utilizes the coherence optimal sampler which
% can be downloaded at www.github.com/CU-UQ

clear; clc;
addpath(genpath('./'))
delete(gcp('nocreate'))
pool_data = parpool;
rng('shuffle'); % Shuffles on local worker only
%Shuffle on each parallel worker
seed_offset = randi(floor(intmax/10),1);
parfor kk = 1:pool_data.NumWorkers
  rng(kk + seed_offset);
end

%%
number_of_design_points = 250
number_of_candidate_points = 1000
dimension_of_points = 2 % This is simply the dimension that the sample points exist in. Typically denoted by d. 
order = 20 % Total order of PCE p
polytype = 'h' % 'L' for Legendre and 'h' for Hermite
sampling_type = 'C' % 'C' for coherence optimal 'S' for MC 
number_of_design_samples = 10^3 

sample_data.M = number_of_candidate_points;
sample_data.dim = dimension_of_points;
sample_data.polytype = polytype;
sample_data.sampling_type = sampling_type;
sample_data.order = order;
sample_data.NumWorkers = pool_data.NumWorkers;

tic;
design_scores = zeros(number_of_design_samples,1);
for k = 1:number_of_design_samples
    [~, design_scores(k)] = design_select_algorithm(number_of_design_points, sample_data); % This function doesn't exist, but is drawing two samples from candidate distribution and computing the score of the associated design
    if mod(k,100) == 0
        fprintf('completed k = %i \n', k) 
        time = toc
    end
end


save_name = ['doe_results_',num2str(number_of_design_points),'N_',num2str(number_of_candidate_points),'M_',num2str(dimension_of_points),'d_',num2str(order),'p_', polytype,'_',sampling_type];   %_',num2str(K),'K'];
time = toc
save(save_name)
