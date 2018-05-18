%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Paul Diaz 
% Date: 5/1/2018 
% Email: pdiaz303@gmail.com or paul.diaz@colorado.edu
% Please cite our paper:
% Diaz P., Doostan A., & Hampton J. "Sparse polynomial chaos expansions via
% compressed sensing and D-optimal design." Computer Methods in Applied
% Mechanics and Engineering. Volume 336, 1 July 2018, Pages 640-666
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Testing  D-Optimal SP (DSP) vs. (SP)
% To run this script in the background use the linux command:
% matlab -nodesktop -nodisplay <initial_ishigami_test.m &> file.out & 
% This script uses Modified subspace pursuit to approximate solutions of
% the ishigami function for varying initial design sizes.

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

order = 12;
load_name = ['ishigami_recon_data_p',num2str(order)];
load(load_name)
load_name = ['ishigami_validation_data_p',num2str(order)];
load(load_name)
sample_sizes = round(linspace(20,130,10))
R = 1000
Num_samps = length(sample_sizes);
tic
N0_percents = [0.5:0.1:.9]
Run_data = cell(length(N0_percents),1);


for m = 1:length(N0_percents)
    alpha = N0_percents(m)
    DSP_stats = zeros(Num_samps,2);
    dsp_err = zeros(R,Num_samps);
    for n = 1:Num_samps %Iterate to next sample size
      N = sample_sizes(n);
      parfor r = 1:R   % Repeat test R times
            data_inds = datasample(1:M,round(0.5*M),'replace',false);
            psi_cands = psi(data_inds,:); u_cands = U(data_inds);
            %Solve with modified subspace pursuit
            [c_hat] = DSP(psi_cands,0,N,u_cands,round(alpha*N),false);
            e2 = norm(psi_validation*c_hat - u_validation)/norm(u_validation);
            dsp_err(r,n) = e2;
            if mod(r,100)==0
                fprintf('r = %i  \n', r)
            end
      end  
       DSP_stats(n,1) =  mean(dsp_err(:,n));
       DSP_stats(n,2) =  std(dsp_err(:,n));
       Run_data{m} = DSP_stats;
       clear Psi 
       save_name = ['initial_ishigami',num2str(dim),'d_',num2str(order),'p_',num2str(R),'R'];   %_',num2str(K),'K'];
       %save(save_name)
       fprintf('Completed sample size N = %i with alpha = %f \n', N,alpha)
       disp(table(sample_sizes(1:n)',DSP_stats(1:n,1),DSP_stats(1:n,2),'VariableNames',{'N' 'Rel_err' 'Std_err'}))
       time = toc
    end
    fprintf('Completed alpha = %f \n', alpha)
    fprintf('Columns: 1 =  Avg relative error, 2 = Std Dev Rel Err, 3 = Avg Cond #, 4 = Std Dev Cond #\n')
    disp(table(sample_sizes',DSP_stats(:,1),DSP_stats(:,2),'VariableNames',{'N' 'Rel_err' 'Std_err'}))
    
    save_name = ['initial_ishigami',num2str(dim),'d_',num2str(order),'p_',num2str(R),'R'];   %_',num2str(K),'K'];
    %save(save_name)
end







