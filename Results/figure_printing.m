%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Paul Diaz 
% Date: 5/1/2018 
% Email: pdiaz303@gmail.com or paul.diaz@colorado.edu
% Please cite our paper:
% Diaz P., Doostan A., & Hampton J. "Sparse polynomial chaos expansions via
% compressed sensing and D-optimal design." Computer Methods in Applied
% Mechanics and Engineering. Volume 336, 1 July 2018, Pages 640-666
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script loads .mat files and generates figures
%% This is for manufacture sparse functions (relative error)
% load('noise_03_60s_20d_2p_240R_4620M_h_S')
% load('noise_03_60s_20d_2p_240R_4620M_h_C')
% load('noise_03_60s_2d_20p_240R_4620M_h_S')
 load('noise_03_60s_2d_20p_240R_4620M_h_C')

close all; 
figure
semilogy(sample_sizes,SP_stats(:,1),'*-',sample_sizes,GDSP_stats(:,1),'o-', sample_sizes, DSP_stats(:,1),'d-',sample_sizes,ORACLE_stats(:,1),'h-','markersize',12,'linewidth',3);
if sampling_type == 'C'
    L = legend('Coh-Opt','D-Coh-Opt','Seq-D-Coh-Op','Oracle','location','best')
else 
    L = legend('MC','D-MC','Seq-D-MC','Oracle','location','best')
end
grid on
set(L,'interpreter','latex','fontsize',18);
xlim([sample_sizes(1),sample_sizes(end)])
ylim([min([SP_stats(:,1); DSP_stats(:,1);GDSP_stats(:,1);ORACLE_stats(:,1)]),max([SP_stats(:,1); DSP_stats(:,1);GDSP_stats(:,1);ORACLE_stats(:,1)])])
xlabel('Sample Size $N$','interpreter','latex','fontsize',24)
ylabel('Rel. Error','interpreter','latex','fontsize',24)
title(['$(d,p) = $','(',num2str(dim),',', num2str(order),')'],'interpreter','latex','fontsize',24);
set(gca,'FontSize',20)

%%
print_name = ['MFS_rel_err_p',num2str(order),'_d',num2str(dim)]
%print_name = ['MFS_standard_rel_err_p',num2str(order),'_d',num2str(dim)]
print(print_name,'-depsc','-r300')

%% This is for manufactured sparse functions (standard dev of relative error)
close all
figure
semilogy(sample_sizes,SP_stats(:,2),'*-',sample_sizes,GDSP_stats(:,2),'o-', sample_sizes, DSP_stats(:,2),'d-','markersize',12,'linewidth',3);
if sampling_type == 'C'
    L = legend('Coh-Opt','D-Coh-Opt','Seq-D-Coh-Op','Oracle','location','best')
else 
    L = legend('MC','D-MC','Seq-D-MC','Oracle','location','best')
end
grid on
set(L,'interpreter','latex','fontsize',18);
xlim([sample_sizes(1),sample_sizes(end)])
ylim([min([SP_stats(:,2); DSP_stats(:,2);GDSP_stats(:,2)]),max([SP_stats(:,2); DSP_stats(:,2);GDSP_stats(:,2)])])
xlabel('Sample Size $N$','interpreter','latex','fontsize',24)
ylabel('Std of rel. error','interpreter','latex','fontsize',24)
title(['$(d,p) = $','(',num2str(dim),',', num2str(order),')'],'interpreter','latex','fontsize',24);
set(gca,'FontSize',20)

%%
print_name = ['MFS_std_p',num2str(order),'_d',num2str(dim)]
%print_name = ['MFS_standard_std_p',num2str(order),'_d',num2str(dim)]
print(print_name,'-depsc','-r300')

%% This is for manufacture sparse functions percent support estimated correctly
close all
figure
plot(sample_sizes,100*SP_stats(:,3),'*-',sample_sizes,100*GDSP_stats(:,3),'o-', sample_sizes, 100*DSP_stats(:,3),'d-','markersize',12,'linewidth',3);
if sampling_type == 'C'
    L = legend('Coh-Opt','D-Coh-Opt','Seq-D-Coh-Op','location','best')
else 
    L = legend('MC','D-MC','Seq-D-MC','location','best')
end
grid on
set(L,'interpreter','latex','fontsize',18);
xlim([sample_sizes(1),sample_sizes(end)])
ylim([100*min([SP_stats(:,3); DSP_stats(:,3);GDSP_stats(:,3)]),100])
xlabel('Sample Size $N$','interpreter','latex','fontsize',24)
ylabel('\% of $\mathcal{S}$ identified on average','interpreter','latex','fontsize',24)
title(['$(d,p) = $','(',num2str(dim),',', num2str(order),')'],'interpreter','latex','fontsize',24);
set(gca,'FontSize',20)

%%
print_name = ['MFS_supp_p',num2str(order),'_d',num2str(dim)]
%print_name = ['MFS_standard_supp_p',num2str(order),'_d',num2str(dim)]
print(print_name,'-depsc','-r300')

%%  This is for the duffing, ishigami, and wingweight (relative error)
%load('ishigami3d_7p_1000R.mat')
%load('ishigami3d_9p_1000R.mat')
load('ishigami3d_12p_1000R.mat')
%load('coh_opt_duffing3d_9p_1000R.mat')
%load('coh_opt_duffing3d_12p_1000R.mat')
%load('wingweight10d_3p_1000R.mat')
%close all;
figure;
semilogy(sample_sizes,SP_stats(:,1),'*-',sample_sizes,GSP_stats(:,1),'o-', sample_sizes, DSP_stats(:,1),'d-','markersize',12,'linewidth',3);
L = legend('Coh-Opt','D-Coh-Opt','Seq-D-Coh-Opt','location','best')
ylim([min([SP_stats(:,1); DSP_stats(:,1);GSP_stats(:,1)]),max([SP_stats(:,1); DSP_stats(:,1);GSP_stats(:,1)])]);
grid on;
set(L,'interpreter','latex','fontsize',18);
xlim([sample_sizes(1),sample_sizes(end)]);
%xlim([sample_sizes(1),303])
xlabel('Sample Size $N$','interpreter','latex','fontsize',24);
ylabel('Relative Validation Error','interpreter','latex','fontsize',24);
%title(['$(d,p) = $','(',num2str(dim),',', num2str(order),')'],'interpreter','latex','fontsize',24);
set(gca,'FontSize',20);
%%
save_name = ['R0_rel_err_p',num2str(order)]
print(save_name,'-depsc','-r300')



%%  This is for the duffing, ishigami, and wingweight (standard dev of relative  error)
close all
figure
semilogy(sample_sizes,SP_stats(:,2),'*-',sample_sizes,GSP_stats(:,2),'o-', sample_sizes, DSP_stats(:,2),'d-','markersize',12,'linewidth',3);
L = legend('Coh-Opt','D-Coh-Opt','Seq-D-Coh-Op','location','best')
grid on
xlim([sample_sizes(1),sample_sizes(end)])
ylim([min([SP_stats(:,2); DSP_stats(:,2);GSP_stats(:,2)]),max([SP_stats(:,2); DSP_stats(:,2);GSP_stats(:,2)])])
set(L,'interpreter','latex','fontsize',18);

%xlim([sample_sizes(1),303])
xlabel('Sample Size $N$','interpreter','latex','fontsize',24)
ylabel('Std of rel. error','interpreter','latex','fontsize',24)
%title(['$(d,p) = $','(',num2str(dim),',', num2str(order),')'],'interpreter','latex','fontsize',24);
set(gca,'FontSize',20);
%%
save_name = ['ww_std_p',num2str(order)]
print(save_name,'-depsc','-r300')


%% This creates histograms of the optimality criterion
close all;
figure
load('doe_results_250N_1000M_2d_20p_h_S')
histogram(design_scores,'normalization','cdf','LineWidth',1)
grid on
ylim([0,1])
hold on
load('doe_results_250N_1000M_2d_20p_h_C')
histogram(design_scores,7,'normalization','cdf','LineStyle','-.','LineWidth',1)
set(gca,'Xscale','log')
xlabel('$\vert \tilde{\bf{M}} \vert^{1/P}$','interpreter','latex','fontsize',24)
ylabel('Cumulative Probability','interpreter','latex','fontsize',24)
title(['$(d,p) = $','(',num2str(dimension_of_points),',', num2str(order),')'],'interpreter','latex','fontsize',24);
L = legend('MC','Coh-Opt','location','NorthWest')
set(L,'interpreter','latex','fontsize',18);
set(gca,'FontSize',20)

%%
save_name = ['doe_test_p',num2str(order),'_d',num2str(dimension_of_points ),'_',polytype]
print(save_name,'-depsc','-r300')

%% Initial Ishigami rel error
load('initial_ishigami3d_12p_1000R.mat')
figure;
hold on;
box on;
DSP_stats = Run_data{1};
the_min = min(DSP_stats(:,1));
the_max = max(DSP_stats(:,1));
linestyles = {'-','--',':','-.','+-'};

for m = 1:length(N0_percents)
    DSP_stats = Run_data{m};
    line_style = linestyles{m};
    plot(sample_sizes, DSP_stats(:,1),line_style,'markersize',12,'linewidth',3);
    if min(DSP_stats(:,1)) < the_min
        the_min = min(DSP_stats(:,1));
    end
    if max(DSP_stats(:,1)) > the_max
        the_max = max(DSP_stats(:,1));
    end
    
end
set(gca, 'YScale', 'log')
L = legend('0.5','0.6','0.7','0.8','0.9','location','best')
Ltitle = get(L,'Title');
set(Ltitle,'String','$\alpha$','interpreter','latex')
ylim([the_min,the_max]);
grid on;
set(L,'interpreter','latex','fontsize',18);
xlim([sample_sizes(1),sample_sizes(end)]);
%xlim([sample_sizes(1),303])
xlabel('Sample Size $N$','interpreter','latex','fontsize',24);
ylabel('Rel. Error','interpreter','latex','fontsize',24);
%title(['$(d,p) = $','(',num2str(dim),',', num2str(order),')'],'interpreter','latex','fontsize',24);
set(gca,'FontSize',20);
%%
save_name = ['initial_ishigami_rel_err_p',num2str(order)]
print(save_name,'-depsc','-r300')



%% Initial Ishigami Std of Rel Error
load('initial_ishigami3d_12p_1000R.mat')
close all;
figure;
hold on;
box on;

DSP_stats = Run_data{1};
the_min = min(DSP_stats(:,2));
the_max = max(DSP_stats(:,2));
linestyles = {'-','--',':','-.','+-'};

for m = 1:length(N0_percents)
    DSP_stats = Run_data{m};
    line_style = linestyles{m};
    semilogy(sample_sizes, DSP_stats(:,2),line_style,'markersize',12,'linewidth',3);
    if min(DSP_stats(:,2)) < the_min
        the_min = min(DSP_stats(:,2));
    end
    if max(DSP_stats(:,2)) > the_max
        the_max = max(DSP_stats(:,2));
    end
end
set(gca, 'YScale', 'log')
L = legend('0.5','0.6','0.7','0.8','0.9','location','best')
Ltitle = get(L,'Title');
set(Ltitle,'String','$\alpha$','interpreter','latex')
ylim([the_min,the_max]);
grid on;
set(L,'interpreter','latex','fontsize',18);
xlim([sample_sizes(1),sample_sizes(end)]);
%xlim([sample_sizes(1),303])
xlabel('Sample Size $N$','interpreter','latex','fontsize',24);
ylabel('Std of rel. error','interpreter','latex','fontsize',24);
%title(['$(d,p) = $','(',num2str(dim),',', num2str(order),')'],'interpreter','latex','fontsize',24);
set(gca,'FontSize',20);
%%
save_name = ['initial_ishigami_std_p',num2str(order)]
print(save_name,'-depsc','-r300')






