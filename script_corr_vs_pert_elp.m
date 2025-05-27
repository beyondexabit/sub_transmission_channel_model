%  close all;
clear all;
clc
addpath("sub_mmf_fibre_channel/")
addpath("sub_common_admin_func/")
%%
%% Fibre design parameters1e-3
dMax          = 0.03;       % maximum displacement - this is more like a surface parameter
surface_type  = '2d_coarse';   % 'fixRho' '2d_fine' '2d_coarse'
delta_core = 0.01;  % relative refractive index diference at the fibre core centre

%% fibre span parameters
Len        = 10; %[10 20 30 50 70 100 200 500 1e3 2e3 4e3 7e3 10e3]; %1e3; % meters
dz       = 10;  % meters
crosstalk= -40; % -40;  % average fibre XT dB/km
method = 'expm';     %  'sa' semi-analytical, 'expm' exp matrix sol.,...

lambda0   = 1550e-9;
f0        = physconst('LightSpeed')/lambda0;

Modes   = [6]; % 10 15 21]; % 15 21]; %[6 10 15 21]; % 10;
nPols   = 2; % 2;


xt_metric        = 'fixed2_pw'; %'maxk4';

coup_step_method = 'surf';
drhomethod       = 'mean';
xtmetric_SVD = xt_metric; %'deg';
statmethod='median';

theta_modes = 'different';
rho =[1e-5:1e-5:1e-4 1e-4:5e-5:2e-3 2e-3:5e-4:2e-2 2e-2:5e-3:1e-1];
num_phi = 10000;
num_theta = 1;

kdmethod = 'man_min'; % 0mat_fminsearch
kd_avg_stat='mean'; %mean
naverage_kd=10000; %100;
polrotmethod='add_oldkuv';
%%
elp_range=[0 1e-3 1.1e-2]; 
for nM=1:length(elp_range)
    for xx=1:length(Len)
        L = Len(xx); %%kdaverage
        nModes=Modes;
        elp=elp_range(nM);
        if nModes ==3
            ndeg=2;
        elseif nModes ==6
            ndeg=4;
        elseif nModes==10
            ndeg=6;
        elseif nModes==15
            ndeg=9;
        elseif nModes==21
            ndeg=12;
        else
            printf('specify non-degenerate modes')
        end

        lzMax = 1;
        %% Load files for N modes
         [~,~,~,~,elp,~,~,~,~,~,~,~,~,~,~,~,pc_elp_x,pc_elp_y,KuvSurf] = mmf_elliptical_fibre_characteristics(nModes,nPols,elp,delta_core,dMax,surface_type);
        pc=[];

        if nPols==2
            for i=1:length(pc_elp_x)
                pc=[pc pc_elp_x(i) pc_elp_y(i)];
            end
        else
            pc=pc_elp_x;
        end

        [rho1,XT_rho] = disp_nModes_xt_CMT2d_DEN_step_elp(rho,num_phi,dz,method,xt_metric,nModes,nPols,pc,KuvSurf,theta_modes,num_theta,statmethod,elp,polrotmethod);
    [drho] = find_drho(rho1,XT_rho,crosstalk,dz,drhomethod);
        [~,~,corr_kuv(:,xx),corr_beta(:,xx),kdrange_kuv(:,xx),kdrange_beta(:,xx)] = find_kd_control_elp_manual(L,dz,nModes,nPols,naverage_kd,coup_step_method,drho,theta_modes,xt_metric,kd_avg_stat,kdmethod,1,0,elp,delta_core,dMax,surface_type,polrotmethod);
    end
    kdrange_kuv_nM(:,nM)=kdrange_kuv;
    corr_kuv_nM(:,nM)=corr_kuv;
    kdrange_beta_nM(:,nM)=kdrange_beta;
    corr_beta_nM(:,nM)=corr_beta;
end

figpathsave='M:\Matlab_saved_figures\ECOC_presentation_Figures';

fs=10;
figure1=figure();
ax2 = axes();
hold(ax2);
box(ax2,'on');

plot(kdrange_kuv_nM(1:3:end,:),corr_kuv_nM(1:3:end,:),'x--','LineWidth',1,'HandleVisibility','off'); hold on
hold on;
set(gca,'ColorOrderIndex',1)
plot(kdrange_beta_nM(1:3:end,:),corr_beta_nM(1:3:end,:),'o-','LineWidth',1,'HandleVisibility','off');
set(gca,'ColorOrderIndex',1)
plot([NaN*ones(4,4)],'-','LineWidth',1); hold on
set(gca,'ColorOrderIndex',1)

plot(kdrange_beta_nM(1:3:end,1),exp(-kdrange_beta_nM(1:3:end,1)/exp(1)),':k','LineWidth',1.5);
plot(kdrange_beta_nM(1:3:end,1),1/exp(2)*ones(length(kdrange_beta_nM(1:3:end,1)),1),'--k','LineWidth',1.5,'HandleVisibility','off');


ylabel('Correlation','FontName','Times New Roman','FontSize',fs);
xlabel('Perturbation Variance - {\itx}','FontName','Times New Roman','FontSize',fs);
lgd = legend(num2str(repmat([Modes],1,1).'),'NumColumns',1,'Location','NorthWest');
lgd.Title.String = '{\itM}-modes';
set (lgd,'FontName','Times New Roman','FontSize',fs)
axis([0 max(max(kdrange_kuv)) -0.05 1.02])
grid on;
set(gca,'FontSize',fs,'FontName','Times New Roman','XScale','linear','YScale','linear');


ax = axes('Visible', 'off');
hold(ax);
plot(kdrange_kuv_nM(1:3:end-150,1),10000*corr_kuv_nM(1:3:end-150,1),'x--','LineWidth',1,'Color','k','HandleVisibility', 'on');
plot(kdrange_beta_nM(1:3:end-150,1),10000*corr_beta_nM(1:3:end-150,1),'o-','LineWidth',1,'Color','k','HandleVisibility', 'on');grid on;
lgd2 = legend((repmat([{'HoP-drift','\Delta{\it\beta}-drift'}],1,1).'),'NumColumns',1,'Location','best','interpreter','tex');
axis([0 max(max(kdrange_kuv)) -0.05 1.02])
lgd2.Title.String =  'Drift models';
set (lgd2,'FontName','Times New Roman','FontSize',fs)
set(gca,'FontSize',fs,'FontName','Times New Roman','XScale','linear','YScale','linear');


% set (gca,'FontName','Times New Roman','FontSize',fs)
set(gca,'XMinorGrid','off')
annotation(figure1,'textbox', 'String',{'exp({\it-x/e})'},...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
annotation(figure1,'textbox', 'String',{'{\ite}^{-2}'},...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Times New Roman',...
    'FitBoxToText','off');
grid on;



set(gcf, 'Position', [0, -100, 480, 400]);
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode', 'manual');
set(gcf,'PaperSize', [480, 400]/100);
set(gcf,'PaperPosition',[0, 0, 480, 400]/100)

saveas(gcf,fullfile(figpathsave,'Fig2_corr_vs_pert'),'fig')

return
