% close all;
% clear all;
clc

addpath("sub_mmf_fibre_channel\")
addpath("sub_common_admin_func\")


%% Fibre design parameters
dMax          = 0.03;       % maximum displacement - this is more like a surface parameter
surface_type  = '2d_coarse';   % 'fixRho' '2d_fine' '2d_coarse'
delta_core = 0.01;  % relative refractive index diference at the fibre core centre

%% fibre span parameters
L        = 1e3; %1e3; % meters
dz       = 10; %1;  % meters
crosstalk1 = -40; %  -30  -40;  % average fibre XT dB/km
method = 'expm';     %  'sa' semi-analytical, 'expm' exp matrix sol.,...

lambda0   = 1550e-9;
f0        = physconst('LightSpeed')/lambda0;

Modes   = 6; %[6 45 78 105 136];
nPols   = 2; % 2;


xt_metric        ='fixed2_pw'; % 'max_1';% 'fixed2_pw'; %'only_max';

coup_step_method = 'surf';
drhomethod       = 'mean';

naverageSVD      = 100; %100
dtRtenv  =10; %[1e-16 1e-12 1e-9 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1 10 50 100];

xtmetric_SVD ='intramodegroup'; %xt_metric; %'intramodegroup'; % xt_metric; %'intramodegroup'; % xt_metric; %'deg';
statmethod='median';

theta_modes = 'different';
rho =[1e-5:1e-5:1e-4 1e-4:5e-5:2e-3 2e-3:5e-4:2e-2 2e-2:5e-3:1e-1]; %[1e-4:5e-5:2e-3 1e-3:1e-6:2e-3 2e-3:5e-4:2e-2 1e-2:5e-3:1e-1];
rho = sort(unique(rho));

num_phi = 10000; %50;
num_theta = 1; %50;

elp_range=[0 1e-4:1e-4:15e-4 16e-4:5e-4:30e-4 35e-4:10e-4:100e-4 2e-2:5e-3:0.2]; %0.2]; % 0.1]; %[0 1e-3 0.1]; % [0:1e-4:15e-4 16e-4:5e-4:30e-4 35e-4:10e-4:100e-4]; %[0:1e-4:1.3e-3 15e-4:5e-4:1e-2]; %[0 6e-4]; %[-0.1 0 0.1]; %[-0.1 -0.05 -0.01 0 0.01 0.05 0.1];

kdmethod = 'mat_fminsearch_betterguess'; %'mat_fminsearch'; %'man_min'; %
kd_avg_stat='mean';
kd_error_tol=1e-2;
naverage_kd=1000;
polrotmethod='add_oldkuv';
%%
for jj=1:length(elp_range)
    elp=elp_range(jj);
    crosstalk =crosstalk1; % -40;  % average fibre XT dB/km)
    %         XT_res_diff_avg_modes = zeros(length(dtRtenv),length(Modes));
    for nM=1:length(Modes)
        nModes=Modes(nM);

        %%

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
        XT_mean=pow2db(mean(db2pow(XT_rho),2));
        % figure()
        % semilogx(rho1,XT_rho)
        % grid on;
        % xlabel('Radial Displacement')
        % ylabel('XT (dB)')

        [drho] = find_drho(rho1,XT_rho,crosstalk,dz,drhomethod);
        drho_elp(jj)=drho;

        [kd_modified_new_dz10,kd_modified_new_mod_drift_dz10] = find_kd_control_elp(1*dz,dz,nModes,nPols,naverage_kd,coup_step_method,drho,theta_modes,xt_metric,kd_avg_stat,kdmethod,kd_error_tol,1,0,elp,delta_core,dMax,surface_type,polrotmethod);
        kd_modified_new=kd_modified_new_dz10/(L/10);
        kd_modified_new_mod_drift=kd_modified_new_mod_drift_dz10/(L/10);

        kd_hop(jj)=kd_modified_new_dz10;
        kd_modified_drift(jj)=kd_modified_new_mod_drift_dz10;

    end

end




for i=1:length(elp_range)
    e=elp_range(i);
    n=1+e - (1/(1+e));
    d=1+e + (1/(1+e));
    oval(i)=n/d/2 *100;
end

fs=10;
figure()
% ax2 = axes();
% hold(ax2);
% box(ax2,'on');
plot(oval,kd_hop,'x--','LineWidth',1,'HandleVisibility', 'on')
% set(gca,'ColorOrderIndex',1)
hold on
plot(oval,kd_modified_drift,'o-','LineWidth',1,'HandleVisibility', 'on')
% set(gca,'ColorOrderIndex',1)
% plot([NaN*ones(3,3)],'-','LineWidth',1);
% set(gca,'ColorOrderIndex',1)
ylabel('Kd','FontName','Times New Roman','FontSize',fs);
xlabel('Ovality (%)','FontName','Times New Roman','FontSize',fs);
lgd = legend(['HoP','\Delta\beta'],'NumColumns',2,'Position',[0.722222222222222 0.479583339293797 0.135416665797432 0.223749994039536]);
lgd.Title.String = ['Drift Model']; %, sprintf('\n'), '{\itXT}=-40 dB/km', '    {\itXT}=-30 dB/km'];
set (lgd,'FontName','Times New Roman','FontSize',fs)
set(gca,'FontSize',fs,'FontName','Times New Roman','XScale','linear','YScale','linear');
grid on;
set(gca,'XMinorGrid','off')
set(gcf, 'Position', [0, -100, 480, 400]);
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode', 'manual');
set(gcf,'PaperSize', [480, 400]/100);
set(gcf,'PaperPosition',[0, 0, 480, 400]/100)

fs=10;
figure()
% ax2 = axes();
% hold(ax2);
% box(ax2,'on');
plot(oval,drho_elp,'x--','LineWidth',1,'HandleVisibility', 'on')

ylabel('drho','FontName','Times New Roman','FontSize',fs);
xlabel('Ovality (%)','FontName','Times New Roman','FontSize',fs);
set(gca,'FontSize',fs,'FontName','Times New Roman','XScale','linear','YScale','linear');
grid on;
set(gca,'XMinorGrid','off')
set(gcf, 'Position', [0, -100, 480, 400]);
set(gcf,'PaperUnits','inches')
set(gcf,'PaperPositionMode', 'manual');
set(gcf,'PaperSize', [480, 400]/100);
set(gcf,'PaperPosition',[0, 0, 480, 400]/100)