function [kd_modified_new,kd_modified_new_mod_drift,corr_kuv,corr_beta,kdrange_kuv,kdrange_beta] = ...
    find_kd_control_elp_manual(L,dz,nModes,nPols,naverage,method,drho,theta_modes,xt_metric,statmethod,kdmethod,incbeta,plotFlag,elp,delta_core,dMax,surface_type,polrotmethod)


Nsec=L/dz;
kdrange_kuv= [0:0.05/Nsec:8/Nsec];  %[5:0.05/Nsec:15/Nsec];% [0:0.05/Nsec:2.5/Nsec];
kdrange_beta=[0:0.05/Nsec:8/Nsec]; %[5:0.05/Nsec:15/Nsec]; %[0:1/Nsec:10/Nsec];

%% try to load
switch kdmethod
    case 'man_min'
        fileName1  = [method,kdmethod,num2str(sum(kdrange_kuv)*length(kdrange_kuv),'%.12f'),num2str(sum(kdrange_beta)*length(kdrange_beta),'%.12f'),statmethod,xt_metric,num2str(theta_modes),'_xt=',num2str(drho),'_',num2str(L),'_',num2str(dz),'_',num2str(nModes),'_',num2str(nPols),'_',num2str(naverage),'_',num2str(nPols)];
       

    case 'mat_fminsearch'
        fileName1  = [method,kdmethod,statmethod,xt_metric,num2str(theta_modes),'_xt=',num2str(drho),'_',num2str(L),'_',num2str(dz),'_',num2str(nModes),'_',num2str(nPols),'_',num2str(naverage),'_',num2str(nPols)];
       
end



tag1       = ['findKd\NModes_NPols_elp4\manual_final\'];

kd_modified_new = loadComputedData(tag1,['kd_modified_new_',fileName1]);
kd_modified_new_mod_drift = loadComputedData(tag1,['kd_modified_new_mod_drift_',fileName1]);
corr_kuv = loadComputedData(tag1,['corr_kuv_',fileName1]);
corr_beta = loadComputedData(tag1,['corr_beta_',fileName1]);
kdrange_kuv1 =loadComputedData(tag1,['kdrange_kuv_',fileName1]);
kdrange_beta1 =loadComputedData(tag1,['kdrange_beta_',fileName1]);


switch kdmethod
    case 'man_min'
        if ~isempty(kdrange_beta)
            if length(kdrange_kuv)==length(kdrange_kuv1) && length(kdrange_beta)==length(kdrange_beta1)
                return
            else
                corr_kuv=[];
                corr_beta=[];
                kd_modified_new=[];
                kd_modified_new_mod_drift=[];
            end
        end
    case 'mat_fminsearch'
        if ~isempty(kd_modified_new_mod_drift)
            return
        end
end

fprintf('Calculating Kd...\n')
tic
%%
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


%% Load files for 6 modes
if strcmp(method,'surf')
    [~,~,~,~,elp,~,~,~,~,~,~,~,~,~,~,~,pc_elp_x,pc_elp_y,KuvSurf] = mmf_elliptical_fibre_characteristics(nModes,nPols,elp,delta_core,dMax,surface_type);
    pc=[];

    if nPols==2
        for i=1:length(pc_elp_x)
            pc=[pc pc_elp_x(i) pc_elp_y(i)];
        end
    else
        pc=pc_elp_x;
    end

    if incbeta==0
        pc=pc.*0;
    end
else
    pc = [];
    KuvSurf = [];
    drho = [];
end

%% Kd minsearch
switch kdmethod
    case 'mat_fminsearch'
        f = @(x)abs(-1/exp(2)+real(find_corr_general_elp(x,KuvSurf,drho,pc,nModes,nPols,naverage,L,dz,method,theta_modes,statmethod,polrotmethod)));%correlationTest_khan_sections(x,naverage,nModes,nPols,rho,XT_rho,crosstalk,B,KuvSurf,dz,L/dz)));
        if plotFlag
            options = optimset('PlotFcns',@optimplotfval);
        else
            options = {};
        end
        [x,~,~,~] = fminsearch(f,0,options);
        kd_modified_new = x;

        f = @(x)abs(-1/exp(2)+real(find_corr_general_mod_drift_elp(x,KuvSurf,drho,pc,nModes,nPols,naverage,L,dz,method,theta_modes,statmethod,polrotmethod)));%correlationTest_khan_sections(x,naverage,nModes,nPols,rho,XT_rho,crosstalk,B,KuvSurf,dz,L/dz)));
        if plotFlag
            options = optimset('PlotFcns',@optimplotfval);
        else
            options = {};
        end
        [x,~,~,~] = fminsearch(f,0,options);
        kd_modified_new_mod_drift = x;
        corr_beta=[];
        corr_kuv=[];
        kdrange_kuv=[];
        kdrange_beta=[];

    case 'man_min'
        corr_kuv=zeros(length(kdrange_kuv),1);
        corr_beta=zeros(length(kdrange_beta),1);
        for h=1:length(kdrange_kuv)
            corr_kuv(h) = real(find_corr_general_elp(kdrange_kuv(h),KuvSurf,drho,pc,nModes,nPols,naverage,L,dz,method,theta_modes,statmethod,polrotmethod));%correlationTest_khan_sections(x,naverage,nModes,nPols,rho,XT_rho,crosstalk,B,KuvSurf,dz,L/dz)));
        end
        for h=1:length(kdrange_beta)
            corr_beta(h) =real(find_corr_general_mod_drift_elp(kdrange_beta(h),KuvSurf,drho,pc,nModes,nPols,naverage,L,dz,method,theta_modes,statmethod,polrotmethod));%correlationTest_khan_sections(x,naverage,nModes,nPols,rho,XT_rho,crosstalk,B,KuvSurf,dz,L/dz)));
        end
        [~,loc_kuv]= min(abs(-1/exp(2)+corr_kuv),[],1);
        kd_modified_new =kdrange_kuv(loc_kuv);
        [~,loc_beta] = min(abs(-1/exp(2)+corr_beta),[],1);
end
hhh = toc;
fprintf(['Finished calculating Kd, time elapsed ',num2str(hhh/60),'mins \n'])


saveComputedData(kdrange_kuv,tag1,['kdrange_kuv_', fileName1]);
saveComputedData(kdrange_beta,tag1,['kdrange_beta_', fileName1]);
saveComputedData(corr_kuv,tag1,['corr_kuv_', fileName1]);
saveComputedData(corr_beta,tag1,['corr_beta_', fileName1]);
saveComputedData(kd_modified_new_mod_drift,tag1,['kd_modified_new_mod_drift_', fileName1]);
saveComputedData(kd_modified_new,tag1,['kd_modified_new_', fileName1]);
end


