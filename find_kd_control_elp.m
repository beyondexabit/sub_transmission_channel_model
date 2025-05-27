function [kd_modified_new,kd_modified_new_mod_drift] = ...
    find_kd_control_elp(L,dz,nModes,nPols,naverage,method,drho,theta_modes,xt_metric,statmethod,kdmethod1,errortol,incbeta,plotFlag,elp,delta_core,dMax,surface_type,polrotmethod)



%% try to load

fileName1  = [method,kdmethod1,num2str(errortol),statmethod,xt_metric,num2str(theta_modes),'_',num2str(drho),'_',num2str(L),'_',num2str(dz),'_',num2str(nModes),'_',num2str(nPols),'_',num2str(naverage),'_',num2str(nPols),num2str(incbeta),num2str(plotFlag),num2str(elp),num2str(delta_core),num2str(dMax),num2str(surface_type),polrotmethod];
tag1       = ['findKd\NModes_NPols_elp4\'];

kd_modified_new = loadComputedData(tag1,['kd_modified_new_',fileName1],0);
kd_modified_new_mod_drift = loadComputedData(tag1,['kd_modified_new_mod_drift_',fileName1],0);

fprintf(['kdmethod =' ,num2str(kdmethod1),'\n'])

if ~isempty(kd_modified_new_mod_drift)
    return
end

fprintf('Calculating Kd...\n')
tic

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
if strcmp(kdmethod1,'mat_fminsearch')
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

elseif contains(kdmethod1,'mat_fminsearch_betterguess')
    tolx = errortol;
    options = optimset('TolX',tolx,'TolFun',tolx);

    [kd_temp,kd_mod_temp] = find_kd_control_elp(10*dz,dz,nModes,nPols,naverage,method,drho,theta_modes,xt_metric,statmethod,'mat_fminsearch',0,incbeta,0,elp,delta_core,dMax,surface_type,polrotmethod);
    f = @(x)abs(-1/exp(2)+real(find_corr_general_elp(x,KuvSurf,drho,pc,nModes,nPols,naverage,L,dz,method,theta_modes,statmethod,polrotmethod)));%correlationTest_khan_sections(x,naverage,nModes,nPols,rho,XT_rho,crosstalk,B,KuvSurf,dz,L/dz)));
    if plotFlag
        options = optimset('PlotFcns',@optimplotfval);
    else
        options = {};
    end
    [x,~,~,~] = fminsearch(f,kd_temp*10/(L/dz),options);
    kd_modified_new = x;

    f = @(x)abs(-1/exp(2)+real(find_corr_general_mod_drift_elp(x,KuvSurf,drho,pc,nModes,nPols,naverage,L,dz,method,theta_modes,statmethod,polrotmethod)));%correlationTest_khan_sections(x,naverage,nModes,nPols,rho,XT_rho,crosstalk,B,KuvSurf,dz,L/dz)));
    
    if plotFlag
        options = optimset('PlotFcns',@optimplotfval);
    else
        options = {};
    end
    [x,~,~,~] = fminsearch(f,kd_mod_temp*10/(L/dz),options);
    kd_modified_new_mod_drift = x;

end
hhh = toc;
fprintf(['Finished calculating Kd, time elapsed ',num2str(hhh/60),'mins \n'])

saveComputedData(kd_modified_new_mod_drift,tag1,['kd_modified_new_mod_drift_', fileName1],0);
saveComputedData(kd_modified_new,tag1,['kd_modified_new_', fileName1],0);
end

