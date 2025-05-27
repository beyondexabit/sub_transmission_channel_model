function [XT]= calc_xt_after_svd(k_d,k_d_md,drho,L,dz,nModes,nPols,coup_step_method,xt_metric,KuvSurf,pc,theta_modes,rep,elp)


fileName1  = [num2str(theta_modes),'_',num2str(k_d),'_',num2str(k_d_md),'_',xt_metric,coup_step_method,'_xt=',num2str(drho,'%.12f'),'_',num2str(L),'_',num2str(dz),'_',num2str(nModes),'_',num2str(nPols),'_',num2str(rep),'_',num2str(elp)];
tag1       = 'findXT_afterSVD\nModesnPols\';

XT=[];
% XT = loadComputedData(tag1,['XT_',fileName1],0);
if isempty(XT)
XT.XT1d_res_drift_modified=[];
end
if ~isempty(XT.XT1d_res_drift_modified)
    return
end

% XT.XT1d_res_khan  = zeros(nModes*nPols,1);
% XT.XT1d_res_drift_modified = zeros(nModes*nPols,1);

lambda0  = 1550e-9;
w0       = 2*pi*physconst('LightSpeed')/lambda0;
B2 = pc/physconst('LightSpeed')*w0;

D= nPols*nModes;
nsteps = ceil(L/dz);
txSig=1;
navg=1;

    %% check accumulated XT
    coupMat = eye(D);
    coupMat_drift = eye(D);
    coupMat_drift_modified = eye(D);

    for kx = 1:nsteps

        dphi = 2*pi * rand(1);
        dx = drho .* cos(dphi);
        dy = drho .* sin(dphi);

        kuvX=zeros(nModes,nModes);
        for k1 = 1:nModes
            for k2 = 1:nModes
                kuvX(k1,k2) = KuvSurf(k1,k2).s(dx , dy);
            end
        end

        if nPols==2          
            kuvX_2D = zeros(nModes*nPols,nModes*nPols);
            kuvX_2D(1:2:end,1:2:end)=kuvX;
            kuvX_2D(1:2:end,2:2:end)=0;
            kuvX_2D(2:2:end,1:2:end)=0;
            kuvX_2D(2:2:end,2:2:end)=kuvX;
           
                [kuvX_modified]=polrot_Kuv_2D(kuvX_2D,nModes,theta_modes);
           

        else
            kuvX_modified=kuvX;
        end
        argumento = 1i*diag(B2)*dz +1i*kuvX_modified*dz;

        coupMat_itr(:,:) = expm(argumento);


        if sum(abs(coupMat_itr(:,:)).^2) > D+0.001
            error('not unitary');
        end

        %%
        rand_mat_real=randn(D);
        rand_mat_imag=randn(D);
        rand_skew_symm_mat_deltash(:,:)=triu(rand_mat_real.',1) - tril(rand_mat_real,-1)+1j*(triu(rand_mat_imag.') + tril(rand_mat_imag,-1));


        coupMat_itr_drift =  expm(argumento + (sqrt(k_d)*(rand_skew_symm_mat_deltash)));
        coupMat_itr_md_drift =  expm(argumento + (sqrt(k_d_md)*diag(diag(rand_skew_symm_mat_deltash))));

        coupMat = coupMat_itr*coupMat;
        coupMat_drift= coupMat_itr_drift*coupMat_drift;
        coupMat_drift_modified= coupMat_itr_md_drift*coupMat_drift_modified;


    end

    [U,~,V]=svd(coupMat);
  
    tx_Sig_preprocess=V*txSig;
    rx_Sig_preprocess= coupMat*tx_Sig_preprocess;
    rx_Sig_postprocess=U'*rx_Sig_preprocess;
    rx_Sig_preprocess_drift= coupMat_drift*tx_Sig_preprocess;
    rx_Sig_postprocess_drift=U'*rx_Sig_preprocess_drift;
    rx_Sig_preprocess_drift_modified= coupMat_drift_modified*tx_Sig_preprocess;
    rx_Sig_postprocess_drift_modified=U'*rx_Sig_preprocess_drift_modified;

%  [rx_Sig_postprocess_drift]=mode_seq(rx_Sig_postprocess_drift);
% [rx_Sig_postprocess_drift_modified]=mode_seq(rx_Sig_postprocess_drift_modified);

    [XT.XT1d_res_khan(:,navg),~]=calc_xt(rx_Sig_postprocess_drift,xt_metric,nModes,nPols,pc);
    [XT.XT1d_res_drift_modified(:,navg),~] = calc_xt(rx_Sig_postprocess_drift_modified,xt_metric,nModes,nPols,pc);

    [XT.XT1d_khan(:,navg),~] = calc_xt(coupMat_drift,xt_metric,nModes,nPols,pc);
    [XT.XT1d_drift_modified(:,navg),~] = calc_xt(coupMat_drift_modified,xt_metric,nModes,nPols,pc);

   saveComputedData(XT,tag1,['XT_',fileName1],0);

end