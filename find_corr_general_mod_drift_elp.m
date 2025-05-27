function [corr_kd_avg_modified_drift] = find_corr_general_mod_drift_elp(k_d,KuvSurf,drho,pc,nModes,nPols,naverage,L,dz,method,theta_modes,statmethod,polrotmethod)

lambda0  = 1550e-9;
f0       = physconst('LightSpeed')/lambda0;
w0     = 2*pi*f0;
B = pc/physconst('LightSpeed')*w0;
D=nModes*nPols;
nsteps = ceil(L/dz);

%%

corr_kd_modified_drift=zeros(naverage,1);
disp(['averaging kd... x=',num2str(k_d)])
parfor jj_var = 1:naverage
    if mod(jj_var,100) == 0
        fprintf([num2str(jj_var),' '])
    end

    coupMat_itr=eye(D);
    coupMat=eye(D);
    coupMat_modified_drift=eye(D);
    for kx = 1:nsteps
        rand_skew_symm_mat_deltash=zeros(D);
        switch method
            case 'surf'
                % Random Fiber Core Displacement
                dphi = 2*pi * rand(1);
                dx = drho .* cos(dphi);
                dy = drho .* sin(dphi);

                kuvX=zeros(nModes,nModes);
                for k1 = 1:nModes*nPols
                    for k2 = 1:nModes*nPols
                        kuvX(k1,k2) = KuvSurf(k1,k2).s(dx , dy);
                    end
                end

                if nPols==2

        [ kuvX_modified,~]=polrot_Kuv_2D(kuvX,nModes,theta_modes,polrotmethod);
   
                else
                    kuvX_modified=kuvX;
                end
                argumento = 1i*diag(B)*dz +1i*kuvX_modified*dz;
            case 'skewkhan'
                rand_mat_real=randn(D);
                rand_mat_imag=randn(D);
                rand_skew_symm_mat_deltash(:,:)=triu(rand_mat_real.',1) - tril(rand_mat_real,-1)+1j*(triu(rand_mat_imag.') + tril(rand_mat_imag,-1));

                argumento = rand_skew_symm_mat_deltash;
        end

        coupMat_itr(:,:) = expm(argumento);


        if sum(abs(coupMat_itr(:,:)).^2) > D+0.001
            error('not unitary');
        end

        %%
        rand_mat_real=randn(D);
        rand_mat_imag=randn(D);
        rand_skew_symm_mat_deltash(:,:)=triu(rand_mat_real.',1) - tril(rand_mat_real,-1)+1j*(triu(rand_mat_imag.') + tril(rand_mat_imag,-1));


        coupMat_itr_modified_drift =  expm(argumento + (sqrt(k_d)*diag(diag(rand_skew_symm_mat_deltash))));

        coupMat=coupMat_itr*coupMat;
        coupMat_modified_drift= coupMat_itr_modified_drift*coupMat_modified_drift;
    end
    corr_kd_modified_drift(jj_var)=1/D*trace(coupMat*coupMat_modified_drift');
end
fprintf(['\n'])

switch statmethod
    case 'median'
        corr_kd_avg_modified_drift=median(corr_kd_modified_drift,1);
    case 'mean'
        corr_kd_avg_modified_drift=mean(corr_kd_modified_drift,1);
end

end