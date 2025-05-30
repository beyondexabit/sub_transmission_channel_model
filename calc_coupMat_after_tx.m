function [Coup]= calc_coupMat_after_tx(drho,L,dz,nModes,nPols,KuvSurf,pc,theta_modes,polrotmethod,dz_fine)

if nargin<10
    dz_fine=dz;
end

lambda0  = 1550e-9;
w0       = 2*pi*physconst('LightSpeed')/lambda0;
B = pc/physconst('LightSpeed')*w0;

D= nPols*nModes;
nsteps = ceil(L/dz);

    %% check accumulated XT
   
  
    for kx = 1:nsteps

        dphi = 2*pi *rand(1);
        dx = drho .* cos(dphi);
        dy = drho .* sin(dphi);

        kuvX=zeros(nModes*nPols,nModes*nPols);
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

        argumento = (1i*diag(B) +1i*kuvX_modified)*dz_fine;

        coupMat_itr(:,:) = expm(argumento);


        if sum(abs(coupMat_itr(:,:)).^2) > D+0.001
            error('not unitary');
        end

    Coup(kx).coupMat=coupMat_itr;
    
    end
    
end