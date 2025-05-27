function [CoupMatLPabxy] = calculateCoupMat_2pol(P,f0,Lspan,dz,XTavg,method,seed)

s = RandStream('swb2712','Seed',seed); RandStream.setGlobalStream(s);
nsteps = ceil(Lspan/dz);
nModes = P.tnom;
nPols  = 2;

if contains(method,'unitaryLc')
    Lc = 10^(-XTavg/10)*1000;
    zc = round(Lc/dz)*dz;
    
    for z = 1:nsteps
        polRot1 = eye(nPols*nModes); polRot2 = eye(nPols*nModes); polRot3 = eye(nPols*nModes);
        polRot4 = eye(nPols*nModes); polRot5 = eye(nPols*nModes); polRot6 = eye(nPols*nModes);
        [polRot1( 1:1: 2, 1:1: 2),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot2( 3:1: 4, 3:1: 4),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot3( 5:1: 6, 5:1: 6),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot4( 7:1: 8, 7:1: 8),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot5( 9:1: 10,9:1:10),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot6(11:1:12,11:1:12),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        polRot = polRot1*polRot2*polRot3*polRot4*polRot5*polRot6;
        
        CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = polRot;
        
        realZ = z*dz;
        if rem(realZ,zc) == 0
            auxCoup = sqrt(1/2) * (randn(nPols*nModes,nPols*nModes) + 1i*randn(nPols*nModes,nPols*nModes));
            [Q,R] = qr(auxCoup);
            CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = Q;
            disp('.')
        end
    end
elseif contains(method,'justPolRot')
    for z = 1:nsteps
        polRot1 = eye(nPols*nModes); polRot2 = eye(nPols*nModes); polRot3 = eye(nPols*nModes);
        polRot4 = eye(nPols*nModes); polRot5 = eye(nPols*nModes); polRot6 = eye(nPols*nModes);
        [polRot1( 1:1: 2, 1:1: 2),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot2( 3:1: 4, 3:1: 4),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot3( 5:1: 6, 5:1: 6),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot4( 7:1: 8, 7:1: 8),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot5( 9:1: 10,9:1:10),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        [polRot6(11:1:12,11:1:12),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
        polRot = polRot1*polRot2*polRot3*polRot4*polRot5*polRot6;
        
        CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = polRot;
    end
    return;
elseif contains(method,'fullCoupAll') || XTavg == Inf
    for z = 1:nsteps
        auxCoup = sqrt(1/2) * (randn(nPols*nModes,nPols*nModes) + 1i*randn(nPols*nModes,nPols*nModes));
        [Q,~] = qr(auxCoup);
        CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = Q;
    end
elseif contains(method,'noCoup') || XTavg == -Inf
    for z = 1:nsteps
        CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = eye(nPols*nModes);
    end
else
    if XTavg ~= 0
    if contains(method,'expm')
        [CoupMatLPab       ] = mmf_CoupMat_expm(P,f0,Lspan,dz,XTavg,method);
    else
        error('not done yet')
        [CoupMatLPab       ] = fmf6CoupMat_DEN_MULT(f0,Lspan,dz,XTavg);
    end
    
    
    
    for p1 = 1:length(CoupMatLPab)
        d1(p1) = sum(abs(CoupMatLPab(p1).c(:)).^2);
    end
    if max(d1) >length(CoupMatLPab(1).c)+.1
        warning('fibre matrix is not be unitary')
    end
    end
    
    CoupMatLPabxyTemp = (zeros(nPols*nModes));
    for z = 1:nsteps
        fprintf('.')
        %% Modal coupling LP??ab,xy
        polRot = (eye(nPols*nModes)); 
        for k1 = 1:nModes
%             polRot1 = gpuArray(eye(nPols*nModes)); 
            [polRot( 1+2*(k1-1):1: 2+2*(k1-1), 1+2*(k1-1):1: 2+2*(k1-1)),~] = qr(sqrt(1/2) * (randn(nPols,nPols) + 1i*randn(nPols,nPols)));
%             polRot = polRot1*polRot;
%             polRot(1:4,1:4)
        end
        
        if XTavg == 0
            Q = eye(nModes);
        else
            [Q]             = (CoupMatLPab(z).c(1:1:nModes,1:1:nModes));
        end
        
        CoupMatLPabxyTemp(1:2:nPols*nModes,1:2:nPols*nModes) = Q;%CoupMatLPab(z).c(1:1:nModes,1:1:nModes);
        CoupMatLPabxyTemp(2:2:nPols*nModes,2:2:nPols*nModes) = Q;%CoupMatLPab(z).c(1:1:nModes,1:1:nModes);
        
        CoupMatLPabxy(z).c(1:1:nPols*nModes,1:1:nPols*nModes) = CoupMatLPabxyTemp*polRot;
    end
    fprintf('\n')
end