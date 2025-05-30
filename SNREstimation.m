function SNR = SNREstimation(x,y,P)

    for j = 1:P.NTrib

        % Transmitted and received symbols:
        xt = x(:,j);
        yr = y(:,j);

        % Constellation points:
        C1 = unique(xt(:,1));

        % SNR Calculation
        V        = zeros(numel(C1), 1);
        CentoidC = zeros(numel(C1), 1);
    
        % Noise variance:
        for i=1:numel(C1)
            pntX = find(xt(:,1)==C1(i));
    
            % Centroids:
            CentoidC(i,1) = mean(yr(pntX));
    
            xt(pntX,1) = (xt(pntX)-mean(xt(pntX)))+CentoidC(i);
    
            % Noise variance:
            V(i) = var(yr(pntX)-xt(pntX));
        end
    
        SNRdB = zeros(1,numel(C1));
        for i=1:numel(C1)
            SNRdB(1, i) = 10*log10(var(xt(:,1))/V(i)); % SNR of the i-th symbol XPOL (dB)
        end
        
        auxInfRemoval = find(SNRdB == inf);
        if numel(auxInfRemoval) ~= 0
            while numel(auxInfRemoval) ~= 0 
                SNRdB(auxInfRemoval(1)) = [];
                auxInfRemoval = find(SNRdB == inf);
            end
        end
    
        % Mean SNR (dB)
        SNRlin = mean(10.^(SNRdB./10));
    
        % SNR and mean SNR:
        SNR(j) = 10*log10(SNRlin);
    
    end

end