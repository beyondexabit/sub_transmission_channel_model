function M= transmittedsignal(M,nModes,nPols,Nsymb) 

k = log2(M.M); % Bits/symbol
n = k*M.Nsymb; %M.symbolrate;  % Transmitted bits
x      = randi([0 1],n,nModes*nPols); %*nPols
modSig = qammod(x,M.M,'InputType','bit');

if M.MSS==1 %% Modes same sequence but circularly shifted
    for k1=2:nModes
        modSig(:,k1)=circshift(modSig(:,k1-1),50);
    end
end


if M.PS==1
    M.txfilter  = comm.RaisedCosineTransmitFilter('RolloffFactor',M.rolloff, 'FilterSpanInSymbols',16,'OutputSamplesPerSymbol',M.nSamp);
    M.rxfilter  = comm.RaisedCosineReceiveFilter('RolloffFactor',M.rolloff, 'FilterSpanInSymbols',16,'InputSamplesPerSymbol',M.nSamp,'DecimationFactor',M.nSamp);
    filtDelay = k*M.M;

    if nPols==2
        for k1 = 1:nModes
            txSig(:,2*k1-1) = 6*1e-2*M.txfilter(modSig(:,k1)); %% ?? Why 6*1e-2 %modSig(:,k1); %
            txSig(:,2*k1  ) = 6*1e-2*M.txfilter(modSig(:,k1));  %% modSig(:,k1); Dual polarization
        end
    else
        txSig = modSig; %txfilter(modSig); % modSig; %;
    end
else
%     if nPols==2
%         for k1 = 1:nModes
%             txSig(:,2*k1-1) = modSig(:,k1);
%             txSig(:,2*k1  ) = modSig(:,k1);
%         end
%     else
        txSig = modSig; %txfilter(modSig); % modSig; %;
%     end
    filtDelay=0;
end

%  if nPols==2
%     M.Tx_sig(:,1:2:nModes*nPols)=modSig;
%     M.Tx_sig(:,2:2:nModes*nPols)=modSig;
% else
    M.Tx_sig=modSig;
% end

M.filtDelay=filtDelay;
M.txSig=txSig;
M.modSig=modSig;


end
