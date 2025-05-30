function [SNReq,SNR]=transmission_fibre_MGDM(navg_txn,drho,F,nModes,nPols,KuvSurf,pc,FiberParameters,nlInd,minStep,maxStep,delta_tol,fineStep,Sin,Sin_pre,P,M,DSP,NOverlap,elp,deltacore,mgscase)

if nargin<22
fileName1  = [num2str(P.OSNRdB),'_',num2str(navg_txn),'_',F.xt_metric,'_xt=',num2str(drho),'_',num2str(F.L),'_',num2str(F.dz),'_',num2str(nModes),'_',num2str(nPols),'_',num2str(M.symbolrate),num2str(F.XTavg),num2str(NOverlap)];
tag1       = 'find_snreq_afterSVD_single\';
else
 fileName1  = [num2str(P.OSNRdB),'_',num2str(navg_txn),'_',F.xt_metric,'_xt=',num2str(drho),'_',num2str(F.L),'_',num2str(F.dz),'_',num2str(nModes),'_',num2str(nPols),'_',num2str(M.symbolrate),num2str(F.XTavg),num2str(NOverlap),num2str(elp),num2str(deltacore),mgscase];
tag1       =  'find_snreq_afterSVD_single\';
end   
SNReq=[];


SNR = loadComputedData(tag1,['SNR_',fileName1],0);
SNReq = loadComputedData(tag1,['SNReq_',fileName1],0);

if isempty(SNReq)
    SNReq.nodrift=[];
end
if ~isempty(SNReq.nodrift)
    return
end



%% get step-by-step coupling matrices
[Coup]= calc_coupMat_after_tx(drho,F.L,F.dz,nModes,nPols,KuvSurf,pc,F.theta_modes,F.polrotmethod);
% Restructuring Coup Matrix Vector
for k1 = 1:F.L/F.dz
    CoupMatLPabVect_nodrift(k1).c=Coup(k1).coupMat;
end


%% Initial Transmission to calculate H without drift
[Sout] = FMF_transmission_NModes_2pol('stochastic',FiberParameters,CoupMatLPabVect_nodrift,F.L,F.dz,nlInd,minStep,maxStep,delta_tol,fineStep,nModes,nPols,Sin);

        Sout.E = NoiseInsertion(Sout.E,P.OSNRdB,P);


disp_param=Sout.disp_param;

if DSP.disp_comp==1
    H_coupMat =disp_comp_simulation(disp_param,Sout.H); %Txmatrix; %
else
    H_coupMat =Sout.H;
end

Aout =Sout.E; % fftshift(ifft( ifftshift( rx_Sig_postprocess,2), [], 2 ),2);

if M.PS==1
    rxfilter= M.rxfilter;
    for k1 = 1:nModes*nPols
        rxSig(:,k1) = rxfilter(Aout(k1,:).');
           end
else
    for k1 = 1:nModes*nPols
        rxSig(:,k1) = Aout(k1,:).';
           end
end


%% snr calculation before equalization
rxSig          = bsxfun(@rdivide,rxSig,sqrt(mean(abs(rxSig).^2)));
P.NTrib=nModes*nPols;
modSigAlt      = bsxfun(@rdivide,M.Tx_sig,sqrt(mean(abs(M.Tx_sig).^2)));%M.Tx_sig;
SNR.nodrift=SNREstimation(modSigAlt,rxSig,P);

%%
% Overlap-save parameters:
DFTsize  = size(rxSig,1); %2^6;
%     NOverlap = 64; %DFTsize/2; %  2^6; % % 1024;

x=rxSig;
P.GPU=0;
P.NTrib=nModes*nPols;


Westim = EqualiserCalc(H_coupMat,inf,P);

for i = 1:P.NTrib
    for j = 1:P.NTrib
        W(:,i,j) = squeeze(Westim(i,j,:));  
    end
end

%  Waux(:,Modes,Modes) = W;
% clearvars W
% W = Waux;

for mm=1:nModes*nPols  %%MISO Size
    clearvars P.EqRequired.v
    
    %             if contains(F.xt_metric,'maxk')
    kkkk =  mm; %str2num(F.xt_metric(5:end));
    for nM = 1:nModes*nPols
        PoutX          = abs(H_coupMat(:,nM)).^2;
        [~,P.EqRequired{nM}.v] =maxk(PoutX,kkkk,"ComparisonMethod",'abs');
    end

    for i = 1:nModes*nPols
        if numel(P.EqRequired{i}.v) > 0 %1
            EqType = 'MISO';
            yMode  = Equalisation(x(:,(P.EqRequired{i}.v)),W(:,i,P.EqRequired{i}.v),DFTsize,NOverlap,EqType,P);
           
            rxSig_eq(:,i) = yMode;
           
        end
    end

    %% snr calculation after equalization
    rxSig_equalized    = bsxfun(@rdivide,rxSig_eq,sqrt(mean(abs(rxSig).^2)));
    P.NTrib=nModes*nPols;
    SNReq.nodrift(:,mm)=SNREstimation(modSigAlt,rxSig_equalized,P);
end


saveComputedData(SNR,tag1,['SNR_',fileName1],0);
saveComputedData(SNReq,tag1,['SNReq_',fileName1],0);

end