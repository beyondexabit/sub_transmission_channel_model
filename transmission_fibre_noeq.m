function [SNR]=transmission_fibre_noeq(kd_modified_new,kd_modified_new_mod_drift,D,navg_txn,drho,F,nModes,nPols,KuvSurf,pc,FiberParameters,nlInd,minStep,maxStep,delta_tol,fineStep,Sin,Sin_pre,P,M,DSP)

fileName1  = [num2str(kd_modified_new),'_',num2str(kd_modified_new_mod_drift),num2str(D.dtRtenv),num2str(P.OSNRdB),'_',num2str(navg_txn),'_',F.xt_metric,'_xt=',num2str(drho),'_',num2str(F.L),'_',num2str(F.dz),'_',num2str(nModes),'_',num2str(nPols),'_',num2str(M.symbolrate),num2str(F.XTavg)];
tag1       = 'find_snr_afterSVD\';
SNR=[];

SNR = loadComputedData(tag1,['SNR_',fileName1]);

if isempty(SNR)
    SNR.md_drift=[];
end
if ~isempty(SNR.md_drift)
    return
end

navg_txn=1;
%%
k_d = kd_modified_new*D.dtRtenv;  %%
k_d_md = kd_modified_new_mod_drift*D.dtRtenv;
for kk=1
    fprintf([num2str(kk),' transmission out of ',num2str(navg_txn)]);
    %% get step-by-step coupling matrices
    [Coup]= calc_coupMat_after_tx(k_d,k_d_md,drho,F.L,F.dz,nModes,nPols,KuvSurf,pc,F.theta_modes);
    % Restructuring Coup Matrix Vector
    for k1 = 1:F.L/F.dz
        CoupMatLPabVect_nodrift(k1).c=Coup(k1).coupMat;
        CoupMatLPabVect_drift(k1).c=Coup(k1).coupMat_drift;
        CoupMatLPabVect_md_drift(k1).c=Coup(k1).coupMat_drift_modified;
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

    %% CaLculating SVD of transmission matrix without drift
    [U,S,V]=pagesvd(H_coupMat); % 6x6x512
    txSig_fft=(fftshift(fft(ifftshift(Sin.E,2),[],2),2));  %6x512;

    %% Transmitter - processed transmitted signal transmission over H with/without drift
    for i = 1:size(V,3)
        tx_Sig_preprocess(:,i) = V(:,:,i)*txSig_fft(:,i);
    end
    Sin_pre.E= fftshift(ifft( ifftshift( tx_Sig_preprocess,2), [], 2 ),2); %tx_Sig_preprocess

    Sout_nodrift = FMF_transmission_NModes_2pol('stochastic',FiberParameters,CoupMatLPabVect_nodrift,F.L,F.dz,nlInd,minStep,maxStep,delta_tol,fineStep,nModes,nPols,Sin_pre);
    disp_param=Sout_nodrift.disp_param;

    Sout_drift = FMF_transmission_NModes_2pol('stochastic',FiberParameters,CoupMatLPabVect_drift,F.L,F.dz,nlInd,minStep,maxStep,delta_tol,fineStep,nModes,nPols,Sin_pre);

    Sout_md_drift = FMF_transmission_NModes_2pol('stochastic',FiberParameters,CoupMatLPabVect_md_drift,F.L,F.dz,nlInd,minStep,maxStep,delta_tol,fineStep,nModes,nPols,Sin_pre);


    Sout_nodrift.E = NoiseInsertion(Sout_nodrift.E,P.OSNRdB,P);
    Sout_drift.E = NoiseInsertion(Sout_drift.E,P.OSNRdB,P);
    Sout_md_drift.E = NoiseInsertion(Sout_md_drift.E,P.OSNRdB,P);
    %% Rx svd processing and dispersion comp

    rx_Sig_preprocess = (fftshift(fft(ifftshift(Sout_nodrift.E,2),[],2),2));
    rx_Sig_preprocess_drift = (fftshift(fft(ifftshift(Sout_drift.E,2),[],2),2)); % Sout_drift.E;
    rx_Sig_preprocess_md_drift = (fftshift(fft(ifftshift(Sout_md_drift.E,2),[],2),2));

    for i = 1:size(U,3)
        rx_Sig_postprocess(:,i) = U(:,:,i)'*rx_Sig_preprocess(:,i);
        rx_Sig_postprocess_drift(:,i) = U(:,:,i)'*rx_Sig_preprocess_drift(:,i);
        rx_Sig_postprocess_md_drift(:,i) = U(:,:,i)'*rx_Sig_preprocess_md_drift(:,i);
    end

    if DSP.disp_comp==1
        rx_Sig_postprocess =disp_comp_simulation(disp_param,rx_Sig_postprocess);
        rx_Sig_postprocess_drift =disp_comp_simulation(disp_param,rx_Sig_postprocess_drift);
        rx_Sig_postprocess_md_drift =disp_comp_simulation(disp_param,rx_Sig_postprocess_md_drift);
    end

    Aout = fftshift(ifft( ifftshift( rx_Sig_postprocess,2), [], 2 ),2);
    Aout_drift = fftshift(ifft( ifftshift( rx_Sig_postprocess_drift,2), [], 2 ),2);
    Aout_md_drift = fftshift(ifft( ifftshift( rx_Sig_postprocess_md_drift,2), [], 2 ),2);

    if M.PS==1
        rxfilter= M.rxfilter;
        for k1 = 1:nModes*nPols
            rxSig(:,k1) = rxfilter(Aout(k1,:).');
            rxSig_drift(:,k1) = rxfilter(Aout_drift(k1,:).');
            rxSig_md_drift(:,k1) = rxfilter(Aout_md_drift(k1,:).');
        end
    else
        for k1 = 1:nModes*nPols
            rxSig(:,k1) = Aout(k1,:).';
            rxSig_drift(:,k1) = Aout_drift(k1,:).';
            rxSig_md_drift(:,k1) = Aout_md_drift(k1,:).';
        end
    end


    %% snr calculation before equalization
    rxSig          = bsxfun(@rdivide,rxSig,sqrt(mean(abs(rxSig).^2)));
    rxSig_drift    = bsxfun(@rdivide,rxSig_drift,sqrt(mean(abs(rxSig_drift).^2)));
    rxSig_md_drift = bsxfun(@rdivide,rxSig_md_drift,sqrt(mean(abs(rxSig_md_drift).^2)));
    P.NTrib=nModes*nPols;
    modSigAlt      = bsxfun(@rdivide,M.Tx_sig,sqrt(mean(abs(M.Tx_sig).^2)));%M.Tx_sig;
    SNR.nodrift(:,kk)=SNREstimation(modSigAlt,rxSig,P);
    SNR.khan_drift(:,kk)=SNREstimation(modSigAlt,rxSig_drift,P);
    SNR.md_drift(:,kk)=SNREstimation(modSigAlt,rxSig_md_drift,P);

  if mod(kk,10) == 0
  fileName1  = [num2str(kd_modified_new),'_',num2str(kd_modified_new_mod_drift),num2str(D.dtRtenv),num2str(P.OSNRdB),'_',num2str(kk),'_',F.xt_metric,'_xt=',num2str(drho),'_',num2str(F.L),'_',num2str(F.dz),'_',num2str(nModes),'_',num2str(nPols),'_',num2str(M.symbolrate),num2str(F.XTavg)];

        saveComputedData(SNR,tag1,['SNR_',fileName1]);
    end

end

saveComputedData(SNR,tag1,['SNR_',fileName1]);
% saveComputedData(SNReq,tag1,['SNReq_',fileName1]);