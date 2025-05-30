function [SNReq]=transmission_fibre_flexible_v2_test(navg_txn,drho,F,nModes,nPols,FiberParameters,nlInd,minStep,maxStep,delta_tol,fineStep,Sin,P,M,DSP,NOverlap,elp,deltacore,mgscase,txtype,Ntrib,best_modes_loc,CoupMatLPabVect_nodrift,modSigAlt1)

fileName1  = [num2str(P.OSNRdB),'_',num2str(navg_txn),'_',F.xt_metric,'_xt=',num2str(drho),'_',num2str(F.L),'_',num2str(F.dz),'_',num2str(nModes),'_',num2str(nPols),'_',num2str(M.symbolrate),num2str(F.XTavg),num2str(NOverlap),num2str(elp),num2str(deltacore),mgscase,'_',txtype,'_',num2str(Ntrib)];
tag1       = 'find_snr_lesstrib_leastinterferes_v1\';

seed=rng(navg_txn);
SNReq=[];

SNReq = loadComputedData(tag1,['SNReq_',fileName1],0);

if isempty(SNReq)
    SNReq.Ntrib=[];
end
if ~isempty(SNReq.Ntrib)
    return
end


SNReq(length(Ntrib)).Ntrib = [];

for in= 1:length(Ntrib)
    Ntrib(in)
    clearvars H_coupMat2 H_coupMat
    P.NTrib=Ntrib(in); %% set this correctly before noise addition
    Sin.E=zeros(length(M.txSig),nModes*nPols);

    chosen_modes= [best_modes_loc(1:Ntrib(in))]; % [worst_modes_loc(1:Ntrib)]; %
    chosen_modes=unique(sort(chosen_modes));
    Sin.E(:,chosen_modes)=M.txSig(:,chosen_modes);
    Sin.E= Sin.E.';

    %% Initial Transmission to calculate H without drift
    [Sout] = FMF_transmission_NModes_2pol('stochastic',FiberParameters,CoupMatLPabVect_nodrift,F.L,F.dz,nlInd,minStep,maxStep,delta_tol,fineStep,nModes,nPols,Sin);
    Sout.E = NoiseInsertion_samelevel_after_coupling(Sout.E,P.OSNRdB,P);
    disp_param=Sout.disp_param;

    if DSP.disp_comp==1
        H_coupMat1 =disp_comp_simulation(disp_param,Sout.H); %Txmatrix; %
    else
        H_coupMat1 =Sout.H;
    end

    rxSig =Sout.E.';
    %%
    % Overlap-save parameters:
    DFTsize  = size(rxSig,1); %2^6;
    P.GPU=0;

    x=rxSig(:,chosen_modes);
    modSigAlt=modSigAlt1(:,chosen_modes);

    H_coupMat2=pagetranspose(H_coupMat1(:,chosen_modes,:));
    H_coupMat=pagetranspose(H_coupMat2(:,chosen_modes,:));

    Westim = EqualiserCalc(H_coupMat,inf,P);
    W = zeros(size(Westim, 3), P.NTrib, P.NTrib);
    for i = 1:P.NTrib
        for j = 1:P.NTrib
            W(:,i,j) = squeeze(Westim(i,j,:));
        end
    end
    rxSig_eq = zeros(DFTsize, P.NTrib);
    EqType = 'MISO';
    for mm=1:P.NTrib  %%MISO Size
        clearvars EqRequired
        for nM = 1:P.NTrib
            PoutX          = abs(H_coupMat(:,nM)).^2;
            [~,EqRequired] =maxk(PoutX,mm,"ComparisonMethod",'abs');
            yMode  = Equalisation(x(:,EqRequired),W(:,nM,EqRequired),DFTsize,NOverlap,EqType,P);

            rxSig_eq(:,nM) = yMode;
        end
        %% snr calculation after equalization
        rxSig_equalized    = bsxfun(@rdivide,rxSig_eq,sqrt(mean(abs(rxSig_eq).^2)));
        SNReq(in).Ntrib(:,mm)=SNREstimation(modSigAlt,rxSig_equalized,P);
    end
end
saveComputedData(SNReq,tag1,['SNReq_',fileName1],0);
end