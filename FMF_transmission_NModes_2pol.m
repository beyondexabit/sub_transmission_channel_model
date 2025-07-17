function [varargout] = FMF_transmission_NModes_2pol(methodSSMF,FiberParameters,CoupMatLPab,L,dz,nlInd,minStep,maxStep,delta_tol,fineStep,nModes,nPols,varargin)
% REVISION
% created by Filipe Ferreira @ 2022

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

%% Source parameters
lambda0  = 1550e-9;
w0       = 2*pi*c0/lambda0;
f0       = w0/(2*pi);

% Calling from Matlab. Extract signals and create matrix Ain
Ain   = varargin{1}.E;


if nPols==2 && length(FiberParameters.lossCoef)==nModes
DMD(1:2:nModes*nPols)     = FiberParameters.ModeDelay-FiberParameters.dgd_pol/2;
DMD(2:2:nModes*nPols)     = FiberParameters.ModeDelay+FiberParameters.dgd_pol/2;
Disp(1:2:nModes*nPols) = FiberParameters.D; 
S(1:2:nModes*nPols)    = FiberParameters.S;
lossCoef(1:2:nModes*nPols)        = FiberParameters.lossCoef;
Disp(2:2:nModes*nPols) = FiberParameters.D; 
S(2:2:nModes*nPols)    = FiberParameters.S;
lossCoef(2:2:nModes*nPols)        = FiberParameters.lossCoef;
else
DMD(1:nModes*nPols)  = FiberParameters.ModeDelay;
Disp(1:nModes*nPols) = FiberParameters.D; 
S(1:nModes*nPols)    = FiberParameters.S;
lossCoef(1:nModes*nPols)        = FiberParameters.lossCoef;
end


% if nPols==2 && length(FiberParameters.lossCoef)==nModes
%     lossCoef(1:2:nModes*nPols)        = FiberParameters.lossCoef;
% lossCoef(2:2:nModes*nPols)        = FiberParameters.lossCoef;
% else
% lossCoef(1:nModes*nPols)        = FiberParameters.lossCoef;
% end

% Define frequency vector f
TimeWindow     	= varargin{1}.T*varargin{1}.dt;             % Duration of signal in seconds
Fs            	= 1/TimeWindow;                             % Frequency spacing
f             	= (-varargin{1}.T/2:varargin{1}.T/2-1).'*Fs;

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

%% Fiber Model General Parameters
nlCoefOrig = FiberParameters.nlCoef;

if nPols==2
nlCoefOrig2pols(1:2:nModes*nPols,1:2:nModes*nPols) = nlCoefOrig;
nlCoefOrig2pols(2:2:nModes*nPols,2:2:nModes*nPols) = nlCoefOrig;
nlCoefOrig2pols(1:2:nModes*nPols,2:2:nModes*nPols) = nlCoefOrig;
nlCoefOrig2pols(2:2:nModes*nPols,1:2:nModes*nPols) = nlCoefOrig;
else
    nlCoefOrig2pols=nlCoefOrig;
end

switch methodSSMF
    case 'stochastic'
        for k1 = 1:nModes*nPols/2 %%  ??? Check with Filipe
            nlCoefWeigthed(2*k1-1,1:2:nModes*nPols) = 2/1;
            nlCoefWeigthed(2*k1  ,2:2:nModes*nPols) = 2/1;
            
            nlCoefWeigthed(2*k1-1,2:2:nModes*nPols) = 2/3;
            nlCoefWeigthed(2*k1  ,1:2:nModes*nPols) = 2/3;
            
            nlCoefWeigthed(2*k1-1,2*k1-1) = 1/1;
            nlCoefWeigthed(2*k1  ,2*k1  ) = 1/1;
            try; nlCoefWeigthed(2*k1-1,2*k1  ) = 2/3; end
            try; nlCoefWeigthed(2*k1  ,2*k1-1) = 2/3; end
        end
        nlCoef2pols = nlCoefWeigthed.*nlCoefOrig2pols;
    case 'wcoupled'
        for k1 = 1:nModes
            nlCoefWeigthed(2*k1-1,1:2:nModes*nPols) = 4/3;
            nlCoefWeigthed(2*k1  ,1:2:nModes*nPols) = 4/3;
            
            nlCoefWeigthed(2*k1-1,2:2:nModes*nPols) = 4/3;
            nlCoefWeigthed(2*k1  ,2:2:nModes*nPols) = 4/3;
            
            nlCoefWeigthed(2*k1-1,2*k1-1) = 8/9;
            nlCoefWeigthed(2*k1  ,2*k1  ) = 8/9;
            try; nlCoefWeigthed(2*k1-1,2*k1  ) = 8/9; end
            try; nlCoefWeigthed(2*k1  ,2*k1-1) = 8/9; end
        end
        nlCoef2pols = nlCoefWeigthed.*nlCoefOrig2pols;
    case 'scoupled'
        nlCoef2pols = 4/3*12/13*mean(mean(nlCoefOrig2pols))*ones(nModes*nPols,nModes*nPols);
    case 'unitaryLc'
        for k1 = 1:nModes
            nlCoefWeigthed(2*k1-1,1:2:nModes*nPols) = 2/1;
            nlCoefWeigthed(2*k1  ,2:2:nModes*nPols) = 2/1;
            
            nlCoefWeigthed(2*k1-1,2:2:nModes*nPols) = 2/3;
            
            nlCoefWeigthed(2*k1  ,1:2:nModes*nPols) = 2/3;
            
            nlCoefWeigthed(2*k1-1,2*k1-1) = 1/1;
            nlCoefWeigthed(2*k1  ,2*k1  ) = 1/1;
            try; nlCoefWeigthed(2*k1-1,2*k1  ) = 2/3; end
            try; nlCoefWeigthed(2*k1  ,2*k1-1) = 2/3; end
        end
        nlCoef2pols = nlCoefWeigthed.*nlCoefOrig2pols;
end
nlCoef2pols = nlInd*nlCoef2pols;

%% Transmission
[Aout,~,H] = mmf_NL_xM_2pol(Ain,f,f0,DMD, Disp,  S, nlCoef2pols,lossCoef,L,dz,dz,CoupMatLPab,minStep,maxStep,delta_tol,fineStep);

% [H]=disp_comp_simulation_parameters(f,lossCoef,DMD,Disp,S,lambda0,dz,H);

disp_param.f=f;
disp_param.lossCoef=lossCoef;
disp_param.DMD=DMD;
disp_param.Disp=Disp;
disp_param.DispS=S;
disp_param.lambda0=lambda0;
disp_param.dz=dz;

varargout{1}       = varargin{1};
varargout{1}.E     = Aout;
varargout{1}.H     = H;
varargout{1}.disp_param = disp_param;
varargout{1}.freq_range     = f;



