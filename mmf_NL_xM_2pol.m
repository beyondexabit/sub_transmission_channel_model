function [A0,A0_fft,H] = mmf_NL_xM_2pol(A,f,f0,DMD,Disp,DispS,nlCoef,lossCoef,L,dz,XT_RefL,CoupMatLPab,minStep,maxStep,delta_tol,fineStep,GPU)
%Six Mode Fiber Transmission
%
% INPUTS
%    A              - signal Nmodes x Nsamples
%    f              - frequency vector 
%    f0             - center frequency
%    DMD            - differential mode delay
%    Disp           - mode dispersion
%    DispS          - mode dipersion slope
%    nlCoef         - nonlinear coeffs
%    lossCoef       - Fiber attenuation (dB/km)
%    L              - Fiber length (m)
%    dz             - Step length (m)
%    XT_RefL        - xt step important when doing NL tx with adapt step
%    CoupMatLPab    - step-by-step coupling matrix
%    minStep        - minimum step size
%    maxStep        - maximum step size
%    delta_tol      - local error bound for adaptative step
%    finestep       - override flag for adaptative step
%
% OUTPUTS
%    A              - signal Nmodes x Nsamples
%
% REFERENCES
%
%
% REVISION
%	26-09-2016 10:08:29 - created by Filipe Ferreira
%

%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

% if GPU
%     A = gpuArray(A);
% end

%% Propagations constants load or calculation
nModes = size(A,1);

%% parameters and global variables
att = lossCoef/(10*log10(exp(1)))/1e3;	% Np/m

w0      = 2*pi*f0;
lambda0 = 2*pi*c0/w0;

%% Dispersion operator for LPuv,ab,xy
beta2 = -Disp*lambda0.^2/(2*pi*c0);		%beta2 (ps^2/m)
beta3 = lambda0^2/(2*pi*c0)^2*...	    %beta3 (ps^3/m)
    (lambda0^2*DispS+2*lambda0*Disp);

D = zeros(nModes,length(f));
D2d = zeros(nModes,nModes,length(f));
fullStep2D = zeros(nModes,nModes,length(f));

for k1 = 1:nModes
    ATT = att(k1); 
    B1  = DMD(k1);    
    B2  = beta2(k1);
    B3  = beta3(k1);
    
    D(k1, :) = -ATT/2 -1i*B1*(2*pi*f) -1i*B2/2*(2*pi*f).^2 -1i*B3/6*(2*pi*f).^3;
    D2d(k1,k1 , :) = -ATT/2 -1i*B1*(2*pi*f) -1i*B2/2*(2*pi*f).^2 -1i*B3/6*(2*pi*f).^3;
    fullStep2D(k1,k1,:) = exp((-ATT/2 -1i*B1*(2*pi*f) -1i*B2/2*(2*pi*f).^2 -1i*B3/6*(2*pi*f).^3)*dz);
end

H = zeros(nModes,nModes,length(f));
for i = 1:nModes
    H(i,i,:) = 1;
end

z  = 0;
% if GPU
%     Q = zeros(nModes,nModes,L/dz,'gpuArray');
% else
    Q = (zeros(nModes,nModes,L/dz));
% end
for n = 0:1e6
    if z == L
        break;
    elseif z+dz > L
        dz = L-z;
    end
    Q(:,:,z/XT_RefL+1) = (CoupMatLPab(z/XT_RefL+1).c(1:1:nModes,1:1:nModes));
    z = z + dz;
end

%% Linear Transmission
if max(max(abs(nlCoef))) == 0
    z  = 0;
    
%     fprintf('START')
    A0       = A;
    A0_fft   = (fftshift(fft(ifftshift(A0,2),[],2),2));
%     if GPU
%         fullStep = gpuArray(exp(D*dz));
%     else
        fullStep = (exp(D*dz));
% end

%     for i = 1:nModes
%         fullStep2D(i,i,:) = exp(D2d(i,i,:)*dz);
%     end
    
    dz0      = dz;
    
    tic
    for n = 0:1e6
        %fprintf('.')
        if rem(n,1000) == 0 && n > 0
            fprintf([' ',num2str(n),' '])
        end
        if rem(n,1000) == 0
            tp1000 = toc;
            time2f = (L/dz-n)/1000*tp1000/60;
            fprintf(['\n',num2str(time2f),' mins \n'])
            tic
        end
        
        %% Step Correction and Segment Type
        if z == L
            A = A0;
            break;
        elseif z+dz > L
            dz = L-z;
        end
        
        %% Mode Coupling
        if rem(z,XT_RefL) == 0
            A0_fft = Q(:,:,z/XT_RefL+1) * A0_fft;
            %tic
            H = pagemtimes(Q(:,:,z/XT_RefL+1),H);
            %toc
        end
        
        %% Dispersion Step
        A0_fft = fullStep .* A0_fft;
%         for i = 1:length(f)
%             H(:,:,i) = fullStep2D(:,:,i)*H(:,:,i); 
%         end
%         tic
        H = pagemtimes(fullStep2D,H);
%         toc
        
        z = z + dz;
    end
    A0     = fftshift(ifft( ifftshift( A0_fft,2), [], 2 ),2);
%     if GPU
%         A0     = gather(A0);
%         A0_fft = gather(A0_fft);
%         H      = gather(H);
%     end
    return;
end

%% Nonlinear Transmission
z  = 0;
quarStep = exp(D*dz/4);
halfStep = exp(D*dz/2);

fprintf('START')
A0       = A;
A0_fft   = fftshift(fft(A0,[],2),2);
dz0      = dz;

for n = 0:1e6
    fprintf('.')
    if rem(n,100) == 0
        fprintf(['\n'])
    end
    
    %% Step Correction and Segment Type
    if z == L
        A = A0;
        break;
    elseif z+dz > L
        dz = L-z;
    end
    
    % in any case, check if segType needs update
    if z+dz == L || dz ~= dz0
        quarStep = exp(D*dz/4);
        halfStep = exp(D*dz/2);
        dz0      = dz;
    end
    
    %% Mode Coupling
    if rem(z,XT_RefL) == 0
        [Q,~] = qr(CoupMatLPab(z/XT_RefL+1).c(1:1:nModes,1:1:nModes));
        A0    = Q * A0;
        A0_fft   = fftshift(fft(A0,[],2),2);
    end
    
    %% Coarse Step
    Ac = ifft( ifftshift(halfStep .*  A0_fft,2), [], 2 );
    Ac = Ac.*exp( -1i * nlCoef*abs(A0).^2*dz);
    if fineStep == 1
        Ac_fft = halfStep .* fftshift(fft(Ac,[],2),2);			            %simple Split Step method
        Ac     = ifft( ifftshift( Ac_fft,2), [], 2 );
    else
        disp(['z = ',num2str(round(z+dz)),' - dz = ',num2str(dz)])
        A0_fft = halfStep .* fftshift(fft(Ac,[],2),2);			            %simple Split Step method
        A0     = ifft( ifftshift( A0_fft,2), [], 2 );
        
        z = z + dz;
    end
    
    %% Fine Step
    while fineStep == 1
        %% Update dispersion step if step length changed
        if dz ~= dz0
            quarStep = exp(D*dz/4);
            halfStep = exp(D*dz/2);
            dz0          = dz;
        end
        
        %% Fine Step - 1st part
        Af = ifft( ifftshift(quarStep .* A0_fft ,2) , [], 2 );
        Af = Af.*exp( -1i * nlCoef*abs(A0).^2*dz/2);
        Af_fft_1   = quarStep .* fftshift( fft(Af,[],2) ,2);
        Af_1       = ifft( ifftshift( Af_fft_1 ,2) , [], 2 );
        
        %% Fine Step - 2nd part
        Af = ifft( ifftshift(quarStep .* Af_fft_1 ,2) , [], 2 );
        Af = Af.*exp( -1i * nlCoef*abs(Af_1).^2*dz/2);
        Af_fft   = quarStep .* fftshift( fft(Af,[],2) ,2);
        Af       = ifft( ifftshift( Af_fft,2) , [], 2 );
        
        %% Relative Local Error (delta)
        delta = norm(Af-Ac,2) / norm(Af,2);
        
        %% Decision based on Precision vs Error
        disp(['z = ',num2str(round(z+dz)),' - dz = ',num2str(dz),' - delta = ',num2str(delta)])
        if delta < 0.5*delta_tol
            z = z + dz;
            if dz * 2^(1/3) < maxStep;
                dz = dz * 2^(1/3);
                if abs(round((z+dz)/XT_RefL)*XT_RefL - (z+dz)) < 1e-12
                    dz = round((z+dz)/XT_RefL)*XT_RefL - z;
                end
            end;
            
            A0     = 4/3*Af - 1/3*Ac;
            A0_fft = 4/3*Af_fft - 1/3*Ac_fft;
            break;
        elseif delta > 0.5*delta_tol && delta < 1*delta_tol
            z = z + dz;
            A0     = 4/3*Af - 1/3*Ac;
            A0_fft = 4/3*Af_fft - 1/3*Ac_fft;
            break;
        elseif delta > 1*delta_tol && delta < 2*delta_tol
            z = z + dz;
            dz = dz / 2^(1/3); if z+10*dz < L && dz < minStep; error('lowest step reached!!!'); end
            A0     = 4/3*Af - 1/3*Ac;
            A0_fft = 4/3*Af_fft - 1/3*Ac_fft;
            break;
        elseif dz <= minStep
            z = z + dz;
            A0     = 4/3*Af - 1/3*Ac;
            A0_fft = 4/3*Af_fft - 1/3*Ac_fft;
            break;
        elseif delta > 2*delta_tol
            dz     = dz/2;
            Ac     = Af_1;
            Ac_fft = Af_fft_1;
            % do it again!;
        end
    end
    
    
end

A0     = gather(A0); 
A0_fft = gather(A0_fft); 
H      = gather(H); 

fprintf('\n mmf DONE!\n')
