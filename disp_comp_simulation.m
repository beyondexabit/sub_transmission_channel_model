function [H]=disp_comp_simulation(P,H)

f=P.f;
dz=P.dz;
lambda0=P.lambda0;

modes_Disp_stat ='mean';

switch modes_Disp_stat
    case 'mean'
lossCoef=repelem(mean(P.lossCoef),size(P.lossCoef,2));
DMD= repelem(mean(P.DMD),size(P.lossCoef,2)) ;
Disp= repelem(mean(P.Disp),size(P.lossCoef,2)) ;
DispS=repelem(mean(P.DispS),size(P.lossCoef,2));

    case '1to1'
lossCoef=P.lossCoef;
DMD= P.DMD;
Disp= P.Disp;
DispS= P.DispS;
end

e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

att   = lossCoef/(10*log10(exp(1)))/1e3;	% Np/m
beta2 = -Disp*lambda0.^2/(2*pi*c0);		%beta2 (ps^2/m)
beta3 = lambda0^2/(2*pi*c0)^2*...	    %beta3 (ps^3/m)
    (lambda0^2*DispS+2*lambda0*Disp);
nModes = size(att,2);

for k1 = 1:nModes
    ATT = att(k1); 
    B1  = DMD(k1);    
    B2  = beta2(k1);
    B3  = beta3(k1);
    
    if size(H,3)==1
 fullStep2D(k1,:) = exp(-(-ATT/2 -1i*B1*(2*pi*f) -1i*B2/2*(2*pi*f).^2 -1i*B3/6*(2*pi*f).^3)*dz);

    else
    fullStep2D(k1,k1,:) = exp(-(-ATT/2 -1i*B1*(2*pi*f) -1i*B2/2*(2*pi*f).^2 -1i*B3/6*(2*pi*f).^3)*dz);

    end
end
 if size(H,3)==1
     H=fullStep2D.*H;
 else
 H = pagemtimes(fullStep2D,H);
 end
end