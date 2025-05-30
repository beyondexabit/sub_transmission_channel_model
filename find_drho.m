function [drho] = find_drho(rho,XT_rho,XTavg,dz,drhomethod,npts)

switch drhomethod
    case 'mean'
aux(:) = pow2db(mean(db2pow(XT_rho),2,'omitnan'));
    case 'max'
[aux(:)] = max(XT_rho,[],2);
    case 'median'
[aux(:)] = median(XT_rho,2);
case 'min'
[aux(:)] = min(XT_rho,[],2);
case 'LP02'
    if size(XT_rho,2)==12
[aux(:)] = pow2db(mean(db2pow(XT_rho(:,7:8)),2,'omitnan'));
    elseif size(XT_rho,2)==6
[aux(:)] = XT_rho(:,4);
    end
    case 'LP11'
        [aux(:)] = pow2db(mean(db2pow(XT_rho(:,3:6)),2,'omitnan'));
end

if nargin<6
    npts=1e-6;
end

rhoExt = rho(1):npts:rho(end);
auxExt = interp1(pow2db(rho(1:end)),aux(1:end),pow2db(rhoExt),'linear'); 
XTavg_dz   = XTavg + 10*log10(dz/1000);
Ind    = 1;

while auxExt(Ind) < XTavg_dz && Ind < length(auxExt)
    Ind = Ind + 1;
end

if Ind==1
    disp('error - increase rho lower value')
end

XT_min  =[auxExt(Ind-1) auxExt(Ind)];
rho_min =[rhoExt(Ind-1) rhoExt(Ind)];

[~,loc] = min(abs(XT_min-XTavg_dz));
drho=rho_min(loc);

% figure()
% semilogx(rho,aux,'-x')
% hold on;
% semilogx(rhoExt,auxExt,'-x')

end