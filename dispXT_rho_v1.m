function [XT1d_mean,XT1d_rho_mean,rho,phi,xDisp,yDisp] = dispXT_rho_v1(P,dMax,surface_type,wiRad,nopF,dz,xtmethod,nPols,theta_modes,rho1,num_phi)
%% dMax:
%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = physconst('LightSpeed');%1/sqrt(e0*u0);

%% input parameter parsing
if nargin < 5
    nopF        = 50;   % number of points somehow
end
if nargin < 4
    wiRad = 0; % no radiant modes
    warning('no radiant modes included')
end
if nargin < 3
    surface_type = 'coarse';
    warning('going for coarse surface')
end
if nargin < 2
    dMax = 0.01;
    warning('dMax assumed to be 0.01')
end
%% Surface definition
switch surface_type
    case '2d_fine'
        rho = [0 1e-6:1e-6:1e-4 2e-4:1e-4:0.01 0.011:0.001:dMax];
        rho = unique(rho); rho = sort(rho);
        phi = -pi:2*pi/(360):pi-2*pi/(360);%
    case '2d_coarse'
        rho = [0 1e-4:4e-4:0.01 0.011:0.004:dMax];
        rho = unique(rho); rho = sort(rho);
        phi = -pi:2*pi/(90):pi-2*pi/(90);%
    case 'input'
        rho=rho1;
        phi=linspace(0,2*pi,num_phi);
end

rDisp = ones(length(phi),1)*rho;
pDisp = phi.'*ones(1,length(rho));
xDisp = rDisp.*cos(pDisp);
yDisp = rDisp.*sin(pDisp);

%% parameters
lambda0  = P.lambda0;
w        = 2*pi*c0/P.lambda0;
w0       = 2*pi*c0/P.lambda0;

%% try to load
P.freq0    = physconst('LightSpeed')/P.lambda0;
P.freq     = physconst('LightSpeed')/P.lambda;
elp = P.elliptical;
fibreTag = getFibreTag(P);

fibreTag = [fibreTag,'_ds=',surface_type,'_dMax=',num2str(dMax),'_wiRad=',num2str(wiRad),num2str(elp),num2str(dz),num2str(xtmethod),num2str(nPols),num2str(length(rho1)),num2str(num_phi1)];

tag1     = ['dispXT'];
XT1d_rho_mean=[];
XT1d_mean =[];
% XT1d_rho_mean = loadComputedData([tag1,'/',fibreTag],['XT1d_rho_mean']);
% XT1d_mean = loadComputedData([tag1,'/',fibreTag],['XT1d_mean']);

if ~isempty(XT1d_rho_mean) && ~isempty(XT1d_mean)
    return;
end

%% Fiber Parameters
disp('loading mode field distribution...')
P1 = P; %P1.getField = 1;
temp   = femFiberAnalyzerLp(P1,1,1);
if wiRad == 1
    nModes = length(temp.neff_arb);
else
    nModes = temp.tnom;
end



pc = temp.neff_arb(1:nModes);

[KuvSurf] = slImpKuvSurfLp(P,dMax,surface_type,wiRad,nopF);


for k1 = 1:nModes
    for k2 = 1:nModes
        deltaBeta1(k1,k2) = ( pc(k1) - pc(k2) ) * w0 / c0;
    end
end

for k1 = 1:nModes
    B(k1) = ( pc(k1) ) * w0 / c0;
end

%% Coupling Coefficients
coupMat = zeros(nModes,nModes,numel(xDisp));
Aout    = zeros(nModes,nModes,numel(xDisp));

XT1d_temp = zeros(nModes*nPols,size(rho,2),size(phi,2));

%%

disp('interpolating kuv...')
kuvX = zeros([nModes,nModes,size(xDisp)]);
for k1 = 1:nModes
    fprintf([num2str((nModes) - k1 + 1),' ']);
    if ~rem(k1,20); fprintf('\n'); end
    for k2 = 1:nModes
        kuvX(k1,k2,:,:) = KuvSurf(k1,k2).s(xDisp , yDisp);
    end
end

disp('calculating disp xt...')
XT1d_mean = zeros(nModes*nPols,size(rho,2));

for kk1 = 1:size(rho,2)
    for kk2 = 1:size(phi,2)

        if nPols==2
            kuvX_2D = zeros(nModes*nPols,nModes*nPols);
            kuvX_2D(1:2:end,1:2:end)=kuvX(:,:,kk2,kk1);
            kuvX_2D(1:2:end,2:2:end)=0;
            kuvX_2D(2:2:end,1:2:end)=0;
            kuvX_2D(2:2:end,2:2:end)=kuvX(:,:,kk2,kk1);
            B_2pol=diag(repelem(B,2));
            [kuvX_modified]=polrot_Kuv_2D(kuvX_2D,nModes,theta_modes);
            coupMat = expm(1i*(B_2pol+kuvX_modified)*dz);
            [XT1d_temp(:,kk1,kk2),~] = calc_xt(coupMat,xtmethod,nModes,nPols,pc);
        else
            kuvX_modified=kuvX;
            coupMat = expm(1i*(diag(B)+kuvX_modified)*dz);
            [XT1d_temp(:,kk1,kk2),~] = calc_xt(coupMat,xtmethod,nModes,nPols,pc);
       end
    end
    XT1d_mean(:,kk1)=pow2db(mean(squeeze(db2pow(XT1d_temp(:,kk1,kk2))),2,'omitnan'));
end
XT1d_rho_mean=XT1d_temp;
fprintf('\n');
% XT1d_rho_mean = squeeze(10*log10(mean(XT1d_mean)));


saveComputedData(xDisp, [tag1,'/',fibreTag],['xDisp']);%tag1,['xDisp_',fileName1]);
saveComputedData(yDisp, [tag1,'/',fibreTag],['yDisp']);%tag1,['yDisp_',fileName1]);
saveComputedData(rho, [tag1,'/',fibreTag],['rho']);%tag1,['rho_',fileName1]);
saveComputedData(phi, [tag1,'/',fibreTag],['phi']);%tag1,['phi_',fileName1]);
saveComputedData(XT1d_rho_mean,    [tag1,'/',fibreTag],['XT1d_rho_mean']);
saveComputedData(XT1d_mean,    [tag1,'/',fibreTag],['XT1d_mean']);
end


%% function generate fibre tag
function fibreTag = getFibreTag(profPar)
if profPar.profType == 1
    fibreTag = [...
        'fem_',num2str(profPar.profType),...
        '_d1',num2str(profPar.d1_0*1000),'a1',num2str(profPar.a1_0),'w1',num2str(profPar.w1),'\',...
        'w2',num2str(profPar.w2),'\',...
        'd3',num2str(profPar.d3_0*1000),'w3',num2str(profPar.w3),'\',...
        'n',num2str(profPar.NoP),'_c',num2str(profPar.cladD),'\',...
        'f0_',num2str(profPar.freq0),'f_',num2str(profPar.freq)];
    if isfield(profPar,'solverTol')
        fibreTag = [fibreTag,'_mSolT',num2str(P.solverTol)];
    end
    if isfield(profPar,'solverMaxIter')
        fibreTag = [fibreTag,'_mSolMI',num2str(P.solverMaxIter)];
    end
    if isfield(profPar,'getField')
        fibreTag = [fibreTag,'_getField',num2str(profPar.getField)];
    end
else
    error('unknow fibre profile')
end
end






