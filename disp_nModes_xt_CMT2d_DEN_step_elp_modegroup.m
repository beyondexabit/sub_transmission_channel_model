function [rho,XT1d_rho_max] = disp_nModes_xt_CMT2d_DEN_step_elp_modegroup(rho,num_phi,step,method,xtmethod,nModes,nPols,pc,KuvSurf,theta_modes,num_theta,statmethod,elp,polrotmethod,deltacore)
%% Condider xtmethod as degenerate - works better

if num_phi < 1000
    num_phi = 1000;
elseif mod(num_phi,1000)
    num_phi = ceil(num_phi/1000)*1000;%10^(round((log10(num_phi))));
end

rho = sort(unique(rho));


%% Physical constans
e0      = 8.854187817e-12;
u0      = 1.25663706e-6;
c0      = 1/sqrt(e0*u0);

%% Source parameters
lambda0  = 1550e-9;
w0       = 2*pi*c0/lambda0;

%% try to load

fileName1  = [num2str(max(rho),'%.12f'),num2str(sum(rho),'%.12f'),num2str(length(rho),'%.12f'),num2str(num_phi),num2str(theta_modes),'_',num2str(step),'_',method,'_DEN_',num2str(xtmethod),'_nModes_',num2str(nModes),'_nPols_',num2str(nPols),num2str(num_theta),num2str(statmethod),num2str(elp),polrotmethod,num2str(deltacore)];

% tag1       = ['Xdisp_xt_CMT2d\den_step_elp_test2\'];
tag1       = ['Xdisp_xt_CMT2d\den_step_elp_new\'];

XT1d_rho_max=[];
rho1=[];
rho1 = loadComputedData(tag1,['rho_',fileName1],0);
XT1d_rho_max = loadComputedData(tag1,['XT1d_rho_max_', fileName1],0);

if  ~isempty(rho1) && ~isempty(XT1d_rho_max)
    rho = rho1;
    return
end
clear rho1 XT1d_rho_max

%% Fiber Parameters
B = diag(pc) * w0 / c0;

%% Mode Field

phi = linspace(0,2*pi,num_phi);


%% Coupling Coefficients
coupMat = zeros(nModes*nPols,nModes*nPols);
XT1d_max = zeros(length(rho),length(phi),nModes*nPols);

%%

nStart = 1;
for kk1 = 1000:1000:num_phi
    fileName1  = [num2str(max(rho),'%.12f'),num2str(sum(rho),'%.12f'),num2str(length(rho),'%.12f'),num2str(kk1),num2str(theta_modes),'_',num2str(step),'_',method,'_DEN_',num2str(xtmethod),'_nModes_',num2str(nModes),'_nPols_',num2str(nPols),num2str(num_theta),num2str(statmethod),num2str(elp),polrotmethod,num2str(deltacore)];
    XT1d_max_temp =[];
    XT1d_max_temp = loadComputedData(tag1,['XT1d_max_', fileName1],0);

    if ~isempty(XT1d_max_temp)
        XT1d_max(:,1:kk1,:) = XT1d_max_temp(:,1:kk1,:);
        nStart = kk1+1;
    else
        nStart = kk1-1000+1;
        break;
    end
end

h=tic;
h2=tic;

for kk2 = nStart:num_phi
    if mod(kk2,100) == 0
        if kk2 > 1
            tt = toc(h); h = tic; fprintf([num2str(kk2),' out of ',num2str(num_phi),' - time to finish =',num2str(tt*(num_phi-kk2)/60/100),' min \n'])
        end
    end
    for kk1 = 1:length(rho)
        rDisp = rho(kk1);
        pDisp =2*pi*rand(1);% phi(kk2);


        xDisp = rDisp.*cos(pDisp);
        yDisp = rDisp.*sin(pDisp);

        kuvX=zeros(nPols*nModes,nPols*nModes);

        for k1=1:nPols*nModes 
            for k2=1:nPols*nModes
                kuvX(k1,k2)= KuvSurf(k1,k2).s(xDisp,yDisp);
            end
        end

        [ kuvX_modified,~]=polrot_Kuv_2D(kuvX,nModes,theta_modes,polrotmethod);

        coupMat = expm(1i*(B+kuvX_modified)*step);
        [XT1d_temp,~] = calc_xt(coupMat,xtmethod,nModes,nPols,pc);

        switch xtmethod
            case  'intermodegroup'
                XT1d_temp = pow2db(mean(db2pow(XT1d_temp),1,'omitnan')); %% mEAN XT OF MODEGROUPS
            otherwise
                XT1d_temp=XT1d_temp;
                XT1d_max(kk1,kk2,1:length(XT1d_temp))=XT1d_temp;
        end

%         switch statmethod
%             case 'median'
%                 XT1d_max(kk1,kk2,:)   = pow2db(median(db2pow(XT1d_temp),2,'omitnan'));
%             case 'mean'
%                 XT1d_max(kk1,kk2,:)   = pow2db(mean(db2pow(XT1d_temp),2,'omitnan'));
%         end


    end

    if mod(kk2,1000) == 0
        fileName1  = [num2str(max(rho),'%.12f'),num2str(sum(rho),'%.12f'),num2str(length(rho),'%.12f'),num2str(kk2),num2str(theta_modes),'_',num2str(step),'_',method,'_DEN_',num2str(xtmethod),'_nModes_',num2str(nModes),'_nPols_',num2str(nPols),num2str(num_theta),num2str(statmethod),num2str(elp),polrotmethod,num2str(deltacore)];

        saveComputedData(rho, tag1,['rho_',fileName1],0);
        saveComputedData(XT1d_max,tag1,['XT1d_max_', fileName1],0);
    end
end


switch statmethod
    case 'median'
        XT1d_rho_max(:,:)  = pow2db(median(db2pow(XT1d_max),2,'omitnan'));
    case 'mean'
        XT1d_rho_max(:,:)  = pow2db(mean(db2pow(XT1d_max),2,'omitnan'));
end

tt2=toc(h2);

saveComputedData(rho, tag1,['rho_',fileName1],0);
saveComputedData(XT1d_rho_max,tag1,['XT1d_rho_max_', fileName1],0);

123;%

