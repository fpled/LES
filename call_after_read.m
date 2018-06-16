clc
clearvars
close all
% rng('default');
myparallel('start');

pathname = fileparts(mfilename('fullpath'));

sigma = 0.45; % surface tension (N/m)
mu_w = 0.1; % dynamic viscosity of water (Pa.s)
mu_o = 0.1; % dynamic viscosity of oil (Pa.s)
rho_w = 1e3; % mass density of water (kg/m3)
rho_o = 900; % mass density of oil (kg/m3)

nn = 13; % number of post-processed variables

% for g=2.^(4:8)
for g=2^4
    tic
    gridname = ['Grid' num2str(g)];
    disp(gridname)
    pathnamegrid = fullfile(pathname,gridname);
    load(fullfile(pathnamegrid,'data.mat'),'N','n','m','p');
    
    fprintf('\nn = %d variables',n);
    fprintf('\n  = %d post-processed variables',nn);
    fprintf('\nN = %d samples',N);
    fprintf('\nm = %d spatial points',m);
    fprintf('\np+1 = %d time steps',p+1);
    fprintf('\n');
    
    % Spatial scheme
    dx = 1/(g+1);
    % Implicit central-difference spatial scheme (second-order accurate, unconditionlly stable)
    Dx = spdiags(repmat([1 -1],g+1,1),[1 -1],g+1,g+1)/(2*dx);
    
    % Time scheme
    dt = 1/(p+1);
    % Explicit backward Euler time scheme (first-order accurate, conditionally stable)
    % Dt = spdiags(repmat([1 -1],p+1,1),[0 -1],p+1,p+1)/dt;
    % Implicit forward Euler time scheme (first-order accurate, unconditionally stable)
    % Dt = spdiags(repmat([1 -1],p+1,1),[1 0],p+1,p+1)/dt;
    % Implicit central-difference time scheme (second-order accurate, unconditionally stable)
    % Dt = spdiags(repmat([-1 2 -1],p+1,1),[1 0 -1],p+1,p+1)/(2*dt);
    Dt = toeplitz([2 -1 zeros(1,p-1)])/(2*dt);
    
    if g<2^7
        load(fullfile(pathnamegrid,'data.mat'),'Y');
        u = Y(:,1:3,:,:);
        C = Y(:,4,:,:);
        clear Y
        YY = zeros(N,nn,m,p+1);
    end
    
    for t=0:p
        time = ['Time ' num2str(t)];
        disp(time)
        
        ut_old = zeros(N,3,m);
        Ct_old = zeros(N,1,m);
        ut_new = zeros(N,3,m);
        Ct_new = zeros(N,1,m);
        if g<2^7
            ut = u(:,:,:,t+1);
            Ct = C(:,:,:,t+1);
            if t>0
                ut_old = u(:,:,:,t);
                Ct_old = C(:,:,:,t);
            end
            if t<p
                ut_new = u(:,:,:,t+2);
                Ct_new = C(:,:,:,t+2);
            end
        else
            load(fullfile(pathnamegrid,['data_t' num2str(t) '.mat']),'Yt');
            ut = Yt(:,1:3,:);
            Ct = Yt(:,4,:);
            if t>0
                load(fullfile(pathnamegrid,['data_t' num2str(t-1) '.mat']),'Yt');
                ut_old = Yt(:,1:3,:);
                Ct_old = Yt(:,4,:);
            end
            if t<p
                load(fullfile(pathnamegrid,['data_t' num2str(t+1) '.mat']),'Yt');
                ut_new = Yt(:,1:3,:);
                Ct_new = Yt(:,4,:);
            end
            clear Yt
        end
        
        YYt = zeros(N,nn,m);
        parfor l=1:N
            samplename = ['Sample ' num2str(l)];
            disp(samplename)
            
            ul = reshape(ut(l,:,:),[3,m]);
            ul_old = reshape(ut_old(l,:,:),[3,m]);
            ul_new = reshape(ut_new(l,:,:),[3,m]);
            Cl = Ct(l,:);
            Cl_old = Ct_old(l,:);
            Cl_new = Ct_new(l,:);
            rhol = Cl.*rho_o + (1-Cl).*rho_w;
            rhol_old = Cl_old.*rho_o + (1-Cl_old).*rho_w;
            rhol_new = Cl_new.*rho_o + (1-Cl_new).*rho_w;
            rhoul = repmat(rhol,3,1).*ul;
            rhoul_old = repmat(rhol_old,3,1).*ul_old;
            rhoul_new = repmat(rhol_new,3,1).*ul_new;
            tauTime = (-rhoul_old+2*rhoul-rhoul_new)/(2*dt);
            
            uijk = zeros(3,g+1,g+1,g+1);
            Cijk = zeros(g+1,g+1,g+1);
            for i=1:g+1
                for k=1:g+1
                    for j=1:g+1
                        ind = (g+1)^2*(i-1)+(g+1)*(k-1)+j;
                        uijk(:,i,j,k) = ul(:,ind);
                        Cijk(i,j,k) = Cl(ind);
                    end
                end
            end
            
            graduijk = zeros(3,3,g+1,g+1,g+1);
            gradCijk = zeros(3,g+1,g+1,g+1);
            for k=1:g+1
                for j=1:g+1
                    for c=1:3
                        u1 = uijk(c,:,j,k);
                        graduijk(1,c,:,j,k) = Dx*u1(:);
                    end
                    C1 = Cijk(:,j,k);
                    gradCijk(1,:,j,k) = Dx*C1(:);
                end
            end
            for k=1:g+1
                for i=1:g+1
                    for c=1:3
                        u2 = uijk(c,i,:,k);
                        graduijk(2,c,i,:,k) = Dx*u2(:);
                    end
                    C2 = Cijk(i,:,k);
                    gradCijk(2,i,:,k) = Dx*C2(:);
                end
            end
            for j=1:g+1
                for i=1:g+1
                    for c=1:3
                        u3 = uijk(c,i,j,:);
                        graduijk(3,c,i,j,:) = Dx*u3(:);
                    end
                    C3 = Cijk(i,j,:);
                    gradCijk(3,i,j,:) = Dx*C3(:);
                end
            end
            
            Sijk = (graduijk+permute(graduijk,[2,1,3,4,5]))/2; 
            
            normalijk = zeros(3,g+1,g+1,g+1);
            for k=1:g+1
                for j=1:g+1
                    for i=1:g+1
                        gradC = gradCijk(:,i,j,k);
                        if all(gradC==0)
                            normalijk(:,i,j,k) = zeros(3,1);
                        else
                            normalijk(:,i,j,k) = gradC./norm(gradC);
                        end
                    end
                end
            end
            
            divSijk = zeros(3,g+1,g+1,g+1);
            divnormalijk = zeros(g+1,g+1,g+1);
            for k=1:g+1
                for j=1:g+1
                    S1 = sum(Sijk(1,:,:,j,k),2);
                    n1 = normalijk(1,:,j,k);
                    divSijk(1,:,j,k) = Dx*S1(:);
                    divnormalijk(:,j,k) = Dx*n1(:);
                end
            end
            for k=1:g+1
                for i=1:g+1
                    S2 = sum(Sijk(2,:,i,:,k),2);
                    n2 = normalijk(2,i,:,k);
                    divn2 = divnormalijk(i,:,k);
                    divSijk(2,i,:,k) = Dx*S2(:);
                    divnormalijk(i,:,k) = divn2(:) + Dx*n2(:);
                end
            end
            for j=1:g+1
                for i=1:g+1
                    S3 = sum(Sijk(3,:,i,j,:),2);
                    n3 = normalijk(3,i,j,:);
                    divn3 = divnormalijk(i,j,:);
                    divSijk(3,i,j,:) = Dx*S3(:);
                    divnormalijk(i,j,:) = divn3(:)+ Dx*n3(:);
                end
            end
            
            tauConv = zeros(3,m);
            tauDiff = zeros(3,m);
            tauSurf = zeros(3,m);
            tauInterf = zeros(1,m);
            for k=1:g+1
                for j=1:g+1
                    for i=1:g+1
                        mu  = Cijk(i,j,k)*mu_o  + (1-Cijk(i,j,k))*mu_w;
                        rho = Cijk(i,j,k)*rho_o + (1-Cijk(i,j,k))*rho_w;
                        
                        uu = uijk(:,i,j,k);
                        gradu = graduijk(:,:,i,j,k);
                        gradC = gradCijk(:,i,j,k);
                        divS  = divSijk(:,i,j,k);
                        kappa = divnormalijk(i,j,k);
                        
                        tConv = rho*gradu'*uu;
                        tDiff = 2*mu*divS;
                        tSurf = sigma*kappa*gradC;
                        tInterf = dot(uu,gradC,1);
                        
                        ind = (g+1)^2*(i-1)+(g+1)*(k-1)+j;
                        tauConv(:,ind) = tConv;
                        tauDiff(:,ind) = tDiff;
                        tauSurf(:,ind) = tSurf;
                        tauInterf(1,ind) = tInterf;
                    end
                end
            end
            
            YYl = cat(1,tauTime,tauConv,tauDiff,tauSurf,tauInterf);
            YYt(l,:,:) = YYl;
        end
        if g>=2^7
            save(fullfile(pathname,gridname,['data_post_t' num2str(t) '.mat']),'YYt');
        else
            YY(:,:,:,t+1) = YYt;
        end
    end
    fprintf('\n');
    if g<2^7
        save(fullfile(pathname,gridname,'data_post.mat'),'YY','nn');
    else
        save(fullfile(pathname,gridname,'data_post.mat'),'nn');
    end
    toc
end

myparallel('stop');
