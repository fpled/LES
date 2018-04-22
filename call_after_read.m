clc
clearvars
close all

pathname = fileparts(mfilename('fullpath'));

sigma = 0.45; % surface tension (N/m)
mu_w = 0.1; % dynamic viscosity of water (Pa.s)
mu_o = 0.1; % dynamic viscosity of oil (Pa.s)
rho_w = 1e3; % mass density of water (kg/m3)
rho_o = 900; % mass density of oil (kg/m3)

nn = 10; % number of post-processed variables

for g=2.^(4:7)
% for g=2^4
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
    Dt = spdiags(repmat([1 -1],p+1,1),[0 -1],p+1,p+1)/dt;
    % Implicit forward Euler time scheme (first-order accurate, unconditionally stable)
    % Dt = spdiags(repmat([1 -1],p+1,1),[1 0],p+1,p+1)/dt;
    
    %D = spdiags(repmat([-1 2 -1],g+1,1),[1 0 -1],g+1,g+1)/(2*dx);
    %D = (2*diag(ones(g+1,1),0) - diag(ones(g,1),-1) - diag(ones(g,1),1))/(2*dx);
    %D = toeplitz([2 -1 zeros(1,g-2)])/(2*dx);
    
    if g<2^8
        load(fullfile(pathnamegrid,'data.mat'),'Y');
        u = Y(:,1:3,:,:);
        C = Y(:,4,:,:);
        YY = zeros(N,nn,m,p+1);
    end
    
    for t=0:p
        time = ['Time ' num2str(t)];
        disp(time)
        
        if g<2^8
            ut = u(:,:,:,t+1);
            Ct = C(:,:,:,t+1);
        else
            load(fullfile(pathnamegrid,['data_t' num2str(t) '.mat']),'Yt');
            ut = Yt(:,1:3,:);
            Ct = Yt(:,4,:);
            YYt = zeros(N,nn,m);
        end
        
        for l=1:N
            % samplename = ['Sample ' num2str(l)];
            % disp(samplename)
            
            ul = reshape(ut(l,:,:),[3,m]);
            Cl = Ct(l,:);
            
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
            clear ul Cl
            
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
            
            tauTime = zeros(3,m);
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
                        tDiff = mu*divS;
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
            
            if g<2^8
                YY(l,1:3,:,t+1) = tauConv;
                YY(l,4:6,:,t+1) = tauDiff;
                YY(l,7:9,:,t+1) = tauSurf;
                YY(l,10,:,t+1) = tauInterf;
            else
                YYt(l,1:3,:) = tauConv;
                YYt(l,4:6,:) = tauDiff;
                YYt(l,7:9,:) = tauSurf;
                YYt(l,10,:) = tauInterf;
            end
        end
        if g>=2^8
            save(fullfile(pathname,gridname,['data_post_t' num2str(t) '.mat']),'YYt');
        end
    end
    fprintf('\n');
    if g<2^8
        save(fullfile(pathname,gridname,'data_post.mat'),'YY','nn');
    else
        save(fullfile(pathname,gridname,'data_post.mat'),'nn');
    end
end
