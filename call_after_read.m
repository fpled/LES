clc
clearvars
close all
% rng('default');
% myparallel('start');

pathname = fileparts(mfilename('fullpath'));
% pathname = '/mnt/tcm13/SV_FP/';

L = 1; % domain size (m)
sigma = 0.45; % surface tension (N/m)
mu = [0.1,0.1]; % dynamic viscosity of [water,oil] (Pa.s)
rho = [1000,900]; % mass density of [water,oil] (kg/m3)

% Time scheme
dt = 5e-3; % time step (s)
dt = 100*dt; % physical time step stored every 100 time iterations

nn = 13; % number of post-processed variables

order = [3 1 2]; % dimension order

% for g=2.^(4:8)
for g=2^4
    tic
    gridname = ['Grid' num2str(g)];
    disp(gridname)
    gridpathname = fullfile(pathname,gridname);
    load(fullfile(gridpathname,'data.mat'),'N','n','m','p');
    
    fprintf('\nn = %d variables',n);
    fprintf('\n  = %d post-processed variables',nn);
    fprintf('\nN = %d samples',N);
    fprintf('\nm = %d spatial points',m);
    fprintf('\np+1 = %d time steps',p+1);
    fprintf('\n');
    
    s = [g+1,g+1,g+1]; % spatial dimensions
    
    % Spatial scheme
    dx = L/g; % spatial step (m)
    % Implicit central-difference spatial scheme (second-order accurate, unconditionlly stable)
    Dx = spdiags(repmat([1 -1],g+1,1),[1 -1],g+1,g+1)/(2*dx);
    
    % Explicit forward Euler time scheme (first-order accurate, conditionally stable)
    % Dt = spdiags(repmat([1 -1],p+1,1),[0 -1],p+1,p+1)/dt;
    % Implicit backward Euler time scheme (first-order accurate, unconditionally stable)
    % Dt = spdiags(repmat([1 -1],p+1,1),[1 0],p+1,p+1)/dt;
    % Implicit central-difference time scheme (second-order accurate, unconditionally stable)
    Dt = spdiags(repmat([1 -1],p+1,1),[1 -1],p+1,p+1)/(2*dt);
    
    if g<2^7
        load(fullfile(gridpathname,'data.mat'),'Y');
        YY = zeros(N,nn,m,p+1);
    end
    
    fprintf('\nPost-processing data');
    fprintf('\n');
    for t=0:p
        time = ['Time ' num2str(t)];
        disp(time)
        
        Yt_old = zeros(N,n,m);
        Yt_new = zeros(N,n,m);
        if g<2^7
            if t>0
                Yt_old = Y(:,:,:,t);
            end
            if t<p
                Yt_new = Y(:,:,:,t+2);
            end
            Yt = Y(:,:,:,t+1);
        else
            if t>0
                load(fullfile(gridpathname,['data_t' num2str(t-1) '.mat']),'Yt');
                Yt_old = Yt;
            end
            if t<p
                load(fullfile(gridpathname,['data_t' num2str(t+1) '.mat']),'Yt');
                Yt_new = Yt;
            end
            load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
        end
        
        ut = Yt(:,1:3,:);
        ut_old = Yt_old(:,1:3,:);
        ut_new = Yt_new(:,1:3,:);
        
        Ct = Yt(:,4,:);
        Ct_old = Yt_old(:,4,:);
        Ct_new = Yt_new(:,4,:);
        
        % rhot = Ct*rho(2) + (1-Ct)*rho(1);
        rhot_old = Ct_old*rho(2) + (1-Ct_old)*rho(1);
        rhot_new = Ct_new*rho(2) + (1-Ct_new)*rho(1);
        % rhout = rhot.*ut;
        rhout_old = rhot_old.*ut_old;
        rhout_new = rhot_new.*ut_new;
        tauTime = (rhout_new-rhout_old)/(2*dt);
        
        ut = permute(reshape(ut,[N,3,s]),[2,order+2,1]);
        Ct = permute(reshape(Ct,[N,s]),[order+1,1]);
        rhot = Ct*rho(2) + (1-Ct)*rho(1);
        mut  = Ct*mu(2) + (1-Ct)*mu(1);
        gradut = grad(ut,Dx);
        gradrhout = grad(shiftdim(rhot,-1).*ut,Dx);
        tauConv = compute_tauConv(ut,gradrhout);
        tauDiff = compute_tauDiff(2*shiftdim(mut,-2).*gradut,Dx);
        gradCt = grad(Ct,Dx);
        ngradCt = normal(gradCt);
        kappa = div(ngradCt,Dx);
        tauSurf = sigma*shiftdim(kappa,-1).*gradCt;
        tauInterf = compute_tauInterf(ut,gradCt);
        tauConv = ipermute(tauConv,[2,order+2,1]);
        tauDiff = ipermute(tauDiff,[2,order+2,1]);
        tauSurf = ipermute(tauSurf,[2,order+2,1]);
        tauInterf = ipermute(tauInterf,[2,order+2,1]);
        tauConv = tauConv(:,:,:);
        tauDiff = tauDiff(:,:,:);
        tauSurf = tauSurf(:,:,:);
        tauInterf = tauInterf(:,:,:);
        YYt = cat(2,tauTime,tauConv,tauDiff,tauSurf,tauInterf);
        
%         YYt = zeros(N,nn,m);
%         parfor l=1:N
%             % samplename = ['Sample ' num2str(l)];
%             % disp(samplename)
%             
%             ul = permute(reshape(ut(l,:,:),[3,s]),[1,order+1]);
%             Cl = permute(reshape(Ct(l,:),s),order);
%             
%             gradul = grad(ul,Dx);
%             rhol = Cl*rho(2) + (1-Cl)*rho(1);
%             mul  = Cl*mu(2) + (1-Cl)*mu(1);
%             tauConv = shiftdim(rhol,-1).*compute_tauConv(ul,gradul);
%             tauDiff = 2*shiftdim(mul,-1).*compute_tauDiff(gradul,Dx);
%             
%             gradCl = grad(Cl,Dx);
%             ngradCl = normal(gradCl);
%             kappa = div(ngradCl,Dx);
%             tauSurf = sigma*shiftdim(kappa,-1).*gradCl;
%             tauInterf = compute_tauInterf(ul,gradCl);
%             
%             tauConv = ipermute(tauConv,[1,order+1]);
%             tauDiff = ipermute(tauDiff,[1,order+1]);
%             tauSurf = ipermute(tauSurf,[1,order+1]);
%             tauInterf = ipermute(tauInterf,[1,order+1]);
%             
%             tauConv = tauConv(:,:);
%             tauDiff = tauDiff(:,:);
%             tauSurf = tauSurf(:,:);
%             tauInterf = tauInterf(:,:);
%             tauTimel = shiftdim(tauTime(l,:,:));
%                         
%             YYl = cat(1,tauTimel,tauConv,tauDiff,tauSurf,tauInterf);
%             YYt(l,:,:) = YYl;
%         end

        if g<2^7
            YY(:,:,:,t+1) = YYt;
        else
            save(fullfile(gridpathname,['data_post_t' num2str(t) '.mat']),'YYt');
        end
    end
    fprintf('\n');
    if g<2^7
        save(fullfile(gridpathname,'data_post.mat'),'YY','nn');
    else
        save(fullfile(gridpathname,'data_post.mat'),'nn');
    end
    toc
end

% myparallel('stop');
