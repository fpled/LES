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
gravity = 9.81; % gravity (m/s2)

% Time scheme
dt = 5e-3; % time step (s)
dt = 100*dt; % physical time step stored every 100 time iterations

nt = 13; % number of post-processed variables
ne = 9; % number of energy variables

order = [3 1 2]; % dimension order

% for g=2.^(4:8)
for g=2^4
    tic
    gridname = ['Grid' num2str(g)];
    disp(gridname)
    gridpathname = fullfile(pathname,gridname);
    load(fullfile(gridpathname,'data.mat'),'N','n','m','p');
    
    fprintf('\nn = %d variables',n);
    fprintf('\n  = %d post-processed variables',nt);
    fprintf('\n  = %d energy variables',ne);
    fprintf('\nN = %d samples',N);
    fprintf('\nm = %d spatial points',m);
    fprintf('\np+1 = %d time steps',p+1);
    fprintf('\n');
    
    sx = [g+1,g+1,g+1]; % spatial dimensions
    
    % Spatial scheme
    dx = L/g; % spatial step (m)
    x = linspace(0,L,g+1); % spatial discretization in each spatial dimension
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
        Tau = zeros(N,nt,m,p+1);
        E = zeros(N,ne,m,p+1);
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
        % Ek = 1/2*rhot.*dot(ut,ut,2);
        Ek_old = 1/2*rhot_old.*dot(ut_old,ut_old,2);
        Ek_new = 1/2*rhot_new.*dot(ut_new,ut_new,2);
        energyKinTime = (Ek_new-Ek_old)/(2*dt);
        
%         Taut = zeros(N,nt,m);
%         Et = zeros(N,ne,m);
%         parfor l=1:N
%             % samplename = ['Sample ' num2str(l)];
%             % disp(samplename)
%             
%             ul = permute(reshape(ut(l,:,:),[3,sx]),[1,order+1]);
%             Cl = permute(reshape(Ct(l,:),sx),order);
%             
%             rhol = Cl*rho(2) + (1-Cl)*rho(1);
%             mul = Cl*mu(2) + (1-Cl)*mu(1);
%             gradul = grad(ul,Dx);
%             rhoul = shiftdim(rhol,-1).*ul;
%             gradrhoul = grad(rhoul,Dx);
%             Sl = (gradul+permute(gradul,[2,1,3:ndims(gradul)]))/2;
%             muSl = shiftdim(mul,-2).*Sl;
%             gradCl = grad(Cl,Dx);
%             ngradCl = normal(gradCl);
%             kappa = div(ngradCl,Dx);
%             u2l = dot(ul,ul,1);
%             rhou2l = shiftdim(rhol,-1).*u2l;
%             % tauConv = shiftdim(rhol,-1).*squeeze(sum(gradul.*shiftdim(ul,-1),2));
%             tauConv = squeeze(sum(gradrhoul.*shiftdim(ul,-1),2));
%             % tauDiff = 2*shiftdim(mul,-1).*div(Sl,Dx);
%             tauDiff = div(2*muSl,Dx);
%             tauSurf = sigma*shiftdim(kappa,-1).*gradCl;
%             tauInterf = dot(ul,gradCl,1);
%             energyConv = shiftdim(div(rhou2l.*ul,Dx),-1);
%             energyGrav = gravity.*rhoul(3,:,:,:);
%             energyPres = zeros(1,g+1,g+1,g+1);
%             energyPresDil = zeros(1,g+1,g+1,g+1);
%             energyKinSpace = dot(rhoul,grad(shiftdim(u2l/2),Dx),1);
%             energyDiff = shiftdim(div(squeeze(dot(2*muSl,repmat(shiftdim(ul,-1),[3,1]),2)),Dx),-1);
%             energyVisc = shiftdim(sum(sum(2*muSl.*gradul,1),2),1);
%             energySurf = dot(tauSurf,ul,1);
%             
%             Taul = cat(1,tauConv,tauDiff,tauSurf,tauInterf);
%             El = cat(1,energyConv,energyGrav,energyPres,energyPresDil,energyKinSpace,energyDiff,energyVisc,energySurf);
%             Taul = ipermute(Taul,[1,order+1]);
%             El = ipermute(El,[1,order+1]);
%             Taul = Taul(:,:);
%             El = El(:,:);
%             tauTimel = shiftdim(tauTime(l,:,:),1);
%             energyKinTimel = shiftdim(energyKinTime(l,:,:),1);
%             Taul = cat(1,tauTimel,Taul);
%             El = cat(1,energyKinTimel,El);
%             
%             Taut(l,:,:) = Taul;
%             Et(l,:,:) = El;
%         end
        
        ut = permute(reshape(ut,[N,3,sx]),[2,order+2,1]);
        Ct = permute(reshape(Ct,[N,sx]),[order+1,1]);
        
        rhot = Ct*rho(2) + (1-Ct)*rho(1);
        mut = Ct*mu(2) + (1-Ct)*mu(1);
        gradut = grad(ut,Dx);
        rhout = shiftdim(rhot,-1).*ut;
        gradrhout = grad(rhout,Dx);
        St = (gradut+permute(gradut,[2,1,3:ndims(gradut)]))/2;
        muSt = shiftdim(mut,-2).*St;
        gradCt = grad(Ct,Dx);
        ngradCt = normal(gradCt);
        kappa = div(ngradCt,Dx);
        u2t = dot(ut,ut,1);
        rhou2t = shiftdim(rhot,-1).*u2t;
        
        % tauConv = shiftdim(rhot,-1).*squeeze(sum(gradut.*shiftdim(ut,-1),2));
        tauConv = squeeze(sum(gradrhout.*shiftdim(ut,-1),2));
        % tauDiff = 2*shiftdim(mut,-1).*div(St,Dx);
        tauDiff = div(2*muSt,Dx);
        tauSurf = sigma*shiftdim(kappa,-1).*gradCt;
        tauInterf = dot(ut,gradCt,1);
        energyConv = shiftdim(div(rhou2t.*ut,Dx),-1);
        energyGrav = gravity.*rhout(3,:,:,:,:);
        energyPres = zeros(1,g+1,g+1,g+1,N);
        energyPresDil = zeros(1,g+1,g+1,g+1,N);
        energyKinSpace = dot(rhout,grad(shiftdim(u2t/2),Dx),1);
        energyDiff = shiftdim(div(squeeze(dot(2*muSt,repmat(shiftdim(ut,-1),[3,1]),2)),Dx),-1);
        energyVisc = shiftdim(sum(sum(2*muSt.*gradut,1),2),1);
        energySurf = dot(tauSurf,ut,1);
        
        Taut = cat(1,tauConv,tauDiff,tauSurf,tauInterf);
        Et = cat(1,energyConv,energyGrav,energyPres,energyPresDil,energyKinSpace,energyDiff,energyVisc,energySurf);
        Taut = ipermute(Taut,[2,order+2,1]);
        Et = ipermute(Et,[2,order+2,1]);
        Taut = Taut(:,:,:);
        Et = Et(:,:,:);
        Taut = cat(2,tauTime,Taut);
        Et = cat(2,energyKinTime,Et);

        if g<2^7
            Tau(:,:,:,t+1) = Taut;
            E(:,:,:,t+1) = Et;
        else
            save(fullfile(gridpathname,['tau_t' num2str(t) '.mat']),'Taut');
            save(fullfile(gridpathname,['energy_t' num2str(t) '.mat']),'Et');
        end
    end
    fprintf('\n');
    if g<2^7
        save(fullfile(gridpathname,'tau.mat'),'Tau');
        save(fullfile(gridpathname,'energy.mat'),'E');
    end
    save(fullfile(gridpathname,'tau.mat'),'nt','-append');
    save(fullfile(gridpathname,'energy.mat'),'ne','-append');
    toc
end

% myparallel('stop');
