clc
clearvars
close all

performPCA = false;
performPCAspace = true;
performPCAtime = true;
computeMean = true;
postProcess = true;
postProcessTau = true;
postProcessEnergy = false;
computeQoI = true;
applyFilter = true;
computeError = true;
constructMesh = true;

displayEigenvalues = false;
displayCovariance = false;
displayQoI = false;
displayError = false;

cmap = 'default';
% cmap = 'gray';
framerate = 5;
fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

pathname = fileparts(mfilename('fullpath'));

sigma = 0.45; % surface tension (N/m)
mu = [0.1,0.1]; % dynamic viscosity of [water,oil] (Pa.s)
rho = [1000,900]; % mass density of [water,oil] (kg/m3)
gravity = 9.81; % gravity (m/s2)

perm = @(u) permute(u,[2,3,4,5,1]);
iperm = @(u) ipermute(u,[2,3,4,5,1]);

tolsvdYc = eps; % relative precision for truncated SVD of Yc
tolsvdZc = eps; % relative precision for truncated SVD of Zc

index = 'coord'; % index for ordering ('coord', 'time')

filterType = 'box'; % 3D filter type ('box' or 'mean' or 'average', 'linear' or 'trapz')

% Spatial grid size
gset = 2.^(4:5); % set of spatial grid sizes
g = gset(end); % current spatial grid size
gref = gset(end); % reference spatial grid size
ng = length(gset); % number of spatial grid sizes

t_Total = tic;
fprintf('Grid %d\n',g)
gridname = ['Grid' num2str(g)];
gridpathname = fullfile(pathname,gridname);
load(fullfile(gridpathname,'data.mat'),'N','n','m','p');
r = n*m;

fprintf('\nn = %d variables',n);
fprintf('\nN = %d samples',N);
fprintf('\nm = %d spatial points',m);
fprintf('\np+1 = %d time steps',p+1);
fprintf('\n');

% Spatial scheme
L = 1; % domain size (m)
sx = [g+1,g+1,g+1]; % sx = nbelem+1; % spatial dimensions
dx = L/g; % spatial step (m)
x = linspace(0,L,g+1); % spatial discretization in each spatial dimension
Dx = spdiags(repmat([1 -1],g+1,1),[1 -1],g+1,g+1)/(2*dx); % Implicit central-difference spatial scheme (second-order accurate, unconditionlly stable)
Dx(1,[1 2]) = [-1 1]/dx; Dx(end,[end-1 end]) = [-1 1]/dx;

% Time scheme
dt = 5e-3; % computing time step (s)
dt = 100*dt; % physical time step stored every 100 computing time steps
% Dt = spdiags(repmat([1 -1],p+1,1),[0 -1],p+1,p+1)/dt; Dt(1,1) = 1/(2*dt); % Explicit backward Euler time scheme (first-order accurate, conditionally stable)
% Dt = spdiags(repmat([1 -1],p+1,1),[1 0],p+1,p+1)/dt; Dt(end,end) = 1/(2*dt); % Implicit forward Euler time scheme (first-order accurate, unconditionally stable)
Dt = spdiags(repmat([1 -1],p+1,1),[1 -1],p+1,p+1)/(2*dt); % Implicit central-difference time scheme (second-order accurate, unconditionally stable)
Dt(1,[1 2]) = [-1 1]/dt; Dt(end,[end-1 end]) = [-1 1]/dt;

if g<2^7
    fprintf('\nLoading DNS data');
    t_load = tic;
    load(fullfile(gridpathname,'data.mat'),'Y');
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
end

if performPCA
%% First reduction step in space
if performPCAspace
    fprintf('\nPerforming PCA in space');
    t_PCA_space = tic;
    Rinit = min(r,N);
    if g<2^7
        mY = mean(Y,1);
        Yc = Y - repmat(mY,[N,1,1,1]); % Yc = Y - mY.*ones(N,1,1,1);
        V = zeros(r,Rinit,p+1);
        clear Y
    else
        mY = zeros(1,n,m,p+1);
    end
    sig = zeros(Rinit,p+1);
    Z = zeros(N,Rinit,p+1);
    errsvdYc = zeros(Rinit,p+1);
    err2Yc = zeros(1,p+1);
    norm2Yc = zeros(1,p+1);
    R = zeros(1,p+1);
    for t=0:p
        t_PCA_spacet = tic;
        if g<2^7
            Yct = Yc(:,:,:,t+1);
        else
            load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
            mYt = mean(Yt,1);
            mY(1,:,:,t+1) = mYt;
            Yct = Yt - repmat(mYt,[N,1,1]); % Yct = Yt - mYt.*ones(N,1,1);
            clear Yt
        end
        Yct = Yct(:,:)';
        % if t==0
        %     [Vt,Sigt,Zt,errsvdYct] = svdtruncate(Yct,Rinit-1);
        % else
        [Vt,Sigt,Zt,errsvdYct] = svdtruncate(Yct,tolsvdYc);
        % end
        Sigt = Sigt/sqrt(N-1);
        sigt = diag(Sigt);
        Zt = Zt*sqrt(N-1);
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            Yct_approx = Vt*Sigt*Zt';
        else
            Yct_approx = Vt*(sigt.*Zt');
        end
        % errYct = norm(Yct_approx-Yct)/norm(Yct);
        errYct = errsvdYct(end);
        Rt = length(sigt);
        time_PCA_spacet = toc(t_PCA_spacet);
        fprintf('\nTime step %2d/%2d : rank R = %d, error = %.3e for Y, elapsed time = %f s',t,p,Rt,errYct,time_PCA_spacet);
        
        norm2Yct = sum(var(Yct,0,2));
        % norm2Yct_approx = sum(var(Yct_approx,0,2));
        err2Yct = errYct^2*norm2Yct;
        % err2Yct = norm2Yct-norm2Yct_approx;
        % sigtf = svdtruncate(Yct,eps);
        % sigtf = sigtf/sqrt(N-1);
        % err2Yct = sum(sigtf.^2)-sum(sigt.^2);
        
        % mZt = mean(Zt,1)';
        % CZt = cov(Zt); % CZt = 1/(N-1)*Zt'*Zt;
        % norm(mZt)
        % norm(CZt-eye(Rt))
        % norm(Vt'*Vt-eye(Rt))
        
        %if verLessThan('matlab','9.1') % compatibility (<R2016b)
        %    CYt_approx = Vt*Sigt.^2*Vt';
        %else
        %    CYt_approx = Vt*(sigt.^2.*Vt');
        %end
        %CYt = cov(Yct'); % CYt = 1/(N-1)*Yct*Yct';
        %errCYt = norm(CYt_approx-CYt)/norm(CYt);
        %fprintf('\n                                           error = %.3e for CY',errCYt);
        
        sig(1:Rt,t+1) = sigt;
        if g<2^7
            V(:,1:Rt,t+1) = Vt;
        else
            save(fullfile(gridpathname,['PCA_space_t' num2str(t) '.mat']),'Vt');
        end
        Z(:,1:Rt,t+1) = Zt;
        R(t+1) = Rt;
        errsvdYc(1:Rt,t+1) = errsvdYct;
        err2Yc(t+1) = err2Yct;
        norm2Yc(t+1) = norm2Yct;
    end
    
    ts = (0:p)*dt;
    errL2 = trapz(ts,err2Yc,2)/trapz(ts,norm2Yc,2);
    fprintf('\nL2-error = %.3e for Y',errL2);
    fprintf('\n');
    
    Rmax = max(R);
    sig = sig(1:Rmax,:);
    if g<2^7
        V = V(:,1:Rmax,:);
    end
    Z = Z(:,1:Rmax,:);
    errsvdYc = errsvdYc(1:Rmax,:);
    
    time_PCA_space = toc(t_PCA_space);
    fprintf('\nelapsed time = %f s',time_PCA_space);
    fprintf('\n');
    
    save(fullfile(gridpathname,'PCA_space.mat'),'mY','sig','Z','R','errsvdYc','time_PCA_space');
    if g<2^7
        save(fullfile(gridpathname,'PCA_space.mat'),'V','-append');
    end
else
    fprintf('\nLoading DNS data from PCA in space');
    t_load = tic;
    load(fullfile(gridpathname,'PCA_space.mat'),'mY','sig','Z','R','errsvdYc','time_PCA_space');
    if g<2^7
        load(fullfile(gridpathname,'PCA_space.mat'),'V');
    end
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
end

%% Second reduction step in time for each coordinate
% if performPCAtime
%     fprintf('\nPCA in time');
%     t_PCA_time = tic;
%     Q = min(p+1,N);
%     Rmax = max(R);
%     s = zeros(Q,Rmax);
%     W = zeros(p+1,Q,Rmax);
%     X = zeros(N,Q,Rmax);
%     errsvdZc = zeros(Q,Rmax);
%     err2Zc = zeros(Rmax,1);
%     Q = zeros(1,Rmax);
%     for a=1:Rmax
%         Zca = Z(:,a,:);
%         Zca = Zca(:,:)';
%         [Wa,Sa,Xa,errsvdZca] = svdtruncate(Zca,tolsvdZc);
%         Sa = Sa/sqrt(N-1);
%         sa = diag(Sa);
%         Xa = Xa*sqrt(N-1);
%         if verLessThan('matlab','9.1') % compatibility (<R2016b)
%             Zca_approx = Wa*Sa*Xa';
%         else
%             Zca_approx = Wa*(sa.*Xa');
%         end
%         % errZca = norm(Zca_approx-Zca)/norm(Zca);
%         errZca = errsvdZca(end);
%         Qa = length(sa);
%         fprintf('\nCoordinate alpha = %2d/%2d : rank Q = %d, error = %.3e for Z',a,Rmax,Qa,errZca);
%         
%         % norm2Zca = sum(var(Zca,0,2));
%         % norm2Zca_approx = sum(var(Zca_approx,0,2));
%         % err2Zca = errZca^2*norm2Zca;
%         % err2Zca = norm2Zca-norm2Zca_approx;
%         % saf = svdtruncate(Zca,eps);
%         % saf = saf/sqrt(N-1);
%         % err2Zca = sum(saf.^2)-sum(sa.^2);
%         
%         % mXa = mean(Xa,1)';
%         % CXa = cov(Xa); % CXa = 1/(N-1)*Xa'*Xa;
%         % norm(mXa)
%         % norm(CXa-eye(Qa))
%         % norm(Wa'*Wa-eye(Qa))
%         
%         if verLessThan('matlab','9.1') % compatibility (<R2016b)
%             CZa_approx = Wa*Sa.^2*Wa';
%         else
%             CZa_approx = Wa*(sa.^2.*Wa');
%         end
%         CZa = cov(Zca'); % CZa = 1/(N-1)*Zca*Zca';
%         errCZa = norm(CZa_approx-CZa)/norm(CZa);
%         fprintf('\n                                        error = %.3e for CZ',errCZa);
%         
%         if displayCovariance
%             figure('Name','Covariance matrix')
%             clf
%             imagesc(CZa)
%             colorbar
%             axis image
%             set(gca,'FontSize',fontsize)
%             xlabel('$k$','Interpreter','latex')
%             ylabel('$k''$','Interpreter','latex')
%             title(['Covariance matrix $[C_{\zeta}(t^k,t^{k''})]_{\alpha,\alpha} = [C_{Z_{\alpha}}]_{k,k''}$ for $\alpha=$' num2str(a)],'Interpreter','latex')
%             mysaveas(gridpathname,['covariance_CZ_a' num2str(a)],formats,renderer);
%             mymatlab2tikz(gridpathname,['covariance_CZ_a' num2str(a) '.tex']);
%         end
%         
%         s(1:Qa,a) = sa;
%         W(:,1:Qa,a) = Wa;
%         X(:,1:Qa,a) = Xa;
%         Q(a) = Qa;
%         errsvdZc(1:Qa,a) = errsvdZca;
%     end
%     
%     Q = max(Q);
%     s = s(1:Q,:);
%     W = W(:,1:Q,:);
%     X = X(:,1:Q,:);
%     errsvdZc = errsvdZc(1:Q,:);
%     
%     time_PCA_time = toc(t_PCA_time);
%     fprintf('\nelapsed time = %f s',time_PCA_time);
%     fprintf('\n');
%     
%     save(fullfile(gridpathname,'PCA_time.mat'),'s','W','Q','errsvdZc','time_PCA_time');
% else
%     fprintf('\nLoading DNS data from PCA in time');
%     t_load = tic;
%     load(fullfile(gridpathname,'PCA_time.mat'),'s','W','Q','errsvdZc','time_PCA_time');
%     time_load = toc(t_load);
%     fprintf('\nelapsed time = %f s',time_load);
%     fprintf('\n');
% end
    
%% Second reduction step in time
if performPCAtime
    fprintf('\nPerforming PCA in time');
    t_PCA_time = tic;
    Rmax = max(R);
    q = (p+1)*Rmax;
    switch index
        case 'coord'
            Zc = Z(:,:)';
        case 'time'
            Zc = permute(Z,[1,3,2]);
            Zc = Zc(:,:)';
    end
    [W,S,X,errsvdZc] = svdtruncate(Zc,tolsvdZc);
    S = S/sqrt(N-1);
    s = diag(S);
    X = X*sqrt(N-1);
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        Zc_approx = W*S*X';
    else
        Zc_approx = W*(s.*X');
    end
    % errZc = norm(Zc_approx-Zc)/norm(Zc);
    errZc = errsvdZc(end);
    Q = length(s);
    fprintf('\nrank R = %d, rank Q = %d, error = %.3e for Z',Rmax,Q,errZc);
    
    % norm2Zc = sum(var(Zc,0,2));
    % norm2Zc_approx = sum(var(Zc_approx,0,2));
    % err2Zc = errZc^2*norm2Zc;
    % err2Zc = norm2Zc-norm2Zc_approx;
    % sf = svdtruncate(Zc,eps);
    % sf = sf/sqrt(N-1);
    % err2Zc = sum(sf.^2)-sum(s.^2);
    
    % mX = mean(X,1)';
    % CX = cov(X); % CX = 1/(N-1)*X'*X;
    % norm(mX)
    % norm(CX-eye(Q))
    % norm(W'*W-eye(Q))
    % norm(CZ_approx*W-(p+1)*W)
    % norm(CZ*W-(p+1)*W)
    % norm(abs(s.^2-(p+1)))
    
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        CZ_approx = W*S.^2*W';
    else
        CZ_approx = W*(s.^2.*W');
    end
    CZ = cov(Zc'); % CZ = 1/(N-1)*Zc*Zc';
    errCZ = norm(CZ_approx-CZ)/norm(CZ);
    fprintf('\n                          error = %.3e for CZ',errCZ);
    
    if displayCovariance
        figure('Name','Covariance matrix')
        clf
        imagesc(CZ)
        colorbar
        axis image
        set(gca,'FontSize',fontsize)
        switch index
            case 'coord'
                xlabel('$K=(k-1)R+\alpha$','Interpreter','latex')
                ylabel('$K''=(k''-1)R+\alpha''$','Interpreter','latex')
            case 'time'
                xlabel('$K=\alpha p+k$','Interpreter','latex')
                ylabel('$K''=\alpha'' p+k''$','Interpreter','latex')
        end
        title('Covariance matrix $[C_{\zeta}(t^k,t^{k''})]_{\alpha,\alpha''} = [C_Z]_{K,K''}$','Interpreter','latex')
        mysaveas(gridpathname,'covariance_CZ',formats,renderer);
        mymatlab2tikz(gridpathname,'covariance_CZ.tex');
        
        %     for t=0:p
        %         switch index
        %             case 'coord'
        %                 ind = t*Rmax+(1:Rmax);
        %             case 'time'
        %                 ind = (0:Rmax-1)*(p+1)+t+1;
        %         end
        %         CZt_approx = CZ_approx(ind,ind);
        %         CZt = CZ(ind,ind);
        %
        %         figure('Name','Covariance matrix')
        %         clf
        %         imagesc(CZt)
        %         colorbar
        %         axis image
        %         set(gca,'FontSize',fontsize)
        %         xlabel('$\alpha$','Interpreter','latex')
        %         ylabel('$\alpha''$','Interpreter','latex')
        %         title(['Covariance matrix $[C_{\zeta}(t^k,t^k)]_{\alpha,\alpha''} = [C_{Z_k}]_{\alpha,\alpha''}$ for $k=$' num2str(t)],'Interpreter','latex')
        %         mysaveas(gridpathname,['covariance_CZ_t' num2str(t*100)],formats,renderer);
        %         mymatlab2tikz(gridpathname,['covariance_CZ_t' num2str(t*100) '.tex']);
        %     end
    end
    
    time_PCA_time = toc(t_PCA_time);
    fprintf('\nelapsed time = %f s',time_PCA_time);
    fprintf('\n');
    
    save(fullfile(gridpathname,'PCA_time.mat'),'s','W','Q','errsvdZc','time_PCA_time');
else
    fprintf('\nLoading DNS data from PCA time');
    t_load = tic;
    load(fullfile(gridpathname,'PCA_time.mat'),'s','W','Q','errsvdZc','time_PCA_time');
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
end
else
    if computeMean
        fprintf('\nComputing mean DNS data');
        t_Mean = tic;
        
        if g<2^7
            mY = mean(Y,1);
        else
            mY = zeros(1,n,m,p+1);
            for t=0:p
                t_Meant = tic;
                load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
                mYt = mean(Yt,1);
                mY(1,:,:,t+1) = mYt;
                time_Meant = toc(t_Meant);
                fprintf('\nTime step %2d/%2d : elapsed time = %f s',t,p,time_Meant);
            end
        end
        
        time_Mean = toc(t_Mean);
        fprintf('\nelapsed time = %f s',time_Mean);
        fprintf('\n');
        
        save(fullfile(gridpathname,'mean_data.mat'),'mY','time_Mean');
    else
        fprintf('\nLoading mean DNS data');
        t_load = tic;
        load(fullfile(gridpathname,'mean_data.mat'),'mY','time_Mean');
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
    end
end

if postProcess
%% Post-processing DNS data
fprintf('\nPost-processing DNS data');
t_PostProcess = tic;
if postProcessTau
    ntau = 13; % number of tau variables
    fprintf('\nn = %d tau variables',ntau);
    mTau = zeros(1,ntau,m,p+1);
    if g<2^7
        Tau = zeros(N,ntau,m,p+1);
    end
end
if postProcessEnergy
    ne = 9; % number of energy variables
    fprintf('\n  = %d energy variables',ne);
    mE = zeros(1,ne,m,p+1);
    if g<2^7
        E = zeros(N,ne,m,p+1);
    end
end

for t=0:p
    t_PostProcesst = tic;
    
    if performPCA
        if g<2^7
            if t>0
                Yt_old = get_reduced_order_representation_at_step(mY,sig,V,s,W,X,R,t,'index',index,'gridpathname',gridpathname);
            end
            if t<p
                Yt_new = get_reduced_order_representation_at_step(mY,sig,V,s,W,X,R,t+2,'index',index,'gridpathname',gridpathname);
            end
            Yt = get_reduced_order_representation_at_step(mY,sig,V,s,W,X,R,t+1,'index',index,'gridpathname',gridpathname);
        else
            if t>0
                Yt_old = get_reduced_order_representation_at_step(mY,sig,[],s,W,X,R,t,'index',index,'gridpathname',gridpathname);
            end
            if t<p
                Yt_new = get_reduced_order_representation_at_step(mY,sig,[],s,W,X,R,t+2,'index',index,'gridpathname',gridpathname);
            end
            Yt = get_reduced_order_representation_at_step(mY,sig,[],s,W,X,R,t+1,'index',index,'gridpathname',gridpathname);
        end
    else
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
    end
    
    Yt = perm(reshape(Yt,[N,n,sx]));
    ut = Yt(1:3,:,:,:,:);
    Ct = Yt(4,:,:,:,:);
    clear Yt
    rhot = Ct*rho(2) + (1-Ct)*rho(1);
    rhout = repmat(rhot,[3,ones(1,4)]).*ut;
    if ~postProcessEnergy
        clear rhot
    end
    if postProcessEnergy
        u2t = dot(ut,ut,1);
        if t==0 || t==p
            Ek = 1/2*rhot.*u2t;
        end
    end
    if t>0
        Yt_old = perm(reshape(Yt_old,[N,n,sx]));
        ut_old = Yt_old(1:3,:,:,:,:);
        Ct_old = Yt_old(4,:,:,:,:);
        clear Yt_old
        rhot_old = Ct_old*rho(2) + (1-Ct_old)*rho(1);
        if postProcessTau
            rhout_old = repmat(rhot_old,[3,ones(1,4)]).*ut_old;
        end
        if postProcessEnergy
            Ek_old = 1/2*rhot_old.*dot(ut_old,ut_old,1);
        end
        clear ut_old Ct_old rhot_old
    end
    if t<p
        Yt_new = perm(reshape(Yt_new,[N,n,sx]));
        ut_new = Yt_new(1:3,:,:,:,:);
        Ct_new = Yt_new(4,:,:,:,:);
        clear Yt_new
        rhot_new = Ct_new*rho(2) + (1-Ct_new)*rho(1);
        if postProcessTau
            rhout_new = repmat(rhot_new,[3,ones(1,4)]).*ut_new;
        end
        if postProcessEnergy
            Ek_new = 1/2*rhot_new.*dot(ut_new,ut_new,1);
        end
        clear ut_new Ct_new rhot_new
    end
    
    gradut = grad(ut,Dx);
    St = (gradut+permute(gradut,[2,1,3:6]))/2;
    if ~postProcessEnergy
        clear gradut
    end
    mut = Ct*mu(2) + (1-Ct)*mu(1);
    muSt = repmat(shiftdim(mut,-1),[3,3,ones(1,4)]).*St;
    clear mut St
    gradCt = grad(Ct,Dx);
    clear Ct
    ngradCt = normal(gradCt);
    kappa = div(ngradCt,Dx);
    clear ngradCt
    
    tauSurft = sigma*repmat(shiftdim(kappa,-1),[3,ones(1,4)]).*gradCt;
    clear kappa
    if ~postProcessTau
        clear gradCt
    end
    if postProcessTau
        if t==0
            tauTimet = (rhout_new-rhout)/dt;
            clear rhout_new
        elseif t==p
            tauTimet = (rhout-rhout_old)/dt;
            clear rhout_old
        else
            tauTimet = (rhout_new-rhout_old)/(2*dt);
            clear rhout_new rhout_old
        end
        divtauConvt = squeeze(sum(grad(rhout,Dx).*repmat(shiftdim(ut,-1),[3,ones(1,5)]),2));
        % divtauConvt = div(permute(repmat(shiftdim(rhout,-1),[3,ones(1,5)]),[2,1,3:6]).*repmat(shiftdim(ut,-1),[3,ones(1,5)]),Dx);
        tauInterft = dot(ut,gradCt,1);
        clear gradCt
        divtauDifft = div(2*muSt,Dx);
        Taut = cat(1,tauTimet,divtauConvt,divtauDifft,tauSurft,tauInterft);
        clear tauTimet divtauConvt divtauDifft tauInterft
        if ~postProcessEnergy
            clear ut rhot rhout muSt tauSurft
        end
        Taut = iperm(Taut);
        Taut = Taut(:,:,:);
        mTaut = mean(Taut,1);
        mTau(1,:,:,t+1) = mTaut;
        clear mTaut
        if g<2^7
            Tau(:,:,:,t+1) = Taut;
        else
            save(fullfile(gridpathname,['data_tau_t' num2str(t) '.mat']),'Taut');
        end
        clear Taut
    end
    if postProcessEnergy
        if t==0
            energyKinTimet = (Ek_new-Ek)/dt;
            clear Ek Ek_new
        elseif t==p
            energyKinTimet = (Ek-Ek_old)/dt;
            clear Ek Ek_old
        else
            energyKinTimet = (Ek_new-Ek_old)/(2*dt);
            clear Ek_new Ek_old
        end
        energyConvt = shiftdim(div(repmat(rhot.*u2t,[3,ones(1,4)]).*ut,Dx),-1);
        clear rhot
        energyKinSpacet = dot(rhout,grad(u2t/2,Dx),1);
        clear u2t
        energyGravt = gravity.*rhout(2,:,:,:,:);
        clear rhout
        energyPrest = zeros(1,g+1,g+1,g+1,N);
        energyPresDilt = zeros(1,g+1,g+1,g+1,N);
        energyDifft = shiftdim(div(squeeze(dot(2*muSt,repmat(shiftdim(ut,-1),[3,ones(1,5)]),2)),Dx),-1);
        energyVisct = shiftdim(sum(sum(2*muSt.*gradut,1),2),1);
        clear gradut muSt
        energySurft = dot(tauSurft,ut,1);
        clear ut tauSurft
        Et = cat(1,energyKinTimet,energyConvt,energyGravt,energyPrest,energyPresDilt,energyKinSpacet,energyDifft,energyVisct,energySurft);
        clear energyKinTimet energyConvt energyGravt energyPrest energyPresDilt energyKinSpacet energyDifft energyVisct energySurft
        Et = iperm(Et);
        Et = Et(:,:,:);
        mEt = mean(Et,1);
        mE(:,:,:,t+1) = mEt;
        clear mEt
        if g<2^7
            E(:,:,:,t+1) = Et;
        else
            save(fullfile(gridpathname,['data_energy_t' num2str(t) '.mat']),'Et');
        end
        clear Et
    end
    time_PostProcesst = toc(t_PostProcesst);
    fprintf('\nTime step %2d/%2d : elapsed time = %f s',t,p,time_PostProcesst);
end

time_PostProcess = toc(t_PostProcess);
fprintf('\nelapsed time = %f s',time_PostProcess);
fprintf('\n');

if postProcessTau
    save(fullfile(gridpathname,'mean_data_tau.mat'),'mTau','ntau');
    if g<2^7
        save(fullfile(gridpathname,'data_tau.mat'),'Tau');
    end
end
if postProcessEnergy
    save(fullfile(gridpathname,'mean_data_energy.mat'),'mE','ne');
    if g<2^7
        save(fullfile(gridpathname,'data_energy.mat'),'E');
    end
end
else
    if postProcessTau
        if g<2^7
            fprintf('\nLoading tau data');
        else
            fprintf('\nLoading mean tau data');
        end
        t_load = tic;
        load(fullfile(gridpathname,'mean_data_tau.mat'),'mTau','ntau');
        if g<2^7
            load(fullfile(gridpathname,'data_tau.mat'),'Tau');
        end
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
    end
    if postProcessEnergy
        if g<2^7
            fprintf('\nLoading energy data');
        else
            fprintf('\nLoading mean energy data');
        end
        t_load = tic;
        load(fullfile(gridpathname,'mean_data_energy.mat'),'mE','ne');
        if g<2^7
            load(fullfile(gridpathname,'data_energy.mat'),'E');
        end
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
    end
end

%% Compute quantities of interest: spatial average in each phase
if computeQoI
    fprintf('\nComputing quantities of interest');
    t_QoI = tic;
    
    Qu = zeros(N,3,2,p+1);
    if postProcessTau
        Qtau = zeros(N,ntau,2,p+1);
    end
    if postProcessEnergy
        Qe = zeros(N,ne,2,p+1);
    end
    
    for t=0:p
        t_QoIt = tic;
        
        if performPCA
            if g<2^7
                Yt = get_reduced_order_representation_at_step(mY,sig,V,s,W,X,R,t+1,'index',index,'gridpathname',gridpathname);
            else
                Yt = get_reduced_order_representation_at_step(mY,sig,[],s,W,X,R,t+1,'index',index,'gridpathname',gridpathname);
            end
        else
            if g<2^7
                Yt = Y(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
            end
        end
        
        Yt = reshape(Yt,[N,n,sx]);
        ut = Yt(:,1:3,:,:,:);
        Ct = Yt(:,4,:,:,:);
        clear Yt
        Qut = cat(3,trapz(x,trapz(x,trapz(x,repmat(1-Ct,[1,3,1,1,1]).*ut,3),4),5),...
            trapz(x,trapz(x,trapz(x,repmat(Ct,[1,3,1,1,1]).*ut,3),4),5));
        Qu(:,:,:,t+1) = Qut;
        clear ut Qut
        
        if postProcessTau
            if g<2^7
                Taut = Tau(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_tau_t' num2str(t) '.mat']),'Taut');
            end
            Taut = reshape(Taut,[N,ntau,sx]);
            Qtaut = cat(3,trapz(x,trapz(x,trapz(x,repmat(1-Ct,[1,ntau,1,1,1]).*Taut,3),4),5),...
                trapz(x,trapz(x,trapz(x,repmat(Ct,[1,ntau,1,1,1]).*Taut,3),4),5));
            Qtau(:,:,:,t+1) = Qtaut;
            clear Taut Qtaut
        end
        
        if postProcessEnergy
            if g<2^7
                Et = E(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_energy_t' num2str(t) '.mat']),'Et');
            end
            Et = reshape(Et,[N,ne,sx]);
            Qet = cat(3,trapz(x,trapz(x,trapz(x,repmat(1-Ct,[1,ne,1,1,1]).*Et,3),4),5),...
                trapz(x,trapz(x,trapz(x,repmat(Ct,[1,ne,1,1,1]).*Et,3),4),5));
            Qe(:,:,:,t+1) = Qet;
            clear Et Qet
        end
        clear Ct
        time_QoIt = toc(t_QoIt);
        fprintf('\nTime step %2d/%2d : elapsed time = %f s',t,p,time_QoIt);
    end
    
    [mQu,stdQu,RQu,IQu] = compute_stats(Qu,dt);
    if postProcessTau
        [mQtau,stdQtau,RQtau,IQtau] = compute_stats(Qtau,dt);
    end
    if postProcessEnergy
        [mQe,stdQe,RQe,IQe] = compute_stats(Qe,dt);
    end
    
    time_QoI = toc(t_QoI);
    fprintf('\nelapsed time = %f s',time_QoI);
    fprintf('\n');
    
    save(fullfile(gridpathname,'data_qoi.mat'),'Qu','mQu','IQu');
    if postProcessTau
        save(fullfile(gridpathname,'data_qoi_tau.mat'),'Qtau','mQtau','IQtau');
    end
    if postProcessEnergy
        save(fullfile(gridpathname,'data_qoi_energy.mat'),'Qe','mQe','IQe');
    end
else
    fprintf('\nLoading quantities of interest');
    t_load = tic;
    load(fullfile(gridpathname,'data_qoi.mat'),'Qu','mQu','IQu');
    if postProcessTau
        load(fullfile(gridpathname,'data_qoi_tau.mat'),'Qtau','mQtau','IQtau');
    end
    if postProcessEnergy
        load(fullfile(gridpathname,'data_qoi_energy.mat'),'Qe','mQe','IQe');
    end
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
end

%% Applying filter
if g==gref
    if applyFilter
        fprintf('\nConstructing filter');
        hset = cell(1,ng-1);
        for ig=1:ng-1
            switch filterType
                case {'box','mean','average'}
                    filterSize = [2^ig+1 2^ig+1 2^ig+1];
                    h = ones(filterSize)/prod(filterSize);
                case {'linear','trapz'}
                    filterSize = [3 3 3];
                    c = 4;
                    h = zeros(filterSize);
                    h(:,:,1) = [1 c 1; c c^2 c; 1 c 1];
                    h(:,:,2) = [c c^2 c; c^2 c^3 c^2; c c^2 c];
                    h(:,:,3) = [1 c 1; c c^2 c; 1 c 1];
                    h = h/prod(2*filterSize);
            end
            hset{ig} = h;
        end
        
        fprintf('\nApplying filter for DNS data');
        t_Filter = tic;
        
        mYrefbar = zeros(ng-1,n,m,p+1);
        for ig=1:ng-1
            gbar = gset(ng-ig);
            gridbarname = ['Grid' num2str(gbar)];
            gridbarpathname = fullfile(pathname,gridbarname);
            mbar = (gbar+1)^3;
            k = g/gbar;
            h = hset{ig};
            H = shiftdim(h,-2);
            
            fprintf('\nFiltering DNS data on coarse grid %d',gbar);
            t_Filterg = tic;
            
            mYbar = zeros(1,n,mbar,p+1);
            if g<2^7
                Ybar = zeros(N,n,mbar,p+1);
            end
            for t=0:p
                t_Filtergt = tic;
                if performPCA
                    if g<2^7
                        Yt = get_reduced_order_representation_at_step(mY,sig,V,s,W,X,R,t+1,'index',index,'gridpathname',gridpathname);
                    else
                        Yt = get_reduced_order_representation_at_step(mY,sig,[],s,W,X,R,t+1,'index',index,'gridpathname',gridpathname);
                    end
                else
                    if g<2^7
                        Yt = Y(:,:,:,t+1);
                    else
                        load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
                    end
                end
                
                Ybart = reshape(Yt,[N,n,sx]);
                clear Yt
                % Ybart = apply_filter(Ybart,filterType,h);
                Ybart = imfilter(Ybart,H,'replicate');
                
                time_Filtergt = toc(t_Filtergt);
                fprintf('\nTime %2d/%2d : elapsed time = %f s',t,p,time_Filtergt);
                
                Yrefbart = Ybart(:,:,:);
                mYrefbar(ig,:,:,t+1) = mean(Yrefbart,1);
                clear Yrefbart
                
                Ybart = reshape(Ybart(:,:,1:k:end,1:k:end,1:k:end),[N,n,mbar]);
                mYbar(1,:,:,t+1) = mean(Ybart,1);
                if gbar<2^7
                    Ybar(:,:,:,t+1) = Ybart;
                else
                    save(fullfile(gridbarpathname,['data_filtered_t' num2str(t) '.mat']),'Ybart');
                end
                clear Ybart
            end
            
            time_Filterg = toc(t_Filterg);
            fprintf('\nelapsed time = %f s for coarse grid %d',time_Filterg,gbar);
            
            % save(fullfile(gridpathname,['mean_data_filtered_grid' num2str(gbar) '.mat']),'mYrefbar');
            % clear mYrefbar
            save(fullfile(gridbarpathname,'mean_data_filtered.mat'),'mYbar');
            clear mYbar
            if gbar<2^7
                save(fullfile(gridbarpathname,'data_filtered.mat'),'Ybar');
                clear Ybar
            end
        end
        time_Filter = toc(t_Filter);
        fprintf('\nelapsed time = %f s for all coarse grids',time_Filter);
        fprintf('\n');
        
        save(fullfile(gridpathname,'mean_data_filtered.mat'),'mYrefbar');
        
        if postProcessTau
            fprintf('\nApplying filter for tau data');
            t_Filter = tic;
            
            mTaurefbar = zeros(ng-1,ntau,m,p+1);
            for ig=1:ng-1
                gbar = gset(ng-ig);
                gridbarname = ['Grid' num2str(gbar)];
                gridbarpathname = fullfile(pathname,gridbarname);
                mbar = (gbar+1)^3;
                k = g/gbar;
                h = hset{ig};
                H = shiftdim(h,-2);
                
                fprintf('\nFiltering tau data on coarse grid %d',gbar);
                t_Filterg = tic;
                
                mTaubar = zeros(1,ntau,mbar,p+1);
                if g<2^7
                    Taubar = zeros(N,ntau,mbar,p+1);
                end
                for t=0:p
                    t_Filtergt = tic;
                    if g<2^7
                        Taut = Tau(:,:,:,t+1);
                    else
                        load(fullfile(gridpathname,['data_tau_t' num2str(t) '.mat']),'Taut');
                    end
                    
                    Taubart = reshape(Taut,[N,ntau,sx]);
                    clear Taut
                    % Taubart = apply_filter(Taubart,filterType,h);
                    Taubart = imfilter(Taubart,H,'replicate');
                    
                    time_Filtergt = toc(t_Filtergt);
                    fprintf('\nTime %2d/%2d : elapsed time = %f s',t,p,time_Filtergt);
                    
                    Taurefbart = Taubart(:,:,:);
                    mTaurefbar(ig,:,:,t+1) = mean(Taurefbart,1);
                    clear Taurefbart
                    
                    Taubart = reshape(Taubart(:,:,1:k:end,1:k:end,1:k:end),[N,ntau,mbar]);
                    mTaubar(1,:,:,t+1) = mean(Taubart,1);
                    if gbar<2^7
                        Taubar(:,:,:,t+1) = Taubart;
                    else
                        save(fullfile(gridbarpathname,['data_tau_filtered_t' num2str(t) '.mat']),'Taubart');
                    end
                    clear Taubart
                end
                
                time_Filterg = toc(t_Filterg);
                fprintf('\nelapsed time = %f s for coarse grid %d',time_Filterg,gbar);
                
                % save(fullfile(gridpathname,['mean_data_tau_filtered_grid' num2str(gbar) '.mat']),'mTaurefbar');
                % clear mTaurefbar
                save(fullfile(gridbarpathname,'mean_data_tau_filtered.mat'),'mTaubar');
                clear mTaubar
                if gbar<2^7
                    save(fullfile(gridbarpathname,'data_tau_filtered.mat'),'Taubar');
                    clear Taubar
                end
            end
            time_Filter = toc(t_Filter);
            fprintf('\nelapsed time = %f s for all coarse grids',time_Filter);
            fprintf('\n');
            
            save(fullfile(gridpathname,'mean_data_tau_filtered.mat'),'mTaurefbar');
        end
    else
        fprintf('\nLoading mean filtered DNS data');
        t_load = tic;
        load(fullfile(gridpathname,'mean_data_filtered.mat'),'mYrefbar');
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
        if postProcessTau
            fprintf('\nLoading mean filtered tau data');
            t_load = tic;
            load(fullfile(gridpathname,'mean_data_tau_filtered.mat'),'mTaurefbar');
            time_load = toc(t_load);
            fprintf('\nelapsed time = %f s',time_load);
            fprintf('\n');
        end
    end
else
    %% Loading filtered DNS data
    if g<2^7
        fprintf('\nLoading filtered DNS data');
    else
        fprintf('\nLoading mean filtered DNS data');
    end
    t_load = tic;
    load(fullfile(gridpathname,'mean_data_filtered.mat'),'mYbar');
    if g<2^7
        load(fullfile(gridpathname,'data_filtered.mat'),'Ybar');
    end
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
    
    %% Loading filtered tau data
    if postProcessTau
        if g<2^7
            fprintf('\nLoading filtered tau data');
        else
            fprintf('\nLoading mean filtered tau data');
        end
        t_load = tic;
        load(fullfile(gridpathname,'mean_data_tau_filtered.mat'),'mTaubar');
        if g<2^7
            load(fullfile(gridpathname,'data_tau_filtered.mat'),'Taubar');
        end
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
    end
        
    %% Computing error
    if computeError
        fprintf('\nComputing error between DNS and filtered data');
        t_Error = tic;
        
        erroru = zeros(1,p+1);
        errorC = zeros(1,p+1);
        normubar = zeros(1,p+1);
        normCbar = zeros(1,p+1);
        
        if postProcessTau
            errortauTime = zeros(1,p+1);
            errordivtauConv = zeros(1,p+1);
            errordivtauDiff = zeros(1,p+1);
            errortauSurf = zeros(1,p+1);
            errortauInterf = zeros(1,p+1);
            normtauTimebar = zeros(1,p+1);
            normdivtauConvbar = zeros(1,p+1);
            normdivtauDiffbar = zeros(1,p+1);
            normtauSurfbar = zeros(1,p+1);
            normtauInterfbar = zeros(1,p+1);
        end
        
        for t=0:p
            t_Errort = tic;
            
            if g<2^7
                Yt = Y(:,:,:,t+1);
                Ybart = Ybar(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
                load(fullfile(gridpathname,['data_filtered_t' num2str(t) '.mat']),'Ybart');
            end
            ut = Yt(:,1:3,:);
            Ct = Yt(:,4,:);
            clear Yt
            ubart = Ybart(:,1:3,:);
            Cbart = Ybart(:,4,:);
            clear Ybart
            
            errorut = sqrt(mean(trapz(trapz(trapz(x,reshape(sum((ubart-ut).^2,2),[N,sx]),2),3),4),1));
            errorCt = sqrt(mean(trapz(trapz(trapz(x,reshape((Cbart-Ct).^2,[N,sx]),2),3),4),1));
            clear ut Ct
            
            normubart = sqrt(mean(trapz(trapz(trapz(x,reshape(sum(ubart.^2,2),[N,sx]),2),3),4),1));
            normCbart = sqrt(mean(trapz(trapz(trapz(x,reshape(Cbart.^2,[N,sx]),2),3),4),1));
            clear ubart Cbart
            
            erroru(t+1) = errorut;
            errorC(t+1) = errorCt;
            clear errorut errorCt
            
            normubar(t+1) = normubart;
            normCbar(t+1) = normCbart;
            clear normubart normCbart
            
            if postProcessTau
                if g<2^7
                    Taut = Tau(:,:,:,t+1);
                    Taubart = Taubar(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,['data_tau_t' num2str(t) '.mat']),'Taut');
                    load(fullfile(gridpathname,['data_tau_filtered_t' num2str(t) '.mat']),'Taubart');
                end
                tauTimet = Taut(:,1:3,:);
                divtauConvt = Taut(:,4:6,:);
                divtauDifft = Taut(:,7:9,:);
                tauSurft = Taut(:,10:12,:);
                tauInterft = Taut(:,13,:);
                clear Taut
                tauTimebart = Taubart(:,1:3,:);
                divtauConvbart = Taubart(:,4:6,:);
                divtauDiffbart = Taubart(:,7:9,:);
                tauSurfbart = Taubart(:,10:12,:);
                tauInterfbart = Taubart(:,13,:);
                clear Taubart
                
                errortauTimet = sqrt(mean(trapz(trapz(trapz(x,reshape(sum((tauTimebart-tauTimet).^2,2),[N,sx]),2),3),4),1));
                errordivtauConvt = sqrt(mean(trapz(trapz(trapz(x,reshape(sum((divtauConvbart-divtauConvt).^2,2),[N,sx]),2),3),4),1));
                errordivtauDifft = sqrt(mean(trapz(trapz(trapz(x,reshape(sum((divtauDiffbart-divtauDifft).^2,2),[N,sx]),2),3),4),1));
                errortauSurft = sqrt(mean(trapz(trapz(trapz(x,reshape(sum((tauSurfbart-tauSurft).^2,2),[N,sx]),2),3),4),1));
                errortauInterft = sqrt(mean(trapz(trapz(trapz(x,reshape((tauInterfbart-tauInterft).^2,[N,sx]),2),3),4),1));
                clear tauTimet divtauConvt divtauDifft tauSurft tauInterft
                
                normtauTimebart = sqrt(mean(trapz(trapz(trapz(x,reshape(sum(tauTimebart.^2,2),[N,sx]),2),3),4),1));
                normdivtauConvbart = sqrt(mean(trapz(trapz(trapz(x,reshape(sum(divtauConvbart.^2,2),[N,sx]),2),3),4),1));
                normdivtauDiffbart = sqrt(mean(trapz(trapz(trapz(x,reshape(sum(divtauDiffbart.^2,2),[N,sx]),2),3),4),1));
                normtauSurfbart = sqrt(mean(trapz(trapz(trapz(x,reshape(sum(tauSurfbart.^2,2),[N,sx]),2),3),4),1));
                normtauInterfbart = sqrt(mean(trapz(trapz(trapz(x,reshape(tauInterfbart.^2,[N,sx]),2),3),4),1));
                clear tauTimebart divtauConvbart divtauDiffbart tauSurfbart tauInterfbart
                
                errortauTime(t+1) = errortauTimet;
                errordivtauConv(t+1) = errordivtauConvt;
                errordivtauDiff(t+1) = errordivtauDifft;
                errortauSurf(t+1) = errortauSurft;
                errortauInterf(t+1) = errortauInterft;
                clear errortauTimet errordivtauConvt errordivtauDifft errortauSurft errortauInterft
                
                normtauTimebar(t+1) = normtauTimebart;
                normdivtauConvbar(t+1) = normdivtauConvbart;
                normdivtauDiffbar(t+1) = normdivtauDiffbart;
                normtauSurfbar(t+1) = normtauSurfbart;
                normtauInterfbar(t+1) = normtauInterfbart;
                clear normtauTimebart normdivtauConvbart normdivtauDiffbart normtauSurfbart normtauInterfbart
            end
            time_Errort = toc(t_Errort);
            fprintf('\nTime %2d/%2d : elapsed time = %f s',t,p,time_Errort);
        end
        
        time_Error = toc(t_Error);
        fprintf('\nelapsed time = %f s',time_Error);
        fprintf('\n');
        
        save(fullfile(gridpathname,'error_filter.mat'),'erroru','errorC','normubar','normCbar','normubar');
        if postProcessTau
            save(fullfile(gridpathname,'error_filter_tau.mat'),'errortauTime','errordivtauConv','errordivtauDiff','errortauSurf','errortauInterf',...
                'normtauTimebar','normdivtauConvbar','normdivtauDiffbar','normtauSurfbar','normtauInterfbar');
        end
    else
        fprintf('\nLoading error between DNS and filtered data');
        t_load = tic;
        load(fullfile(gridpathname,'error_filter.mat'),'erroru','errorC','normubar','normCbar','normubar');
        if postProcessTau
            load(fullfile(gridpathname,'error_filter_tau.mat'),'errortauTime','errordivtauConv','errordivtauDiff','errortauSurf','errortauInterf',...
            'normtauTimebar','normdivtauConvbar','normdivtauDiffbar','normtauSurfbar','normtauInterfbar');
        end
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
    end
end

%% Outputs
% Display eigenvalues
if displayEigenvalues && performPCA
    figure('Name','Evolution of eigenvalues w.r.t order at each time')
    clf
    nCols = 5;
    leg = cell(1,p+1);
    hdl = zeros(1,p+1);
    color = distinguishable_colors(p+1);
    c = 0;
    for t=0:p
        c = c+1;
        Rt = R(t+1);
        hdl(c) = semilogy(1:Rt,sig(1:Rt,t+1).^2,'LineStyle','-','Color',color(c,:),'LineWidth',1);
        leg{c} = ['$t = ' num2str(t*dt) '$ s'];
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$\alpha$','Interpreter','latex')
    ylabel('$\lambda_{\alpha}(t)$','Interpreter','latex')
    gridLegend(hdl,nCols,leg,'Location','BestOutside','Interpreter','latex');
    % legend(leg{:},'Location','NorthEastOutside')
    mysaveas(gridpathname,'eigenvalues_CY_order',formats,renderer);
    mymatlab2tikz(gridpathname,'eigenvalues_CY_order.tex');
    
    figure('Name','Evolution of eigenvalues w.r.t time for each order')
    clf
    nCols = 5;
    leg = cell(1,Rmax);
    hdl = zeros(1,Rmax);
    color = distinguishable_colors(Rmax);
    c = 0;
    t = 0:p;
    for r=1:Rmax
        c = c+1;
        hdl(c) = semilogy(t*dt,sig(r,t+1).^2,'LineStyle','-','Color',color(c,:),'LineWidth',1);
        leg{c} = ['$\alpha = ' num2str(r) '$'];
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$t$ [s]','Interpreter','latex')
    ylabel('$\lambda_{\alpha}(t)$','Interpreter','latex')
    gridLegend(hdl,nCols,leg,'Location','BestOutside','Interpreter','latex');
    % legend(leg{:},'Location','NorthEastOutside')
    mysaveas(gridpathname,'eigenvalues_CY_time',formats,renderer);
    mymatlab2tikz(gridpathname,'eigenvalues_CY_time.tex');
    
    figure('Name','Evolution of errors')
    clf
    nCols = 5;
    leg = cell(1,p+1);
    hdl = zeros(1,p+1);
    color = distinguishable_colors(p+1);
    c = 0;
    for t=0:p
        c = c+1;
        Rt = R(t+1);
        hdl(c) = semilogy(1:Rt,errsvdYc(1:Rt,t+1).^2,'LineStyle','-','Color',color(c,:),'LineWidth',1);
        leg{c} = ['$t = ' num2str(t*dt) '$ s'];
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$R$','Interpreter','latex')
    ylabel('$\varepsilon_{Y}(R;t)$','Interpreter','latex')
    gridLegend(hdl,nCols,leg,'Location','BestOutside','Interpreter','latex');
    % legend(leg{:},'Location','NorthEastOutside')
    mysaveas(gridpathname,'error_svdYc',formats,renderer);
    mymatlab2tikz(gridpathname,'error_svdYc.tex');
    
    figure('Name','Evolution of eigenvalues')
    clf
    semilogy(1:Q,s(:).^2,'LineStyle','-','Color','b','LineWidth',1);
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$\beta$','Interpreter','latex')
    ylabel('$\Lambda_{\beta}$','Interpreter','latex')
    mysaveas(gridpathname,'eigenvalues_CZ',formats,renderer);
    mymatlab2tikz(gridpathname,'eigenvalues_CZ.tex');
    
    figure('Name','Evolution of errors')
    clf
    semilogy(1:Q,errsvdZc(:).^2,'LineStyle','-','Color','b','LineWidth',1);
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$Q$','Interpreter','latex')
    ylabel('$\varepsilon_{Z}(Q)$','Interpreter','latex')
    mysaveas(gridpathname,'error_svdZc',formats,renderer);
    mymatlab2tikz(gridpathname,'error_svdZc.tex');
end

% Mean solution
mU = reshape(mY(1,1:3,:,:),[3*m,p+1]);
mC = reshape(mY(1,4,:,:),[m,p+1]);
if postProcessTau
    mtauTime = reshape(mTau(1,1:3,:,:),[3*m,p+1]);
    mdivtauConv = reshape(mTau(1,4:6,:,:),[3*m,p+1]);
    mdivtauDiff = reshape(mTau(1,7:9,:,:),[3*m,p+1]);
    mtauSurf = reshape(mTau(1,10:12,:,:),[3*m,p+1]);
    mtauInterf = reshape(mTau(1,13,:,:),[m,p+1]);
end

% Mean filtered solution
if g==gref
    mUbar = zeros(ng-1,3*m,p+1);
    mCbar = zeros(ng-1,m,p+1);
    for ig=1:ng-1
        gbar = gset(ng-ig);
        mUbar(ig,:,:) = reshape(mYrefbar(ig,1:3,:,:),[3*m,p+1]);
        mCbar(ig,:,:) = reshape(mYrefbar(ig,4,:,:),[m,p+1]);
    end
    clear mYrefbar
    if postProcessTau
        mtauTimebar = zeros(ng-1,3*m,p+1);
        mdivtauConvbar = zeros(ng-1,3*m,p+1);
        mdivtauDiffbar = zeros(ng-1,3*m,p+1);
        mtauSurfbar = zeros(ng-1,3*m,p+1);
        mtauInterfbar = zeros(ng-1,m,p+1);
        for ig=1:ng-1
            mtauTimebar(ig,:,:) = reshape(mTaurefbar(ig,1:3,:,:),[3*m,p+1]);
            mdivtauConvbar(ig,:,:) = reshape(mTaurefbar(ig,4:6,:,:),[3*m,p+1]);
            mdivtauDiffbar(ig,:,:) = reshape(mTaurefbar(ig,7:9,:,:),[3*m,p+1]);
            mtauSurfbar(ig,:,:) = reshape(mTaurefbar(ig,10:12,:,:),[3*m,p+1]);
            mtauInterfbar(ig,:,:) = reshape(mTaurefbar(ig,13,:,:),[m,p+1]);
        end
        clear mTaurefbar
    end
else
    mUbar = reshape(mYbar(1,1:3,:,:),[3*m,p+1]);
    mCbar = reshape(mYbar(1,4,:,:),[m,p+1]);
    clear mYbar
    if postProcessTau
        mtauTimebar = reshape(mTaubar(1,1:3,:,:),[3*m,p+1]);
        mdivtauConvbar = reshape(mTaubar(1,4:6,:,:),[3*m,p+1]);
        mdivtauDiffbar = reshape(mTaubar(1,7:9,:,:),[3*m,p+1]);
        mtauSurfbar = reshape(mTaubar(1,10:12,:,:),[3*m,p+1]);
        mtauInterfbar = reshape(mTaubar(1,13,:,:),[m,p+1]);
        clear mTaubar
    end
end

% Variance of solution
vU = zeros(3*m,p+1);
vC = zeros(m,p+1);
if postProcessTau
    vtauTime = zeros(3*m,p+1);
    vdivtauConv = zeros(3*m,p+1);
    vdivtauDiff = zeros(3*m,p+1);
    vtauSurf = zeros(3*m,p+1);
    vtauInterf = zeros(m,p+1);
end

if performPCA
switch index
    case 'coord'
        W = permute(reshape(W',[Q,Rmax,p+1]),[2,1,3]);
        % Z_approx = reshape(Zc_approx',[N,Rmax,p+1]);
    case 'time'
        W = permute(reshape(W',[Q,p+1,Rmax]),[3,1,2]);
        % Z_approx = permute(reshape(Zc_approx',[N,p+1,Rmax]),[1,3,2]);
end
end
for t=0:p
    if performPCA
        Rt = R(t+1);
        if g<2^7
            Vt = V(:,1:Rt,t+1);
        else
            load(fullfile(gridpathname,['PCA_space_t' num2str(t) '.mat']),'Vt');
        end
        sigt = sig(1:Rt,t+1);
        Wt = W(1:Rt,:,t+1);
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            % Zct_approx = Wt*diag(s)*X'; % Zct_approx = Z_approx(:,1:Rt,t+1)';
            % Yct_approx = Vt*diag(sigt)*Zct_approx;
            Uct_approx = Vt*diag(sigt)*Wt;
        else
            % Zct_approx = Wt*(s.*X'); % Zct_approx = Z_approx(:,1:Rt,t+1)';
            % Yct_approx = Vt*(sigt.*Zct_approx);
            Uct_approx = Vt*(sigt.*Wt);
        end
        
        % CYt_approx = cov(Yct_approx'); % CYt_approx = 1/(N-1)*Yct_approx*Yct_approx';
        % if verLessThan('matlab','9.1') % compatibility (<R2016b)
        %     CYt_approx = Uct_approx*diag(s).^2*Uct_approx';
        % else
        %     CYt_approx = Uct_approx*(s.^2.*Uct_approx');
        % end
        % vYt_approx = diag(CYt_approx);
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            vYt_approx = sum((diag(s)*Uct_approx').^2)';
        else
            vYt_approx = sum((s.*Uct_approx').^2)';
        end
        
        indut = repmat((0:m-1)*n,[3,1])+repmat((1:3)',[1,m]);
        indCt = (0:m-1)*n+4;
        
        vUt = vYt_approx(indut(:));
        vCt = vYt_approx(indCt(:));
    else
        mYt = mY(:,:,:,t+1);
        if g<2^7
            Yt = Y(:,:,:,t+1);
        else
            load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
        end
        Yct = Yt - repmat(mYt,[N,1,1]); % Yct = Yt - mYt.*ones(N,1,1);
        clear Yt
        
        ut = Yct(:,1:3,:);
        Ct = Yct(:,4,:);
        clear Yct
        
        vUt = 1/(N-1)*sum(ut(:,:).^2)';
        vCt = 1/(N-1)*sum(Ct(:,:).^2)';
        clear ut Ct
    end
    
    if postProcessTau
        mTaut = mTau(:,:,:,t+1);
        if g<2^7
            Taut = Tau(:,:,:,t+1);
        else
            load(fullfile(gridpathname,['data_tau_t' num2str(t) '.mat']),'Taut');
        end
        Tauct = Taut - repmat(mTaut,[N,1,1]); % Tauct = Taut - mTaut.*ones(N,1,1);
        clear Taut
        
        tauTimet = Tauct(:,1:3,:);
        divtauConvt = Tauct(:,4:6,:);
        divtauDifft = Tauct(:,7:9,:);
        tauSurft = Tauct(:,10:12,:);
        tauInterft = Tauct(:,13,:);
        clear Tauct
        
        vtauTimet = 1/(N-1)*sum(tauTimet(:,:).^2)';
        vdivtauConvt = 1/(N-1)*sum(divtauConvt(:,:).^2)';
        vdivtauDifft = 1/(N-1)*sum(divtauDifft(:,:).^2)';
        vtauSurft = 1/(N-1)*sum(tauSurft(:,:).^2)';
        vtauInterft = 1/(N-1)*sum(tauInterft(:,:).^2)';
        clear tauTimet divtauConvt divtauDifft tauSurft tauInterft
    end
    
    vU(:,t+1) = vUt;
    vC(:,t+1) = vCt;
    if postProcessTau
        vtauTime(:,t+1) = vtauTimet;
        vdivtauConv(:,t+1) = vdivtauConvt;
        vdivtauDiff(:,t+1) = vdivtauDifft;
        vtauSurf(:,t+1) = vtauSurft;
        vtauInterf(:,t+1) = vtauInterft;
    end
end

% Display error between data and filtered data
if g~=gref && displayError
    ts = (0:p)*dt;
    
    figure('Name','Error between data and filtered data')
    clf
    hdl(1) = plot(ts,erroru./normubar,'LineStyle','-','Color','b','LineWidth',1);
    hold on
    hdl(2) = plot(ts,errorC./normCbar,'LineStyle','-','Color','r','LineWidth',1);
    if postProcessTau
        hdl(3) = plot(ts,errortauTime./normtauTimebar,'LineStyle','-','Color','g','LineWidth',1);
        hdl(4) = plot(ts,errordivtauConv./normdivtauConvbar,'LineStyle','-','Color','m','LineWidth',1);
        hdl(5) = plot(ts,errordivtauDiff./normdivtauDiffbar,'LineStyle','-','Color','c','LineWidth',1);
        hdl(6) = plot(ts,errortauSurf./normtauSurfbar,'LineStyle','-','Color','y','LineWidth',1);
        hdl(7) = plot(ts,errortauInterf./normtauInterfbar,'LineStyle','-','Color','k','LineWidth',1);
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$t$ [s]','Interpreter','latex')
    ylabel('Normalized error','Interpreter','latex')
    leg = {'$||u||$','$\chi$'};
    if postProcessTau
        leg = [leg,'$||\tau_{\mathrm{time}}||$','$||\nabla \cdot \tau_{\mathrm{conv}}||$','$||\nabla \cdot \tau_{\mathrm{diff}}||$','$||\tau_{\mathrm{surf}}||$','$\tau_{\mathrm{interf}}$'];
    end
    l = legend(leg{:},'Location','NorthEast');
    set(l,'Interpreter','latex')
    mysaveas(gridpathname,'error_filter',formats,renderer);
    mymatlab2tikz(gridpathname,'error_filter.tex');
end

% Display quantities of interest
if displayQoI
    % Mean
    mQu = reshape(mQu(1,1:3,:,:),[3,2,p+1]);
    if postProcessTau
        mQtauTime = reshape(mQtau(1,1:3,:,:),[3,2,p+1]);
        mQdivtauConv = reshape(mQtau(1,4:6,:,:),[3,2,p+1]);
        mQdivtauDiff = reshape(mQtau(1,7:9,:,:),[3,2,p+1]);
        mQtauSurf = reshape(mQtau(1,10:12,:,:),[3,2,p+1]);
        mQtauInterf = reshape(mQtau(1,13,:,:),[1,2,p+1]);
        % Correlation
        IQtauTime = IQtau(1:3,:,1:3,:,:);
        IQdivtauConv = IQtau(4:6,:,4:6,:,:);
        IQdivtauDiff = IQtau(7:9,:,7:9,:,:);
        IQtauSurf = IQtau(10:12,:,10:12,:,:);
        IQtauInterf = IQtau(13,:,13,:,:);
    end
    if postProcessEnergy
        mQenergyKinTime = reshape(mQe(1,1,:,:),[1,2,p+1]);
        mQenergyConv = reshape(mQe(1,2,:,:),[1,2,p+1]);
        mQenergyGrav = reshape(mQe(1,3,:,:),[1,2,p+1]);
        mQenergyPres = reshape(mQe(1,4,:,:),[1,2,p+1]);
        mQenergyPresDil = reshape(mQe(1,5,:,:),[1,2,p+1]);
        mQenergyKinSpace = reshape(mQe(1,6,:,:),[1,2,p+1]);
        mQenergyDiff = reshape(mQe(1,7,:,:),[1,2,p+1]);
        mQenergyVisc = reshape(mQe(1,8,:,:),[1,2,p+1]);
        mQenergySurf = reshape(mQe(1,9,:,:),[1,2,p+1]);
        % Correlation
        IQenergyKinTime = IQe(1,:,1,:,:);
        IQenergyConv = IQe(2,:,2,:,:);
        IQenergyGrav = IQe(3,:,3,:,:);
        IQenergyPres = IQe(4,:,4,:,:);
        IQenergyPresDil = IQe(5,:,5,:,:);
        IQenergyKinSpace = IQe(6,:,6,:,:);
        IQenergyDiff = IQe(7,:,7,:,:);
        IQenergyVisc = IQe(8,:,8,:,:);
        IQenergySurf = IQe(9,:,9,:,:);
    end
    
    ts = (0:p)*dt;
    
    figure('Name','Mean of spatial average of velocity u')
    clf
    hdl(1) = plot(ts,squeeze(mQu(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
    hold on
    hdl(2) = plot(ts,squeeze(mQu(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
    hdl(3) = plot(ts,squeeze(mQu(2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
    hdl(4) = plot(ts,squeeze(mQu(2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
    hdl(5) = plot(ts,squeeze(mQu(3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
    hdl(6) = plot(ts,squeeze(mQu(3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$t$ [s]','Interpreter','latex')
    ylabel('$u$','Interpreter','latex')
    leg = {'component 1 in phase 1','component 1 in phase 2','component 2 in phase 1','component 2 in phase 2','component 3 in phase 1','component 3 in phase 2'};
    legend(leg{:},'Location','NorthEast')
    mysaveas(gridpathname,'mean_u',formats,renderer);
    mymatlab2tikz(gridpathname,'mean_u.tex');
    
    if postProcessTau
        figure('Name','Mean of spatial average of tauTime')
        clf
        hdl(1) = plot(ts,squeeze(mQtauTime(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(mQtauTime(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hdl(3) = plot(ts,squeeze(mQtauTime(2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
        hdl(4) = plot(ts,squeeze(mQtauTime(2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
        hdl(5) = plot(ts,squeeze(mQtauTime(3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
        hdl(6) = plot(ts,squeeze(mQtauTime(3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$t$ [s]','Interpreter','latex')
        ylabel('$\tau_{\mathrm{time}}$','Interpreter','latex')
        leg = {'component 1 in phase 1','component 1 in phase 2','component 2 in phase 1','component 2 in phase 2','component 3 in phase 1','component 3 in phase 2'};
        legend(leg{:},'Location','NorthEast')
        mysaveas(gridpathname,'mean_tauTime',formats,renderer);
        mymatlab2tikz(gridpathname,'mean_tauTime.tex');
        
        figure('Name','Mean of spatial average of div(tauConv)')
        clf
        hdl(1) = plot(ts,squeeze(mQdivtauConv(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(mQdivtauConv(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hdl(3) = plot(ts,squeeze(mQdivtauConv(2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
        hdl(4) = plot(ts,squeeze(mQdivtauConv(2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
        hdl(5) = plot(ts,squeeze(mQdivtauConv(3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
        hdl(6) = plot(ts,squeeze(mQdivtauConv(3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$t$ [s]','Interpreter','latex')
        ylabel('$\nabla \cdot \tau_{\mathrm{conv}}$','Interpreter','latex')
        leg = {'component 1 in phase 1','component 1 in phase 2','component 2 in phase 1','component 2 in phase 2','component 3 in phase 1','component 3 in phase 2'};
        legend(leg{:},'Location','NorthEast')
        mysaveas(gridpathname,'mean_divtauConv',formats,renderer);
        mymatlab2tikz(gridpathname,'mean_divtauConv.tex');
        
        figure('Name','Mean of spatial average of div(tauDiff)')
        clf
        hdl(1) = plot(ts,squeeze(mQdivtauDiff(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(mQdivtauDiff(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hdl(3) = plot(ts,squeeze(mQdivtauDiff(2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
        hdl(4) = plot(ts,squeeze(mQdivtauDiff(2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
        hdl(5) = plot(ts,squeeze(mQdivtauDiff(3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
        hdl(6) = plot(ts,squeeze(mQdivtauDiff(3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$t$ [s]','Interpreter','latex')
        ylabel('$\nabla \cdot \tau_{\mathrm{diff}}$','Interpreter','latex')
        leg = {'component 1 in phase 1','component 1 in phase 2','component 2 in phase 1','component 2 in phase 2','component 3 in phase 1','component 3 in phase 2'};
        legend(leg{:},'Location','NorthEast')
        mysaveas(gridpathname,'mean_divtauDiff',formats,renderer);
        mymatlab2tikz(gridpathname,'mean_divtauDiff.tex');
        
        figure('Name','Mean of spatial average of tauSurf')
        clf
        hdl(1) = plot(ts,squeeze(mQtauSurf(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(mQtauSurf(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hdl(3) = plot(ts,squeeze(mQtauSurf(2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
        hdl(4) = plot(ts,squeeze(mQtauSurf(2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
        hdl(5) = plot(ts,squeeze(mQtauSurf(3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
        hdl(6) = plot(ts,squeeze(mQtauSurf(3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$t$ [s]','Interpreter','latex')
        ylabel('$\tau_{\mathrm{surf}}$','Interpreter','latex')
        leg = {'component 1 in phase 1','component 1 in phase 2','component 2 in phase 1','component 2 in phase 2','component 3 in phase 1','component 3 in phase 2'};
        legend(leg{:},'Location','NorthEast')
        mysaveas(gridpathname,'mean_tauSurf',formats,renderer);
        mymatlab2tikz(gridpathname,'mean_tauSurf.tex');
        
        figure('Name','Mean of spatial average of tauInterf')
        clf
        hdl(1) = plot(ts,squeeze(mQtauInterf(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(mQtauInterf(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$t$ [s]','Interpreter','latex')
        ylabel('$\tau_{\mathrm{interf}}$','Interpreter','latex')
        leg = {'phase 1','phase 2'};
        legend(leg{:},'Location','NorthEast')
        mysaveas(gridpathname,'mean_tauInterf',formats,renderer);
        mymatlab2tikz(gridpathname,'mean_tauInterf.tex');
    end
    
    figure('Name','Power of u')
    clf
    hdl(1) = plot(ts,squeeze(IQu(1,1,1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
    hold on
    hdl(2) = plot(ts,squeeze(IQu(1,2,1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
    hdl(3) = plot(ts,squeeze(IQu(2,1,2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
    hdl(4) = plot(ts,squeeze(IQu(2,2,2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
    hdl(5) = plot(ts,squeeze(IQu(3,1,3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
    hdl(6) = plot(ts,squeeze(IQu(3,2,3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$\tau$ [s]','Interpreter','latex')
    ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter','latex')
    leg = {'u_1 in phase 1 - u_1 in phase 1','u_1 in phase 2 - u_1 in phase 2',...
        'u_2 in phase 1 - u_2 in phase 1','u_2 in phase 2 - u_2 in phase 2',...
        'u_3 in phase 1 - u_3 in phase 1','u_3 in phase 2 - u_3 in phase 2'};
    legend(leg{:},'Location','NorthEast')
    mysaveas(gridpathname,'power_u',formats,renderer);
    mymatlab2tikz(gridpathname,'power_u.tex');
    
    if postProcessTau
        figure('Name','Power of tauTime')
        clf
        hdl(1) = plot(ts,squeeze(IQtauTime(1,1,1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(IQtauTime(1,2,1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hdl(3) = plot(ts,squeeze(IQtauTime(2,1,2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
        hdl(4) = plot(ts,squeeze(IQtauTime(2,2,2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
        hdl(5) = plot(ts,squeeze(IQtauTime(3,1,3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
        hdl(6) = plot(ts,squeeze(IQtauTime(3,2,3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$\tau$ [s]','Interpreter','latex')
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter','latex')
        leg = {'$\tau_{\mathrm{time}\,1}$ in phase 1','$\tau_{\mathrm{time}\,1}$ in phase 2',...
            '$\tau_{\mathrm{time}\,2}$ in phase 1','$\tau_{\mathrm{time}\,2}$ in phase 2',...
            '$\tau_{\mathrm{time}\,3}$ in phase 1','$\tau_{\mathrm{time}\,3}$ in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter','latex')
        mysaveas(gridpathname,'power_tauTime',formats,renderer);
        mymatlab2tikz(gridpathname,'power_tauTime.tex');
        
        figure('Name','Power of div(tauConv)')
        clf
        hdl(1) = plot(ts,squeeze(IQdivtauConv(1,1,1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(IQdivtauConv(1,2,1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hdl(3) = plot(ts,squeeze(IQdivtauConv(2,1,2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
        hdl(4) = plot(ts,squeeze(IQdivtauConv(2,2,2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
        hdl(5) = plot(ts,squeeze(IQdivtauConv(3,1,3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
        hdl(6) = plot(ts,squeeze(IQdivtauConv(3,2,3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$\tau$ [s]','Interpreter','latex')
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter','latex')
        leg = {'$(\nabla \cdot \tau_{\mathrm{conv}})_1$ in phase 1','$(\nabla \cdot \tau_{\mathrm{conv}})_1$ in phase 2',...
            '$(\nabla \cdot \tau_{\mathrm{conv}})_2$ in phase 1','$(\nabla \cdot \tau_{\mathrm{conv}})_2$ in phase 2',...
            '$(\nabla \cdot \tau_{\mathrm{conv}})_3$ in phase 1','$(\nabla \cdot \tau_{\mathrm{conv}})_3$ in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter','latex')
        mysaveas(gridpathname,'power_divtauConv',formats,renderer);
        mymatlab2tikz(gridpathname,'power_divtauConv.tex');
        
        figure('Name','Power of div(tauDiff)')
        clf
        hdl(1) = plot(ts,squeeze(IQdivtauDiff(1,1,1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(IQdivtauDiff(1,2,1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hdl(3) = plot(ts,squeeze(IQdivtauDiff(2,1,2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
        hdl(4) = plot(ts,squeeze(IQdivtauDiff(2,2,2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
        hdl(5) = plot(ts,squeeze(IQdivtauDiff(3,1,3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
        hdl(6) = plot(ts,squeeze(IQdivtauDiff(3,2,3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$\tau$ [s]','Interpreter','latex')
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter','latex')
        leg = {'$(\nabla \cdot \tau_{\mathrm{diff}})_1$ in phase 1','$(\nabla \cdot \tau_{\mathrm{diff}})_1$ in phase 2',...
            '$(\nabla \cdot \tau_{\mathrm{diff}})_2$ in phase 1','$(\nabla \cdot \tau_{\mathrm{diff}})_2$ in phase 2',...
            '$(\nabla \cdot \tau_{\mathrm{diff}})_3$ in phase 1','$(\nabla \cdot \tau_{\mathrm{diff}})_3$ in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter','latex')
        mysaveas(gridpathname,'power_divtauDiff',formats,renderer);
        mymatlab2tikz(gridpathname,'power_divtauDiff.tex');
        
        figure('Name','Power of tauSurf')
        clf
        hdl(1) = plot(ts,squeeze(IQtauSurf(1,1,1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(IQtauSurf(1,2,1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hdl(3) = plot(ts,squeeze(IQtauSurf(2,1,2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
        hdl(4) = plot(ts,squeeze(IQtauSurf(2,2,2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
        hdl(5) = plot(ts,squeeze(IQtauSurf(3,1,3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
        hdl(6) = plot(ts,squeeze(IQtauSurf(3,2,3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$\tau$ [s]','Interpreter','latex')
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter','latex')
        leg = {'$\tau_{\mathrm{surf}\,1}$ in phase 1','$\tau_{\mathrm{surf}\,1}$ in phase 2',...
            '$\tau_{\mathrm{surf}\,2}$ in phase 1','$\tau_{\mathrm{surf}\,2}$ in phase 2',...
            '$\tau_{\mathrm{surf}\,3}$ in phase 1','$\tau_{\mathrm{surf}\,3}$ in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter','latex')
        mysaveas(gridpathname,'power_tauSurf',formats,renderer);
        mymatlab2tikz(gridpathname,'power_tauSurf.tex');
        
        figure('Name','Power of tauInterf')
        clf
        hdl(1) = plot(ts,squeeze(IQtauInterf(1,1,1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(IQtauInterf(1,2,1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$\tau$ [s]','Interpreter','latex')
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter','latex')
        leg = {'$\tau_{\mathrm{interf}}$ in phase 1','$\tau_{\mathrm{interf}}$ in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter','latex')
        mysaveas(gridpathname,'power_tauInterf',formats,renderer);
        mymatlab2tikz(gridpathname,'power_tauInterf.tex');
    end
end

%% Spatial mesh
if constructMesh
    fprintf('\nConstructing spatial mesh');
    t_Mesh = tic;
    
    D = DOMAIN(3,[0.0,0.0,0.0],[L,L,L]);
    elemtype = 'CUB8';
    nbelem = [g,g,g];
    M = build_model(D,'nbelem',nbelem,'elemtype',elemtype);
    coord = getcoord(getnode(M));
    M = setnode(M,NODE(coord(:,[2,1,3])));
    M = final(M,DDL(DDLVECT('U',M.syscoord)));
    
    time_Mesh = toc(t_Mesh);
    fprintf('\nelapsed time = %f s',time_Mesh);
    fprintf('\n');
    
    save(fullfile(gridpathname,'mesh.mat'),'M','time_Mesh');
else
    fprintf('\nLoading spatial mesh');
    t_load = tic;
    load(fullfile(gridpathname,'mesh.mat'),'M','time_Mesh');
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
end

%% Mean
for t=0:p
    mUt = mU(:,t+1);
    mCt = mC(:,t+1);
    fields = {mUt,mCt};
    fieldnames = {'velocity','phase'};
    if postProcessTau
        mtauTimet = mtauTime(:,t+1);
        mdivtauConvt = mdivtauConv(:,t+1);
        mdivtauDifft = mdivtauDiff(:,t+1);
        mtauSurft = mtauSurf(:,t+1);
        mtauInterft = mtauInterf(:,t+1);
        fields = [fields,mtauTimet,mdivtauConvt,mdivtauDifft,mtauSurft,mtauInterft];
        fieldnames = [fieldnames,'tauTime','div(tauConv)','div(tauDiff)','tauSurf','tauInterf'];
    end
%     if postProcessEnergy
%         menergyKinTimet = menergyKinTime(:,t+1);
%         menergyConvt = menergyConv(:,t+1);
%         menergyGravt = menergyGrav(:,t+1);
%         menergyPrest = menergyPres(:,t+1);
%         menergyPresDilt = menergyPresDil(:,t+1);
%         menergyKinSpacet = menergyKinSpace(:,t+1);
%         menergyDifft = menergyDiff(:,t+1);
%         menergyVisct = menergyVisc(:,t+1);
%         menergySurft = menergySurf(:,t+1);
%         fields = [fields,menergyKinTimet,menergyConvt,menergyGravt,...
%             menergyPrest,menergyPresDilt,menergyKinSpacet,...
%             menergyDifft,menergyVisct,menergySurft];
%         fieldnames = [fieldnames,'kinetic energy','convection energy','gravity energy',...
%             'power of external pressure forces','pressure-dilatation energy transfer','transport of gradient of kinetic energy',...
%             'energy exchange with kinetic energy','power of external viscous stresses','capillary kinetic energy'];
%     end
    
    if g==gref
        for ig=1:ng-1
            gbar = gset(ng-ig);
            mUbart = mUbar(ig,:,t+1);
            mCbart = mCbar(ig,:,t+1);
            legUbart = ['velocity filtered on grid ' num2str(gbar)];
            legCbart = ['phase filtered on grid ' num2str(gbar)];
            fields = [fields,mUbart,mCbart];
            fieldnames = [fieldnames,legUbart,legCbart];
            if postProcessTau
                mtauTimebart = mtauTimebar(ig,:,t+1);
                mdivtauConvbart = mdivtauConvbar(ig,:,t+1);
                mdivtauDiffbart = mdivtauDiffbar(ig,:,t+1);
                mtauSurfbart = mtauSurfbar(ig,:,t+1);
                mtauInterfbart = mtauInterfbar(ig,:,t+1);
                legtauTimebart = ['tauTime filtered on grid ' num2str(gbar)];
                legdivtauConvbart = ['div(tauConv) filtered on grid ' num2str(gbar)];
                legdivtauDiffbart = ['div(tauDiff) filtered on grid ' num2str(gbar)];
                legtauSurfbart = ['tauSurf filtered on grid ' num2str(gbar)];
                legtauInterfbart = ['tauInterf filtered on grid ' num2str(gbar)];
                fields = [fields,mtauTimebart,mdivtauConvbart,mdivtauDiffbart,mtauSurfbart,mtauInterfbart];
                fieldnames = [fieldnames,legtauTimebart,legdivtauConvbart,legdivtauDiffbart,legtauSurfbart,legtauInterfbart];
            end
%             if postProcessEnergy
%                 menergyKinTimebart = menergyKinTimebar(ig,:,t+1);
%                 menergyConvbart = menergyConvbar(ig,:,t+1);
%                 menergyGravbart = menergyGravbar(ig,:,t+1);
%                 menergyPresbart = menergyPresbar(ig,:,t+1);
%                 menergyPresDilbart = menergyPresDilbar(ig,:,t+1);
%                 menergyKinSpacebart = menergyKinSpacebar(ig,:,t+1);
%                 menergyDiffbart = menergyDiffbar(ig,:,t+1);
%                 menergyViscbart = menergyViscbar(ig,:,t+1);
%                 menergySurfbart = menergySurfbar(ig,:,t+1);
%                 legenergyKinTimebart = ['kinetic energy filtered on grid ' num2str(gbar)];
%                 legenergyConvbart = ['convection energy filtered on grid ' num2str(gbar)];
%                 legenergyGravbart = ['gravity energy filtered on grid ' num2str(gbar)];
%                 legenergyPresbart = ['power of external pressure forces filtered on grid ' num2str(gbar)];
%                 legenergyPresDilbart = ['pressure-dilatation energy transfer filtered on grid ' num2str(gbar)];
%                 legenergyKinSpacebart = ['transport of gradient of kinetic energy filtered on grid ' num2str(gbar)];
%                 legenergyDiffbart = ['energy exchange with kinetic energy filtered on grid ' num2str(gbar)];
%                 legenergyViscbart = ['power of external viscous stresses filtered on grid ' num2str(gbar)];
%                 legenergySurfbart = ['capillary kinetic energy filtered on grid ' num2str(gbar)];
%                 fields = [fields,menergyKinTimebart,menergyConvbart,menergyGravbart,...
%                     menergyPresbart,menergyPresDilbart,menergyKinSpacebart,...
%                     menergyDiffbart,menergyViscbart,menergySurfbart];
%                 fieldnames = [fieldnames,legenergyKinTimebart,legenergyConvbart,legenergyGravbart,...
%                     legenergyPresbart,legenergyPresDilbart,legenergyKinSpacebart,...
%                     legenergyDiffbart,legenergyViscbart,legenergySurfbart];
%             end
        end
        
    else
        mUbart = mUbar(:,t+1);
        mCbart = mCbar(:,t+1);
        fields = [fields,mUbart,mCbart];
        fieldnames = [fieldnames,'velocity filtered','phase filtered'];
        if postProcessTau
            mtauTimebart = mtauTimebar(:,t+1);
            mdivtauConvbart = mdivtauConvbar(:,t+1);
            mdivtauDiffbart = mdivtauDiffbar(:,t+1);
            mtauSurfbart = mtauSurfbar(:,t+1);
            mtauInterfbart = mtauInterfbar(:,t+1);
            fields = [fields,mtauTimebart,mdivtauConvbart,mdivtauDiffbart,mtauSurfbart,mtauInterfbart];
            fieldnames = [fieldnames,'tauTime filtered','div(tauConv) filtered','div(tauDiff) filtered','tauSurf filtered','tauInterf filtered'];
        end
%         if postProcessEnergy
%             menergyKinTimebart = menergyKinTimebar(:,t+1);
%             menergyConvbart = menergyConvbar(:,t+1);
%             menergyGravbart = menergyGravbar(:,t+1);
%             menergyPresbart = menergyPresbar(:,t+1);
%             menergyPresDilbart = menergyPresDilbar(:,t+1);
%             menergyKinSpacebart = menergyKinSpacebar(:,t+1);
%             menergyDiffbart = menergyDiffbar(:,t+1);
%             menergyViscbart = menergyViscbar(:,t+1);
%             menergySurfbart = menergySurfbar(:,t+1);
%             fields = [fields,menergyKinTimebart,menergyConvbart,menergyGravbart,...
%                 menergyPresbart,menergyPresDilbart,menergyKinSpacebart,...
%                 menergyDiffbart,menergyViscbart,menergySurfbart];
%             fieldnames = [fieldnames,'kinetic energy filtered','convection energy filtered','gravity energy filtered',...
%                 'power of external pressure forces filtered','pressure-dilatation energy transfer filtered','transport of gradient of kinetic energy filtered',...
%                 'energy exchange with kinetic energy filtered','power of external viscous stresses filtered','capillary kinetic energy filtered'];
%         end
    end
    write_vtk_mesh(M,fields,[],fieldnames,[],gridpathname,['diphasic_fluids_grid' num2str(g) '_mean'],1,t);
end
make_pvd_file(gridpathname,['diphasic_fluids_grid' num2str(g) '_mean'],1,p+1);

%% Variance
% for t=0:p
%     vUt = vU(:,t+1);
%     vCt = vC(:,t+1);
%     fields = {vUt,vCt};
%     fieldnames = {'velocity','phase'};
%     if postProcessTau
%         vtauTimet = vtauTime(:,t+1);
%         vdivtauConvt = vdivtauConv(:,t+1);
%         vdivtauDifft = vdivtauDiff(:,t+1);
%         vtauSurft = vtauSurf(:,t+1);
%         vtauInterft = vtauInterf(:,t+1);
%         fields = [fields,vtauTimet,vdivtauConvt,vdivtauDifft,vtauSurft,vtauInterft];
%         fieldnames = [fieldnames,'tauTime','div(tauConv)','div(tauDiff)','tauSurf','tauInterf'];
%     end
% %     if postProcessEnergy
% %         venergyKinTimet = venergyKinTime(:,t+1);
% %         venergyConvt = venergyConv(:,t+1);
% %         venergyGravt = venergyGrav(:,t+1);
% %         venergyPrest = venergyPres(:,t+1);
% %         venergyPresDilt = venergyPresDil(:,t+1);
% %         venergyKinSpacet = venergyKinSpace(:,t+1);
% %         venergyDifft = venergyDiff(:,t+1);
% %         venergyVisct = venergyVisc(:,t+1);
% %         venergySurft = venergySurf(:,t+1);
% %         fields = [fields,venergyKinTimet,venergyConvt,venergyGravt,...
% %             venergyPrest,venergyPresDilt,venergyKinSpacet,...
% %             venergyDifft,venergyVisct,venergySurft];
% %         fieldnames = [fieldnames,'kinetic energy','convection energy','gravity energy',...
% %             'power of external pressure forces','pressure-dilatation energy transfer','transport of gradient of kinetic energy',...
% %             'energy exchange with kinetic energy','power of external viscous stresses','capillary kinetic energy'];
% %     end
%     
%     if g==gref
%         for ig=1:ng-1
%             gbar = gset(ng-ig);
%             vUbart = vUbar(ig,:,t+1);
%             vCbart = vCbar(ig,:,t+1);
%             legUbart = ['velocity filtered on grid ' num2str(gbar)];
%             legCbart = ['phase filtered on grid ' num2str(gbar)];
%             fields = [fields,vUbart,vCbart];
%             fieldnames = [fieldnames,legUbart,legCbart];
%             if postProcessTau
%                 vtauTimebart = vtauTimebar(ig,:,t+1);
%                 mdivtauConvbart = mdivtauConvbar(ig,:,t+1);
%                 mdivtauDiffbart = mdivtauDiffbar(ig,:,t+1);
%                 vtauSurfbart = vtauSurfbar(ig,:,t+1);
%                 vtauInterfbart = vtauInterfbar(ig,:,t+1);
%                 legtauTimebart = ['tauTime filtered on grid ' num2str(gbar)];
%                 legdivtauConvbart = ['div(tauConv) filtered on grid ' num2str(gbar)];
%                 legdivtauDiffbart = ['div(tauDiff) filtered on grid ' num2str(gbar)];
%                 legtauSurfbart = ['tauSurf filtered on grid ' num2str(gbar)];
%                 legtauInterfbart = ['tauInterf filtered on grid ' num2str(gbar)];
%                 fields = [fields,vtauTimebart,mdivtauConvbart,mdivtauDiffbart,vtauSurfbart,vtauInterfbart];
%                 fieldnames = [fieldnames,legtauTimebart,legdivtauConvbart,legdivtauDiffbart,legtauSurfbart,legtauInterfbart];
%             end
% %             if postProcessEnergy
% %                 venergyKinTimebart = venergyKinTimebar(ig,:,t+1);
% %                 venergyConvbart = venergyConvbar(ig,:,t+1);
% %                 venergyGravbart = venergyGravbar(ig,:,t+1);
% %                 venergyPresbart = venergyPresbar(ig,:,t+1);
% %                 venergyPresDilbart = venergyPresDilbar(ig,:,t+1);
% %                 venergyKinSpacebart = venergyKinSpacebar(ig,:,t+1);
% %                 venergyDiffbart = venergyDiffbar(ig,:,t+1);
% %                 venergyViscbart = venergyViscbar(ig,:,t+1);
% %                 venergySurfbart = venergySurfbar(ig,:,t+1);
% %                 legenergyKinTimebart = ['kinetic energy filtered on grid ' num2str(gbar)];
% %                 legenergyConvbart = ['convection energy filtered on grid ' num2str(gbar)];
% %                 legenergyGravbart = ['gravity energy filtered on grid ' num2str(gbar)];
% %                 legenergyPresbart = ['power of external pressure forces filtered on grid ' num2str(gbar)];
% %                 legenergyPresDilbart = ['pressure-dilatation energy transfer filtered on grid ' num2str(gbar)];
% %                 legenergyKinSpacebart = ['transport of gradient of kinetic energy filtered on grid ' num2str(gbar)];
% %                 legenergyDiffbart = ['energy exchange with kinetic energy filtered on grid ' num2str(gbar)];
% %                 legenergyViscbart = ['power of external viscous stresses filtered on grid ' num2str(gbar)];
% %                 legenergySurfbart = ['capillary kinetic energy filtered on grid ' num2str(gbar)];
% %                 fields = [fields,venergyKinTimebart,venergyConvbart,venergyGravbart,...
% %                     venergyPresbart,venergyPresDilbart,venergyKinSpacebart,...
% %                     venergyDiffbart,venergyViscbart,venergySurfbart];
% %                 fieldnames = [fieldnames,legenergyKinTimebart,legenergyConvbart,legenergyGravbart,...
% %                     legenergyPresbart,legenergyPresDilbart,legenergyKinSpacebart,...
% %                     legenergyDiffbart,legenergyViscbart,legenergySurfbart];
% %             end
%         end
%         
%     else
%         vUbart = vUbar(:,t+1);
%         vCbart = vCbar(:,t+1);
%         fields = [fields,mUbart,mCbart];
%         fieldnames = [fieldnames,'velocity filtered','phase filtered'];
%         if postProcessTau
%             vtauTimebart = vtauTimebar(:,t+1);
%             vdivtauConvbart = vdivtauConvbar(:,t+1);
%             vdivtauDiffbart = vdivtauDiffbar(:,t+1);
%             vtauSurfbart = vtauSurfbar(:,t+1);
%             vtauInterfbart = vtauInterfbar(:,t+1);
%             fields = [fields,vtauTimebart,vdivtauConvbart,vdivtauDiffbart,vtauSurfbart,vtauInterfbart];
%             fieldnames = [fieldnames,'tauTime filtered','div(tauConv) filtered','div(tauDiff) filtered','tauSurf filtered','tauInterf filtered'];
%         end
% %         if postProcessEnergy
% %             venergyKinTimebart = venergyKinTimebar(:,t+1);
% %             venergyConvbart = venergyConvbar(:,t+1);
% %             venergyGravbart = venergyGravbar(:,t+1);
% %             venergyPresbart = venergyPresbar(:,t+1);
% %             venergyPresDilbart = venergyPresDilbar(:,t+1);
% %             venergyKinSpacebart = venergyKinSpacebar(:,t+1);
% %             venergyDiffbart = venergyDiffbar(:,t+1);
% %             venergyViscbart = venergyViscbar(:,t+1);
% %             venergySurfbart = venergySurfbar(:,t+1);
% %             fields = [fields,venergyKinTimebart,venergyConvbart,venergyGravbart,...
% %                 venergyPresbart,venergyPresDilbart,venergyKinSpacebart,...
% %                 venergyDiffbart,venergyViscbart,venergySurfbart];
% %             fieldnames = [fieldnames,'kinetic energy filtered','convection energy filtered','gravity energy filtered',...
% %                 'power of external pressure forces filtered','pressure-dilatation energy transfer filtered','transport of gradient of kinetic energy filtered',...
% %                 'energy exchange with kinetic energy filtered','power of external viscous stresses filtered','capillary kinetic energy filtered'];
% %         end
%     end
%     write_vtk_mesh(M,fields,[],fieldnames,[],gridpathname,['diphasic_fluids_grid' num2str(g) '_variance'],1,t);
% end
% make_pvd_file(gridpathname,['diphasic_fluids_grid' num2str(g) '_variance'],1,p+1);

time_Total = toc(t_Total);
fprintf('\nElapsed time = %f s\n',time_Total);
