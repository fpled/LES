clc
clearvars
close all

computePCAspace = true;
computePCAtime = true;
postProcess = true;
applyFilter = true;

displaySolution = false;
displayEigenvalues = false;
displayCovariance = false;
displayQoI = false;

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
gset = 2.^(4:6); % set of spatial grid sizes
g = gset(3); % current spatial grid size
gref = gset(end); % reference spatial grid size
ng = length(gset); % number of spatial grid sizes

t_Total = tic;
gridname = ['Grid' num2str(g)];
fprintf([gridname '\n'])
gridpathname = fullfile(pathname,gridname);
load(fullfile(gridpathname,'data.mat'),'N','n','m','p');
r = n*m;

fprintf('\nn = %d variables',n);
fprintf('\nN = %d samples',N);
fprintf('\nm = %d spatial points',m);
fprintf('\np+1 = %d time steps',p+1);
fprintf('\n');

% Spatial mesh
L = 1; % domain size (m)
D = DOMAIN(3,[0.0,0.0,0.0],[L,L,L]);
elemtype = 'CUB8';
nbelem = [g,g,g];
M = build_model(D,'nbelem',nbelem,'elemtype',elemtype);
coord = getcoord(getnode(M));
M = setnode(M,NODE(coord(:,[2,1,3])));
M = final(M,DDL(DDLVECT('U',M.syscoord)));
Mscal = final(M,DDL('C'));

% Spatial scheme
sx = [g+1,g+1,g+1]; % sx = nbelem+1; % spatial dimensions
dx = L/g; % spatial step (m)
x = linspace(0,L,g+1); % spatial discretization in each spatial dimension
Dx = spdiags(repmat([1 -1],g+1,1),[1 -1],g+1,g+1)/(2*dx); % Implicit central-difference spatial scheme (second-order accurate, unconditionlly stable)

% Time scheme
dt = 5e-3; % computing time step (s)
dt = 100*dt; % physical time step stored every 100 computing time steps
% Dt = spdiags(repmat([1 -1],p+1,1),[0 -1],p+1,p+1)/dt; % Explicit backward Euler time scheme (first-order accurate, conditionally stable)
% Dt = spdiags(repmat([1 -1],p+1,1),[1 0],p+1,p+1)/dt; % Implicit forward Euler time scheme (first-order accurate, unconditionally stable)
Dt = spdiags(repmat([1 -1],p+1,1),[1 -1],p+1,p+1)/(2*dt); % Implicit central-difference time scheme (second-order accurate, unconditionally stable)
tf = p*dt;
T = TIMEMODEL(0,tf,p);

%% First reduction step in space
if computePCAspace
    fprintf('\nPCA in space');
    t_PCA_space = tic;
    Rinit = min(r,N);
    if g<2^7
        load(fullfile(gridpathname,'data.mat'),'Y');
        mY = mean(Y,1);
        Yc = Y - repmat(mY,[N,1,1,1]); % Yc = Y - mY.*ones(N,1,1,1);
        clear Y
        V = zeros(r,Rinit,p+1);
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
        fprintf('\nTime step %2d/%2d : rank R = %d, error = %.3e for Y',t,p,Rt,errYct);
        
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
    fprintf('\n');
    
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
    load(fullfile(gridpathname,'PCA_space.mat'),'mY','sig','Z','R','errsvdYc','time_PCA_space');
    if g<2^7
        load(fullfile(gridpathname,'PCA_space.mat'),'V');
    end
end

%% Second reduction step in time for each coordinate
% if computePCAtime
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
%     fprintf('\n');
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
%     load(fullfile(gridpathname,'PCA_time.mat'),'s','W','Q','errsvdZc','time_PCA_time');
% end
    
%% Second reduction step in time
if computePCAtime
    fprintf('\nPCA in time');
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
    fprintf('\n');
    
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
    load(fullfile(gridpathname,'PCA_time.mat'),'s','W','Q','errsvdZc','time_PCA_time');
end

%% Post-processing data
if postProcess
    fprintf('\nPost-processing data');
    t_PostProcess = tic;
    ntau = 13; % number of post-processed variables
    ne = 9; % number of energy variables
    fprintf('\nn = %d post-processed variables',ntau);
    fprintf('\n  = %d energy variables',ne);
    
    mTau = zeros(1,ntau,m,p+1);
    if g<2^7
        Tau = zeros(N,ntau,m,p+1);
    end
    Qu = zeros(N,3,2,p+1);
    Qtau = zeros(N,ntau,2,p+1);
    Qe = zeros(N,ne,2,p+1);
    
    for t=0:p
        fprintf('\nTime step %2d/%2d',t,p);
        
        if g<2^7
            Yt = get_reduced_order_representation_at_step(mY,sig,V,s,W,X,R,t+1,'index',index,'gridpathname',gridpathname);
            if t==0
                Yt_old = get_reduced_order_representation_at_step(mY,sig,V,s,W,X,R,1,'index',index,'gridpathname',gridpathname);
            else
                Yt_old = get_reduced_order_representation_at_step(mY,sig,V,s,W,X,R,t,'index',index,'gridpathname',gridpathname);
            end
            if t==p
                Yt_new = get_reduced_order_representation_at_step(mY,sig,V,s,W,X,R,p+1,'index',index,'gridpathname',gridpathname);
            else
                Yt_new = get_reduced_order_representation_at_step(mY,sig,V,s,W,X,R,t+2,'index',index,'gridpathname',gridpathname);
            end
        else
            Yt = get_reduced_order_representation_at_step(mY,sig,[],s,W,X,R,t+1,'index',index,'gridpathname',gridpathname);
            if t==0
                Yt_old = get_reduced_order_representation_at_step(mY,sig,[],s,W,X,R,1,'index',index,'gridpathname',gridpathname);
            else
                Yt_old = get_reduced_order_representation_at_step(mY,sig,[],s,W,X,R,t,'index',index,'gridpathname',gridpathname);
            end
            if t==p
                Yt_new = get_reduced_order_representation_at_step(mY,sig,[],s,W,X,R,p+1,'index',index,'gridpathname',gridpathname);
            else
                Yt_new = get_reduced_order_representation_at_step(mY,sig,[],s,W,X,R,t+2,'index',index,'gridpathname',gridpathname);
            end
        end
        
        Yt = perm(reshape(Yt,[N,n,sx]));
        Yt_old = perm(reshape(Yt_old,[N,n,sx]));
        Yt_new = perm(reshape(Yt_new,[N,n,sx]));
        
        ut = Yt(1:3,:,:,:,:);
        ut_old = Yt_old(1:3,:,:,:,:);
        ut_new = Yt_new(1:3,:,:,:,:);
        
        Ct = Yt(4,:,:,:,:);
        Ct_old = Yt_old(4,:,:,:,:);
        Ct_new = Yt_new(4,:,:,:,:);
        
        % rhot = Ct*rho(2) + (1-Ct)*rho(1);
        rhot_old = Ct_old*rho(2) + (1-Ct_old)*rho(1);
        rhot_new = Ct_new*rho(2) + (1-Ct_new)*rho(1);
        % rhout = repmat(rhot,[3,ones(1,4)]).*ut;
        rhout_old = repmat(rhot_old,[3,ones(1,4)]).*ut_old;
        rhout_new = repmat(rhot_new,[3,ones(1,4)]).*ut_new;
        tauTimet = (rhout_new-rhout_old)/(2*dt);
        % u2t = dot(ut,ut,1);
        u2t_old = dot(ut_old,ut_old,1);
        u2t_new = dot(ut_new,ut_new,1);
        % Ek = 1/2*rhot.*u2t;
        Ek_old = 1/2*rhot_old.*u2t_old;
        Ek_new = 1/2*rhot_new.*u2t_new;
        energyKinTimet = (Ek_new-Ek_old)/(2*dt);
        
        rhot = Ct*rho(2) + (1-Ct)*rho(1);
        mut = Ct*mu(2) + (1-Ct)*mu(1);
        gradut = grad(ut,Dx);
        rhout = repmat(rhot,[3,ones(1,4)]).*ut;
        gradrhout = grad(rhout,Dx);
        St = (gradut+permute(gradut,[2,1,3:6]))/2;
        muSt = repmat(shiftdim(mut,-1),[3,3,ones(1,4)]).*St;
        gradCt = grad(Ct,Dx);
        ngradCt = normal(gradCt);
        kappa = div(ngradCt,Dx);
        u2t = dot(ut,ut,1);
        rhou2t = rhot.*u2t;
        
        % Post-processed variables
        divtauConvt = squeeze(sum(gradrhout.*repmat(shiftdim(ut,-1),[3,ones(1,5)]),2));
        % divtauConvt = div(permute(repmat(shiftdim(rhout,-1),[3,ones(1,5)]),[2,1,3:6]).*repmat(shiftdim(ut,-1),[3,ones(1,5)]),Dx);
        divtauDifft = div(2*muSt,Dx);
        tauSurft = sigma*repmat(shiftdim(kappa,-1),[3,ones(1,4)]).*gradCt;
        tauInterft = dot(ut,gradCt,1);
        energyConvt = shiftdim(div(repmat(rhou2t,[3,ones(1,4)]).*ut,Dx),-1);
        energyGravt = gravity.*rhout(2,:,:,:,:);
        energyPrest = zeros(1,g+1,g+1,g+1,N);
        energyPresDilt = zeros(1,g+1,g+1,g+1,N);
        energyKinSpacet = dot(rhout,grad(u2t/2,Dx),1);
        energyDifft = shiftdim(div(squeeze(dot(2*muSt,repmat(shiftdim(ut,-1),[3,ones(1,5)]),2)),Dx),-1);
        energyVisct = shiftdim(sum(sum(2*muSt.*gradut,1),2),1);
        energySurft = dot(tauSurft,ut,1);
        
        Taut = cat(1,tauTimet,divtauConvt,divtauDifft,tauSurft,tauInterft);
        Et = cat(1,energyKinTimet,energyConvt,energyGravt,energyPrest,energyPresDilt,energyKinSpacet,energyDifft,energyVisct,energySurft);
        
        % Quantities of interest: spatial average in each phase
        Qut = cat(2,reshape(trapz(x,trapz(x,trapz(x,repmat(1-Ct,[3,ones(1,4)]).*ut,2),3),4),[3,1,N]),...
            reshape(trapz(x,trapz(x,trapz(x,repmat(Ct,[3,ones(1,4)]).*ut,2),3),4),[3,1,N]));
        Qtaut = cat(2,reshape(trapz(x,trapz(x,trapz(x,repmat(1-Ct,[ntau,ones(1,4)]).*Taut,2),3),4),[ntau,1,N]),...
            reshape(trapz(x,trapz(x,trapz(x,repmat(Ct,[ntau,ones(1,4)]).*Taut,2),3),4),[ntau,1,N]));
        Qet = cat(2,reshape(trapz(x,trapz(x,trapz(x,repmat(1-Ct,[ne,ones(1,4)]).*Et,2),3),4),[ne,1,N]),...
            reshape(trapz(x,trapz(x,trapz(x,repmat(Ct,[ne,ones(1,4)]).*Et,2),3),4),[ne,1,N]));
        
        Qut = shiftdim(Qut,2);
        Qtaut = shiftdim(Qtaut,2);
        Qet = shiftdim(Qet,2);
        
        Qu(:,:,:,t+1) = Qut;
        Qtau(:,:,:,t+1) = Qtaut;
        Qe(:,:,:,t+1) = Qet;
        
        Taut = iperm(Taut);
        Taut = Taut(:,:,:);
        
        mTaut = mean(Taut,1);
        mTau(:,:,:,t+1) = mTaut;
        
        if g<2^7
            Tau(:,:,:,t+1) = Taut;
        else
            save(fullfile(gridpathname,['data_post_t' num2str(t) '.mat']),'Taut');
        end
    end
    fprintf('\n');
    
    mQu = mean(Qu,1);
    mQtau = mean(Qtau,1);
    mQe = mean(Qe,1);
    
    Quc = Qu - repmat(mQu,N,1,1,1);
    Qtauc = Qtau - repmat(mQtau,N,1,1,1);
    Qec = Qe - repmat(mQe,N,1,1,1);
    
    stdQu = zeros(3*2,p+1);
    stdQtau = zeros(ntau*2,p+1);
    stdQe = zeros(ne*2,p+1);
    for t=0:p
        Qut = Qu(:,:,:,t+1);
        Qtaut = Qtau(:,:,:,t+1);
        Qet = Qe(:,:,:,t+1);
        stdQu(:,t+1) = std(Qut(:,:));
        stdQtau(:,t+1) = std(Qtaut(:,:));
        stdQe(:,t+1) = std(Qet(:,:));
        % Quct = Quc(:,:,:,t+1);
        % Qtauct = Qtauc(:,:,:,t+1);
        % Qect = Qec(:,:,:,t+1);
        % stdQu(:,t+1) = sqrt(1/(N-1)*sum(Quct(:,:).^2));
        % stdQtau(:,t+1) = sqrt(1/(N-1)*sum(Qtauct(:,:).^2));
        % stdQe(:,t+1) = sqrt(1/(N-1)*sum(Qect(:,:).^2));
    end
    
    RQu = zeros(3*2,3*2,p+1,p+1);
    RQtau = zeros(ntau*2,ntau*2,p+1,p+1);
    RQe = zeros(ne*2,ne*2,p+1,p+1);
    for t=0:p
        Quct = Quc(:,:,:,t+1);
        Qtauct = Qtauc(:,:,:,t+1);
        Qect = Qec(:,:,:,t+1);
        Quct = Quct(:,:)./repmat(stdQu(:,t+1)',[N,1]);
        Qtauct = Qtauct(:,:)./repmat(stdQtau(:,t+1)',[N,1]);
        Qect = Qect(:,:)./repmat(stdQe(:,t+1)',[N,1]);
        for tt=0:p
            Quctt = Quc(:,:,:,tt+1);
            Qtauctt = Qtauc(:,:,:,tt+1);
            Qectt = Qec(:,:,:,tt+1);
            Quctt = Quctt(:,:)./repmat(stdQu(:,tt+1)',[N,1]);
            Qtauctt = Qtauctt(:,:)./repmat(stdQtau(:,tt+1)',[N,1]);
            Qectt = Qectt(:,:)./repmat(stdQe(:,tt+1)',[N,1]);
            RQu(:,:,t+1,tt+1) = 1/(N-1)*Quct(:,:)'*Quctt(:,:);
            RQtau(:,:,t+1,tt+1) = 1/(N-1)*Qtauct(:,:)'*Qtauctt(:,:);
            RQe(:,:,t+1,tt+1) = 1/(N-1)*Qect(:,:)'*Qectt(:,:);
        end
    end
    
    IQu = zeros(3*2,3*2,p+1);
    IQtau = zeros(ntau*2,ntau*2,p+1);
    IQe = zeros(ne*2,ne*2,p+1);
    for t=0:p-1
        funu = zeros(3*2,3*2,p-t+1);
        funtau = zeros(ntau*2,ntau*2,p-t+1);
        fune = zeros(ne*2,ne*2,p-t+1);
        for tt=0:p-t
            funu(:,:,tt+1) = RQu(:,:,tt+t+1,tt+1);
            funtau(:,:,tt+1) = RQtau(:,:,tt+t+1,tt+1);
            fune(:,:,tt+1) = RQe(:,:,tt+t+1,tt+1);
        end
        IQut = 1/((p-t)*dt)*trapz((0:p-t)*dt,funu,3);
        IQtaut = 1/((p-t)*dt)*trapz((0:p-t)*dt,funtau,3);
        IQet = 1/((p-t)*dt)*trapz((0:p-t)*dt,fune,3);
        IQu(:,:,t+1) = IQut;
        IQtau(:,:,t+1) = IQtaut;
        IQe(:,:,t+1) = IQet;
    end
    IQu = reshape(IQu,3,2,3,2,p+1);
    IQtau = reshape(IQtau,ntau,2,ntau,2,p+1);
    IQe = reshape(IQe,ne,2,ne,2,p+1);
    
    time_PostProcess = toc(t_PostProcess);
    fprintf('\nelapsed time = %f s',time_PostProcess);
    fprintf('\n');
    
    save(fullfile(gridpathname,'data_post.mat'),'mTau','ntau','ne',...
        'Qu','Qtau','Qe','mQu','mQtau','mQe','IQu','IQtau','IQe','time_PostProcess');
    if g<2^7
        save(fullfile(gridpathname,'data_post.mat'),'Tau','-append');
    end
else
    load(fullfile(gridpathname,'data_post.mat'),'mTau','ntau','ne',...
        'Qu','Qtau','Qe','mQu','mQtau','mQe','IQu','IQtau','IQe','time_PostProcess');
    if g<2^7
        load(fullfile(gridpathname,'data_post.mat'),'Tau');
    end
end

%% Applying filter
if applyFilter
    if g==gref
        fprintf('\nApplying filter');
        t_Filter = tic;
        
        if g<2^7
            load(fullfile(gridpathname,'data.mat'),'Y');
            load(fullfile(gridpathname,'data_post.mat'),'Tau');
        end
        
        for ig=1:ng-1
            gbar = gset(ng-ig);
            fprintf(['\nFiltered Grid' num2str(gbar)]);
            
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
            
            mYbar = zeros(1,n,m,p+1);
            mTaubar = zeros(1,ntau,m,p+1);
            if g<2^7
                Ybar = zeros(N,n,m,p+1);
                Taubar = zeros(N,ntau,m,p+1);
            end
            
            for t=0:p
                fprintf('\nTime %2d/%2d',t,p);
                
                if g<2^7
                    Yt = Y(:,:,:,t+1);
                    Taut = Tau(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
                    load(fullfile(gridpathname,['data_post_t' num2str(t) '.mat']),'Taut');
                end
                YTaut = cat(2,Yt,Taut);
                
                YTaut = reshape(YTaut,[N,n+ntau,sx]);
                
                YTaubart = YTaut;
                for l=1:N
                    for i=1:(n+ntau)
                        YTaubartl = squeeze(YTaubart(l,i,:,:,:));
                        switch filterType
                            case {'box','mean','average'}
                                YTaubartl = imboxfilt3(YTaubartl,filterSize,'Padding','replicate','NormalizationFactor',1/prod(filterSize));
                                % YTaubartl = imfilter(YTaubartl,h,'replicate');
                            case {'linear','trapz'}
                                YTaubartl = imfilter(YTaubartl,h,'replicate');
                        end
                        YTaubart(l,i,:,:,:) = YTaubartl;
                    end
                end
                YTaubart = YTaubart(:,:,:);
                
                Ybart = YTaubart(:,1:4,:);
                Taubart = YTaubart(:,5:end,:);
                
                mYbart = mean(Ybart,1);
                mTaubart = mean(Taubart,1);
                
                mYbar(:,:,:,t+1) = mYbart;
                mTaubar(:,:,:,t+1) = mTaubart;
                
                if g<2^7
                    Ybar(:,:,:,t+1) = Ybart;
                    Taubar(:,:,:,t+1) = Taubart;
                else
                    save(fullfile(gridpathname,['data_filtered_grid' num2str(gbar) '_t' num2str(t) '.mat']),'Ybart');
                    save(fullfile(gridpathname,['data_post_filtered_grid' num2str(gbar) '_t' num2str(t) '.mat']),'Taubart');
                end
            end
            fprintf('\n');
            
            time_Filter = toc(t_Filter);
            fprintf('\nelapsed time = %f s',time_Filter);
            fprintf('\n');
            
            save(fullfile(gridpathname,['data_filtered_grid' num2str(gbar) '.mat']),'mYbar');
            save(fullfile(gridpathname,['data_post_filtered_grid' num2str(gbar) '.mat']),'mTaubar');
            if g<2^7
                save(fullfile(gridpathname,['data_filtered_grid' num2str(gbar) '.mat']),'Ybar','-append');
                save(fullfile(gridpathname,['data_post_filtered_grid' num2str(gbar) '.mat']),'Taubar','-append');
            end
        end
    else
        fprintf('\nExtracting filtered data');
        gridrefname = ['Grid' num2str(gref)];
        fprintf(['\nReference ' gridrefname]);
        gridrefpathname = fullfile(pathname,gridrefname);
        
        nbelemref = [gref,gref,gref];
        Mref = build_model(D,'nbelem',nbelemref,'elemtype',elemtype);
        coordref = getcoord(getnode(Mref));
        Mref = setnode(Mref,NODE(coordref(:,[2,1,3])));
        Mref = final(Mref,DDL(DDLVECT('U',Mref.syscoord)));
        Mrefscal = final(Mref,DDL('C'));
        
        P = calc_P(Mref,M);
        Pscal = calc_P(Mrefscal,Mscal);
        
        ng = length(gset);
        kg = find(gset==g);
        ig = ng-kg;   
        load(fullfile(gridrefpathname,['data_filtered_grid' num2str(g) '.mat']),'mYbar');
        load(fullfile(gridrefpathname,['data_post_filtered_grid' num2str(g) '.mat']),'mTaubar');
        if gref<2^7
            load(fullfile(gridrefpathname,['data_filtered_grid' num2str(g) '.mat']),'Ybar');
            load(fullfile(gridrefpathname,['data_post_filtered_grid' num2str(g) '.mat']),'Taubar');
        end
        
        mYbar_new = zeros(1,n,m,p+1);
        mTaubar_new = zeros(1,ntau,m,p+1);
        if g<2^7
            Ybar_new = zeros(N,n,m,p+1);
            Taubar_new = zeros(N,ntau,m,p+1);
        end
        for t=0:p
            fprintf('\nTime %2d/%2d',t,p);
            
            mYbart = mYbar(:,:,:,t+1);
            mTaubart = mTaubar(:,:,:,t+1);
            if gref<2^7
                Ybart = Ybar(:,:,:,t+1);
                Taubart = Taubar(:,:,:,t+1);
            else
                load(fullfile(gridrefpathname,['data_filtered_grid' num2str(g) '_t' num2str(t) '.mat']),'Ybart');
                load(fullfile(gridrefpathname,['data_post_filtered_grid' num2str(g) '_t' num2str(t) '.mat']),'Taubart');
            end
            
            mUbart = mYbart(1,1:3,:);
            mCbart = mYbart(1,4,:);
            mtauTimebart = mTaubart(1,1:3,:);
            mdivtauConvbart = mTaubart(1,4:6,:);
            mdivtauDiffbart = mTaubart(1,7:9,:);
            mtauSurfbart = mTaubart(1,10:12,:);
            mtauInterfbart = mTaubart(1,13,:);
            
            Ubart = Ybart(:,1:3,:);
            Cbart = Ybart(:,4,:);
            tauTimebart = Taubart(:,1:3,:);
            divtauConvbart = Taubart(:,4:6,:);
            divtauDiffbart = Taubart(:,7:9,:);
            tauSurfbart = Taubart(:,10:12,:);
            tauInterfbart = Taubart(:,13,:);
            
            mUbart = P*mUbart(:);
            mCbart = Pscal*mCbart(:);
            mtauTimebart = P*mtauTimebart(:);
            mdivtauConvbart = P*mdivtauConvbart(:);
            mdivtauDiffbart = P*mdivtauDiffbart(:);
            mtauSurfbart = P*mtauSurfbart(:);
            mtauInterfbart = Pscal*mtauInterfbart(:);
            
            Ubart = P*Ubart(:,:)';
            Cbart = Pscal*Cbart(:,:)';
            tauTimebart = P*tauTimebart(:,:)';
            divtauConvbart = P*divtauConvbart(:,:)';
            divtauDiffbart = P*divtauDiffbart(:,:)';
            tauSurfbart = P*tauSurfbart(:,:)';
            tauInterfbart = Pscal*tauInterfbart(:,:)';
            
            mUbart = reshape(mUbart,[3,m]);
            mCbart = reshape(mCbart,[1,m]);
            mtauTimebart = reshape(mtauTimebart,[3,m]);
            mdivtauConvbart = reshape(mdivtauConvbart,[3,m]);
            mdivtauDiffbart = reshape(mdivtauDiffbart,[3,m]);
            mtauSurfbart = reshape(mtauSurfbart,[3,m]);
            mtauInterfbart = reshape(mtauInterfbart,[1,m]);
            
            Ubart = reshape(Ubart',[N,3,m]);
            Cbart = reshape(Cbart',[N,1,m]);
            tauTimebart = reshape(tauTimebart',[N,3,m]);
            divtauConvbart = reshape(divtauConvbart',[N,3,m]);
            divtauDiffbart = reshape(divtauDiffbart',[N,3,m]);
            tauSurfbart = reshape(tauSurfbart',[N,3,m]);
            tauInterfbart = reshape(tauInterfbart',[N,1,m]);
            
            mYbart = cat(1,mUbart,mCbart);
            mTaubart = cat(1,mtauTimebart,mdivtauConvbart,mdivtauDiffbart,mtauSurfbart,mtauInterfbart);
            
            Ybart = cat(2,Ubart,Cbart);
            Taubart = cat(2,tauTimebart,divtauConvbart,divtauDiffbart,tauSurfbart,tauInterfbart);
            
            mYbar_new(1,:,:,t+1) = mYbart;
            mTaubar_new(1,:,:,t+1) = mTaubart;
            if g<2^7
                Ybar_new(:,:,:,t+1) = Ybart;
                Taubar_new(:,:,:,t+1) = Taubart;
            else
                save(fullfile(gridpathname,['data_filtered_t' num2str(t) '.mat']),'Ybart');
                save(fullfile(gridpathname,['data_post_filtered_t' num2str(t) '.mat']),'Taubart');
            end
        end
        
        mYbar = mYbar_new;
        mTaubar = mTaubar_new;
        if g<2^7
            Ybar = Ybar_new;
            Taubar = Taubar_new;
        end
        
        save(fullfile(gridpathname,'data_filtered.mat'),'mYbar');
        save(fullfile(gridpathname,'data_post_filtered.mat'),'mTaubar');
        if g<2^7
            save(fullfile(gridpathname,'data_filtered.mat'),'Ybar','-append');
            save(fullfile(gridpathname,'data_post_filtered.mat'),'Taubar','-append');
        end
    end
end

%% Outputs
% Display eigenvalues
if displayEigenvalues
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
mtauTime = reshape(mTau(1,1:3,:,:),[3*m,p+1]);
mdivtauConv = reshape(mTau(1,4:6,:,:),[3*m,p+1]);
mdivtauDiff = reshape(mTau(1,7:9,:,:),[3*m,p+1]);
mtauSurf = reshape(mTau(1,10:12,:,:),[3*m,p+1]);
mtauInterf = reshape(mTau(1,13,:,:),[m,p+1]);

mU = TIMEMATRIX(mU,T);
mC = TIMEMATRIX(mC,T);
mtauTime = TIMEMATRIX(mtauTime,T);
mdivtauConv = TIMEMATRIX(mdivtauConv,T);
mdivtauDiff = TIMEMATRIX(mdivtauDiff,T);
mtauSurf = TIMEMATRIX(mtauSurf,T);
mtauInterf = TIMEMATRIX(mtauInterf,T);

% Mean filtered solution
if g~=gref
    mUbar = reshape(mYbar(1,1:3,:,:),[3*m,p+1]);
    mCbar = reshape(mYbar(1,4,:,:),[m,p+1]);
    mtauTimebar = reshape(mTaubar(1,1:3,:,:),[3*m,p+1]);
    mdivtauConvbar = reshape(mTaubar(1,4:6,:,:),[3*m,p+1]);
    mdivtauDiffbar = reshape(mTaubar(1,7:9,:,:),[3*m,p+1]);
    mtauSurfbar = reshape(mTaubar(1,10:12,:,:),[3*m,p+1]);
    mtauInterfbar = reshape(mTaubar(1,13,:,:),[m,p+1]);
    
    mUbar = TIMEMATRIX(mUbar,T);
    mCbar = TIMEMATRIX(mCbar,T);
    mtauTimebar = TIMEMATRIX(mtauTimebar,T);
    mdivtauConvbar = TIMEMATRIX(mdivtauConvbar,T);
    mdivtauDiffbar = TIMEMATRIX(mdivtauDiffbar,T);
    mtauSurfbar = TIMEMATRIX(mtauSurfbar,T);
    mtauInterfbar = TIMEMATRIX(mtauInterfbar,T);
end

% Variance of solution
vu = zeros(3*m,p+1);
vC = zeros(m,p+1);
vtauTime = zeros(3*m,p+1);
vdivtauConv = zeros(3*m,p+1);
vdivtauDiff = zeros(3*m,p+1);
vtauSurf = zeros(3*m,p+1);
vtauInterf = zeros(m,p+1);

switch index
    case 'coord'
        W = permute(reshape(W',[Q,Rmax,p+1]),[2,1,3]);
        % Z_approx = reshape(Zc_approx',[N,Rmax,p+1]);
    case 'time'
        W = permute(reshape(W',[Q,p+1,Rmax]),[3,1,2]);
        % Z_approx = permute(reshape(Zc_approx',[N,p+1,Rmax]),[1,3,2]);
end

for t=0:p
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
    
    mTaut = mTau(:,:,:,t+1);
    if g<2^7
        Taut = Tau(:,:,:,t+1);
    else
        load(fullfile(gridpathname,['data_post_t' num2str(t) '.mat']),'Taut');
    end
    Tauct = Taut - repmat(mTaut,[N,1,1]); % Tauct = Taut - mTaut.*ones(N,1,1);
    
    tauTimet = Tauct(:,1:3,:);
    divtauConvt = Tauct(:,4:6,:);
    divtauDifft = Tauct(:,7:9,:);
    tauSurft = Tauct(:,10:12,:);
    tauInterft = Tauct(:,13,:);
    
    indut = repmat((0:m-1)*n,[3,1])+repmat((1:3)',[1,m]);
    indCt = (0:m-1)*n+4;
    
    vut = vYt_approx(indut(:));
    vCt = vYt_approx(indCt(:));
    vtauTimet = 1/(N-1)*sum(tauTimet(:,:).^2)';
    vdivtauConvt = 1/(N-1)*sum(divtauConvt(:,:).^2)';
    vdivtauDifft = 1/(N-1)*sum(divtauDifft(:,:).^2)';
    vtauSurft = 1/(N-1)*sum(tauSurft(:,:).^2)';
    vtauInterft = 1/(N-1)*sum(tauInterft(:,:).^2)';
    
    vu(:,t+1) = vut;
    vC(:,t+1) = vCt;
    vtauTime(:,t+1) = vtauTimet;
    vdivtauConv(:,t+1) = vdivtauConvt;
    vdivtauDiff(:,t+1) = vdivtauDifft;
    vtauSurf(:,t+1) = vtauSurft;
    vtauInterf(:,t+1) = vtauInterft;
end
fprintf('\n');

vu = TIMEMATRIX(vu,T);
vC = TIMEMATRIX(vC,T);
vtauTime = TIMEMATRIX(vtauTime,T);
vdivtauConv = TIMEMATRIX(vdivtauConv,T);
vdivtauDiff = TIMEMATRIX(vdivtauDiff,T);
vtauSurf = TIMEMATRIX(vtauSurf,T);
vtauInterf = TIMEMATRIX(vtauInterf,T);

% Display quantities of interest
if displayQoI
    % Mean
    mQu = reshape(mQu(1,1:3,:,:),[3,2,p+1]);
    mQtauTime = reshape(mQtau(1,1:3,:,:),[3,2,p+1]);
    mQdivtauConv = reshape(mQtau(1,4:6,:,:),[3,2,p+1]);
    mQdivtauDiff = reshape(mQtau(1,7:9,:,:),[3,2,p+1]);
    mQtauSurf = reshape(mQtau(1,10:12,:,:),[3,2,p+1]);
    mQtauInterf = reshape(mQtau(1,13,:,:),[1,2,p+1]);
    mQenergyKinTime = reshape(mQe(1,1,:,:),[1,2,p+1]);
    mQenergyConv = reshape(mQe(1,2,:,:),[1,2,p+1]);
    mQenergyGrav = reshape(mQe(1,3,:,:),[1,2,p+1]);
    mQenergyPres = reshape(mQe(1,4,:,:),[1,2,p+1]);
    mQenergyPresDil = reshape(mQe(1,5,:,:),[1,2,p+1]);
    mQenergyKinSpace = reshape(mQe(1,6,:,:),[1,2,p+1]);
    mQenergyDiff = reshape(mQe(1,7,:,:),[1,2,p+1]);
    mQenergyVisc = reshape(mQe(1,8,:,:),[1,2,p+1]);
    mQenergySurf = reshape(mQe(1,9,:,:),[1,2,p+1]);
    
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
    
    % Correlation
    IQtauTime = IQtau(1:3,:,1:3,:,:);
    IQdivtauConv = IQtau(4:6,:,4:6,:,:);
    IQdivtauDiff = IQtau(7:9,:,7:9,:,:);
    IQtauSurf = IQtau(10:12,:,10:12,:,:);
    IQtauInterf = IQtau(13,:,13,:,:);
    IQenergyKinTime = IQe(1,:,1,:,:);
    IQenergyConv = IQe(2,:,2,:,:);
    IQenergyGrav = IQe(3,:,3,:,:);
    IQenergyPres = IQe(4,:,4,:,:);
    IQenergyPresDil = IQe(5,:,5,:,:);
    IQenergyKinSpace = IQe(6,:,6,:,:);
    IQenergyDiff = IQe(7,:,7,:,:);
    IQenergyVisc = IQe(8,:,8,:,:);
    IQenergySurf = IQe(9,:,9,:,:);
    
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

% Display solution
if displaySolution
    t = 11;
    ampl = 0;
    % ampl = getsize(M)/max(abs(mut))/5;
    az_z = -37.5; % azimuth for the default 3D view with vertical z-axis
    el_z = 30; % elevation for the default 3D view with vertical z-axis
    az = az_z-90; % azimuth for the default 3D view with vertical y-axis
    el = -el_z; % elevation for the default 3D view with vertical y-axis
    for i=1:3
        evolSolution(M,mU,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_mean_u' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
            'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        evolSolution(M,mtauTime,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_mean_tauTime' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
            'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        evolSolution(M,mdivtauConv,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_mean_divtauConv' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
            'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        evolSolution(M,mdivtauDiff,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_mean_divtauDiff' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
            'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        evolSolution(M,mtauSurf,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_mean_tauSurf' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
            'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        
%         evolSolution(M,vu,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_var_u' num2str(i)],'pathname',gridpathname,'FrameRate',framerate);,...
%             'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%         evolSolution(M,vtauTime,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_var_tauTime' num2str(i)],'pathname',gridpathname,'FrameRate',framerate);,...
%             'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%         evolSolution(M,vdivtauConv,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_var_divtauConv' num2str(i)],'pathname',gridpathname,'FrameRate',framerate);,...
%             'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%         evolSolution(M,vdivtauDiff,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_var_divtauDiff' num2str(i)],'pathname',gridpathname,'FrameRate',framerate);,...
%             'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%         evolSolution(M,vtauSurf,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_var_tauSurf' num2str(i)],'pathname',gridpathname,'FrameRate',framerate);,...
%             'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        
        figure('Name',['Mean of velocity u_' num2str(i)])
        clf
        plot_sol(M,getmatrixatstep(mU,t+1),'displ',i,'ampl',ampl);
        title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        view(az,el)
        camup([0 1 0])
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(gridpathname,['mean_u' num2str(i) '_t' num2str(t*100)],formats,renderer);
        
        figure('Name',['Mean of tauTime_' num2str(i)])
        clf
        plot_sol(M,getmatrixatstep(mtauTime,t+1),'displ',i,'ampl',ampl);
        title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        view(az,el)
        camup([0 1 0])
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(gridpathname,['mean_tauTime' num2str(i) '_t' num2str(t*100)],formats,renderer);
        
        figure('Name',['Mean of div(tauConv)_' num2str(i)])
        clf
        plot_sol(M,getmatrixatstep(mdivtauConv,t+1),'displ',i,'ampl',ampl);
        title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        view(az,el)
        camup([0 1 0])
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(gridpathname,['mean_divtauConv' num2str(i) '_t' num2str(t*100)],formats,renderer);
        
        figure('Name',['Mean of div(tauDiff)_' num2str(i)])
        clf
        plot_sol(M,getmatrixatstep(mdivtauDiff,t+1),'displ',i,'ampl',ampl);
        title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        view(az,el)
        camup([0 1 0])
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(gridpathname,['mean_divtauDiff' num2str(i) '_t' num2str(t*100)],formats,renderer);
        
        figure('Name',['Mean of tauSurf ' num2str(i)])
        clf
        plot_sol(M,getmatrixatstep(mtauSurf,t+1),'displ',i,'ampl',ampl);
        title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        view(az,el)
        camup([0 1 0])
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(gridpathname,['mean_tauSurf' num2str(i) '_t' num2str(t*100)],formats,renderer);
        
%         figure('Name',['Variance of velocity u_' num2str(i)])
%         clf
%         plot_sol(M,getmatrixatstep(vu,t+1),'displ',i,'ampl',ampl);
%         title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
%         colormap(cmap)
%         colorbar
%         axis on
%         box on
%         view(az,el)
%         camup([0 1 0])
%         set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%         mysaveas(gridpathname,['var_u' num2str(i) '_t' num2str(t*100)],formats,renderer);
%         
%         figure('Name',['Variance of tauTime_' num2str(i)])
%         clf
%         plot_sol(M,getmatrixatstep(vtauTime,t+1),'displ',i,'ampl',ampl);
%         title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
%         colormap(cmap)
%         colorbar
%         axis on
%         box on
%         view(az,el)
%         camup([0 1 0])
%         set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%         mysaveas(gridpathname,['var_tauTime' num2str(i) '_t' num2str(t*100)],formats,renderer);
%         
%         figure('Name',['Variance of div(tauConv)_' num2str(i)])
%         clf
%         plot_sol(M,getmatrixatstep(vdivtauConv,t+1),'displ',i,'ampl',ampl);
%         title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
%         colormap(cmap)
%         colorbar
%         axis on
%         box on
%         view(az,el)
%         camup([0 1 0])
%         set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%         mysaveas(gridpathname,['var_divtauConv' num2str(i) '_t' num2str(t*100)],formats,renderer);
%         
%         figure('Name',['Variance of div(tauDiff)_' num2str(i)])
%         clf
%         plot_sol(M,getmatrixatstep(vdivtauDiff,t+1),'displ',i,'ampl',ampl);
%         title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
%         colormap(cmap)
%         colorbar
%         axis on
%         box on
%         view(az,el)
%         camup([0 1 0])
%         set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%         mysaveas(gridpathname,['var_divtauDiff' num2str(i) '_t' num2str(t*100)],formats,renderer);
%         
%         figure('Name',['Variance of tauSurf_' num2str(i)])
%         clf
%         plot_sol(M,getmatrixatstep(vtauSurf,t+1),'displ',i,'ampl',ampl);
%         title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
%         colormap(cmap)
%         colorbar
%         axis on
%         box on
%         view(az,el)
%         camup([0 1 0])
%         set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%         mysaveas(gridpathname,['var_tauSurf' num2str(i) '_t' num2str(t*100)],formats,renderer);
    end
    
    evolSolution(Mscal,mC,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename','evol_mean_C','pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
    evolSolution(Mscal,mtauInterf,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename','evol_mean_tauInterf','pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
    
%     evolSolution(Mscal,vC,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename','evol_var_C','pathname',gridpathname,'FrameRate',framerate,...
%         'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%     evolSolution(Mscal,vtauInterf,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename','evol_var_tauInterf','pathname',gridpathname,'FrameRate',framerate,...
%         'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
    
    figure('Name','Mean of indicator function C')
    clf
    plot_sol(Mscal,getmatrixatstep(mC,t+1));
    title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['mean_C_t' num2str(t*100)],formats,renderer);
    
    figure('Name','Mean of tauInterf')
    clf
    plot_sol(Mscal,getmatrixatstep(mtauInterf,t+1));
    title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['mean_tauInterf_t' num2str(t*100)],formats,renderer);
    
%     figure('Name','Variance of indicator function C')
%     clf
%     plot_sol(Mscal,getmatrixatstep(vC,t+1));
%     title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
%     colormap(cmap)
%     colorbar
%     axis on
%     box on
%     view(az,el)
%     camup([0 1 0])
%     set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%     mysaveas(gridname,['var_C_t' num2str(t*100)],formats,renderer);
%     
%     figure('Name','Variance of tauInterf')
%     clf
%     plot_sol(Mscal,getmatrixatstep(vtauInterf,t+1));
%     title(['time ' num2str(t*dt,'%.2f') ' s'],'FontSize',fontsize)
%     colormap(cmap)
%     colorbar
%     axis on
%     box on
%     view(az,el)
%     camup([0 1 0])
%     set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%     mysaveas(gridname,['var_tauInterf_t' num2str(t*100)],formats,renderer);
end

for t=0:p
    mUt = getmatrixatstep(mU,t+1);
    mCt = getmatrixatstep(mC,t+1);
    mtauTimet = getmatrixatstep(mtauTime,t+1);
    mdivtauConvt = getmatrixatstep(mdivtauConv,t+1);
    mdivtauDifft = getmatrixatstep(mdivtauDiff,t+1);
    mtauSurft = getmatrixatstep(mtauSurf,t+1);
    mtauInterft = getmatrixatstep(mtauInterf,t+1);
    
    vut = getmatrixatstep(vu,t+1);
    vCt = getmatrixatstep(vC,t+1);
    vtauTimet = getmatrixatstep(vtauTime,t+1);
    vdivtauConvt = getmatrixatstep(vdivtauConv,t+1);
    vdivtauDifft = getmatrixatstep(vdivtauDiff,t+1);
    vtauSurft = getmatrixatstep(vtauSurf,t+1);
    vtauInterft = getmatrixatstep(vtauInterf,t+1);
    
    if g==gref
        write_vtk_mesh(M,{mUt,mCt,mtauTimet,mdivtauConvt,mdivtauDifft,mtauSurft,mtauInterft,...
            %menergyKinTimet,menergyConvt,menergyGravt,menergyPrest,menergyPresDilt,...
            %menergyKinSpacet,menergyDifft,menergyVisct,menergySurft...
            },[],...
            {'velocity','phase','tauTime','div(tauConv)','div(tauDiff)','tauSurf','tauInterf',...
            %'kinetic energy','convection energy','gravity energy','power of external pressure forces','pressure-dilatation energy transfer',...
            %'transport of gradient of kinetic energy','energy exchange with kinetic energy','power of external viscous stresses','capillary kinetic energy'...
            },[],...
            gridpathname,'diphasic_fluids_mean',1,t);
        write_vtk_mesh(M,{vut,vCt,vtauTimet,vdivtauConvt,vdivtauDifft,vtauSurft,vtauInterft...
            %venergyKinTimet,venergyConvt,venergyGravt,venergyPrest,venergyPresDilt,...
            %venergyKinSpacet,venergyDifft,venergyVisct,venergySurft...
            },[],...
            {'velocity','phase','tauTime','div(tauConv)','div(tauDiff)','tauSurf','tauInterf'...
            %'kinetic energy','convection energy','gravity energy','power of external pressure forces','pressure-dilatation energy transfer',...
            %'transport of gradient of kinetic energy','energy exchange with kinetic energy','power of external viscous stresses','capillary kinetic energy'...
            },[],...
            gridpathname,'diphasic_fluids_variance',1,t);
        
    else
        mUbart = getmatrixatstep(mUbar,t+1);
        mCbart = getmatrixatstep(mCbar,t+1);
        mtauTimebart = getmatrixatstep(mtauTimebar,t+1);
        mdivtauConvbart = getmatrixatstep(mdivtauConvbar,t+1);
        mdivtauDiffbart = getmatrixatstep(mdivtauDiffbar,t+1);
        mtauSurfbart = getmatrixatstep(mtauSurfbar,t+1);
        mtauInterfbart = getmatrixatstep(mtauInterfbar,t+1);
        
        write_vtk_mesh(M,{mUt,mCt,mtauTimet,mdivtauConvt,mdivtauDifft,mtauSurft,mtauInterft,...
            mUbart,mCbart,mtauTimebart,mdivtauConvbart,mdivtauDiffbart,mtauSurfbart,mtauInterfbart...
            %menergyKinTimet,menergyConvt,menergyGravt,menergyPrest,menergyPresDilt,...
            %menergyKinSpacet,menergyDifft,menergyVisct,menergySurft...
            },[],...
            {'velocity','phase','tauTime','div(tauConv)','div(tauDiff)','tauSurf','tauInterf',...
            'velocity filtered','phase filtered','tauTime filtered','div(tauConv) filtered','div(tauDiff) filtered','tauSurf filtered','tauInterf filtered'...
            %'kinetic energy','convection energy','gravity energy','power of external pressure forces','pressure-dilatation energy transfer',...
            %'transport of gradient of kinetic energy','energy exchange with kinetic energy','power of external viscous stresses','capillary kinetic energy'...
            },[],...
            gridpathname,'diphasic_fluids_mean',1,t);
        write_vtk_mesh(M,{vut,vCt,vtauTimet,vdivtauConvt,vdivtauDifft,vtauSurft,vtauInterft...
            %venergyKinTimet,venergyConvt,venergyGravt,venergyPrest,venergyPresDilt,...
            %venergyKinSpacet,venergyDifft,venergyVisct,venergySurft...
            },[],...
            {'velocity','phase','tauTime','div(tauConv)','div(tauDiff)','tauSurf','tauInterf'...
            %'kinetic energy','convection energy','gravity energy','power of external pressure forces','pressure-dilatation energy transfer',...
            %'transport of gradient of kinetic energy','energy exchange with kinetic energy','power of external viscous stresses','capillary kinetic energy'...
            },[],...
            gridpathname,'diphasic_fluids_variance',1,t);
    end
end
make_pvd_file(gridpathname,'diphasic_fluids_mean',1,p+1);
make_pvd_file(gridpathname,'diphasic_fluids_variance',1,p+1);

time_Total = toc(t_Total);
fprintf('Elapsed time = %f s\n',time_Total);
