clc
clearvars
close all

performPCA = true;
performPCAspace = true;
performPCAtime = true;
constructMesh = false;

displayEigenvalues = true;
displayCovariance = false;

cmap = 'default';
% cmap = flipud(gray);
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

filterType = 'box'; % 3D filter type ('box' or 'mean' or 'average', 'linear' or 'trapz')

% Spatial grid size
gset = 2.^(4:8); % set of spatial grid sizes
g = gset(1); % current spatial grid size
gref = gset(end); % reference spatial grid size
ng = length(gset); % number of spatial grid sizes

t_Total = tic;
fprintf('Grid %d\n',g);
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
dim = 3; % spatial dimension
L = 1; % domain size (m)
sx = repmat(g+1,[1,dim]); % spatial dimensions
dx = L/g; % spatial step (m)
x = linspace(0,L,g+1); % spatial discretization in each spatial dimension
e = ones(g+1,1);
% Dxb = spdiags([-e e],[-1 0],g+1,g+1)/dx; % Explicit backward-difference spatial scheme (first-order accurate, conditionally stable)
% Dxf = spdiags([-e e],[0 1],g+1,g+1)/dx; % Implicit forward-difference spatial scheme (first-order accurate, unconditionally stable)
Dx = spdiags([-e e],[-1 1],g+1,g+1)/(2*dx); % Implicit central-difference spatial scheme (second-order accurate, unconditionally stable)
Dx(1,[1 2]) = [-1 1]/dx; Dx(end,[end-1 end]) = [-1 1]/dx; % Single-sided differences at boundary points
DxN = Dx; DxN(1,:) = 0; DxN(end,:) = 0; % homogeneous Neumann BCs
Dxx = spdiags([e -2*e e],-1:1,g+1,g+1)/(dx^2); % Implicit central-difference spatial scheme (second-order accurate, unconditionally stable)
DxxN = Dxx + sparse([1 g+1],[2 g],1,g+1,g+1)/(dx^2);

% Time scheme
dt = 5e-3; % computing time step (s)
dt = 100*dt; % physical time step stored every 100 computing time steps
% Dtb = spdiags(repmat([1 -1],p+1,1),[0 -1],p+1,p+1)/dt; Dt(1,1) = 1/(2*dt); % Explicit backward Euler time scheme (first-order accurate, conditionally stable)
% Dtf = spdiags(repmat([1 -1],p+1,1),[1 0],p+1,p+1)/dt; Dt(end,end) = 1/(2*dt); % Implicit forward Euler time scheme (first-order accurate, unconditionally stable)
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

%% Spatial mesh
if constructMesh
    fprintf('\nConstructing spatial mesh');
    t_Mesh = tic;
    
    D = DOMAIN(dim,zeros(1,dim),L*ones(1,dim));
    elemtype = 'CUB8';
    nbelem = repmat(g,[1,dim]);
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

%% Samples of random solution
fprintf('\nSaving samples');
t_save = tic;
for t=0:p
    t_savet = tic;
    if g<2^7
        Yt = Y(:,:,:,t+1);
    else
        load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
    end
    Ut = reshape(Yt(:,1:dim,:),[N,dim*m]);
    Ct = reshape(Yt(:,dim+1,:),[N,m]);
    parfor l=1:N
        Ult = Ut(l,:);
        Clt = Ct(l,:);
        fields = {Ult,Clt};
        fieldnames = {'velocity','phase'};
        write_vtk_mesh(M,fields,[],fieldnames,[],gridpathname,['diphasic_fluids_grid' num2str(g) '_sample_' num2str(l)],1,t);
    end
    time_savet = toc(t_savet);
    fprintf('\nTime step %2d/%2d : elapsed time = %f s',t,p,time_savet);
end
for l=1:N
    make_pvd_file(gridpathname,['diphasic_fluids_grid' num2str(g) '_sample_' num2str(l)],1,p+1);
end
time_save = toc(t_save);
fprintf('\nelapsed time = %f s',time_save);
fprintf('\n');

%% Single PCA
if performPCA
    fprintf('\nPerforming PCA');
    t_PCA = tic;
    r = n*m*(p+1);
    % Rinit = min(r,N);
    if g>=2^7
        Y = zeros(N,n,m,p+1);
        for t=0:p
            load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
            Y(:,:,:,t+1) = Yt;
            clear Yt
        end
    end
    [Y,Ya,Yb] = scaling(Y);
    mY = mean(Y,1);
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        Yc = bsxfun(@minus,Y,mY);
    else
        Yc = Y - mY;
    end
    % Yc = Y - repmat(mY,[N,1,1,1]); % Yc = Y - mY.*ones(N,1,1,1);
%     clear Y
    
    Yc = Yc(:,:)';
    [Phi,Sig,X,errsvdYc] = svdtruncate(Yc,tolsvdYc);
    Sig = Sig/sqrt(N-1);
    sig = diag(Sig);
    X = X*sqrt(N-1);
    % if verLessThan('matlab','9.1') % compatibility (<R2016b)
    %     Yc_approx = Phi*Sig*X';
    % else
    %     Yc_approx = Phi*(sig.*X');
    % end
    
    % [coeff,score,latent] = pca(Yc');
    % Phi = coeff;
    % sig = sqrt(latent);
    % Sig = diag(sig);
    % if verLessThan('matlab','9.1') % compatibility (<R2016b)
    %     X = score*diag(1./sig);
    % else
    %     X = score./sig';
    % end
    % % Yc_approx = coeff*score';
    
    % errYc = norm(Yc_approx-Yc)/norm(Yc);
    errYc = errsvdYc(end);
    R = length(sig);
    fprintf('\nrank R = %d, error = %.3e for Y',R,errYc);
    
    % norm2Yc = sum(var(Yc,0,2));
    % norm2Yc_approx = sum(var(Yc_approx,0,2));
    % err2Yc = errYc^2*norm2Yc;
    % err2Yc = norm2Yc-norm2Yc_approx;
    % sigf = svdtruncate(Yc,eps);
    % sigf = sigf/sqrt(N-1);
    % err2Yc = sum(sigf.^2)-sum(sig.^2);
    
    % mX = mean(X,1)';
    % CX = cov(X); % CX = 1/(N-1)*X'*X;
    % norm(mX)
    % norm(CX-eye(R))
    % norm(Phi'*Phi-eye(R))
    
    %if verLessThan('matlab','9.1') % compatibility (<R2016b)
    %    CY_approx = Phi*Sig.^2*Phi';
    %else
    %    CY_approx = Phi*(sig.^2.*Phi');
    %end
    %CY = cov(Yc'); % CY = 1/(N-1)*Yc*Yc';
    %errCY = norm(CY_approx-CY)/norm(CY);
    %fprintf('\n                          error = %.3e for CY',errCY);
    
    time_PCA = toc(t_PCA);
    fprintf('\nelapsed time = %f s',time_PCA);
    fprintf('\n');
    
    save(fullfile(gridpathname,'scaling.mat'),'Ya','Yb');
    save(fullfile(gridpathname,'data_PCA.mat'),'mY','Phi','sig','X','R','errsvdYc','time_PCA');
else
    fprintf('\nLoading DNS data from PCA');
    t_load = tic;
    load(fullfile(gridpathname,'scaling.mat'),'Ya','Yb');
    load(fullfile(gridpathname,'data_PCA.mat'),'mY','Phi','sig','X','R','errsvdYc','time_PCA');
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
end

%% Display eigenvalues
if displayEigenvalues
    figure('Name','Evolution of eigenvalues w.r.t order')
    clf
    R = length(sig);
    semilogy(1:R,sig(:).^2,'LineStyle','-','Color','b','LineWidth',1);
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\lambda_{\alpha}$','Interpreter',interpreter)
    mysaveas(gridpathname,'eigenvalues_CY_order',formats,renderer);
    mymatlab2tikz(gridpathname,'eigenvalues_CY_order_single_PCA.tex');
    
    figure('Name','Evolution of errors')
    clf
    R = length(sig);
    semilogy(1:R,errsvdYc.^2,'LineStyle','-','Color','b','LineWidth',1);
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$R$','Interpreter',interpreter)
    ylabel('$\varepsilon_{Y}(R)$','Interpreter',interpreter)
    mysaveas(gridpathname,'error_svdYc',formats,renderer);
    mymatlab2tikz(gridpathname,'error_svdYc_single_PCA.tex');
end

%% Samples of single PCA representation
fprintf('\nSaving samples of single PCA representation');
t_save = tic;
for t=0:p
    t_savet = tic;
    Yt = get_single_PCA_at_step(mY,Phi,sig,X,p,t+1,'gridpathname',gridpathname);
    Yat = Ya(1,:,:,t+1);
    Ybt = Yb(1,:,:,t+1);
    Yt = unscaling(Yt,Yat,Ybt);
    clear Yat Ybt
    
    Ut = reshape(Yt(:,1:dim,:),[N,dim*m]);
    Ct = reshape(Yt(:,dim+1,:),[N,m]);
    parfor l=1:N
        Ult = Ut(l,:);
        Clt = Ct(l,:);
        fields = {Ult,Clt};
        fieldnames = {'velocity','phase'};
        write_vtk_mesh(M,fields,[],fieldnames,[],gridpathname,['diphasic_fluids_grid' num2str(g) '_sample_' num2str(l) '_single_PCA'],1,t);
    end
    time_savet = toc(t_savet);
    fprintf('\nTime step %2d/%2d : elapsed time = %f s',t,p,time_savet);
end
for l=1:N
    make_pvd_file(gridpathname,['diphasic_fluids_grid' num2str(g) '_sample_' num2str(l) '_single_PCA'],1,p+1);
end
time_save = toc(t_save);
fprintf('\nelapsed time = %f s',time_save);
fprintf('\n');

%% Double PCA
%% First reduction step in space
if performPCAspace
    fprintf('\nPerforming PCA in space');
    t_PCA_space = tic;
    r = n*m;
    Rinit = min(r,N);
    if g<2^7
        [Y,Ya,Yb] = scaling(Y);
        mY = mean(Y,1);
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            Yc = bsxfun(@minus,Y,mY);
        else
            Yc = Y - mY;
        end
        % Yc = Y - repmat(mY,[N,1,1,1]); % Yc = Y - mY.*ones(N,1,1,1);
%         clear Y
    else
        Ya = zeros(1,n,m,p+1);
        Yb = zeros(1,n,m,p+1);
        mY = zeros(1,n,m,p+1);
    end
    if g<2^7
        V = zeros(r,Rinit,p+1);
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
            [Yt,Yat,Ybt] = scaling(Yt);
            Ya(1,:,:,t+1) = Yat;
            Yb(1,:,:,t+1) = Ybt;
            clear Yat Ybt
            mYt = mean(Yt,1);
            mY(1,:,:,t+1) = mYt;
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                Yct = bsxfun(@minus,Yt,mYt);
            else
                Yct = Yt - mYt;
            end
            % Yct = Yt - repmat(mYt,[N,1,1]); % Yct = Yt - mYt.*ones(N,1,1);
            clear Yt mYt
        end
        Yct = Yct(:,:)';
        [Vt,Sigt,Zt,errsvdYct] = svdtruncate(Yct,tolsvdYc);
        Sigt = Sigt/sqrt(N-1);
        sigt = diag(Sigt);
        Zt = Zt*sqrt(N-1);
        % if verLessThan('matlab','9.1') % compatibility (<R2016b)
        %     Yct_approx = Vt*Sigt*Zt';
        % else
        %     Yct_approx = Vt*(sigt.*Zt');
        % end
        
        % [coefft,scoret,latentt] = pca(Yct');
        % Vt = coefft;
        % sigt = sqrt(latentt);
        % Sigt = diag(sigt);
        % if verLessThan('matlab','9.1') % compatibility (<R2016b)
        %     Zt = scoret*diag(1./sigt);
        % else
        %     Zt = scoret./sigt';
        % end
        % % Yct_approx = coefft*scoret';
        
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
        
        if g<2^7
            V(:,1:Rt,t+1) = Vt;
        else
            save(fullfile(gridpathname,['data_PCA_space_t' num2str(t) '.mat']),'Vt');
        end
        sig(1:Rt,t+1) = sigt;
        Z(:,1:Rt,t+1) = Zt;
        R(t+1) = Rt;
        errsvdYc(1:Rt,t+1) = errsvdYct;
        err2Yc(t+1) = err2Yct;
        norm2Yc(t+1) = norm2Yct;
        clear Yct Vt sigt Zt Rt errsvdYct err2Yct norm2Yct
    end
    
    ts = (0:p)*dt;
    errL2 = trapz(ts,err2Yc,2)/trapz(ts,norm2Yc,2);
    fprintf('\nL2-error = %.3e for Y',errL2);
    fprintf('\n');
    
    Rmax = max(R);
    if g<2^7
        V = V(:,1:Rmax,:);
    end
    sig = sig(1:Rmax,:);
    Z = Z(:,1:Rmax,:);
    errsvdYc = errsvdYc(1:Rmax,:);
    
    time_PCA_space = toc(t_PCA_space);
    fprintf('\nelapsed time = %f s',time_PCA_space);
    fprintf('\n');
    
    save(fullfile(gridpathname,'scaling.mat'),'Ya','Yb');
    save(fullfile(gridpathname,'data_PCA_space.mat'),'mY','sig','Z','R','errsvdYc','time_PCA_space');
    if g<2^7
        save(fullfile(gridpathname,'data_PCA_space.mat'),'V','-append');
    end
else
    fprintf('\nLoading DNS data from PCA in space');
    t_load = tic;
    load(fullfile(gridpathname,'scaling.mat'),'Ya','Yb');
    load(fullfile(gridpathname,'data_PCA_space.mat'),'mY','sig','Z','R','errsvdYc','time_PCA_space');
    Rmax = max(R);
    if g<2^7
        load(fullfile(gridpathname,'data_PCA_space.mat'),'V');
    end
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
end

%% Second reduction step in time
if performPCAtime
    fprintf('\nPerforming PCA in time');
    t_PCA_time = tic;
    Rmax = max(R);
    q = (p+1)*Rmax;
    Zr = Z(:,:)';
    mZ = mean(Zr,2);
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        Zc = bsxfun(@minus,Zr,mZ);
    else
        Zc = Zr - mZ;
    end
    % Zc = Zr - repmat(mZ,[1,N]); % Zc = Zr - mZ.*ones(1,N);
    [W,S,X,errsvdZc] = svdtruncate(Zc,tolsvdZc);
    S = S/sqrt(N-1);
    s = diag(S);
    X = X*sqrt(N-1);
    % if verLessThan('matlab','9.1') % compatibility (<R2016b)
    %     Zc_approx = W*S*X';
    % else
    %     Zc_approx = W*(s.*X');
    % end
    
    % [coeff,score,latent] = pca(Zc');
    % W = coeff;
    % s = sqrt(latent);
    % S = diag(s);
    % if verLessThan('matlab','9.1') % compatibility (<R2016b)
    %     X = score*diag(1./s);
    %     % norm(score'-S*X')
    % else
    %     X = score./s';
    %     % norm(score'-s.*X')
    % end
    % % Zc_approx = coeff*score';
    
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
        xlabel('$K''=(k''-1)R+\alpha''$','Interpreter',interpreter)
        ylabel('$K=(k-1)R+\alpha$','Interpreter',interpreter)
        title('Covariance matrix $[C_{\zeta}(t^k,t^{k''})]_{\alpha,\alpha''} = [C_Z]_{K,K''}$','Interpreter',interpreter)
        mysaveas(gridpathname,'covariance_CZ',formats,renderer);
        mymatlab2tikz(gridpathname,'covariance_CZ.tex');
        
%         for t=0:p
%             ind = t*Rmax+(1:Rmax);
%             CZt_approx = CZ_approx(ind,ind);
%             CZt = CZ(ind,ind);
%             
%             figure('Name','Covariance matrix')
%             clf
%             imagesc(CZt)
%             colorbar
%             axis image
%             set(gca,'FontSize',fontsize)
%             xlabel('$\alpha''$','Interpreter',interpreter)
%             ylabel('$\alpha$','Interpreter',interpreter)
%             title(['Covariance matrix $[C_{\zeta}(t^k,t^k)]_{\alpha,\alpha''} = [C_{Z_k}]_{\alpha,\alpha''}$ for $k=$' num2str(t)],'Interpreter',interpreter)
%             mysaveas(gridpathname,['covariance_CZ_t' num2str(t*100)],formats,renderer);
%             mymatlab2tikz(gridpathname,['covariance_CZ_t' num2str(t*100) '.tex']);
%         end
    end
    
    time_PCA_time = toc(t_PCA_time);
    fprintf('\nelapsed time = %f s',time_PCA_time);
    fprintf('\n');
    
    save(fullfile(gridpathname,'data_PCA_time.mat'),'mZ','W','s','X','Q','errsvdZc','time_PCA_time');
else
    fprintf('\nLoading DNS data from PCA time');
    t_load = tic;
    load(fullfile(gridpathname,'data_PCA_time.mat'),'mZ','W','s','X','Q','errsvdZc','time_PCA_time');
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
end

%% Display eigenvalues
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
    xlabel('$\alpha$','Interpreter',interpreter)
    ylabel('$\lambda_{\alpha}(t)$','Interpreter',interpreter)
    gridLegend(hdl,nCols,leg,'Location','BestOutside','Interpreter',interpreter);
    % legend(leg{:},'Location','NorthEastOutside')
    mysaveas(gridpathname,'eigenvalues_CY_order',formats,renderer);
    mymatlab2tikz(gridpathname,'eigenvalues_CY_order_double_PCA.tex');
    
    figure('Name','Evolution of eigenvalues w.r.t time for each order')
    clf
    nCols = 5;
    Rmax = max(R);
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
    xlabel('$t$ [s]','Interpreter',interpreter)
    ylabel('$\lambda_{\alpha}(t)$','Interpreter',interpreter)
    gridLegend(hdl,nCols,leg,'Location','BestOutside','Interpreter',interpreter);
    % legend(leg{:},'Location','NorthEastOutside')
    mysaveas(gridpathname,'eigenvalues_CY_time',formats,renderer);
    mymatlab2tikz(gridpathname,'eigenvalues_CY_time_double_PCA.tex');
    
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
    xlabel('$R$','Interpreter',interpreter)
    ylabel('$\varepsilon_{Y}(R;t)$','Interpreter',interpreter)
    gridLegend(hdl,nCols,leg,'Location','BestOutside','Interpreter',interpreter);
    % legend(leg{:},'Location','NorthEastOutside')
    mysaveas(gridpathname,'error_svdYc',formats,renderer);
    mymatlab2tikz(gridpathname,'error_svdYc_double_PCA.tex');
    
    figure('Name','Evolution of eigenvalues')
    clf
    semilogy(1:Q,s(:).^2,'LineStyle','-','Color','b','LineWidth',1);
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$\beta$','Interpreter',interpreter)
    ylabel('$\Lambda_{\beta}$','Interpreter',interpreter)
    mysaveas(gridpathname,'eigenvalues_CZ',formats,renderer);
    mymatlab2tikz(gridpathname,'eigenvalues_CZ_double_PCA.tex');
    
    figure('Name','Evolution of errors')
    clf
    semilogy(1:Q,errsvdZc(:).^2,'LineStyle','-','Color','b','LineWidth',1);
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$Q$','Interpreter',interpreter)
    ylabel('$\varepsilon_{Z}(Q)$','Interpreter',interpreter)
    mysaveas(gridpathname,'error_svdZc',formats,renderer);
    mymatlab2tikz(gridpathname,'error_svdZc_double_PCA.tex');
end

%% Samples of double PCA representation
fprintf('\nSaving samples of double PCA representation');
t_save = tic;
for t=0:p
    t_savet = tic;
    if g<2^7
        Yt = get_double_PCA_at_step(mY,V,sig,mZ,W,s,X,R,t+1,'gridpathname',gridpathname);
    else
        Yt = get_double_PCA_at_step(mY,[],sig,mZ,W,s,X,R,t+1,'gridpathname',gridpathname);
    end
    Yat = Ya(1,:,:,t+1);
    Ybt = Yb(1,:,:,t+1);
    Yt = unscaling(Yt,Yat,Ybt);
    clear Yat Ybt
    Ut = reshape(Yt(:,1:dim,:),[N,dim*m]);
    Ct = reshape(Yt(:,dim+1,:),[N,m]);
    parfor l=1:N
        Ult = Ut(l,:);
        Clt = Ct(l,:);
        fields = {Ult,Clt};
        fieldnames = {'velocity','phase'};
        write_vtk_mesh(M,fields,[],fieldnames,[],gridpathname,['diphasic_fluids_grid' num2str(g) '_sample_' num2str(l) '_double_PCA'],1,t);
    end
    time_savet = toc(t_savet);
    fprintf('\nTime step %2d/%2d : elapsed time = %f s',t,p,time_savet);
end
for l=1:N
    make_pvd_file(gridpathname,['diphasic_fluids_grid' num2str(g) '_sample_' num2str(l) '_double_PCA'],1,p+1);
end
time_save = toc(t_save);
fprintf('\nelapsed time = %f s',time_save);
fprintf('\n');

time_Total = toc(t_Total);
fprintf('\nElapsed time = %f s\n',time_Total);
