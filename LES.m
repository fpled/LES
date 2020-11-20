clc
clearvars
close all
rng('default');
myparallel('start');

usePCA = 'single'; % 'no', 'single', 'double'
PostProcessingTau = true;
PostProcessingPressure = true;
PostProcessingEnergy = true;
Filtering = false;
QoI = false;
saveMean = true;
saveStd = true;
saveSamples = true;

performPCA = true;
performPCAspace = true;
performPCAtime = true;
computeMean = true;
postProcess = true;
computeQoI = true;
applyFilter = true;
computeError = true;
constructMesh = true;

displayEigenvalues = true;
displayCovariance = false;
displayQoI = false;
displayError = false;

useGPU = false;

cmap = 'default';
% cmap = flipud(gray);
framerate = 5;
fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

pathname = fileparts(mfilename('fullpath'));
% pathname = '/mnt/usbdisk2/pled/LES';

sigma = 0.45; % surface tension (N/m)
mu = [0.1,0.1]; % dynamic viscosity of [water,oil] (Pa.s)
rho = [1000,900]; % mass density of [water,oil] (kg/m3)
gravity = 9.81; % gravity (m/s2)

perm = @(u) permute(u,[2,3,4,5,1]);
iperm = @(u) ipermute(u,[2,3,4,5,1]);

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

fprintf('\nn = %d variables',n);
fprintf('\nN = %d samples',N);
fprintf('\nm = %d spatial points',m);
fprintf('\np+1 = %d time steps',p+1);
fprintf('\n');

s = PrincipalComponentAnalysis('tol',eps,'maxRank',Inf,'checkOrthonormality',false);
sSpace = PrincipalComponentAnalysis('tol',eps,'maxRank',Inf,'checkOrthonormality',false);
sTime = PrincipalComponentAnalysis('tol',eps,'maxRank',Inf,'checkOrthonormality',false);

filterType = 'box'; % 3D filter type ('box' or 'mean' or 'average', 'linear' or 'trapz')

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
if dim==1
%     Grad = Dx;
%     GradN = DxN;
%     Div = Dx;
%     DivN = DxN;
%     Laplacian = Dxx;
    LaplacianN = DxxN;
elseif dim==2
    I1D = speye(g+1);
%     Ex = sparse(1,1,1,dim,1);
%     Ey = sparse(2,1,1,dim,1);
%     Gradx = kron(I1D,Dx);
%     Grady = kron(Dx,I1D);
%     Grad = kron(Gradx,Ex) + kron(Grady,Ey); % Gradx = Grad(1:dim:end,:); Grady = Grad(2:dim:end,:);
    GradxN = kron(I1D,DxN);
    GradyN = kron(DxN,I1D);
%     GradN = kron(GradxN,Ex) + kron(GradyN,Ey); % GradxN = GradN(1:dim:end,:); GradyN = GradN(2:dim:end,:);
%     Div = kron(Gradx,Ex') + kron(Grady,Ey');
%     DivN = kron(GradxN,Ex') + kron(GradyN,Ey');
%     clear Ex Ey
%     Laplacian = kron(I1D,Dxx) + kron(Dxx,I1D);
    LaplacianN = kron(I1D,DxxN) + kron(DxxN,I1D);
    clear I1D
elseif dim==3
    I1D = speye(g+1);
    I2D = speye((g+1)^2); % I2 = kron(I1D,I1D);
%     Ex = sparse(1,1,1,dim,1);
%     Ey = sparse(2,1,1,dim,1);
%     Ez = sparse(3,1,1,dim,1);
%     Gradx = kron(I1D,kron(I1D,Dx));
%     Grady = kron(I1D,kron(Dx,I1D));
%     Gradz = kron(kron(Dx,I1D),I1D);
%     Grad = kron(Gradx,Ex) + kron(Grady,Ey) + kron(Gradz,Ez); % Gradx = Grad(1:dim:end,:); Grady = Grad(2:dim:end,:); Gradz = Grad(3:dim:end,:);
    GradxN = kron(I1D,kron(I1D,DxN));
    GradyN = kron(I1D,kron(DxN,I1D));
    GradzN = kron(kron(DxN,I1D),I1D);
%     GradN = kron(GradxN,Ex) + kron(GradyN,Ey) + kron(GradzN,Ez); % GradxN = GradN(1:dim:end,:); GradyN = GradN(2:dim:end,:); GradzN = GradN(3:dim:end,:);
%     Div = kron(Gradx,Ex') + kron(Grady,Ey') + kron(Gradz,Ez');
%     DivN = kron(GradxN,Ex') + kron(GradyN,Ey') + kron(GradzN,Ez');
%     clear Ex Ey Ez
%     Laplacian1D = Dxx;
%     Laplacian2D = kron(I1D,Dxx) + kron(Dxx,I1D);
%     Laplacian = kron(I2D,Laplacian1D) + kron(Laplacian2D,I1D);
%     clear Laplacian1 Laplacian2
    Laplacian1DN = DxxN;
    Laplacian2DN = kron(I1D,DxxN) + kron(DxxN,I1D);
    LaplacianN = kron(I2D,Laplacian1DN) + kron(Laplacian2DN,I1D);
    clear I1D I2D Laplacian1DN Laplacian2DN
end

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

%% Statistical reduction of DNS data
switch usePCA
    case 'no'
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
    case 'single'
        if performPCA
            fprintf('\nPerforming PCA');
            if g>=2^7
                Y = zeros(N,n,m,p+1);
                for t=0:p
                    load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
                    Y(:,:,:,t+1) = Yt;
                    clear Yt
                end
            end
            t_PCA = tic;
            Y = Y(:,:)';
            [Y,Ya,Yb] = s.scaling(Y);
            [V,sv,X,errsvdY,mY] = s.principalComponents(Y);
            errY = errsvdY(end);
            R = length(sv);
            time_PCA = toc(t_PCA);
            fprintf('\nrank R = %d, error = %.3e for Y',R,errY);
            fprintf('\nelapsed time = %f s',time_PCA);
            fprintf('\n');
            
            %CY = cov(Y');
            %CY_approx = s.cov(V,sv);
            %errCY = norm(CY_approx-CY)/norm(CY);
            %fprintf('\n                          error = %.3e for CY',errCY);
            
            Ya = reshape(Ya,[n*m,p+1]);
            Yb = reshape(Yb,[n*m,p+1]);
            Y = reshape(Y,[n*m,p+1,N]);
            norm2Y = sum(var(Y,0,3));
            % Y_approx = s.reconstruction(mY,V,sv,X);
            % Y_approx = reshape(Y_approx,[n*m,p+1,N]);
            % norm2Y_approx = sum(var(Y_approx,0,3));
            mY = reshape(mY,[n*m,1,p+1]);
            V = permute(reshape(V',[R,n*m,p+1]),[2 1 3]);
            norm2Y_approx = arrayfun(@(t) sum(sSpace.var(V(:,:,t+1),sv)),0:p);
            err2Y = abs(norm2Y-norm2Y_approx);
            
            ts = (0:p)*dt;
            errL2 = trapz(ts,err2Y,2)/trapz(ts,norm2Y,2);
            fprintf('\nL2-error = %.3e for Y',errL2);
            fprintf('\n');
            
            save(fullfile(gridpathname,'scaling_PCA.mat'),'Ya','Yb');
            save(fullfile(gridpathname,'data_PCA.mat'),'s','mY','V','sv','X','R','errsvdY','errL2','time_PCA');
        else
            fprintf('\nLoading DNS data from PCA');
            t_load = tic;
            load(fullfile(gridpathname,'scaling_PCA.mat'),'Ya','Yb');
            load(fullfile(gridpathname,'data_PCA.mat'),'s','mY','V','sv','X','R','errsvdY','errL2','time_PCA');
            time_load = toc(t_load);
            fprintf('\nelapsed time = %f s',time_load);
            fprintf('\n');
        end
        
    case 'double'
        %% First reduction step in space
        if performPCAspace
            fprintf('\nPerforming PCA in space');
            t_PCA_space = tic;
            r = n*m;
            Rinit = min(r,N);
            Ya = zeros(r,p+1);
            Yb = zeros(r,p+1);
            mY = zeros(r,1,p+1);
            if g<2^7
                V = zeros(r,Rinit,p+1);
            end
            sv = zeros(Rinit,p+1);
            Z = zeros(N,Rinit,p+1);
            errsvdY = zeros(Rinit,p+1);
            err2Y = zeros(1,p+1);
            norm2Y = zeros(1,p+1);
            R = zeros(1,p+1);
            for t=0:p
                if g<2^7
                    Yt = Y(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
                end
                Yt = Yt(:,:)';
                t_PCA_spacet = tic;
                [Yt,Yat,Ybt] = sSpace.scaling(Yt);
                [Vt,svt,Zt,errsvdYt,mYt] = sSpace.principalComponents(Yt);
                errYt = errsvdYt(end);
                Rt = length(svt);
                time_PCA_spacet = toc(t_PCA_spacet);
                fprintf('\nTime step %2d/%2d : rank R = %d, error = %.3e for Y, elapsed time = %f s',t,p,Rt,errYt,time_PCA_spacet);
                
                %CYt = cov(Yt');
                %CYt_approx = sSpace.cov(Vt,svt);
                %errCYt = norm(CYt_approx-CYt)/norm(CYt);
                %fprintf('\n                          error = %.3e for CY',errCYt);
                
                norm2Yt = sum(var(Yt,0,2));
                err2Yt = errYt^2*norm2Yt;
                % Yt_approx = sSpace.reconstruction(mYt,Vt,svt,Zt);
                % norm2Yt_approx = sum(var(Yt_approx,0,2));
                % norm2Yt_approx = sum(sSpace.var(Vt,svt));
                % err2Yt = norm2Yt-norm2Yt_approx;
                % svft = sSpace.singularValues(Yt);
                % err2Yt = sum(svft.^2)-sum(svt.^2);
                
                mY(:,1,t+1) = mYt;
                Ya(:,t+1) = Yat;
                Yb(:,t+1) = Ybt;
                if g<2^7
                    V(:,1:Rt,t+1) = Vt;
                else
                    save(fullfile(gridpathname,['data_PCA_space_t' num2str(t) '.mat']),'Vt');
                end
                sv(1:Rt,t+1) = svt;
                Z(:,1:Rt,t+1) = Zt;
                R(t+1) = Rt;
                errsvdY(1:Rt,t+1) = errsvdYt;
                err2Y(t+1) = err2Yt;
                norm2Y(t+1) = norm2Yt;
                clear Yt mYt Yat Ybt Vt svt Zt Rt errsvdYt err2Yt norm2Yt
            end
            
            ts = (0:p)*dt;
            errL2Y = trapz(ts,err2Y,2)/trapz(ts,norm2Y,2);
            fprintf('\nL2-error = %.3e for Y',errL2Y);
            fprintf('\n');
            
            Rmax = max(R);
            if g<2^7
                V = V(:,1:Rmax,:);
            end
            sv = sv(1:Rmax,:);
            Z = Z(:,1:Rmax,:);
            errsvdY = errsvdY(1:Rmax,:);
            
            time_PCA_space = toc(t_PCA_space);
            fprintf('\nelapsed time = %f s',time_PCA_space);
            fprintf('\n');
            
            save(fullfile(gridpathname,'scaling_PCA_space.mat'),'Ya','Yb');
            save(fullfile(gridpathname,'data_PCA_space.mat'),'sSpace','mY','sv','Z','R','errsvdY','errL2Y','time_PCA_space');
            if g<2^7
                save(fullfile(gridpathname,'data_PCA_space.mat'),'V','-append');
            end
        else
            fprintf('\nLoading DNS data from PCA in space');
            t_load = tic;
            load(fullfile(gridpathname,'scaling_PCA_space.mat'),'Ya','Yb');
            load(fullfile(gridpathname,'data_PCA_space.mat'),'sSpace','mY','sv','Z','R','errsvdY','errL2Y','time_PCA_space');
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
            Rmax = max(R);
            q = (p+1)*Rmax;
            t_PCA_time = tic;
            Z = Z(:,:)';
            [W,sw,X,errsvdZ,mZ] = sTime.principalComponents(Z);
            errZ = errsvdZ(end);
            Q = length(sw);
            time_PCA_time = toc(t_PCA_time);
            fprintf('\nrank R = %d, rank Q = %d, error = %.3e for Z',Rmax,Q,errZ);
            
            CZ = cov(Z');
            CZ_approx = sTime.cov(W,sw);
            errCZ = norm(CZ_approx-CZ)/norm(CZ);
            fprintf('\n                          error = %.3e for CZ',errCZ);
            
            Z = reshape(Z,[Rmax,p+1,N]);
            norm2Z = sum(var(Z,0,3));
            % Z_approx = sTime.reconstruction(mZ,W,sw,X);
            % Z_approx = reshape(Z_approx,[Rmax,p+1,N]);
            % norm2Z_approx = sum(var(Z_approx,0,3));
            W = reshape(W',[Q,Rmax,p+1]);
            norm2Z_approx = arrayfun(@(t) sum(sTime.var(W(:,1:R(t+1),t+1)',sw)),0:p);
            err2Z = abs(norm2Z-norm2Z_approx);
            
            ts = (0:p)*dt;
            errL2Z = trapz(ts,err2Z,2)/trapz(ts,norm2Z,2);
            fprintf('\nL2-error = %.3e for Z',errL2Z);
            clear Z
            
            mZ = reshape(mZ,[1,Rmax,p+1]);
            err2Y = zeros(1,p+1);
            for t=0:p
                Rt = R(t+1);
                if g<2^7
                    Vt = V(:,1:Rt,t+1);
                else
                    save(fullfile(gridpathname,['data_PCA_space_t' num2str(t) '.mat']),'Vt');
                end
                svt = sv(1:Rt,t+1);
                Wt = W(:,1:Rt,t+1);
                norm2Yt_approx = sum(sSpace.varDouble(Vt,svt,Wt,sw));
                % mYt = mY(:,:,t+1);
                % mZt = mZ(1,1:Rt,t+1);
                % Zt = sTime.reconstruction(mZt',Wt',sw,X);
                % Yt_approx = sSpace.reconstruction(mYt,Vt,svt,Zt');
                % norm2Yt_approx = sum(var(Yt_approx,0,2));
                norm2Yt = norm2Y(:,t+1);
                err2Yt = abs(norm2Yt-norm2Yt_approx);
                err2Y(:,t+1) = err2Yt;
            end
            
            ts = (0:p)*dt;
            errL2Y = trapz(ts,err2Y,2)/trapz(ts,norm2Y,2);
            fprintf('\nL2-error = %.3e for Y',errL2Y);
            fprintf('\n');
            
            fprintf('\nelapsed time = %f s',time_PCA_time);
            fprintf('\n');
            
            save(fullfile(gridpathname,'data_PCA_time.mat'),'sTime','mZ','W','sw','X','Q','errsvdZ','errL2Z','errL2Y','time_PCA_time');
        else
            fprintf('\nLoading DNS data from PCA time');
            t_load = tic;
            load(fullfile(gridpathname,'data_PCA_time.mat'),'sTime','mZ','W','sw','X','Q','errsvdZ','errL2Z','time_PCA_time');
            time_load = toc(t_load);
            fprintf('\nelapsed time = %f s',time_load);
            fprintf('\n');
        end
        
    otherwise
        error('Method not implemented.')
end

if postProcess && (PostProcessingTau || PostProcessingPressure || PostProcessingEnergy)
%% Post-processing DNS data
fprintf('\nPost-processing DNS data');
t_PostProcess = tic;
if PostProcessingTau
    ntau = 4*dim+1; % number of tau variables
    fprintf('\nn = %d tau variables',ntau);
    mTau = zeros(1,ntau,m,p+1);
    if g<2^7
        Tau = zeros(N,ntau,m,p+1);
    end
elseif PostProcessingPressure
    load(fullfile(gridpathname,'mean_data_tau.mat'),'mTau','ntau');
    if g<2^7
        load(fullfile(gridpathname,'data_tau.mat'),'Tau');
    end
end
if PostProcessingPressure
    fprintf('\nn = %d pressure variable',1);
    mPres = zeros(1,1,m,p+1);
    if g<2^7
        Pres = zeros(N,1,m,p+1);
    end
end
if PostProcessingEnergy
    if ~PostProcessingPressure
        load(fullfile(gridpathname,'mean_data_pressure.mat'),'mPres');
        if g<2^7
            load(fullfile(gridpathname,'data_pressure.mat'),'Pres');
        end
    end
    ne = 9; % number of energy variables
    fprintf('\nn = %d energy variables',ne);
    mE = zeros(1,ne,m,p+1);
    if g<2^7
        E = zeros(N,ne,m,p+1);
    end
end

for t=0:p
    t_PostProcesst = tic;
    
    switch usePCA
        case 'no'
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
            
        case 'single'
            if t>0
                Yt_old = s.reconstructionAtStep(mY,V,sv,X,t);
                Yt_old = s.unscaling(Yt_old,Ya(:,t),Yb(:,t))';
            end
            if t<p
                Yt_new = s.reconstructionAtStep(mY,V,sv,X,t+2);
                Yt_new = s.unscaling(Yt_new,Ya(:,t+2),Yb(:,t+2))';
            end
            Yt = s.reconstructionAtStep(mY,V,sv,X,t+1);
            Yt = s.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
            
        case 'double'
            if t>0
                if g<2^7
                    Yt_old = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t-1) '.mat']);
                else
                    Yt_old = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t-1) '.mat']);
                end
                Yt_old = sSpace.unscaling(Yt_old,Ya(:,t),Yb(:,t))';
            end
            if t<p
                if g<2^7
                    Yt_new = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+2,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t+1) '.mat']);
                else
                    Yt_new = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+2,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t+1) '.mat']);
                end
                Yt_new = sSpace.unscaling(Yt_new,Ya(:,t+2),Yb(:,t+2))';
            end
            if g<2^7
                Yt = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t) '.mat']);
            else
                Yt = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t) '.mat']);
            end
            Yt = sSpace.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
    end
    
    Yt = perm(reshape(Yt,[N,n,sx]));
    ut = Yt(1:dim,:,:,:,:);
    Ct = Yt(dim+1,:,:,:,:);
    clear Yt
    rhot = Ct*rho(2) + (1-Ct)*rho(1);
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        rhout = bsxfun(@times,rhot,ut);
    else
        rhout = rhot.*ut;
    end
    % rhout = repmat(rhot,[dim,ones(1,dim+1)]).*ut;
    if PostProcessingEnergy
        u2t = dot(ut,ut,1);
        if t==0 || t==p
            Ek = 1/2*rhot.*u2t;
        end
    end
    if t>0
        Yt_old = perm(reshape(Yt_old,[N,n,sx]));
        ut_old = Yt_old(1:dim,:,:,:,:);
        Ct_old = Yt_old(dim+1,:,:,:,:);
        clear Yt_old
        rhot_old = Ct_old*rho(2) + (1-Ct_old)*rho(1);
        if PostProcessingTau
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                rhout_old = bsxfun(@times,rhot_old,ut_old);
            else
                rhout_old = rhot_old.*ut_old;
            end
            % rhout_old = repmat(rhot_old,[dim,ones(1,dim+1)]).*ut_old;
        end
        if PostProcessingEnergy
            Ek_old = 1/2*rhot_old.*dot(ut_old,ut_old,1);
        end
        clear ut_old Ct_old rhot_old
    end
    if t<p
        Yt_new = perm(reshape(Yt_new,[N,n,sx]));
        ut_new = Yt_new(1:dim,:,:,:,:);
        Ct_new = Yt_new(dim+1,:,:,:,:);
        clear Yt_new
        rhot_new = Ct_new*rho(2) + (1-Ct_new)*rho(1);
        if PostProcessingTau
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                rhout_new = bsxfun(@times,rhot_new,ut_new);
            else
                rhout_new = rhot_new.*ut_new;
            end
            % rhout_new = repmat(rhot_new,[dim,ones(1,dim+1)]).*ut_new;
        end
        if PostProcessingEnergy
            Ek_new = 1/2*rhot_new.*dot(ut_new,ut_new,1);
        end
        clear ut_new Ct_new rhot_new
    end
    
    % Check boundary conditions: zero normal velocity at boundaries 
    % ut1 = ut(1,[1,end],:,:,:);
    % ut2 = ut(2,:,[1,end],:,:);
    % ut3 = ut(3,:,:,[1,end],:);
    % norm(ut1(:))
    % norm(ut2(:))
    % norm(ut3(:))
    % Check incompressibility condition
    % divut = div(ut,Dx);
    % mdivut = mean(divut(:,:,:,:,:),5);

    gradut = grad(ut,Dx);
    St = (gradut+permute(gradut,[2,1,3+(0:dim)]))/2;
    if ~PostProcessingEnergy
        clear gradut
    end
    mut = Ct*mu(2) + (1-Ct)*mu(1);
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        muSt = bsxfun(@times,shiftdim(mut,-1),St);
    else
        muSt = shiftdim(mut,-1).*St;
    end
    % muSt = repmat(shiftdim(mut,-1),[dim,dim,ones(1,dim+1)]).*St;
    clear mut St
    gradCt = grad(Ct,Dx);
    clear Ct
    ngradCt = normal(gradCt);
    kappa = div(ngradCt,Dx);
    clear ngradCt
    
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        tauSurft = sigma*bsxfun(@times,kappa,gradCt);
    else
        tauSurft = sigma*kappa.*gradCt;
    end
    % tauSurft = sigma*repmat(kappa,[dim,ones(1,dim+1)]).*gradCt;
    clear kappa
    if ~PostProcessingTau
        clear gradCt
    end
    if PostProcessingTau
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
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            divtauConvt = squeeze(sum(bsxfun(@times,grad(rhout,Dx),shiftdim(ut,-1)),2));
        else
            divtauConvt = squeeze(sum(grad(rhout,Dx).*shiftdim(ut,-1),2));
        end
        % divtauConvt = squeeze(sum(grad(rhout,Dx).*repmat(shiftdim(ut,-1),[dim,ones(1,dim+2)]),2));
        % other formulae
        % divtauConvt = div(permute(repmat(shiftdim(rhout,-1),[dim,ones(1,dim+2)]),[2,1,dim:(dim+3)]).*repmat(shiftdim(ut,-1),[dim,ones(1,dim+2)]),Dx);
        tauInterft = dot(ut,gradCt,1);
        clear gradCt
        divtauDifft = div(2*muSt,Dx);
        Taut = cat(1,tauTimet,divtauConvt,divtauDifft,tauSurft,tauInterft);
        clear tauTimet tauInterft
        if ~PostProcessingPressure
            clear divtauConvt divtauDifft
        end
        if ~PostProcessingEnergy
            clear ut rhout muSt
        end
        if ~PostProcessingPressure
            clear rhot tauSurft
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
    
    if PostProcessingPressure
        if ~PostProcessingTau
            if g<2^7
                Taut = Tau(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_tau_t' num2str(t) '.mat']),'Taut');
            end
            Taut = perm(reshape(Taut,[N,ntau,sx]));
            divtauConvt = Taut(dim+(1:dim),:,:,:,:);
            divtauDifft = Taut(2*dim+(1:dim),:,:,:,:);
            tauSurft = Taut(3*dim+(1:dim),:,:,:,:);
            clear Taut
        end
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            b = div(bsxfun(@rdivide,divtauConvt-divtauDifft-tauSurft,rhot),Dx);
        else
            b = div((divtauConvt-divtauDifft-tauSurft)./rhot,Dx);
        end
        % b = div((divtauConvt-divtauDifft-tauSurft)./repmat(rhot,[dim,ones(1,dim+1)]),Dx);
        clear divtauConvt divtauDifft
        if ~PostProcessingEnergy
            clear tauSurft
        end
        gradinvrhot = grad(1./rhot,Dx);
        % Implicit central-difference spatial scheme
        B = reshape(b,[m,N]);
        Gradinvrhot = reshape(gradinvrhot,[dim,m,N]);
        clear b gradinvrhot
        Rhot = reshape(rhot,[m,N]);
        prest = zeros(m,N);
        parfor l=1:N
            bl = B(:,l);
            Rhotl = Rhot(:,l);
            Gradinvrhotl = Gradinvrhot(:,:,l);
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                Al = bsxfun(@times,Gradinvrhotl(1,:),GradxN) + bsxfun(@times,Gradinvrhotl(2,:),GradyN) + bsxfun(@times,Gradinvrhotl(3,:),GradzN);
                Al = Al + bsxfun(@rdivide,LaplacianN,Rhotl);
            else
                Al = Gradinvrhotl(1,:).*GradxN + Gradinvrhotl(2,:).*GradyN + Gradinvrhotl(3,:).*GradzN;
                Al = Al + LaplacianN./Rhotl;
            end
            % Al = repmat(Gradinvrhotl(1,:),[m,1]).*GradxN + repmat(Gradinvrhotl(2,:),[m,1]).*GradyN + repmat(Gradinvrhotl(3,:),[m,1]).*GradzN;
            % Al = Al + LaplacianN./repmat(Rhotl,[1,m]);
            prestl = Al\bl;
            prest(:,l) = prestl;
        end
        Prest = reshape(prest,[1,sx,N]);
        clear B Rhot Gradinvrhot prest
        if ~PostProcessingEnergy
            clear rhot
        end
        
        Prest = iperm(Prest);
        Prest = Prest(:,:,:);
        mPrest = mean(Prest,1);
        mPres(1,:,:,t+1) = mPrest;
        clear mPrest
        if g<2^7
            Pres(:,:,:,t+1) = Prest;
        else
            save(fullfile(gridpathname,['data_pressure_t' num2str(t) '.mat']),'Prest');
        end
        if ~PostProcessingEnergy
            clear Prest
        end
    end
    
    if PostProcessingEnergy
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
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            energyConvt = div(bsxfun(@times,(rhot.*u2t),ut),Dx);
        else
            energyConvt = div((rhot.*u2t).*ut,Dx);
        end
        % energyConvt = div(repmat(rhot.*u2t,[dim,ones(1,dim+1)]).*ut,Dx);
        clear rhot
        energyKinSpacet = dot(rhout,grad(u2t/2,Dx),1);
        clear u2t
        energyGravt = gravity.*rhout(2,:,:,:,:);
        clear rhout
        if ~PostProcessingPressure
            if g<2^7
                Prest = Pres(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_pressure_t' num2str(t) '.mat']),'Prest');
            end
        end
        prest = perm(reshape(Prest,[N,1,sx]));
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            energyPrest = div(bsxfun(@times,prest,ut),Dx);
        else
            energyPrest = div(prest.*ut,Dx);
        end
        % energyPrest = div(repmat(prest,[dim,ones(1,dim+1)]).*ut,Dx);
        energyPresDilt = prest.*div(ut,Dx);
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            energyDifft = div(squeeze(sum(bsxfun(@times,2*muSt,shiftdim(ut,-1)),2)),Dx);
        else
            energyDifft = div(squeeze(sum(2*muSt.*shiftdim(ut,-1),2)),Dx);
        end
        % energyDifft = div(squeeze(dot(2*muSt,repmat(shiftdim(ut,-1),[dim,ones(1,dim+2)]),2)),Dx);
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

if PostProcessingTau
    save(fullfile(gridpathname,'mean_data_tau.mat'),'mTau','ntau');
    if g<2^7
        save(fullfile(gridpathname,'data_tau.mat'),'Tau');
    end
end
if PostProcessingPressure
    save(fullfile(gridpathname,'mean_data_pressure.mat'),'mPres');
    if g<2^7
        save(fullfile(gridpathname,'data_pressure.mat'),'Pres');
    end
end
if PostProcessingEnergy
    save(fullfile(gridpathname,'mean_data_energy.mat'),'mE','ne');
    if g<2^7
        save(fullfile(gridpathname,'data_energy.mat'),'E');
    end
end
else
    if PostProcessingTau
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
    if PostProcessingPressure
        if g<2^7
            fprintf('\nLoading pressure data');
        else
            fprintf('\nLoading mean pressure data');
        end
        t_load = tic;
        if PostProcessingPressure
            load(fullfile(gridpathname,'mean_data_pressure.mat'),'mPres');
            if g<2^7
                load(fullfile(gridpathname,'data_pressure.mat'),'Pres');
            end
        end
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
    end
    if PostProcessingEnergy
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

if QoI
%% Compute quantities of interest: spatial average in each phase
if computeQoI
    fprintf('\nComputing quantities of interest');
    t_QoI = tic;
    
    Qu = zeros(N,dim,2,p+1);
    if PostProcessingTau
        Qtau = zeros(N,ntau,2,p+1);
    end
    if PostProcessingEnergy
        Qe = zeros(N,ne,2,p+1);
    end
    
    for t=0:p
        t_QoIt = tic;
        
        switch usePCA
            case 'no'
                if g<2^7
                    Yt = Y(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
                end
                
            case 'single'
                Yt = s.reconstructionAtStep(mY,V,sv,X,t+1);
                Yt = s.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
                
            case 'double'
                if g<2^7
                    Yt = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t) '.mat']);
                else
                    Yt = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t) '.mat']);
                end
                Yt = sSpace.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
        end
        
        Yt = reshape(Yt,[N,n,sx]);
        ut = Yt(:,1:dim,:,:,:);
        Ct = Yt(:,dim+1,:,:,:);
        clear Yt
        Qut = int_trapz(x,Ct,ut,dim);
        Qu(:,:,:,t+1) = Qut;
        clear ut Qut
        
        if PostProcessingTau
            if g<2^7
                Taut = Tau(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_tau_t' num2str(t) '.mat']),'Taut');
            end
            Taut = reshape(Taut,[N,ntau,sx]);
            Qtaut = int_trapz(x,Ct,Taut,dim);
            Qtau(:,:,:,t+1) = Qtaut;
            clear Taut Qtaut
        end
        
        if PostProcessingPressure
            if g<2^7
                Prest = Pres(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_pressure_t' num2str(t) '.mat']),'Prest');
            end
            Prest = reshape(Prest,[N,1,sx]);
            Qprest = int_trapz(x,Ct,Prest,dim);
            Qpres(:,:,:,t+1) = Qprest;
            clear Prest Qprest
        end
        
        if PostProcessingEnergy
            if g<2^7
                Et = E(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_energy_t' num2str(t) '.mat']),'Et');
            end
            Et = reshape(Et,[N,ne,sx]);
            Qet = int_trapz(x,Ct,Et,dim);
            Qe(:,:,:,t+1) = Qet;
            clear Et Qet
        end
        clear Ct
        time_QoIt = toc(t_QoIt);
        fprintf('\nTime step %2d/%2d : elapsed time = %f s',t,p,time_QoIt);
    end
    
    [mQu,stdQu,RQu,IQu] = compute_stats(Qu,dt);
    if PostProcessingTau
        [mQtau,stdQtau,RQtau,IQtau] = compute_stats(Qtau,dt);
    end
    if PostProcessingPressure
        [mQpres,stdQpres,RQpres,IQpres] = compute_stats(Qpres,dt);
    end
    if PostProcessingEnergy
        [mQe,stdQe,RQe,IQe] = compute_stats(Qe,dt);
    end
    
    time_QoI = toc(t_QoI);
    fprintf('\nelapsed time = %f s',time_QoI);
    fprintf('\n');
    
    save(fullfile(gridpathname,'data_qoi.mat'),'Qu','mQu','IQu');
    if PostProcessingTau
        save(fullfile(gridpathname,'data_qoi_tau.mat'),'Qtau','mQtau','IQtau');
    end
    if PostProcessingPressure
        save(fullfile(gridpathname,'data_qoi_pressure.mat'),'Qpres','mQpres','IQpres');
    end
    if PostProcessingEnergy
        save(fullfile(gridpathname,'data_qoi_energy.mat'),'Qe','mQe','IQe');
    end
else
    fprintf('\nLoading quantities of interest');
    t_load = tic;
    load(fullfile(gridpathname,'data_qoi.mat'),'Qu','mQu','IQu');
    if PostProcessingTau
        load(fullfile(gridpathname,'data_qoi_tau.mat'),'Qtau','mQtau','IQtau');
    end
    if PostProcessingPressure
        load(fullfile(gridpathname,'data_qoi_pressure.mat'),'Qpres','mQpres','IQpres');
    end
    if PostProcessingEnergy
        load(fullfile(gridpathname,'data_qoi_energy.mat'),'Qe','mQe','IQe');
    end
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
end
end

if Filtering
%% Applying filter
if g==gref && applyFilter
    fprintf('\nConstructing filter');
    hset = cell(1,ng-1);
    for ig=1:ng-1
        switch filterType
            case {'box','mean','average'}
                filterSize = (2^ig+1)*ones(1,dim);
                h = ones(filterSize)/prod(filterSize);
            case {'linear','trapz'}
                filterSize = 3*ones(1,3);
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
    
    %mYrefbar = zeros(ng-1,n,m,p+1);
    for ig=1:ng-1
        gbar = gset(ng-ig);
        gridbarname = ['Grid' num2str(gbar)];
        gridbarpathname = fullfile(pathname,gridbarname);
        mbar = (gbar+1)^dim;
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
            switch usePCA
                case 'no'
                    if g<2^7
                        Yt = Y(:,:,:,t+1);
                    else
                        load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
                    end
                    
                case 'single'
                    Yt = s.reconstructionAtStep(mY,V,sv,X,t+1);
                    Yt = s.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
                    
                case 'double'
                    if g<2^7
                        Yt = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t) '.mat']);
                    else
                        Yt = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t) '.mat']);
                    end
                    Yt = sSpace.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
            end
            
            Ybart = reshape(Yt,[N,n,sx]);
            clear Yt
            if useGPU
                Ybart = gpuArray(Ybart);
            end
            Ybart = apply_filter(Ybart,filterType,h);
            % Ybart = imfilter(Ybart,H,'replicate');
            if useGPU
                Ybart = gather(Ybart);
            end
            
            time_Filtergt = toc(t_Filtergt);
            fprintf('\nTime %2d/%2d : elapsed time = %f s',t,p,time_Filtergt);
            
            %Yrefbart = Ybart(:,:,:);
            %mYrefbar(ig,:,:,t+1) = mean(Yrefbart,1);
            %clear Yrefbart
            
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
    
    %save(fullfile(gridpathname,'mean_data_filtered.mat'),'mYrefbar');
    %clear mYrefbar
    
    if PostProcessingTau
        fprintf('\nApplying filter for tau data');
        t_Filter = tic;
        
        %mTaurefbar = zeros(ng-1,ntau,m,p+1);
        for ig=1:ng-1
            gbar = gset(ng-ig);
            gridbarname = ['Grid' num2str(gbar)];
            gridbarpathname = fullfile(pathname,gridbarname);
            mbar = (gbar+1)^dim;
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
                if useGPU
                    Taubart = gpuArray(Taubart);
                end
                Taubart = apply_filter(Taubart,filterType,h);
                % Taubart = imfilter(Taubart,H,'replicate');
                if useGPU
                    Taubart = gather(Taubart);
                end
                
                time_Filtergt = toc(t_Filtergt);
                fprintf('\nTime %2d/%2d : elapsed time = %f s',t,p,time_Filtergt);
                
                %Taurefbart = Taubart(:,:,:);
                %mTaurefbar(ig,:,:,t+1) = mean(Taurefbart,1);
                %clear Taurefbart
                
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
        
        %save(fullfile(gridpathname,'mean_data_tau_filtered.mat'),'mTaurefbar');
        %clear mTaurefbar
    end
    
    if PostProcessingPressure
        fprintf('\nApplying filter for pressure data');
        t_Filter = tic;
        
        %mPresrefbar = zeros(ng-1,1,m,p+1);
        for ig=1:ng-1
            gbar = gset(ng-ig);
            gridbarname = ['Grid' num2str(gbar)];
            gridbarpathname = fullfile(pathname,gridbarname);
            mbar = (gbar+1)^dim;
            k = g/gbar;
            h = hset{ig};
            H = shiftdim(h,-2);
            
            fprintf('\nFiltering pressure implicit data on coarse grid %d',gbar);
            t_Filterg = tic;
            
            mPresbar = zeros(1,1,mbar,p+1);
            if g<2^7
                Presbar = zeros(N,1,mbar,p+1);
            end
            for t=0:p
                t_Filtergt = tic;
                if g<2^7
                    Prest = Pres(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,['data_pressure_t' num2str(t) '.mat']),'Prest');
                end
                
                Presbart = reshape(Prest,[N,1,sx]);
                clear Prest
                if useGPU
                    Presbart = gpuArray(Presbart);
                end
                Presbart = apply_filter(Presbart,filterType,h);
                % Presbart = imfilter(Presbart,H,'replicate');
                if useGPU
                    Presbart = gather(Presbart);
                end
                
                time_Filtergt = toc(t_Filtergt);
                fprintf('\nTime %2d/%2d : elapsed time = %f s',t,p,time_Filtergt);
                
                %Presrefbart = Presbart(:,:,:);
                %mPresrefbar(ig,:,:,t+1) = mean(Presrefbart,1);
                %clear Presrefbart
                
                Presbart = reshape(Presbart(:,:,1:k:end,1:k:end,1:k:end),[N,1,mbar]);
                mPresbar(1,:,:,t+1) = mean(Presbart,1);
                if gbar<2^7
                    Presbar(:,:,:,t+1) = Presbart;
                else
                    save(fullfile(gridbarpathname,['data_pressure_filtered_t' num2str(t) '.mat']),'Presbart');
                end
                clear Presbart
            end
            
            time_Filterg = toc(t_Filterg);
            fprintf('\nelapsed time = %f s for coarse grid %d',time_Filterg,gbar);
            
            save(fullfile(gridbarpathname,'mean_data_pressure_filtered.mat'),'mPresbar');
            clear mPresbar
            if gbar<2^7
                save(fullfile(gridbarpathname,'data_pressure_filtered.mat'),'Presbar');
                clear Presbar
            end
        end
        time_Filter = toc(t_Filter);
        fprintf('\nelapsed time = %f s for all coarse grids',time_Filter);
        fprintf('\n');
        
        %save(fullfile(gridpathname,'mean_data_pressure_filtered.mat'),'mPresrefbar');
        %clear mPresrefbar
    end
    
    if PostProcessingEnergy
        fprintf('\nApplying filter for energy data');
        t_Filter = tic;
        
        %mErefbar = zeros(ng-1,ne,m,p+1);
        for ig=1:ng-1
            gbar = gset(ng-ig);
            gridbarname = ['Grid' num2str(gbar)];
            gridbarpathname = fullfile(pathname,gridbarname);
            mbar = (gbar+1)^dim;
            k = g/gbar;
            h = hset{ig};
            H = shiftdim(h,-2);
            
            fprintf('\nFiltering tau data on coarse grid %d',gbar);
            t_Filterg = tic;
            
            mEbar = zeros(1,ne,mbar,p+1);
            if g<2^7
                Ebar = zeros(N,ne,mbar,p+1);
            end
            for t=0:p
                t_Filtergt = tic;
                if g<2^7
                    Et = E(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,['data_energy_t' num2str(t) '.mat']),'Et');
                end
                
                Ebart = reshape(Et,[N,ne,sx]);
                clear Et
                if useGPU
                    Ebart = gpuArray(Ebart);
                end
                Ebart = apply_filter(Ebart,filterType,h);
                % Ebart = imfilter(Taubart,H,'replicate');
                if useGPU
                    Ebart = gather(Ebart);
                end
                
                time_Filtergt = toc(t_Filtergt);
                fprintf('\nTime %2d/%2d : elapsed time = %f s',t,p,time_Filtergt);
                
                %Erefbart = Ebart(:,:,:);
                %mErefbar(ig,:,:,t+1) = mean(Erefbart,1);
                %clear Erefbart
                
                Ebart = reshape(Ebart(:,:,1:k:end,1:k:end,1:k:end),[N,ne,mbar]);
                mEbar(1,:,:,t+1) = mean(Ebart,1);
                if gbar<2^7
                    Ebar(:,:,:,t+1) = Ebart;
                else
                    save(fullfile(gridbarpathname,['data_energy_filtered_t' num2str(t) '.mat']),'Ebart');
                end
                clear Ebart
            end
            
            time_Filterg = toc(t_Filterg);
            fprintf('\nelapsed time = %f s for coarse grid %d',time_Filterg,gbar);
            
            save(fullfile(gridbarpathname,'mean_data_energy_filtered.mat'),'mEbar');
            clear mEbar
            if gbar<2^7
                save(fullfile(gridbarpathname,'data_energt_filtered.mat'),'Ebar');
                clear Ebar
            end
        end
        time_Filter = toc(t_Filter);
        fprintf('\nelapsed time = %f s for all coarse grids',time_Filter);
        fprintf('\n');
        
        %save(fullfile(gridpathname,'mean_data_energy_filtered.mat'),'mErefbar');
        %clear mErefbar
    end
elseif g~=gref 
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
    if PostProcessingTau
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
    
    %% Loading filtered pressure data
    if PostProcessingPressure
        if g<2^7
            fprintf('\nLoading filtered pressure data');
        else
            fprintf('\nLoading mean filtered pressure data');
        end
        t_load = tic;
        load(fullfile(gridpathname,'mean_data_pressure_filtered.mat'),'mPresbar');
        if g<2^7
            load(fullfile(gridpathname,'data_pressure_filtered.mat'),'Presbar');
        end
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
    end
    
    %% Loading filtered energy data
    if PostProcessingTau
        if g<2^7
            fprintf('\nLoading filtered energy data');
        else
            fprintf('\nLoading mean filtered energy data');
        end
        t_load = tic;
        load(fullfile(gridpathname,'mean_data_energy_filtered.mat'),'mEbar');
        if g<2^7
            load(fullfile(gridpathname,'data_energy_filtered.mat'),'Ebar');
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
        
        if PostProcessingTau
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
        
        if PostProcessingPressure
            errorpres = zeros(1,p+1);
        end
        
        if PostProcessingTau
            errorenergyKinTime = zeros(1,p+1);
            errorenergyConv = zeros(1,p+1);
            errorenergyGrav = zeros(1,p+1);
            errorenergyPres = zeros(1,p+1);
            errorenergyPresDil = zeros(1,p+1);
            errorenergyKinSpace = zeros(1,p+1);
            errorenergyDiff = zeros(1,p+1);
            errorenergyVisc = zeros(1,p+1);
            errorenergySurf = zeros(1,p+1);
            normenergyKinTime = zeros(1,p+1);
            normenergyConv = zeros(1,p+1);
            normenergyGrav = zeros(1,p+1);
            normenergyPres = zeros(1,p+1);
            normenergyPresDil = zeros(1,p+1);
            normenergyKinSpace = zeros(1,p+1);
            normenergyDiff = zeros(1,p+1);
            normenergyVisc = zeros(1,p+1);
            normenergySurf = zeros(1,p+1);
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
            ut = Yt(:,1:dim,:);
            Ct = Yt(:,dim+1,:);
            clear Yt
            ubart = Ybart(:,1:dim,:);
            Cbart = Ybart(:,dim+1,:);
            clear Ybart
            
            errorut = compute_norm(x,reshape(sum((ubart-ut).^2,2),[N,sx]),dim);
            errorCt = compute_norm(x,reshape((Cbart-Ct).^2,[N,sx]),dim);
            clear ut Ct
            
            normubart = compute_norm(x,reshape(sum(ubart.^2,2),[N,sx]),dim);
            normCbart = compute_norm(x,reshape(Cbart.^2,[N,sx]),dim);
            clear ubart Cbart
            
            erroru(t+1) = errorut;
            errorC(t+1) = errorCt;
            clear errorut errorCt
            
            normubar(t+1) = normubart;
            normCbar(t+1) = normCbart;
            clear normubart normCbart
            
            if PostProcessingTau
                if g<2^7
                    Taut = Tau(:,:,:,t+1);
                    Taubart = Taubar(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,['data_tau_t' num2str(t) '.mat']),'Taut');
                    load(fullfile(gridpathname,['data_tau_filtered_t' num2str(t) '.mat']),'Taubart');
                end
                tauTimet = Taut(:,1:dim,:);
                divtauConvt = Taut(:,dim+(1:dim),:);
                divtauDifft = Taut(:,2*dim+(1:dim),:);
                tauSurft = Taut(:,3*dim+(1:dim),:);
                tauInterft = Taut(:,4*dim+1,:);
                clear Taut
                tauTimebart = Taubart(:,1:dim,:);
                divtauConvbart = Taubart(:,dim+(1:dim),:);
                divtauDiffbart = Taubart(:,2*dim+(1:dim),:);
                tauSurfbart = Taubart(:,3*dim+(1:dim),:);
                tauInterfbart = Taubart(:,4*dim+1,:);
                clear Taubart
                
                errortauTimet = compute_norm(x,reshape(sum((tauTimebart-tauTimet).^2,2),[N,sx]),dim);
                errordivtauConvt = compute_norm(x,reshape(sum((divtauConvbart-divtauConvt).^2,2),[N,sx]),dim);
                errordivtauDifft = compute_norm(x,reshape(sum((divtauDiffbart-divtauDifft).^2,2),[N,sx]),dim);
                errortauSurft = compute_norm(x,reshape(sum((tauSurfbart-tauSurft).^2,2),[N,sx]),dim);
                errortauInterft = compute_norm(x,reshape((tauInterfbart-tauInterft).^2,[N,sx]),dim);
                clear tauTimet divtauConvt divtauDifft tauSurft tauInterft
                
                normtauTimebart = compute_norm(x,reshape(sum(tauTimebart.^2,2),[N,sx]),dim);
                normdivtauConvbart = compute_norm(x,reshape(sum(divtauConvbart.^2,2),[N,sx]),dim);
                normdivtauDiffbart = compute_norm(x,reshape(sum(divtauDiffbart.^2,2),[N,sx]),dim);
                normtauSurfbart = compute_norm(x,reshape(sum(tauSurfbart.^2,2),[N,sx]),dim);
                normtauInterfbart = compute_norm(x,reshape(tauInterfbart.^2,[N,sx]),dim);
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
            
            if PostProcessingPressure
                if g<2^7
                    Prest = Pres(:,:,:,t+1);
                    Presbart = Presbar(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,['data_pressure_t' num2str(t) '.mat']),'Prest');
                    load(fullfile(gridpathname,['data_pressure_filtered_t' num2str(t) '.mat']),'Presbart');
                end
                prest = Prest(:,1,:);
                clear Prest
                presbart = Presbart(:,1,:);
                clear Presbart
                
                errorprest = compute_norm(x,reshape((presbart-prest).^2,[N,sx]),dim);
                clear prest
                
                normpresbart = compute_norm(x,reshape(presbart.^2,[N,sx]),dim);
                clear presbart
                
                errorpres(t+1) = errorprest;
                clear errorprest
                
                normpresbar(t+1) = normpresbart;
                clear normpresbart
            end
            
            if PostProcessingEnergy
                if g<2^7
                    Et = E(:,:,:,t+1);
                    Ebart = Ebar(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,['data_energy_t' num2str(t) '.mat']),'Et');
                    load(fullfile(gridpathname,['data_energy_filtered_t' num2str(t) '.mat']),'Ebart');
                end
                energyKinTimet = Et(:,1,:);
                energyConvt = Et(:,2,:);
                energyGravt = Et(:,3,:);
                energyPrest = Et(:,4,:);
                energyPresDilt = Et(:,5,:);
                energyKinSpacet = Et(:,6,:);
                energyDifft = Et(:,7,:);
                energyVisct = Et(:,8,:);
                energySurft = Et(:,9,:);
                clear Et
                energyKinTimebart = Ebart(:,1,:);
                energyConvbart = Ebart(:,2,:);
                energyGravbart = Ebart(:,3,:);
                energyPresbart = Ebart(:,4,:);
                energyPresDilbart = Ebart(:,5,:);
                energyKinSpacebart = Ebart(:,6,:);
                energyDiffbart = Ebart(:,7,:);
                energyViscbart = Ebart(:,8,:);
                energySurfbart = Ebart(:,9,:);
                clear Ebart
                
                errorenergyKinTimet = compute_norm(x,reshape((energyKinTimebart-energyKinTimet).^2,[N,sx]),dim);
                errorenergyConvt = compute_norm(x,reshape((energyConvbart-energyConvt).^2,[N,sx]),dim);
                errorenergyGravt = compute_norm(x,reshape((energyGravbart-energyGravt).^2,[N,sx]),dim);
                errorenergyPrest = compute_norm(x,reshape((energyPresbart-energyPrest).^2,[N,sx]),dim);
                errorenergyPresDilt = compute_norm(x,reshape((energyPresDilbart-energyPresDilt).^2,[N,sx]),dim);
                errorenergyKinSpacet = compute_norm(x,reshape((energyKinSpacebart-energyKinSpacet).^2,[N,sx]),dim);
                errorenergyDifft = compute_norm(x,reshape((energyDiffbart-energyDifft).^2,[N,sx]),dim);
                errorenergyVisct = compute_norm(x,reshape((energyViscbart-energyVisct).^2,[N,sx]),dim);
                errorenergySurft = compute_norm(x,reshape((energySurfbart-energySurft).^2,[N,sx]),dim);
                clear energyKinTimet energyConvt energyGravt energyPrest energyPresDilt energyKinSpacet energyDifft energyVisct energySurft
                
                normenergyKinTimebart = compute_norm(x,reshape(energyKinTimebart.^2,[N,sx]),dim);
                normenergyConvbart = compute_norm(x,reshape(energyConvbart.^2,[N,sx]),dim);
                normenergyGravbart = compute_norm(x,reshape(energyGravbart.^2,[N,sx]),dim);
                normenergyPresbart = compute_norm(x,reshape(energyPresbart.^2,[N,sx]),dim);
                normenergyPresDilbart = compute_norm(x,reshape(energyPresDilbart.^2,[N,sx]),dim);
                normenergyKinSpacebart = compute_norm(x,reshape(energyKinSpacebart.^2,[N,sx]),dim);
                normenergyDiffbart = compute_norm(x,reshape(energyDiffbart.^2,[N,sx]),dim);
                normenergyViscbart = compute_norm(x,reshape(energyViscbart.^2,[N,sx]),dim);
                normenergySurfbart = compute_norm(x,reshape(energySurfbart.^2,[N,sx]),dim);
                clear energyKinTimebart energyConvbart energyGravbart energyPresbart energyPresDilbart energyKinSpacebart energyDiffbart energyViscbart energySurfbart
                
                errorenergyKinTime(t+1) = errorenergyKinTimet;
                errorenergyConv(t+1) = errorenergyConvt;
                errorenergyGrav(t+1) = errorenergyGravt;
                errorenergyPres(t+1) = errorenergyPrest;
                errorenergyPresDil(t+1) = errorenergyPresDilt;
                errorenergyKinSpace(t+1) = errorenergyKinSpacet;
                errorenergyDiff(t+1) = errorenergyDifft;
                errorenergyViscbar(t+1) = errorenergyViscbart;
                errorenergySurf(t+1) = errorenergySurft;
                clear errorenergyKinTimet errorenergyConvt errorenergyGravt errorenergyPrest errorenergyPresDilt errorenergyKinSpacet errorenergyDifft errorenergyViscbart errorenergySurft
                
                normenergyKinTimebar(t+1) = normenergyKinTimebart;
                normenergyConvbar(t+1) = normenergyConvbart;
                normenergyGravbar(t+1) = normenergyGravbart;
                normenergyPresbar(t+1) = normenergyPresbart;
                normenergyPresDilbar(t+1) = normenergyPresDilbart;
                normenergyKinSpacebar(t+1) = normenergyKinSpacebart;
                normenergyDiffbar(t+1) = normenergyDiffbart;
                normenergyViscbar(t+1) = normenergyViscbart;
                normenergySurfbar(t+1) = normenergySurfbart;
                clear normenergyKinTimebart normenergyConvbart normenergyGravbart normenergyPresbart normenergyPresDilbart normenergyKinSpacebart normenergyDiffbart normenergyViscbart normenergySurfbart
            end
            
            time_Errort = toc(t_Errort);
            fprintf('\nTime %2d/%2d : elapsed time = %f s',t,p,time_Errort);
        end
        
        time_Error = toc(t_Error);
        fprintf('\nelapsed time = %f s',time_Error);
        fprintf('\n');
        
        save(fullfile(gridpathname,'error_filter.mat'),'erroru','errorC','normubar','normCbar');
        if PostProcessingTau
            save(fullfile(gridpathname,'error_filter_tau.mat'),'errortauTime','errordivtauConv','errordivtauDiff','errortauSurf','errortauInterf',...
                'normtauTimebar','normdivtauConvbar','normdivtauDiffbar','normtauSurfbar','normtauInterfbar');
        end
        if PostProcessingPressure
            save(fullfile(gridpathname,'error_filter_pressure.mat'),'errorpres','normpresbar');
        end
        if PostProcessingEnergy
            save(fullfile(gridpathname,'error_filter_energy.mat'),'errorenergyKinTime''errorenergyConv','errorenergyGrav','errorenergyPres','errorenergyPresDil','errorenergyKinSpace','errorenergyViscbar','errorenergySurf',...
                'normenergyKinTimebar','normenergyConvbar','normenergyGravbar','normenergyPresbar','normenergyPresDilbar','normenergyKinSpacebar','normenergyDiffbar','normenergyViscbar','normenergySurfbar');
        end
    else
        fprintf('\nLoading error between DNS and filtered data');
        t_load = tic;
        load(fullfile(gridpathname,'error_filter.mat'),'erroru','errorC','normubar','normCbar');
        if PostProcessingTau
            load(fullfile(gridpathname,'error_filter_tau.mat'),'errortauTime','errordivtauConv','errordivtauDiff','errortauSurf','errortauInterf',...
            'normtauTimebar','normdivtauConvbar','normdivtauDiffbar','normtauSurfbar','normtauInterfbar');
        end
        if PostProcessingPressure
            load(fullfile(gridpathname,'error_filter_pressure.mat'),'errorpres','normpresbar');
        end
        if PostProcessingEnergy
            load(fullfile(gridpathname,'error_filter_energy.mat'),'errorenergyKinTime''errorenergyConv','errorenergyGrav','errorenergyPres','errorenergyPresDil','errorenergyKinSpace','errorenergyViscbar','errorenergySurf',...
                'normenergyKinTimebar','normenergyConvbar','normenergyGravbar','normenergyPresbar','normenergyPresDilbar','normenergyKinSpacebar','normenergyDiffbar','normenergyViscbar','normenergySurfbar');
        end
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
    end
end
end

%% Outputs
%% Display eigenvalues
if displayEigenvalues
    switch usePCA
        case 'single'
            figure('Name','Evolution of eigenvalues w.r.t order')
            clf
            R = length(sv);
            semilogy(1:R,sv(:).^2,'LineStyle','-','Color','b','LineWidth',1);
            grid on
            box on
            set(gca,'FontSize',10)
            xlabel('$\alpha$','Interpreter',interpreter)
            ylabel('$\lambda_{\alpha}$','Interpreter',interpreter)
            mysaveas(gridpathname,'eigenvalues_CY_order',formats,renderer);
            mymatlab2tikz(gridpathname,'eigenvalues_CY_order_single_PCA.tex');
            
            figure('Name','Evolution of errors')
            clf
            R = length(sv);
            semilogy(1:R,errsvdY.^2,'LineStyle','-','Color','b','LineWidth',1);
            grid on
            box on
            set(gca,'FontSize',10)
            xlabel('$R$','Interpreter',interpreter)
            ylabel('$\varepsilon_{Y}(R)$','Interpreter',interpreter)
            mysaveas(gridpathname,'error_svdYc',formats,renderer);
            mymatlab2tikz(gridpathname,'error_svdY_single_PCA.tex');
            
        case 'double'
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
                hdl(c) = semilogy(1:Rt,sv(1:Rt,t+1).^2,'LineStyle','-','Color',color(c,:),'LineWidth',1);
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
                hdl(c) = semilogy(t*dt,sv(r,t+1).^2,'LineStyle','-','Color',color(c,:),'LineWidth',1);
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
                hdl(c) = semilogy(1:Rt,errsvdY(1:Rt,t+1).^2,'LineStyle','-','Color',color(c,:),'LineWidth',1);
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
            mymatlab2tikz(gridpathname,'error_svdY_double_PCA.tex');
            
            figure('Name','Evolution of eigenvalues')
            clf
            semilogy(1:Q,sw(:).^2,'LineStyle','-','Color','b','LineWidth',1);
            grid on
            box on
            set(gca,'FontSize',fontsize)
            xlabel('$\beta$','Interpreter',interpreter)
            ylabel('$\Lambda_{\beta}$','Interpreter',interpreter)
            mysaveas(gridpathname,'eigenvalues_CZ',formats,renderer);
            mymatlab2tikz(gridpathname,'eigenvalues_CZ_double_PCA.tex');
            
            figure('Name','Evolution of errors')
            clf
            semilogy(1:Q,errsvdZ(:).^2,'LineStyle','-','Color','b','LineWidth',1);
            grid on
            box on
            set(gca,'FontSize',fontsize)
            xlabel('$Q$','Interpreter',interpreter)
            ylabel('$\varepsilon_{Z}(Q)$','Interpreter',interpreter)
            mysaveas(gridpathname,'error_svdZc',formats,renderer);
            mymatlab2tikz(gridpathname,'error_svdZ_double_PCA.tex');
            
    end
end

%% Display error between data and filtered data
if g~=gref && displayError && Filtering
    ts = (0:p)*dt;
    
    figure('Name','Error between data and filtered data')
    clf
    hdl(1) = plot(ts,erroru./normubar,'LineStyle','-','Color','b','LineWidth',1);
    hold on
    hdl(2) = plot(ts,errorC./normCbar,'LineStyle','-','Color','r','LineWidth',1);
    cpt = 2;
    if PostProcessingTau
        hdl(cpt+1) = plot(ts,errortauTime./normtauTimebar,'LineStyle','-','Color','g','LineWidth',1);
        hdl(cpt+2) = plot(ts,errordivtauConv./normdivtauConvbar,'LineStyle','-','Color','m','LineWidth',1);
        hdl(cpt+3) = plot(ts,errordivtauDiff./normdivtauDiffbar,'LineStyle','-','Color','c','LineWidth',1);
        hdl(cpt+4) = plot(ts,errortauSurf./normtauSurfbar,'LineStyle','-','Color','y','LineWidth',1);
        hdl(cpt+5) = plot(ts,errortauInterf./normtauInterfbar,'LineStyle','-','Color','k','LineWidth',1);
        cpt = cpt+5;
    end
    if PostProcessingPressure
        hdl(cpt+1) = plot(ts,errorpres./normpresbar,'LineStyle','--','Color','g','LineWidth',1);
        cpt = cpt+1;
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$t$ [s]','Interpreter',interpreter)
    ylabel('Normalized error','Interpreter',interpreter)
    leg = {'$||u||$','$\chi$'};
    if PostProcessingTau
        leg = [leg,'$||\tau_{\mathrm{time}}||$','$||\nabla \cdot \tau_{\mathrm{conv}}||$','$||\nabla \cdot \tau_{\mathrm{diff}}||$','$||\tau_{\mathrm{surf}}||$','$\tau_{\mathrm{interf}}$'];
    end
    if PostProcessingPressure
        leg = [leg,'$p$'];
    end
    l = legend(leg{:},'Location','NorthEast');
    set(l,'Interpreter',interpreter)
    mysaveas(gridpathname,'error_filter',formats,renderer);
    mymatlab2tikz(gridpathname,'error_filter.tex');
end

%% Display quantities of interest
if displayQoI && QoI
    % Mean
    mQu = reshape(mQu(1,1:dim,:,:),[3,2,p+1]);
    if PostProcessingTau
        mQtauTime = reshape(mQtau(1,1:dim,:,:),[3,2,p+1]);
        mQdivtauConv = reshape(mQtau(1,dim+(1:dim),:,:),[3,2,p+1]);
        mQdivtauDiff = reshape(mQtau(1,2*dim+(1:dim),:,:),[3,2,p+1]);
        mQtauSurf = reshape(mQtau(1,3*dim+(1:dim),:,:),[3,2,p+1]);
        mQtauInterf = reshape(mQtau(1,4*dim+1,:,:),[1,2,p+1]);
        % Correlation
        IQtauTime = IQtau(1:dim,:,1:dim,:,:);
        IQdivtauConv = IQtau(dim+(1:dim),:,dim+(1:dim),:,:);
        IQdivtauDiff = IQtau(2*dim+(1:dim),:,2*dim+(1:dim),:,:);
        IQtauSurf = IQtau(3*dim+(1:dim),:,3*dim+(1:dim),:,:);
        IQtauInterf = IQtau(4*dim+1,:,4*dim+1,:,:);
    end
    if PostProcessingPressure
        mQpres = reshape(mQpres(1,1:dim,:,:),[1,2,p+1]);
    end
    if PostProcessingEnergy
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
    xlabel('$t$ [s]','Interpreter',interpreter)
    ylabel('$u$','Interpreter',interpreter)
    leg = {'component 1 in phase 1','component 1 in phase 2','component 2 in phase 1','component 2 in phase 2','component 3 in phase 1','component 3 in phase 2'};
    legend(leg{:},'Location','NorthEast')
    mysaveas(gridpathname,'mean_u',formats,renderer);
    mymatlab2tikz(gridpathname,'mean_u.tex');
    
    if PostProcessingTau
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
        xlabel('$t$ [s]','Interpreter',interpreter)
        ylabel('$\tau_{\mathrm{time}}$','Interpreter',interpreter)
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
        xlabel('$t$ [s]','Interpreter',interpreter)
        ylabel('$\nabla \cdot \tau_{\mathrm{conv}}$','Interpreter',interpreter)
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
        xlabel('$t$ [s]','Interpreter',interpreter)
        ylabel('$\nabla \cdot \tau_{\mathrm{diff}}$','Interpreter',interpreter)
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
        xlabel('$t$ [s]','Interpreter',interpreter)
        ylabel('$\tau_{\mathrm{surf}}$','Interpreter',interpreter)
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
        xlabel('$t$ [s]','Interpreter',interpreter)
        ylabel('$\tau_{\mathrm{interf}}$','Interpreter',interpreter)
        leg = {'phase 1','phase 2'};
        legend(leg{:},'Location','NorthEast')
        mysaveas(gridpathname,'mean_tauInterf',formats,renderer);
        mymatlab2tikz(gridpathname,'mean_tauInterf.tex');
    end
    
    if PostProcessingPressure
        figure('Name','Mean of spatial average of pressure p')
        clf
        hdl(1) = plot(ts,squeeze(mQpres(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(mQpres(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$t$ [s]','Interpreter',interpreter)
        ylabel('$p$','Interpreter',interpreter)
        leg = {'phase 1','phase 2'};
        legend(leg{:},'Location','NorthEast')
        mysaveas(gridpathname,'mean_p',formats,renderer);
        mymatlab2tikz(gridpathname,'mean_p.tex');
    end
    
    figure('Name','Power of velocity u')
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
    xlabel('$\tau$ [s]','Interpreter',interpreter)
    ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter',interpreter)
    leg = {'u_1 in phase 1 - u_1 in phase 1','u_1 in phase 2 - u_1 in phase 2',...
        'u_2 in phase 1 - u_2 in phase 1','u_2 in phase 2 - u_2 in phase 2',...
        'u_3 in phase 1 - u_3 in phase 1','u_3 in phase 2 - u_3 in phase 2'};
    legend(leg{:},'Location','NorthEast')
    mysaveas(gridpathname,'power_u',formats,renderer);
    mymatlab2tikz(gridpathname,'power_u.tex');
    
    if PostProcessingTau
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
        xlabel('$\tau$ [s]','Interpreter',interpreter)
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter',interpreter)
        leg = {'$\tau_{\mathrm{time}\,1}$ in phase 1','$\tau_{\mathrm{time}\,1}$ in phase 2',...
            '$\tau_{\mathrm{time}\,2}$ in phase 1','$\tau_{\mathrm{time}\,2}$ in phase 2',...
            '$\tau_{\mathrm{time}\,3}$ in phase 1','$\tau_{\mathrm{time}\,3}$ in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter',interpreter)
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
        xlabel('$\tau$ [s]','Interpreter',interpreter)
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter',interpreter)
        leg = {'$(\nabla \cdot \tau_{\mathrm{conv}})_1$ in phase 1','$(\nabla \cdot \tau_{\mathrm{conv}})_1$ in phase 2',...
            '$(\nabla \cdot \tau_{\mathrm{conv}})_2$ in phase 1','$(\nabla \cdot \tau_{\mathrm{conv}})_2$ in phase 2',...
            '$(\nabla \cdot \tau_{\mathrm{conv}})_3$ in phase 1','$(\nabla \cdot \tau_{\mathrm{conv}})_3$ in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter',interpreter)
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
        xlabel('$\tau$ [s]','Interpreter',interpreter)
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter',interpreter)
        leg = {'$(\nabla \cdot \tau_{\mathrm{diff}})_1$ in phase 1','$(\nabla \cdot \tau_{\mathrm{diff}})_1$ in phase 2',...
            '$(\nabla \cdot \tau_{\mathrm{diff}})_2$ in phase 1','$(\nabla \cdot \tau_{\mathrm{diff}})_2$ in phase 2',...
            '$(\nabla \cdot \tau_{\mathrm{diff}})_3$ in phase 1','$(\nabla \cdot \tau_{\mathrm{diff}})_3$ in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter',interpreter)
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
        xlabel('$\tau$ [s]','Interpreter',interpreter)
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter',interpreter)
        leg = {'$\tau_{\mathrm{surf}\,1}$ in phase 1','$\tau_{\mathrm{surf}\,1}$ in phase 2',...
            '$\tau_{\mathrm{surf}\,2}$ in phase 1','$\tau_{\mathrm{surf}\,2}$ in phase 2',...
            '$\tau_{\mathrm{surf}\,3}$ in phase 1','$\tau_{\mathrm{surf}\,3}$ in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter',interpreter)
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
        xlabel('$\tau$ [s]','Interpreter',interpreter)
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter',interpreter)
        leg = {'$\tau_{\mathrm{interf}}$ in phase 1','$\tau_{\mathrm{interf}}$ in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter',interpreter)
        mysaveas(gridpathname,'power_tauInterf',formats,renderer);
        mymatlab2tikz(gridpathname,'power_tauInterf.tex');
    end
    
    if PostProcessingPressure
        figure('Name','Power of pressure')
        clf
        hdl(1) = plot(ts,squeeze(IQpres(1,1,1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(IQpres(1,2,1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$\tau$ [s]','Interpreter',interpreter)
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter',interpreter)
        leg = {'$p$ in phase 1','$p$ in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter',interpreter)
        mysaveas(gridpathname,'power_p',formats,renderer);
        mymatlab2tikz(gridpathname,'power_p.tex');
    end
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

%% Mean solution
mYr = reshape(mY,[1,n,m,p+1]);
mU = reshape(mYr(1,1:dim,:,:),[dim*m,p+1]);
mC = reshape(mYr(1,dim+1,:,:),[m,p+1]);
if PostProcessingTau
    mtauTime = reshape(mTau(1,1:dim,:,:),[dim*m,p+1]);
    mdivtauConv = reshape(mTau(1,dim+(1:dim),:,:),[dim*m,p+1]);
    mdivtauDiff = reshape(mTau(1,2*dim+(1:dim),:,:),[dim*m,p+1]);
    mtauSurf = reshape(mTau(1,3*dim+(1:dim),:,:),[dim*m,p+1]);
    mtauInterf = reshape(mTau(1,4*dim+1,:,:),[m,p+1]);
end
if PostProcessingPressure
    mpres = reshape(mPres(1,1,:,:),[m,p+1]);
end

if g~=gref && Filtering
    % Mean filtered solution
    mUbar = reshape(mYbar(1,1:dim,:,:),[dim*m,p+1]);
    mCbar = reshape(mYbar(1,dim+1,:,:),[m,p+1]);
    clear mYbar
    if PostProcessingTau
        mtauTimebar = reshape(mTaubar(1,1:dim,:,:),[dim*m,p+1]);
        mdivtauConvbar = reshape(mTaubar(1,dim+(1:dim),:,:),[dim*m,p+1]);
        mdivtauDiffbar = reshape(mTaubar(1,2*dim+(1:dim),:,:),[dim*m,p+1]);
        mtauSurfbar = reshape(mTaubar(1,3*dim+(1:dim),:,:),[dim*m,p+1]);
        mtauInterfbar = reshape(mTaubar(1,4*dim+1,:,:),[m,p+1]);
        clear mTaubar
    end
    if PostProcessingPressure
        mpresbar = reshape(mPresbar(1,1,:,:),[m,p+1]);
        clear mPresbar
    end
end

if saveMean
fprintf('\nSaving mean solution');
t_save = tic;
for t=0:p
    mUt = mU(:,t+1);
    mCt = mC(:,t+1);
    fields = {mUt,mCt};
    fieldnames = {'velocity','phase'};
    if PostProcessingTau
        mtauTimet = mtauTime(:,t+1);
        mdivtauConvt = mdivtauConv(:,t+1);
        mdivtauDifft = mdivtauDiff(:,t+1);
        mtauSurft = mtauSurf(:,t+1);
        mtauInterft = mtauInterf(:,t+1);
        fields = [fields,mtauTimet,mdivtauConvt,mdivtauDifft,mtauSurft,mtauInterft];
        fieldnames = [fieldnames,'tauTime','div(tauConv)','div(tauDiff)','tauSurf','tauInterf'];
    end
    if PostProcessingPressure
        mprest = mpres(:,t+1);
        fields = [fields,mprest];
        fieldnames = [fieldnames,'pressure'];
    end
%     if PostProcessingEnergy
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
    
    if g~=gref && Filtering
        mUbart = mUbar(:,t+1);
        mCbart = mCbar(:,t+1);
        fields = [fields,mUbart,mCbart];
        fieldnames = [fieldnames,'velocity filtered','phase filtered'];
        if PostProcessingTau
            mtauTimebart = mtauTimebar(:,t+1);
            mdivtauConvbart = mdivtauConvbar(:,t+1);
            mdivtauDiffbart = mdivtauDiffbar(:,t+1);
            mtauSurfbart = mtauSurfbar(:,t+1);
            mtauInterfbart = mtauInterfbar(:,t+1);
            fields = [fields,mtauTimebart,mdivtauConvbart,mdivtauDiffbart,mtauSurfbart,mtauInterfbart];
            fieldnames = [fieldnames,'tauTime filtered','div(tauConv) filtered','div(tauDiff) filtered','tauSurf filtered','tauInterf filtered'];
        end
        if PostProcessingPressure
            mpresbart = mpresbar(:,t+1);
            fields = [fields,mpresbart];
            fieldnames = [fieldnames,'pressure filtered'];
        end
%         if PostProcessingEnergy
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
time_save = toc(t_save);
fprintf('\nelapsed time = %f s',time_save);
fprintf('\n');
end

%% Standard deviations of random solution
stdU = zeros(dim*m,p+1);
stdC = zeros(m,p+1);
if PostProcessingTau
    stdtauTime = zeros(dim*m,p+1);
    stddivtauConv = zeros(dim*m,p+1);
    stddivtauDiff = zeros(dim*m,p+1);
    stdtauSurf = zeros(dim*m,p+1);
    stdtauInterf = zeros(m,p+1);
end
if PostProcessingPressure
    stdpres = zeros(m,p+1);
end

for t=0:p
    switch usePCA
        case 'no'
            if g<2^7
                Yt = Y(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
            end
            ut = Yt(:,1:dim,:);
            Ct = Yt(:,dim+1,:);
            clear Yt
            stdUt = std(ut(:,:))';
            stdCt = std(Ct(:,:))';
            clear ut Ct
            
        case 'single'
            [Yt,mYt,Vt] = s.reconstructionAtStep(mY,V,sv,X,t+1);
            %CYt = cov(Yt');
            %CYt = s.cov(Vt,sv);
            %stdYt = sqrt(diag(CYt));
            % stdYt = std(Yt,0,2);
            stdYt = s.std(Vt,sv);
            stdYt = s.unscaling(stdYt,Ya(:,t+1),zeros(size(Yb(:,t+1))))';
            stdYt = reshape(stdYt,[n,m]);
            stdUt = stdYt(1:dim,:);
            stdCt = stdYt(dim+1,:);
            clear stdYt
            stdUt = stdUt(:);
            stdCt = stdCt(:);
            
        case 'double'
            if g<2^7
                [Yt,mYt,Vt,svt,mZt,Wt,Rt] = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t) '.mat']);
            else
                [Yt,mYt,Vt,svt,mZt,Wt,Rt] = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t) '.mat']);
            end
            %CYt = cov(Yt');
            %stdYt = sqrt(diag(CYt));
            % stdYt = std(Yt,0,2);
            stdYt = sSpace.stdDouble(Vt,svt,Wt,sw);
            stdYt = s.unscaling(stdYt,Ya(:,t+1),zeros(size(Yb(:,t+1))))';
            stdYt = reshape(stdYt,[n,m]);
            clear Yat Ybt
            stdUt = stdYt(1:dim,:);
            stdCt = stdYt(dim+1,:);
            clear stdYt
            stdUt = stdUt(:);
            stdCt = stdCt(:);
    end
    
    if PostProcessingTau
        if g<2^7
            Taut = Tau(:,:,:,t+1);
        else
            load(fullfile(gridpathname,['data_tau_t' num2str(t) '.mat']),'Taut');
        end
        tauTimet = Taut(:,1:dim,:);
        divtauConvt = Taut(:,dim+(1:dim),:);
        divtauDifft = Taut(:,2*dim+(1:dim),:);
        tauSurft = Taut(:,3*dim+(1:dim),:);
        tauInterft = Taut(:,4*dim+1,:);
        clear Taut
        stdtauTimet = std(tauTimet(:,:))';
        stddivtauConvt = std(divtauConvt(:,:))';
        stddivtauDifft = std(divtauDifft(:,:))';
        stdtauSurft = std(tauSurft(:,:))';
        stdtauInterft = std(tauInterft(:,:))';
        clear tauTimet divtauConvt divtauDifft tauSurft tauInterft
    end
    
    if PostProcessingPressure
        if g<2^7
            Prest = Pres(:,:,:,t+1);
        else
            load(fullfile(gridpathname,['data_pressure_t' num2str(t) '.mat']),'Prest');
        end
        stdprest = std(Prest(:,:))';
        clear Prest
    end
    
    stdU(:,t+1) = stdUt;
    stdC(:,t+1) = stdCt;
    if PostProcessingTau
        stdtauTime(:,t+1) = stdtauTimet;
        stddivtauConv(:,t+1) = stddivtauConvt;
        stddivtauDiff(:,t+1) = stddivtauDifft;
        stdtauSurf(:,t+1) = stdtauSurft;
        stdtauInterf(:,t+1) = stdtauInterft;
    end
    if PostProcessingPressure
        stdpres(:,t+1) = stdprest;
    end
end

if saveStd
fprintf('\nSaving standard deviation of random solution');
t_save = tic;
for t=0:p
    stdUt = stdU(:,t+1);
    stdCt = stdC(:,t+1);
    fields = {stdUt,stdCt};
    fieldnames = {'velocity','phase'};
    if PostProcessingTau
        stdtauTimet = stdtauTime(:,t+1);
        stddivtauConvt = stddivtauConv(:,t+1);
        stddivtauDifft = stddivtauDiff(:,t+1);
        stdtauSurft = stdtauSurf(:,t+1);
        stdtauInterft = stdtauInterf(:,t+1);
        fields = [fields,stdtauTimet,stddivtauConvt,stddivtauDifft,stdtauSurft,stdtauInterft];
        fieldnames = [fieldnames,'tauTime','div(tauConv)','div(tauDiff)','tauSurf','tauInterf'];
    end
    if PostProcessingPressure
        stdprest = stdpres(:,t+1);
        fields = [fields,stdprest];
        fieldnames = [fieldnames,'pressure'];
    end
    
    if g~=gref && Filtering
        stdUbart = stdUbar(:,t+1);
        stdCbart = stdCbar(:,t+1);
        fields = [fields,stdUbart,stdCbart];
        fieldnames = [fieldnames,'velocity filtered','phase filtered'];
        if PostProcessingTau
            stdtauTimebart = stdtauTimebar(:,t+1);
            stddivtauConvbart = stddivtauConvbar(:,t+1);
            stddivtauDiffbart = stddivtauDiffbar(:,t+1);
            stdtauSurfbart = stdtauSurfbar(:,t+1);
            stdtauInterfbart = stdtauInterfbar(:,t+1);
            fields = [fields,stdtauTimebart,stddivtauConvbart,stddivtauDiffbart,stdtauSurfbart,stdtauInterfbart];
            fieldnames = [fieldnames,'tauTime filtered','div(tauConv) filtered','div(tauDiff) filtered','tauSurf filtered','tauInterf filtered'];
        end
        if PostProcessingPressure
            stdpresbart = stdpresbar(:,t+1);
            fields = [fields,stdpresbart];
            fieldnames = [fieldnames,'pressure implicit filtered'];
        end
    end
    write_vtk_mesh(M,fields,[],fieldnames,[],gridpathname,['diphasic_fluids_grid' num2str(g) '_std'],1,t);
end
make_pvd_file(gridpathname,['diphasic_fluids_grid' num2str(g) '_std'],1,p+1);
time_save = toc(t_save);
fprintf('\nelapsed time = %f s',time_save);
fprintf('\n');
end

%% Samples of random solution
if saveSamples
fprintf('\nSaving samples of random solution');
t_save = tic;
for t=0:p
    t_savet = tic;
    switch usePCA
        case 'no'
            if g<2^7
                Yt = Y(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
            end
            
        case 'single'
            Yt = s.reconstructionAtStep(mY,V,sv,X,t+1);
            Yt = s.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
            
        case 'double'
            if g<2^7
                Yt = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t) '.mat']);
            else
                Yt = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',['PCAspace_t' num2str(t) '.mat']);
            end
            Yt = sSpace.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
    end
    
    Ut = reshape(Yt(:,1:dim,:),[N,dim*m]);
    Ct = reshape(Yt(:,dim+1,:),[N,m]);
    clear Yt
    if PostProcessingTau
        Taut = Tau(:,:,:,t+1);
        tauTimet = reshape(Taut(:,1:dim,:),[N,dim*m]);
        divtauConvt = reshape(Taut(:,dim+(1:dim),:),[N,dim*m]);
        divtauDifft = reshape(Tau(:,2*dim+(1:dim),:,t+1),[N,dim*m]);
        tauSurft = reshape(Tau(:,3*dim+(1:dim),:,t+1),[N,dim*m]);
        tauInterft = reshape(Tau(:,4*dim+1,:,t+1),[N,m]);
        clear Taut
    end
    if PostProcessingPressure
        Prest = reshape(Pres(:,1,:,t+1),[N,m]);
    end
    
    for l=1:N
        Ult = Ut(l,:);
        Clt = Ct(l,:);
        fields = {Ult,Clt};
        fieldnames = {'velocity','phase'};
        if PostProcessingTau
            tauTimelt = tauTimet(l,:);
            divtauConvlt = divtauConvt(l,:);
            divtauDifflt = divtauDifft(l,:);
            tauSurflt = tauSurft(l,:);
            tauInterflt = tauInterft(l,:);
            fields = [fields,tauTimelt,divtauConvlt,divtauDifflt,tauSurflt,tauInterflt];
            fieldnames = [fieldnames,'tauTime','div(tauConv)','div(tauDiff)','tauSurf','tauInterf'];
        end
        if PostProcessingPressure
            preslt = prest(l,:);
            fields = [fields,preslt];
            fieldnames = [fieldnames,'pressure'];
        end
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
end

time_Total = toc(t_Total);
fprintf('\nElapsed time = %f s\n',time_Total);

myparallel('stop');
