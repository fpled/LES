clc
clearvars
close all
% rng('default');
myparallel('start');

usePCA = 'no'; % 'no', 'single', 'double'
PostProcessingTau = true;
PostProcessingPressure = true;
PostProcessingEnergy = true;
computeStatistics = true;
QoI = true;
Filtering = false;
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

displayEigenvalues = false;
displayCovariance = false;
displayStatistics = true;
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
gset = 2.^(4:5); % set of spatial grid sizes
g = gset(2); % current spatial grid size
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
switch usePCA
    case 'no'
        prefix = '';
    otherwise
        prefix = [usePCA '_PCA_'];
end   

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

%% Loading DNS data
if g<2^7 && (strcmpi(usePCA,'no') || (strcmpi(usePCA,'single') && performPCA) || (strcmpi(usePCA,'double') && performPCAspace))
    fprintf('\nLoading DNS data');
    t_load = tic;
    load(fullfile(gridpathname,'data.mat'),'Y');
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
    Y = Y(1:10,:,:,:);
    N = size(Y,1);
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
                    fprintf('\nTime step %2d/%d : elapsed time = %f s',t,p,time_Meant);
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
            clear Y
            % Y_approx = s.reconstruction(mY,V,sv,X);
            % Y_approx = reshape(Y_approx,[n*m,p+1,N]);
            % norm2Y_approx = sum(var(Y_approx,0,3));
            mY = reshape(mY,[n*m,1,p+1]);
            V = permute(reshape(V',[R,n*m,p+1]),[2 1 3]);
            norm2Y_approx = arrayfun(@(t) sum(s.var(V(:,:,t+1),sv)),0:p);
            err2Y = abs(norm2Y-norm2Y_approx);
            
            ts = (0:p)*dt;
            errL2 = trapz(ts,err2Y,2)/trapz(ts,norm2Y,2);
            fprintf('\nL2-error = %.3e for Y',errL2);
            fprintf('\n');
            
            save(fullfile(gridpathname,[prefix 'scaling.mat']),'Ya','Yb');
            save(fullfile(gridpathname,[prefix 'data.mat']),'s','mY','V','sv','X','R','errsvdY','errL2','time_PCA');
        else
            fprintf('\nLoading DNS data from PCA');
            t_load = tic;
            load(fullfile(gridpathname,[prefix 'scaling.mat']),'Ya','Yb');
            load(fullfile(gridpathname,[prefix 'data.mat']),'s','mY','V','sv','X','R','errsvdY','errL2','time_PCA');
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
                fprintf('\nTime step %2d/%d : rank R = %d, error = %.3e for Y, elapsed time = %f s',t,p,Rt,errYt,time_PCA_spacet);
                
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
                    save(fullfile(gridpathname,[prefix 'space_data_t' num2str(t) '.mat']),'Vt');
                end
                sv(1:Rt,t+1) = svt;
                Z(:,1:Rt,t+1) = Zt;
                R(t+1) = Rt;
                errsvdY(1:Rt,t+1) = errsvdYt;
                err2Y(t+1) = err2Yt;
                norm2Y(t+1) = norm2Yt;
                clear Yt mYt Yat Ybt Vt svt Zt Rt errsvdYt err2Yt norm2Yt
            end
            clear Y
            
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
            
            save(fullfile(gridpathname,[prefix 'scaling.mat']),'Ya','Yb');
            save(fullfile(gridpathname,[prefix 'space_data.mat']),'sSpace','mY','sv','Z','R','errsvdY','errL2Y','time_PCA_space');
            if g<2^7
                save(fullfile(gridpathname,[prefix 'space_data.mat']),'V','-append');
            end
        else
            fprintf('\nLoading DNS data from PCA in space');
            t_load = tic;
            load(fullfile(gridpathname,[prefix 'scaling.mat']),'Ya','Yb');
            load(fullfile(gridpathname,[prefix 'space_data.mat']),'sSpace','mY','sv','Z','R','errsvdY','errL2Y','time_PCA_space');
            Rmax = max(R);
            if g<2^7
                load(fullfile(gridpathname,[prefix 'space_data.mat']),'V');
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
            clear Z
            % Z_approx = sTime.reconstruction(mZ,W,sw,X);
            % Z_approx = reshape(Z_approx,[Rmax,p+1,N]);
            % norm2Z_approx = sum(var(Z_approx,0,3));
            W = reshape(W',[Q,Rmax,p+1]);
            norm2Z_approx = arrayfun(@(t) sum(sTime.var(W(:,1:R(t+1),t+1)',sw)),0:p);
            err2Z = abs(norm2Z-norm2Z_approx);
            
            ts = (0:p)*dt;
            errL2Z = trapz(ts,err2Z,2)/trapz(ts,norm2Z,2);
            fprintf('\nL2-error = %.3e for Z',errL2Z);
            
            mZ = reshape(mZ,[1,Rmax,p+1]);
            err2Y = zeros(1,p+1);
            for t=0:p
                Rt = R(t+1);
                if g<2^7
                    Vt = V(:,1:Rt,t+1);
                else
                    load(fullfile(gridpathname,[prefix 'space_data_t' num2str(t) '.mat']),'Vt');
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
            
            save(fullfile(gridpathname,[prefix 'time_data.mat']),'sTime','mZ','W','sw','X','Q','errsvdZ','errL2Z','errL2Y','time_PCA_time');
        else
            fprintf('\nLoading DNS data from PCA time');
            t_load = tic;
            load(fullfile(gridpathname,[prefix 'time_data.mat']),'sTime','mZ','W','sw','X','Q','errsvdZ','errL2Z','time_PCA_time');
            time_load = toc(t_load);
            fprintf('\nelapsed time = %f s',time_load);
            fprintf('\n');
        end
        
    otherwise
        error('Method not implemented.')
end

if postProcess && (PostProcessingTau || PostProcessingPressure || PostProcessingEnergy)
%% Postprocessing DNS data
fprintf('\nPostprocessing DNS data');
t_PostProcess = tic;
if PostProcessingTau
    ntau = 4*dim+1; % number of tau variables
    fprintf('\nn = %d tau variables',ntau);
    mTau = zeros(1,ntau,m,p+1);
    if g<2^7
        Tau = zeros(N,ntau,m,p+1);
    end
end
if PostProcessingPressure
    if ~PostProcessingTau
        load(fullfile(gridpathname,[prefix 'mean_data_tau.mat']),'mTau','ntau');
        if g<2^7
            load(fullfile(gridpathname,[prefix 'data_tau.mat']),'Tau');
        end
    end
    fprintf('\nn = %d pressure variable',1);
    mpres = zeros(1,1,m,p+1);
    if g<2^7
        pres = zeros(N,1,m,p+1);
    end
end
if PostProcessingEnergy
    if ~PostProcessingPressure
        load(fullfile(gridpathname,[prefix 'mean_data_pressure.mat']),'mpres');
        if g<2^7
            load(fullfile(gridpathname,[prefix 'data_pressure.mat']),'pres');
        end
    end
    ne = 9; % number of energy variables
    fprintf('\nn = %d energy variables',ne);
    mE = zeros(1,ne,m,p+1);
    if g<2^7
        E = zeros(N,ne,m,p+1);
    end
end
fprintf('\n');
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
                    Yt_old = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t-1) '.mat']);
                else
                    Yt_old = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t-1) '.mat']);
                end
                Yt_old = sSpace.unscaling(Yt_old,Ya(:,t),Yb(:,t))';
            end
            if t<p
                if g<2^7
                    Yt_new = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+2,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t+1) '.mat']);
                else
                    Yt_new = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+2,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t+1) '.mat']);
                end
                Yt_new = sSpace.unscaling(Yt_new,Ya(:,t+2),Yb(:,t+2))';
            end
            if g<2^7
                Yt = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
            else
                Yt = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
            end
            Yt = sSpace.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
    end
    Yt = perm(reshape(Yt,[N,n,sx]));
    ut = Yt(1:dim,:,:,:,:);
    pt = Yt(dim+1,:,:,:,:);
    Ct = Yt(dim+2,:,:,:,:);
    clear Yt
    
%     tConv = zeros(3,g+1,g+1,g+1,N);
%     tDiff = zeros(3,g+1,g+1,g+1,N);
%     tSurf = zeros(3,g+1,g+1,g+1,N);
%     tInterf = zeros(1,g+1,g+1,g+1,N);
%     for l=1:N
%         samplename = ['Sample ' num2str(l)];
%         disp(samplename)
%         
%         % ul = reshape(ut(:,:,:,:,l),[dim,m]);
%         % Cl = reshape(Ct(:,:,:,:,l),[1,m]);
%         % uijk = zeros(3,g+1,g+1,g+1);
%         % Cijk = zeros(g+1,g+1,g+1);
%         % for i=1:g+1
%         %     for j=1:g+1
%         %         for k=1:g+1
%         %             ind = (g+1)^2*(k-1)+(g+1)*(j-1)+i;
%         %             uijk(:,i,j,k) = ul(:,ind);
%         %             Cijk(1,i,j,k) = Cl(1,ind);
%         %         end
%         %     end
%         % end
%         uijk = ut(:,:,:,:,l);
%         Cijk = Ct(:,:,:,:,l);
%         
%         graduijk = zeros(3,3,g+1,g+1,g+1);
%         gradCijk = zeros(3,g+1,g+1,g+1);
%         divuijk = zeros(1,g+1,g+1,g+1);
%         for k=1:g+1
%             for j=1:g+1
%                 for c=1:3
%                     u1 = uijk(c,:,j,k);
%                     graduijk(1,c,:,j,k) = Dx*u1(:);
%                 end
%                 C1 = Cijk(1,:,j,k);
%                 gradCijk(1,:,j,k) = Dx*C1(:);
%                 u1 = uijk(1,:,j,k);
%                 divu1 = Dx*u1(:);
%                 divuijk(1,:,j,k) = divuijk(1,:,j,k) + reshape(divu1(:),size(divuijk(1,:,j,k)));
%             end
%         end
%         for k=1:g+1
%             for i=1:g+1
%                 for c=1:3
%                     u2 = uijk(c,i,:,k);
%                     graduijk(2,c,i,:,k) = Dx*u2(:);
%                 end
%                 C2 = Cijk(1,i,:,k);
%                 gradCijk(2,i,:,k) = Dx*C2(:);
%                 u2 = uijk(2,i,:,k);
%                 divu2 = Dx*u2(:);
%                 divuijk(1,i,:,k) = divuijk(1,i,:,k) + reshape(divu2(:),size(divuijk(1,i,:,k)));
%             end
%         end
%         for j=1:g+1
%             for i=1:g+1
%                 for c=1:3
%                     u3 = uijk(c,i,j,:);
%                     graduijk(3,c,i,j,:) = Dx*u3(:);
%                 end
%                 C3 = Cijk(1,i,j,:);
%                 gradCijk(3,i,j,:) = Dx*C3(:);
%                 u3 = uijk(3,i,j,:);
%                 divu3 = Dx*u3(:);
%                 divuijk(1,i,j,:) = divuijk(1,i,j,:) + reshape(divu3(:),size(divuijk(1,i,j,:)));
%             end
%         end
%         
%         Sijk = (graduijk+permute(graduijk,[2,1,3,4,5]))/2;
%         
%         normalijk = zeros(3,g+1,g+1,g+1);
%         for k=1:g+1
%             for j=1:g+1
%                 for i=1:g+1
%                     gradC = gradCijk(:,i,j,k);
%                     if all(gradC==0)
%                         normalijk(:,i,j,k) = zeros(3,1);
%                     else
%                         normalijk(:,i,j,k) = gradC./norm(gradC);
%                     end
%                 end
%             end
%         end
%         
%         divSijk = zeros(3,g+1,g+1,g+1);
%         divnormalijk = zeros(g+1,g+1,g+1);
%         for k=1:g+1
%             for j=1:g+1
%                 S1 = sum(Sijk(1,:,:,j,k),2);
%                 n1 = normalijk(1,:,j,k);
%                 divSijk(1,:,j,k) = Dx*S1(:);
%                 divnormalijk(:,j,k) = Dx*n1(:);
%             end
%         end
%         for k=1:g+1
%             for i=1:g+1
%                 S2 = sum(Sijk(2,:,i,:,k),2);
%                 n2 = normalijk(2,i,:,k);
%                 divn2 = divnormalijk(i,:,k);
%                 divSijk(2,i,:,k) = Dx*S2(:);
%                 divnormalijk(i,:,k) = divn2(:) + Dx*n2(:);
%             end
%         end
%         for j=1:g+1
%             for i=1:g+1
%                 S3 = sum(Sijk(3,:,i,j,:),2);
%                 n3 = normalijk(3,i,j,:);
%                 divn3 = divnormalijk(i,j,:);
%                 divSijk(3,i,j,:) = Dx*S3(:);
%                 divnormalijk(i,j,:) = divn3(:)+ Dx*n3(:);
%             end
%         end
%         
%         for k=1:g+1
%             for j=1:g+1
%                 for i=1:g+1
%                     mmu  = Cijk(1,i,j,k)*mu(2)  + (1-Cijk(1,i,j,k))*mu(1);
%                     rrho = Cijk(1,i,j,k)*rho(2) + (1-Cijk(1,i,j,k))*rho(1);
%                     
%                     uu = uijk(:,i,j,k);
%                     gradu = graduijk(:,:,i,j,k);
%                     gradC = gradCijk(:,i,j,k);
%                     divS  = divSijk(:,i,j,k);
%                     kappa = divnormalijk(i,j,k);
%                     
%                     tConv(:,i,j,k,l) = rrho*gradu*uu;
%                     tDiff(:,i,j,k,l) = 2*mmu*divS;
%                     tSurf(:,i,j,k,l) = sigma*kappa*gradC;
%                     tInterf(1,i,j,k,l) = dot(uu,gradC,1);
%                 end
%             end
%         end
%     end
    
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
        Ct_old = Yt_old(dim+2,:,:,:,:);
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
        Ct_new = Yt_new(dim+2,:,:,:,:);
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
    % divut = reshape(gradut(1,1,:,:,:,:)+gradut(2,2,:,:,:,:)+gradut(3,3,:,:,:,:),[1,sx,N]);
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
            % divtauConvt_approx = rhot.*squeeze(bsxfun(@times,sum(grad(ut,Dx),shiftdim(ut,-1)),2));
        else
            divtauConvt = squeeze(sum(grad(rhout,Dx).*shiftdim(ut,-1),2));
            % divtauConvt_approx = rhot.*squeeze(sum(grad(ut,Dx).*shiftdim(ut,-1),2));
        end
        % divtauConvt = squeeze(sum(grad(rhout,Dx).*repmat(shiftdim(ut,-1),[dim,ones(1,ndims(ut))]),2));
        % other formulae
        % if verLessThan('matlab','9.1') % compatibility (<R2016b)
        %     divtauConvt = div(permute(repmat(shiftdim(rhout,-1),[dim,ones(1,ndims(rhout))]),[2,1,dim:(dim+3)]).*repmat(shiftdim(ut,-1),[dim,ones(1,ndims(ut))]),Dx);
        % else
        %     divtauConvt = div(permute(shiftdim(rhout,-1),[2,1,3:ndims(rhout)+1]).*shiftdim(ut,-1),Dx);
        % end
        % divtauConvt = div(permute(repmat(shiftdim(rhout,-1),[dim,ones(1,ndims(rhout))]),[2,1,3:ndims(rhout)+1]).*repmat(shiftdim(ut,-1),[dim,ones(1,ndims(ut))]),Dx);
        tauInterft = dot(ut,gradCt,1);
        clear gradCt
        divtauDifft = div(2*muSt,Dx); % divtauDifft = 2*mut.*div(St,Dx);
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
            save(fullfile(gridpathname,[prefix 'data_tau_t' num2str(t) '.mat']),'Taut');
        end
        clear Taut
    end
    
    if PostProcessingPressure
        if ~PostProcessingTau
            if g<2^7
                Taut = Tau(:,:,:,t+1);
            else
                load(fullfile(gridpathname,[prefix 'data_tau_t' num2str(t) '.mat']),'Taut');
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
        prest = solve_pressure_problem(B,Rhot,Gradinvrhot,GradxN,GradyN,GradzN,LaplacianN);
        clear B Rhot Gradinvrhot
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            preshydrostatict = gravity.*bsxfun(@times,rhot,shiftdim(x(:)-L/2,-2));
        else
            preshydrostatict = gravity.*rhot.*shiftdim(x(:)-L/2,-2);
        end
        preshydrostatict = reshape(preshydrostatict,[m,N]);
        prest = prest + preshydrostatict;
        if ~PostProcessingEnergy
            clear rhot
        end
        prest = reshape(prest',[N,1,m]);
        mprest = mean(prest,1);
        mpres(1,:,:,t+1) = mprest;
        clear mPrest
        if g<2^7
            pres(:,:,:,t+1) = prest;
        else
            save(fullfile(gridpathname,[prefix 'data_pressure_t' num2str(t) '.mat']),'prest');
        end
        if ~PostProcessingEnergy
            clear prest
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
        energyKinGradt = dot(rhout,grad(u2t/2,Dx),1);
        clear u2t
        energyGravt = gravity.*rhout(2,:,:,:,:);
        clear rhout
%         if ~PostProcessingPressure
%             if g<2^7
%                 prest = pres(:,:,:,t+1);
%             else
%                 load(fullfile(gridpathname,[prefix 'data_pressure_t' num2str(t) '.mat']),'prest');
%             end
%         end
%         prest = perm(reshape(prest,[N,1,sx]));
        prest = pt;
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
        Et = cat(1,energyKinTimet,energyConvt,energyGravt,energyPrest,energyPresDilt,energyKinGradt,energyDifft,energyVisct,energySurft);
        clear energyKinTimet energyConvt energyGravt energyPrest energyPresDilt energyKinGradt energyDifft energyVisct energySurft
        Et = iperm(Et);
        Et = Et(:,:,:);
        mEt = mean(Et,1);
        mE(:,:,:,t+1) = mEt;
        clear mEt
        if g<2^7
            E(:,:,:,t+1) = Et;
        else
            save(fullfile(gridpathname,[prefix 'data_energy_t' num2str(t) '.mat']),'Et');
        end
        clear Et
    end
    time_PostProcesst = toc(t_PostProcesst);
    fprintf('Time step %2d/%d : elapsed time = %f s\n',t,p,time_PostProcesst);
end

time_PostProcess = toc(t_PostProcess);
fprintf('\nelapsed time = %f s',time_PostProcess);
fprintf('\n');

if PostProcessingTau
    save(fullfile(gridpathname,[prefix 'mean_data_tau.mat']),'mTau','ntau');
    if g<2^7
        save(fullfile(gridpathname,[prefix 'data_tau.mat']),'Tau');
    end
end
if PostProcessingPressure
    save(fullfile(gridpathname,[prefix 'mean_data_pressure.mat']),'mpres');
    if g<2^7
        save(fullfile(gridpathname,[prefix 'data_pressure.mat']),'pres');
    end
end
if PostProcessingEnergy
    save(fullfile(gridpathname,[prefix 'mean_data_energy.mat']),'mE','ne');
    if g<2^7
        save(fullfile(gridpathname,[prefix 'data_energy.mat']),'E');
    end
end
else
    if PostProcessingTau
        t_load = tic;
        fprintf('\nLoading mean tau data');
        load(fullfile(gridpathname,[prefix 'mean_data_tau.mat']),'mTau','ntau');
        if g<2^7
            fprintf('\nLoading tau data');
            load(fullfile(gridpathname,[prefix 'data_tau.mat']),'Tau');
        end
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
    end
    if PostProcessingPressure
        t_load = tic;
        fprintf('\nLoading mean pressure data');
        load(fullfile(gridpathname,[prefix 'mean_data_pressure.mat']),'mpres');
        if g<2^7
            fprintf('\nLoading pressure data');
            load(fullfile(gridpathname,[prefix 'data_pressure.mat']),'pres');
        end
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
    end
    if PostProcessingEnergy
        t_load = tic;
        fprintf('\nLoading mean energy data');
        load(fullfile(gridpathname,[prefix 'mean_data_energy.mat']),'mE','ne');
        if g<2^7
            fprintf('\nLoading energy data');
            load(fullfile(gridpathname,[prefix 'data_energy.mat']),'E');
        end
        time_load = toc(t_load);
        fprintf('\nelapsed time = %f s',time_load);
        fprintf('\n');
    end
end

%% Compute statistics of stochastic processes
if computeStatistics
    fprintf('\nComputing statistics of stochastic processes');
    t_Stat = tic;
    norm2Y = zeros(3,p+1);
    norm2mY = zeros(3,p+1);
    if PostProcessingTau
        norm2Tau = zeros(5,p+1);
        norm2mTau = zeros(5,p+1);
    end
    if PostProcessingPressure
        norm2pres = zeros(1,p+1);
        norm2mpres = zeros(1,p+1);
    end
    if PostProcessingEnergy
        norm2E = zeros(ne,p+1);
        norm2mE = zeros(ne,p+1);
    end
    for t=0:p
        t_Statt = tic;
        components = {1:dim,dim+1,dim+2};
        mYt = mY(1,:,:,t+1);
        % norm2mY(t+1) = arrayfun(@(i) compute_norm(x,reshape(sum(mYt(1,components{i},:).^2,2),[1,sx]),dim),1:3);
        mYt = reshape(mYt,[n,m]);
        norm2mY(:,t+1) = arrayfun(@(i) sum(sum(mYt(components{i},:).^2)),1:3); % norm2mY(:,t+1) = arrayfun(@(i) norm(mYt(components{i},:),'fro')^2,1:3);
        clear mYt
        switch usePCA
            case 'no'
                if g<2^7
                    Yt = Y(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
                end
                vYt = var(Yt);
                clear Yt
                
            case 'single'
                % Yt = s.reconstructionAtStep(mY,V,sv,X,t+1);
                % vYt = var(Yt,0,2);
                % vYt = s.unscaling(vYt,Ya(:,t+1).^2,zeros(size(Yb(:,t+1))))';
                Vt = V(:,:,t+1);
                vYt = s.var(Vt,sv);
                vYt = s.unscaling(vYt,Ya(:,t+1).^2,zeros(size(Yb(:,t+1))))';
                clear Vt
                
            case 'double'
                % if g<2^7
                %     Yt = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
                % else
                %     Yt = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
                % end
                % vYt = var(Yt,0,2);
                % vYt = sSpace.unscaling(vYt,Ya(:,t+1).^2,zeros(size(Yb(:,t+1))))';
                if g<2^7
                    [Vt,svt,Wt,Rt] = sSpace.getPrincipalComponentsDoubleAtStep(V,sv,W,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
                else
                    [Vt,svt,Wt,Rt] = sSpace.getPrincipalComponentsDoubleAtStep([],sv,W,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
                end
                vYt = sSpace.varDouble(Vt,svt,Wt,sw);
                vYt = sSpace.unscaling(vYt,Ya(:,t+1).^2,zeros(size(Yb(:,t+1))))';
                clear Vt svt Wt
        end
        vYt = reshape(vYt,[n,m]);
        norm2Y(:,t+1) = arrayfun(@(i) sum(sum(vYt(components{i},:))),1:3);
        
        if PostProcessingTau
            components = {1:dim,dim+(1:dim),2*dim+(1:dim),3*dim+(1:dim),4*dim+1};
            mTaut = reshape(mTau(1,:,:,t+1),[ntau,m]);
            norm2mTau(:,t+1) = arrayfun(@(i) sum(sum(mTaut(components{i},:).^2)),1:5); % norm2mTau(:,t+1) = arrayfun(@(i) norm(mTaut(components{i},:),'fro')^2,1:5);
            clear mTaut
            if g<2^7
                Taut = Tau(:,:,:,t+1);
            else
                load(fullfile(gridpathname,[prefix 'data_tau_t' num2str(t) '.mat']),'Taut');
            end
            vTaut = reshape(var(Taut),[ntau,m]);
            clear Taut
            norm2Tau(:,t+1) = arrayfun(@(i) sum(sum(vTaut(components{i},:))),1:5);
            clear vTaut
        end
        if PostProcessingPressure
            mprest = mpres(1,:,:,t+1);
            norm2mpres(t+1) = sum(mprest(:).^2); % norm2mpres(t+1) = norm(mprest(:))^2;
            clear mprest
            if g<2^7
                prest = pres(:,:,:,t+1);
            else
                load(fullfile(gridpathname,[prefix 'data_pressure_t' num2str(t) '.mat']),'prest');
            end
            vprest = var(prest);
            clear prest
            norm2pres(t+1) = sum(vprest(:));
            clear vprest
        end
        if PostProcessingEnergy
            mEt = reshape(mE(1,:,:,t+1),[ne,m]);
            norm2mE(:,t+1) = arrayfun(@(i) sum(sum(mEt(i,:).^2)),1:ne); % norm2mE(:,t+1) = arrayfun(@(i) norm(mEt(i,:),'fro')^2,1:ne);
            clear mEt
            if g<2^7
                Et = E(:,:,:,t+1);
            else
                load(fullfile(gridpathname,[prefix 'data_energy_t' num2str(t) '.mat']),'Et');
            end
            vEt = reshape(var(Et),[ne,m]);
            clear Et
            norm2E(:,t+1) = arrayfun(@(i) sum(sum(vEt(i,:))),1:ne);
            clear vEt
        end
        time_Statt = toc(t_Statt);
        fprintf('\nTime %2d/%d : elapsed time = %f s',t,p,time_Statt);
    end
    
    time_Stat = toc(t_Stat);
    fprintf('\nelapsed time = %f s',time_Stat);
    fprintf('\n');
    
    save(fullfile(gridpathname,[prefix 'statistics.mat']),'norm2Y','norm2mY');
    if PostProcessingTau
        save(fullfile(gridpathname,[prefix 'statistics_tau.mat']),'norm2Tau','norm2mTau');
    end
    if PostProcessingPressure
        save(fullfile(gridpathname,[prefix 'statistics_pressure.mat']),'norm2pres','norm2mpres');
    end
    if PostProcessingEnergy
        save(fullfile(gridpathname,[prefix 'statistics_energy.mat']),'norm2E','norm2mE');
    end
else
    fprintf('\nLoading statistics of stochastic processes');
    t_load = tic;
    load(fullfile(gridpathname,[prefix 'statistics.mat']),'norm2Y','norm2mY');
    if PostProcessingTau
        load(fullfile(gridpathname,[prefix 'statistics_tau.mat']),'norm2Tau','norm2mTau');
    end
    if PostProcessingPressure
        load(fullfile(gridpathname,[prefix 'statistics_pressure.mat']),'norm2pres','norm2mpres');
    end
    if PostProcessingEnergy
        load(fullfile(gridpathname,[prefix 'statistics_energy.mat']),'norm2E','norm2mE');
    end
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
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
    if PostProcessingPressure
        Qpres = zeros(N,2,2,p+1);
    else
        Qpres = zeros(N,1,2,p+1);
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
                    Yt = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
                else
                    Yt = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
                end
                Yt = sSpace.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
        end
        Yt = reshape(Yt,[N,n,sx]);
        ut = Yt(:,1:dim,:,:,:);
        pt = Yt(:,dim+1,:,:,:);
        Ct = Yt(:,dim+2,:,:,:);
        clear Yt
        Qut = int_trapz(x,Ct,ut,dim);
        Qu(:,:,:,t+1) = Qut;
        clear ut Qut
        Qpt = int_trapz(x,Ct,pt,dim);
        Qpres(:,1,:,t+1) = Qpt;
        clear pt Qpt
        if PostProcessingTau
            if g<2^7
                Taut = Tau(:,:,:,t+1);
            else
                load(fullfile(gridpathname,[prefix 'data_tau_t' num2str(t) '.mat']),'Taut');
            end
            Taut = reshape(Taut,[N,ntau,sx]);
            Qtaut = int_trapz(x,Ct,Taut,dim);
            Qtau(:,:,:,t+1) = Qtaut;
            clear Taut Qtaut
        end
        if PostProcessingPressure
            if g<2^7
                prest = pres(:,:,:,t+1);
            else
                load(fullfile(gridpathname,[prefix 'data_pressure_t' num2str(t) '.mat']),'Prest');
            end
            prest = reshape(prest,[N,1,sx]);
            Qprest = int_trapz(x,Ct,prest,dim);
            Qpres(:,2,:,t+1) = Qprest;
            clear Prest Qprest
        end
        if PostProcessingEnergy
            if g<2^7
                Et = E(:,:,:,t+1);
            else
                load(fullfile(gridpathname,[prefix 'data_energy_t' num2str(t) '.mat']),'Et');
            end
            Et = reshape(Et,[N,ne,sx]);
            Qet = int_trapz(x,Ct,Et,dim);
            Qe(:,:,:,t+1) = Qet;
            clear Et Qet
        end
        clear Ct
        time_QoIt = toc(t_QoIt);
        fprintf('\nTime step %2d/%d : elapsed time = %f s',t,p,time_QoIt);
    end
    
    [mQu,stdQu,RQu,IQu] = compute_stats_qoi(Qu,dt);
    if PostProcessingTau
        [mQtau,stdQtau,RQtau,IQtau] = compute_stats_qoi(Qtau,dt);
    end
    if PostProcessingPressure
        [mQpres,stdQpres,RQpres,IQpres] = compute_stats_qoi(Qpres,dt);
    end
    if PostProcessingEnergy
        [mQe,stdQe,RQe,IQe] = compute_stats_qoi(Qe,dt);
    end
    
    time_QoI = toc(t_QoI);
    fprintf('\nelapsed time = %f s',time_QoI);
    fprintf('\n');
    save(fullfile(gridpathname,[prefix 'QoI.mat']),'Qu','mQu','IQu');
    if PostProcessingTau
        save(fullfile(gridpathname,[prefix 'QoI_tau.mat']),'Qtau','mQtau','IQtau');
    end
    if PostProcessingPressure
        save(fullfile(gridpathname,[prefix 'QoI_pressure.mat']),'Qpres','mQpres','IQpres');
    end
    if PostProcessingEnergy
        save(fullfile(gridpathname,[prefix 'QoI_energy.mat']),'Qe','mQe','IQe');
    end
else
    fprintf('\nLoading quantities of interest');
    t_load = tic;
    load(fullfile(gridpathname,[prefix 'QoI.mat']),'Qu','mQu','IQu');
    if PostProcessingTau
        load(fullfile(gridpathname,[prefix 'QoI_tau.mat']),'Qtau','mQtau','IQtau');
    end
    if PostProcessingPressure
        load(fullfile(gridpathname,[prefix 'QoI_pressure.mat']),'Qpres','mQpres','IQpres');
    end
    if PostProcessingEnergy
        load(fullfile(gridpathname,[prefix 'QoI_energy.mat']),'Qe','mQe','IQe');
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
                        Yt = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
                    else
                        Yt = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
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
            fprintf('\nTime %2d/%d : elapsed time = %f s',t,p,time_Filtergt);
            
            Ybart = reshape(Ybart(:,:,1:k:end,1:k:end,1:k:end),[N,n,mbar]);
            mYbar(1,:,:,t+1) = mean(Ybart,1);
            if gbar<2^7
                Ybar(:,:,:,t+1) = Ybart;
            else
                save(fullfile(gridbarpathname,[prefix 'data_filtered_t' num2str(t) '.mat']),'Ybart');
            end
            clear Ybart
        end
        time_Filterg = toc(t_Filterg);
        fprintf('\nelapsed time = %f s for coarse grid %d',time_Filterg,gbar);
        
        save(fullfile(gridbarpathname,[prefix 'mean_data_filtered.mat']),'mYbar');
        clear mYbar
        if gbar<2^7
            save(fullfile(gridbarpathname,[prefix 'data_filtered.mat']),'Ybar');
            clear Ybar
        end
    end
    time_Filter = toc(t_Filter);
    fprintf('\nelapsed time = %f s for all coarse grids',time_Filter);
    fprintf('\n');
    
    if PostProcessingTau
        fprintf('\nApplying filter for tau data');
        t_Filter = tic;
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
                    load(fullfile(gridpathname,[prefix 'data_tau_t' num2str(t) '.mat']),'Taut');
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
                fprintf('\nTime %2d/%d : elapsed time = %f s',t,p,time_Filtergt);
                
                Taubart = reshape(Taubart(:,:,1:k:end,1:k:end,1:k:end),[N,ntau,mbar]);
                mTaubar(1,:,:,t+1) = mean(Taubart,1);
                if gbar<2^7
                    Taubar(:,:,:,t+1) = Taubart;
                else
                    save(fullfile(gridbarpathname,[prefix 'data_tau_filtered_t' num2str(t) '.mat']),'Taubart');
                end
                clear Taubart
            end
            time_Filterg = toc(t_Filterg);
            fprintf('\nelapsed time = %f s for coarse grid %d',time_Filterg,gbar);
            
            save(fullfile(gridbarpathname,[prefix 'mean_data_tau_filtered.mat']),'mTaubar');
            clear mTaubar
            if gbar<2^7
                save(fullfile(gridbarpathname,[prefix 'data_tau_filtered.mat']),'Taubar');
                clear Taubar
            end
        end
        time_Filter = toc(t_Filter);
        fprintf('\nelapsed time = %f s for all coarse grids',time_Filter);
        fprintf('\n');
    end
    
    if PostProcessingPressure
        fprintf('\nApplying filter for pressure data');
        t_Filter = tic;
        for ig=1:ng-1
            gbar = gset(ng-ig);
            gridbarname = ['Grid' num2str(gbar)];
            gridbarpathname = fullfile(pathname,gridbarname);
            mbar = (gbar+1)^dim;
            k = g/gbar;
            h = hset{ig};
            H = shiftdim(h,-2);
            
            fprintf('\nFiltering pressure data on coarse grid %d',gbar);
            t_Filterg = tic;
            mPresbar = zeros(1,1,mbar,p+1);
            if g<2^7
                Presbar = zeros(N,1,mbar,p+1);
            end
            for t=0:p
                t_Filtergt = tic;
                if g<2^7
                    prest = pres(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,[prefix 'data_pressure_t' num2str(t) '.mat']),'Prest');
                end
                Presbart = reshape(prest,[N,1,sx]);
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
                fprintf('\nTime %2d/%d : elapsed time = %f s',t,p,time_Filtergt);
                
                Presbart = reshape(Presbart(:,:,1:k:end,1:k:end,1:k:end),[N,1,mbar]);
                mPresbar(1,:,:,t+1) = mean(Presbart,1);
                if gbar<2^7
                    Presbar(:,:,:,t+1) = Presbart;
                else
                    save(fullfile(gridbarpathname,[prefix 'data_pressure_filtered_t' num2str(t) '.mat']),'Presbart');
                end
                clear Presbart
            end
            time_Filterg = toc(t_Filterg);
            fprintf('\nelapsed time = %f s for coarse grid %d',time_Filterg,gbar);
            
            save(fullfile(gridbarpathname,[prefix 'mean_data_pressure_filtered.mat']),'mPresbar');
            clear mPresbar
            if gbar<2^7
                save(fullfile(gridbarpathname,[prefix 'data_pressure_filtered.mat']),'Presbar');
                clear Presbar
            end
        end
        time_Filter = toc(t_Filter);
        fprintf('\nelapsed time = %f s for all coarse grids',time_Filter);
        fprintf('\n');
    end
    
    if PostProcessingEnergy
        fprintf('\nApplying filter for energy data');
        t_Filter = tic;
        for ig=1:ng-1
            gbar = gset(ng-ig);
            gridbarname = ['Grid' num2str(gbar)];
            gridbarpathname = fullfile(pathname,gridbarname);
            mbar = (gbar+1)^dim;
            k = g/gbar;
            h = hset{ig};
            H = shiftdim(h,-2);
            
            fprintf('\nFiltering energy data on coarse grid %d',gbar);
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
                    load(fullfile(gridpathname,[prefix 'data_energy_t' num2str(t) '.mat']),'Et');
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
                fprintf('\nTime %2d/%d : elapsed time = %f s',t,p,time_Filtergt);
                
                Ebart = reshape(Ebart(:,:,1:k:end,1:k:end,1:k:end),[N,ne,mbar]);
                mEbar(1,:,:,t+1) = mean(Ebart,1);
                if gbar<2^7
                    Ebar(:,:,:,t+1) = Ebart;
                else
                    save(fullfile(gridbarpathname,[prefix 'data_energy_filtered_t' num2str(t) '.mat']),'Ebart');
                end
                clear Ebart
            end
            time_Filterg = toc(t_Filterg);
            fprintf('\nelapsed time = %f s for coarse grid %d',time_Filterg,gbar);
            
            save(fullfile(gridbarpathname,[prefix 'mean_data_energy_filtered.mat']),'mEbar');
            clear mEbar
            if gbar<2^7
                save(fullfile(gridbarpathname,[prefix 'data_energy_filtered.mat']),'Ebar');
                clear Ebar
            end
        end
        time_Filter = toc(t_Filter);
        fprintf('\nelapsed time = %f s for all coarse grids',time_Filter);
        fprintf('\n');
    end
elseif g~=gref 
    %% Loading filtered data
    if g<2^7
        fprintf('\nLoading filtered data');
    else
        fprintf('\nLoading mean filtered data');
    end
    t_load = tic;
    load(fullfile(gridpathname,[prefix 'mean_data_filtered.mat']),'mYbar');
    if g<2^7
        load(fullfile(gridpathname,[prefix 'data_filtered.mat']),'Ybar');
    end
    if PostProcessingTau
        load(fullfile(gridpathname,[prefix 'mean_data_tau_filtered.mat']),'mTaubar');
        if g<2^7
            load(fullfile(gridpathname,[prefix 'data_tau_filtered.mat']),'Taubar');
        end
    end
    if PostProcessingPressure
        load(fullfile(gridpathname,[prefix 'mean_data_pressure_filtered.mat']),'mPresbar');
        if g<2^7
            load(fullfile(gridpathname,[prefix 'data_pressure_filtered.mat']),'Presbar');
        end
    end
    if PostProcessingEnergy
        load(fullfile(gridpathname,[prefix 'mean_data_energy_filtered.mat']),'mEbar');
        if g<2^7
            load(fullfile(gridpathname,[prefix 'data_energy_filtered.mat']),'Ebar');
        end
    end
    time_load = toc(t_load);
    fprintf('\nelapsed time = %f s',time_load);
    fprintf('\n');
        
    %% Computing error
    if computeError
        fprintf('\nComputing error between DNS and filtered data');
        t_Error = tic;
        errorY = zeros(3,p+1);
        normYbar = zeros(3,p+1);
        if PostProcessingTau
            errorTau = zeros(5,p+1);
            normTaubar = zeros(5,p+1);
        end
        if PostProcessingPressure
            errorPres = zeros(1,p+1);
            normPresbar = zeros(1,t+1);
        end
        if PostProcessingEnergy
            errorE = zeros(9,p+1);
            normEbar = zeros(9,p+1);
        end
        
        for t=0:p
            t_Errort = tic;
            components = {1:dim,dim+1,dim+2};
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
                        Yt = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
                    else
                        Yt = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
                    end
                    Yt = sSpace.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
            end
            Yt = reshape(Yt,[N,n,m]);
            if g<2^7
                Ybart = Ybar(:,:,:,t+1);
            else
                load(fullfile(gridpathname,[prefix 'data_filtered_t' num2str(t) '.mat']),'Ybart');
            end
            normYbar(:,t+1) = arrayfun(@(i) compute_norm(x,reshape(sum(Ybart(:,components{i},:).^2,2),[N,sx]),dim),1:3);
            errorYt = Ybart - Yt;
            clear Yt Ybart
            errorY(:,t+1) = arrayfun(@(i) compute_norm(x,reshape(sum(errorYt(:,components{i},:).^2,2),[N,sx]),dim),1:3);
            clear errorYt
            if PostProcessingTau
                components = {1:dim,dim+(1:dim),2*dim+(1:dim),3*dim+(1:dim),4*dim+1};
                if g<2^7
                    Taut = Tau(:,:,:,t+1);
                    Taubart = Taubar(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,[prefix 'data_tau_t' num2str(t) '.mat']),'Taut');
                    load(fullfile(gridpathname,[prefix 'data_tau_filtered_t' num2str(t) '.mat']),'Taubart');
                end
                normTaubar(:,t+1) = arrayfun(@(i) compute_norm(x,reshape(sum(Taubart(:,components{i},:).^2,2),[N,sx]),dim),1:5);
                errorTaut = Taubart - Taut;
                clear Taut Taubart
                errorTau(:,t+1) = arrayfun(@(i) compute_norm(x,reshape(sum(errorTaut(:,components{i},:).^2,2),[N,sx]),dim),1:5);
                clear errorTaut
            end
            if PostProcessingPressure
                if g<2^7
                    prest = pres(:,:,:,t+1);
                    Presbart = Presbar(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,[prefix 'data_pressure_t' num2str(t) '.mat']),'Prest');
                    load(fullfile(gridpathname,[prefix 'data_pressure_filtered_t' num2str(t) '.mat']),'Presbart');
                end
                normPresbar(t+1) = compute_norm(x,reshape(Presbart(:,1,:).^2,[N,sx]),dim);
                errorPrest = Presbart - prest;
                clear Prest Presbart
                errorPres(t+1) = compute_norm(x,reshape(errorPrest(:,1,:).^2,[N,sx]),dim);
                clear errorPrest
            end
            if PostProcessingEnergy
                if g<2^7
                    Et = E(:,:,:,t+1);
                    Ebart = Ebar(:,:,:,t+1);
                else
                    load(fullfile(gridpathname,[prefix 'data_energy_t' num2str(t) '.mat']),'Et');
                    load(fullfile(gridpathname,[prefix 'data_energy_filtered_t' num2str(t) '.mat']),'Ebart');
                end
                normEbar(:,t+1) = arrayfun(@(i) compute_norm(x,reshape(Ebart(:,i,:).^2,[N,sx]),dim),1:ne);
                errorEt = Ebart - Et;
                clear Et Ebart
                errorE(:,t+1) = arrayfun(@(i) compute_norm(x,reshape(errorEt(:,i,:).^2,[N,sx]),dim),1:ne);
                clear errorEt
            end
            time_Errort = toc(t_Errort);
            fprintf('\nTime %2d/%d : elapsed time = %f s',t,p,time_Errort);
        end
        
        time_Error = toc(t_Error);
        fprintf('\nelapsed time = %f s',time_Error);
        fprintf('\n');
        
        save(fullfile(gridpathname,[prefix 'error_filter.mat']),'errorY','normYbar');
        if PostProcessingTau
            save(fullfile(gridpathname,[prefix 'error_filter_tau.mat']),'errorTau','normTaubar');
        end
        if PostProcessingPressure
            save(fullfile(gridpathname,[prefix 'error_filter_pressure.mat']),'errorPres','normPresbar');
        end
        if PostProcessingEnergy
            save(fullfile(gridpathname,[prefix 'error_filter_energy.mat']),'errorE','normEbar');
        end
    else
        fprintf('\nLoading error between DNS and filtered data');
        t_load = tic;
        load(fullfile(gridpathname,[prefix 'error_filter.mat']),'errorY','normYbar');
        if PostProcessingTau
            load(fullfile(gridpathname,[prefix 'error_filter_tau.mat']),'errorTau','normTaubar');
        end
        if PostProcessingPressure
            load(fullfile(gridpathname,[prefix 'error_filter_pressure.mat']),'errorPres','normPresbar');
        end
        if PostProcessingEnergy
            load(fullfile(gridpathname,[prefix 'error_filter_energy.mat']),'errorE','normEbar');
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
            mysaveas(gridpathname,[prefix 'eigenvalues_CY_order'],formats,renderer);
            mymatlab2tikz(gridpathname,[prefix 'eigenvalues_CY_order.tex']);
            
            figure('Name','Evolution of errors')
            clf
            R = length(sv);
            semilogy(1:R,errsvdY.^2,'LineStyle','-','Color','b','LineWidth',1);
            grid on
            box on
            set(gca,'FontSize',10)
            xlabel('$R$','Interpreter',interpreter)
            ylabel('$\varepsilon_{Y}(R)$','Interpreter',interpreter)
            mysaveas(gridpathname,[prefix 'error_svdY'],formats,renderer);
            mymatlab2tikz(gridpathname,[prefix 'error_svdY.tex']);
            
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
            mysaveas(gridpathname,[prefix 'eigenvalues_CY_order'],formats,renderer);
            mymatlab2tikz(gridpathname,[prefix 'eigenvalues_CY_order.tex']);
            
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
            mysaveas(gridpathname,[prefix 'eigenvalues_CY_time'],formats,renderer);
            mymatlab2tikz(gridpathname,[prefix 'eigenvalues_CY_time.tex']);
            
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
            mysaveas(gridpathname,[prefix 'error_svdY'],formats,renderer);
            mymatlab2tikz(gridpathname,[prefix 'error_svdY.tex']);
            
            figure('Name','Evolution of eigenvalues')
            clf
            semilogy(1:Q,sw(:).^2,'LineStyle','-','Color','b','LineWidth',1);
            grid on
            box on
            set(gca,'FontSize',fontsize)
            xlabel('$\beta$','Interpreter',interpreter)
            ylabel('$\Lambda_{\beta}$','Interpreter',interpreter)
            mysaveas(gridpathname,[prefix 'eigenvalues_CZ'],formats,renderer);
            mymatlab2tikz(gridpathname,[prefix 'PCA_eigenvalues_CZ.tex']);
            
            figure('Name','Evolution of errors')
            clf
            semilogy(1:Q,errsvdZ(:).^2,'LineStyle','-','Color','b','LineWidth',1);
            grid on
            box on
            set(gca,'FontSize',fontsize)
            xlabel('$Q$','Interpreter',interpreter)
            ylabel('$\varepsilon_{Z}(Q)$','Interpreter',interpreter)
            mysaveas(gridpathname,[prefix 'error_svdZ'],formats,renderer);
            mymatlab2tikz(gridpathname,[prefix 'error_svdZ.tex']);
    end
end

%% Dsiplay statistics of stochastic processes
if displayStatistics
    ts = (0:p)*dt;
    errL2Y = trapz(ts,norm2Y,2)./trapz(ts,norm2mY,2);
    fprintf('\nL2-norm = %.3e for velocity U',errL2Y(1));
    fprintf('\nL2-norm = %.3e for phase C',errL2Y(3));
    fprintf('\nL2-norm = %.3e for pressure p',errL2Y(2));
    if PostProcessingPressure
        errL2pres = trapz(ts,norm2pres,2)./trapz(ts,norm2mpres,2);
        fprintf('\nL2-norm = %.3e for pressure post-processed p',errL2pres);
    end
    if PostProcessingTau
        errL2Tau = trapz(ts,norm2Tau,2)./trapz(ts,norm2mTau,2);
        fprintf('\nL2-norm = %.3e for tauTime',errL2Tau(1));
        fprintf('\nL2-norm = %.3e for div(tauConv)',errL2Tau(2));
        fprintf('\nL2-norm = %.3e for div(tauDiff)',errL2Tau(3));
        fprintf('\nL2-norm = %.3e for tauSurf',errL2Tau(4));
        fprintf('\nL2-norm = %.3e for tauInterf',errL2Tau(5));
    end
    if PostProcessingEnergy
        errL2E = trapz(ts,norm2E,2)./trapz(ts,norm2mE,2);
        fprintf('\nL2-norm = %.3e for energyKinTime',errL2E(1));
        fprintf('\nL2-norm = %.3e for energyConv',errL2E(2));
        fprintf('\nL2-norm = %.3e for energyGrav',errL2E(3));
        fprintf('\nL2-norm = %.3e for energyPres',errL2E(4));
        fprintf('\nL2-norm = %.3e for energyPresDil',errL2E(5));
        fprintf('\nL2-norm = %.3e for energyKinGrad',errL2E(6));
        fprintf('\nL2-norm = %.3e for energyDiff',errL2E(7));
        fprintf('\nL2-norm = %.3e for energyVisc',errL2E(8));
        fprintf('\nL2-norm = %.3e for energySurf',errL2E(9));
    end
    fprintf('\n');
    
    figure('Name','Dispersion fot DNS data')
    clf
    hdl(1) = plot(ts(2:end),norm2Y(2,2:end)./norm2mY(2,2:end),'LineStyle','-','Color','r','LineWidth',1);
    hold on
    hdl(2) = plot(ts(2:end),norm2Y(1,2:end)./norm2mY(1,2:end),'LineStyle','-','Color','b','LineWidth',1);
    hdl(3) = plot(ts(2:end),norm2Y(3,2:end)./norm2mY(3,2:end),'LineStyle','-','Color','g','LineWidth',1);
    if PostProcessingPressure
        hdl(4) = plot(ts(2:end),norm2pres(:,2:end)./norm2mpres(:,2:end),'LineStyle','--','Color','g','LineWidth',1);
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$t$ [s]','Interpreter',interpreter)
    ylabel('Dispersion','Interpreter',interpreter)
    leg = {'$\chi$','$||u||$','$p$'};
    if PostProcessingPressure
        leg = [leg,'$p_{\mathrm{post}}$'];
    end
    l = legend(leg{:},'Location','NorthEast');
    set(l,'Interpreter',interpreter)
    mysaveas(gridpathname,[prefix 'dispersion'],formats,renderer);
    mymatlab2tikz(gridpathname,[prefix 'dispersion.tex']);
    
    if PostProcessingTau
        figure('Name','Dispersion for tau data')
        clf
        hdl(1) = plot(ts(2:end),norm2Y(2,2:end)./norm2mY(2,2:end),'LineStyle','-','Color','r','LineWidth',1);
        hold on
        hdl(2) = plot(ts(2:end),norm2Tau(1,2:end)./norm2mTau(1,2:end),'LineStyle','-','Color','b','LineWidth',1);
        hdl(3) = plot(ts(2:end),norm2Tau(2,2:end)./norm2mTau(2,2:end),'LineStyle','-','Color','g','LineWidth',1);
        hdl(4) = plot(ts(2:end),norm2Tau(3,2:end)./norm2mTau(3,2:end),'LineStyle','-','Color','m','LineWidth',1);
        hdl(5) = plot(ts(2:end),norm2Tau(4,2:end)./norm2mTau(4,2:end),'LineStyle','-','Color','c','LineWidth',1);
        hdl(6) = plot(ts(2:end),norm2Tau(5,2:end)./norm2mTau(5,2:end),'LineStyle','-','Color',[1 0.5 0],'LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$t$ [s]','Interpreter',interpreter)
        ylabel('Dispersion','Interpreter',interpreter)
        leg = {'$\chi$','$||\tau_{\mathrm{time}}||$','$||\nabla \cdot \tau_{\mathrm{conv}}||$','$||\nabla \cdot \tau_{\mathrm{diff}}||$','$||\tau_{\mathrm{surf}}||$','$\tau_{\mathrm{interf}}$'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter',interpreter)
        mysaveas(gridpathname,[prefix 'dispersion_tau'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'dispersion_tau.tex']);
    end
    
    if PostProcessingEnergy
        figure('Name','Dispersion for energy data')
        hdl(1) = plot(ts(2:end),norm2Y(2,2:end)./norm2mY(2,2:end),'LineStyle','-','Color','r','LineWidth',1);
        hold on
        hdl(2) = plot(ts(2:end),norm2E(1,2:end)./norm2mE(1,2:end),'LineStyle','-','Color','b','LineWidth',1);
        hdl(3) = plot(ts(2:end),norm2E(2,2:end)./norm2mE(2,2:end),'LineStyle','-','Color','g','LineWidth',1);
        hdl(4) = plot(ts(2:end),norm2E(3,2:end)./norm2mE(3,2:end),'LineStyle','-','Color','m','LineWidth',1);
        hdl(5) = plot(ts(2:end),norm2E(4,2:end)./norm2mE(4,2:end),'LineStyle','-','Color','c','LineWidth',1);
        hdl(6) = plot(ts(2:end),norm2E(5,2:end)./norm2mE(5,2:end),'LineStyle','-','Color',[1 0.5 0],'LineWidth',1);
        hdl(7) = plot(ts(2:end),norm2E(6,2:end)./norm2mE(6,2:end),'LineStyle','-','Color','k','LineWidth',1);
        hdl(8) = plot(ts(2:end),norm2E(7,2:end)./norm2mE(7,2:end),'LineStyle','-','Color','y','LineWidth',1);
        hdl(9) = plot(ts(2:end),norm2E(8,2:end)./norm2mE(8,2:end),'LineStyle','-','Color',([1 0 1]+[1 0 0])/2,'LineWidth',1);
        hdl(10) = plot(ts(2:end),norm2E(9,2:end)./norm2mE(9,2:end),'LineStyle','-','Color',([1 0 1]+[0 0 1])/2,'LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$t$ [s]','Interpreter',interpreter)
        ylabel('Dispersion','Interpreter',interpreter)
        leg = {'$\chi$','$e_{\mathrm{kin time}}$','$e_{\mathrm{conv}}$','$e_{\mathrm{grav}}$','$e_{\mathrm{pres}}$','$e_{\mathrm{pres-dil transfer}}$','$e_{\mathrm{grad kin}}$','$e_{\mathrm{diff}}$','$e_{\mathrm{visc}}$','$e_{\mathrm{surf}}$'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter',interpreter)
        mysaveas(gridpathname,[prefix 'dispersion_energy'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'dispersion_energy.tex']);
    end
end  

%% Display quantities of interest
if displayQoI && QoI
    % Mean
    mQu = reshape(mQu(1,1:dim,:,:),[dim,2,p+1]);
    if PostProcessingTau
        % Mean
        mQtauTime = reshape(mQtau(1,1:dim,:,:),[dim,2,p+1]);
        mQdivtauConv = reshape(mQtau(1,dim+(1:dim),:,:),[dim,2,p+1]);
        mQdivtauDiff = reshape(mQtau(1,2*dim+(1:dim),:,:),[dim,2,p+1]);
        mQtauSurf = reshape(mQtau(1,3*dim+(1:dim),:,:),[dim,2,p+1]);
        mQtauInterf = reshape(mQtau(1,4*dim+1,:,:),[1,2,p+1]);
        % Correlation
        IQtauTime = IQtau(1:dim,:,1:dim,:,:);
        IQdivtauConv = IQtau(dim+(1:dim),:,dim+(1:dim),:,:);
        IQdivtauDiff = IQtau(2*dim+(1:dim),:,2*dim+(1:dim),:,:);
        IQtauSurf = IQtau(3*dim+(1:dim),:,3*dim+(1:dim),:,:);
        IQtauInterf = IQtau(4*dim+1,:,4*dim+1,:,:);
    end
    if PostProcessingPressure
        mQpres = reshape(mQpres(1,:,:,:),[2,2,p+1]);
    end
    if PostProcessingEnergy
        % Mean
        mQenergyKinTime = reshape(mQe(1,1,:,:),[1,2,p+1]);
        mQenergyConv = reshape(mQe(1,2,:,:),[1,2,p+1]);
        mQenergyGrav = reshape(mQe(1,3,:,:),[1,2,p+1]);
        mQenergyPres = reshape(mQe(1,4,:,:),[1,2,p+1]);
        mQenergyPresDil = reshape(mQe(1,5,:,:),[1,2,p+1]);
        mQenergyKinGrad = reshape(mQe(1,6,:,:),[1,2,p+1]);
        mQenergyDiff = reshape(mQe(1,7,:,:),[1,2,p+1]);
        mQenergyVisc = reshape(mQe(1,8,:,:),[1,2,p+1]);
        mQenergySurf = reshape(mQe(1,9,:,:),[1,2,p+1]);
        % Correlation
        IQenergyKinTime = IQe(1,:,1,:,:);
        IQenergyConv = IQe(2,:,2,:,:);
        IQenergyGrav = IQe(3,:,3,:,:);
        IQenergyPres = IQe(4,:,4,:,:);
        IQenergyPresDil = IQe(5,:,5,:,:);
        IQenergyKinGrad = IQe(6,:,6,:,:);
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
    mysaveas(gridpathname,[prefix 'mean_u'],formats,renderer);
    mymatlab2tikz(gridpathname,[prefix 'mean_u.tex']);
    
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
    mysaveas(gridpathname,[prefix 'mean_p'],formats,renderer);
    mymatlab2tikz(gridpathname,[prefix 'mean_p.tex']);
    
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
        mysaveas(gridpathname,[prefix 'mean_tauTime'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'mean_tauTime.tex']);
        
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
        mysaveas(gridpathname,[prefix 'mean_divtauConv'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'mean_divtauConv.tex']);
        
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
        mysaveas(gridpathname,[prefix 'mean_divtauDiff'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'mean_divtauDiff.tex']);
        
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
        mysaveas(gridpathname,[prefix 'mean_tauSurf'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'mean_tauSurf.tex']);
        
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
        mysaveas(gridpathname,[prefix 'mean_tauInterf'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'mean_tauInterf.tex']);
    end
    
    if PostProcessingPressure
        figure('Name','Mean of spatial average of pressure p')
        clf
        hdl(1) = plot(ts,squeeze(mQpres(2,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(mQpres(2,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$t$ [s]','Interpreter',interpreter)
        ylabel('$p$ postprocessed','Interpreter',interpreter)
        leg = {'phase 1','phase 2'};
        legend(leg{:},'Location','NorthEast')
        mysaveas(gridpathname,[prefix 'mean_p_post'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'mean_p_post.tex']);
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
    mysaveas(gridpathname,[prefix 'power_u'],formats,renderer);
    mymatlab2tikz(gridpathname,[prefix 'power_u.tex']);
    
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
    mysaveas(gridpathname,[prefix 'power_p'],formats,renderer);
    mymatlab2tikz(gridpathname,[prefix 'power_p.tex']);
    
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
        mysaveas(gridpathname,[prefix 'power_tauTime'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'power_tauTime.tex']);
        
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
        mysaveas(gridpathname,[prefix 'power_divtauConv'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'power_divtauConv.tex']);
        
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
        mysaveas(gridpathname,[prefix 'power_divtauDiff'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'power_divtauDiff.tex']);
        
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
        mysaveas(gridpathname,[prefix 'power_tauSurf'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'power_tauSurf.tex']);
        
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
        mysaveas(gridpathname,[prefix 'power_tauInterf'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'power_tauInterf.tex']);
    end
    
    if PostProcessingPressure
        figure('Name','Power of pressure')
        clf
        hdl(1) = plot(ts,squeeze(IQpres(2,1,2,1,:)),'LineStyle','-','Color','b','LineWidth',1);
        hold on
        hdl(2) = plot(ts,squeeze(IQpres(2,2,2,2,:)),'LineStyle','--','Color','b','LineWidth',1);
        hold off
        grid on
        box on
        set(gca,'FontSize',10)
        xlabel('$\tau$ [s]','Interpreter',interpreter)
        ylabel(['$\displaystyle\frac{1}{T-\tau} \int_0^{T-\tau} {[R_Q(t+\tau,t)]}_{jj} \, dt$'],'Interpreter',interpreter)
        leg = {'$p$ postprocessed in phase 1','$p$ postprocessed in phase 2'};
        l = legend(leg{:},'Location','NorthEast');
        set(l,'Interpreter',interpreter)
        mysaveas(gridpathname,[prefix 'power_p_post'],formats,renderer);
        mymatlab2tikz(gridpathname,[prefix 'power_p_post.tex']);
    end
end

%% Display error between data and filtered data
if g~=gref && displayError && Filtering
    ts = (0:p)*dt;
    orange = [1 0.5 0];
    figure('Name','Error between data and filtered data')
    clf
    hdl(1) = plot(ts,errorY(1,:)./normYbar(1,:),'LineStyle','-','Color','b','LineWidth',1);
    hold on
    hdl(2) = plot(ts,errorY(2,:)./normYbar(2,:),'LineStyle','-','Color',orange,'LineWidth',1);
    hdl(3) = plot(ts,errorY(3,:)./normYbar(3,:),'LineStyle','-','Color','r','LineWidth',1);
    cpt = 3;
    if PostProcessingTau
        hdl(cpt+1) = plot(ts,errorTau(1,:)./normTaubar(1,:),'LineStyle','-','Color','g','LineWidth',1);
        hdl(cpt+2) = plot(ts,errorTau(2,:)./normTaubar(2,:),'LineStyle','-','Color','m','LineWidth',1);
        hdl(cpt+3) = plot(ts,errorTau(3,:)./normTaubar(3,:),'LineStyle','-','Color','c','LineWidth',1);
        hdl(cpt+4) = plot(ts,errorTau(4,:)./normTaubar(4,:),'LineStyle','-','Color','y','LineWidth',1);
        hdl(cpt+5) = plot(ts,errorTau(5,:)./normTaubar(5,:),'LineStyle','-','Color','k','LineWidth',1);
        cpt = cpt+5;
    end
    if PostProcessingPressure
        hdl(cpt+1) = plot(ts,errorPres./normPresbar,'LineStyle','--','Color',orange,'LineWidth',1);
        cpt = cpt+1;
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$t$ [s]','Interpreter',interpreter)
    ylabel('Normalized error','Interpreter',interpreter)
    leg = {'$||u||$','$p$','$\chi$'};
    if PostProcessingTau
        leg = [leg,'$||\tau_{\mathrm{time}}||$','$||\nabla \cdot \tau_{\mathrm{conv}}||$','$||\nabla \cdot \tau_{\mathrm{diff}}||$','$||\tau_{\mathrm{surf}}||$','$\tau_{\mathrm{interf}}$'];
    end
    if PostProcessingPressure
        leg = [leg,'$p_{\mathrm{post}}$'];
    end
    l = legend(leg{:},'Location','NorthEast');
    set(l,'Interpreter',interpreter)
    mysaveas(gridpathname,[prefix 'error_filter'],formats,renderer);
    mymatlab2tikz(gridpathname,[prefix 'error_filter.tex']);
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
switch usePCA
    case 'no'
        mYr = mY;
    otherwise
        mYr = s.unscaling(mY(:,:),Ya,Yb);
        mYr = reshape(mYr,[1,n,m,p+1]);
end
mU = reshape(mYr(1,1:dim,:,:),[dim*m,p+1]);
mp = reshape(mYr(1,dim+1,:,:),[m,p+1]);
mC = reshape(mYr(1,dim+2,:,:),[m,p+1]);
if PostProcessingTau
    mtauTime = reshape(mTau(1,1:dim,:,:),[dim*m,p+1]);
    mdivtauConv = reshape(mTau(1,dim+(1:dim),:,:),[dim*m,p+1]);
    mdivtauDiff = reshape(mTau(1,2*dim+(1:dim),:,:),[dim*m,p+1]);
    mtauSurf = reshape(mTau(1,3*dim+(1:dim),:,:),[dim*m,p+1]);
    mtauInterf = reshape(mTau(1,4*dim+1,:,:),[m,p+1]);
end
if PostProcessingPressure
    mpres = reshape(mpres(1,1,:,:),[m,p+1]);
end

if g~=gref && Filtering
    % Mean filtered solution
    mUbar = reshape(mYbar(1,1:dim,:,:),[dim*m,p+1]);
    mpbar = reshape(mYbar(1,dim+1,:,:),[m,p+1]);
    mCbar = reshape(mYbar(1,dim+2,:,:),[m,p+1]);
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
    mpt = mp(:,t+1);
    mCt = mC(:,t+1);
    fields = {mUt,mpt,mCt};
    fieldnames = {'velocity','pressure','phase'};
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
        fieldnames = [fieldnames,'pressure postprocessed'];
    end
%     if PostProcessingEnergy
%         menergyKinTimet = menergyKinTime(:,t+1);
%         menergyConvt = menergyConv(:,t+1);
%         menergyGravt = menergyGrav(:,t+1);
%         menergyPrest = menergyPres(:,t+1);
%         menergyPresDilt = menergyPresDil(:,t+1);
%         menergyKinGradt = menergyKinGrad(:,t+1);
%         menergyDifft = menergyDiff(:,t+1);
%         menergyVisct = menergyVisc(:,t+1);
%         menergySurft = menergySurf(:,t+1);
%         fields = [fields,menergyKinTimet,menergyConvt,menergyGravt,...
%             menergyPrest,menergyPresDilt,menergyKinGradt,...
%             menergyDifft,menergyVisct,menergySurft];
%         fieldnames = [fieldnames,'kinetic energy','convection energy','gravity energy',...
%             'power of external pressure forces','pressure-dilatation energy transfer','transport of gradient of kinetic energy',...
%             'energy exchange with kinetic energy','power of external viscous stresses','capillary kinetic energy'];
%     end
    
    if g~=gref && Filtering
        mUbart = mUbar(:,t+1);
        mpbart = mpbar(:,t+1);
        mCbart = mCbar(:,t+1);
        fields = [fields,mUbart,mpbart,mCbart];
        fieldnames = [fieldnames,'velocity filtered','pressure filtered','phase filtered'];
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
            fieldnames = [fieldnames,'pressure postprocessed filtered'];
        end
%         if PostProcessingEnergy
%             menergyKinTimebart = menergyKinTimebar(:,t+1);
%             menergyConvbart = menergyConvbar(:,t+1);
%             menergyGravbart = menergyGravbar(:,t+1);
%             menergyPresbart = menergyPresbar(:,t+1);
%             menergyPresDilbart = menergyPresDilbar(:,t+1);
%             menergyKinGradbart = menergyKinGradbar(:,t+1);
%             menergyDiffbart = menergyDiffbar(:,t+1);
%             menergyViscbart = menergyViscbar(:,t+1);
%             menergySurfbart = menergySurfbar(:,t+1);
%             fields = [fields,menergyKinTimebart,menergyConvbart,menergyGravbart,...
%                 menergyPresbart,menergyPresDilbart,menergyKinGradbart,...
%                 menergyDiffbart,menergyViscbart,menergySurfbart];
%             fieldnames = [fieldnames,'kinetic energy filtered','convection energy filtered','gravity energy filtered',...
%                 'power of external pressure forces filtered','pressure-dilatation energy transfer filtered','transport of gradient of kinetic energy filtered',...
%                 'energy exchange with kinetic energy filtered','power of external viscous stresses filtered','capillary kinetic energy filtered'];
%         end
    end
    write_vtk_mesh(M,fields,[],fieldnames,[],gridpathname,[prefix 'diphasic_fluids_grid' num2str(g) '_mean'],1,t);
end
make_pvd_file(gridpathname,[prefix 'diphasic_fluids_grid' num2str(g) '_mean'],1,p+1);
time_save = toc(t_save);
fprintf('\nelapsed time = %f s',time_save);
fprintf('\n');
end

%% Standard deviation of random solution
stdU = zeros(dim*m,p+1);
stdp = zeros(m,p+1);
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

if g~=gref && Filtering
    % Standard deviation of filtered solution
    stdUbar = zeros(dim*m,p+1);
    stdpbar = zeros(m,p+1);
    stdCbar = zeros(m,p+1);
    if PostProcessingTau
        stdtauTimebar = zeros(dim*m,p+1);
        stddivtauConvbar = zeros(dim*m,p+1);
        stddivtauDiffbar = zeros(dim*m,p+1);
        stdtauSurfbar = zeros(dim*m,p+1);
        stdtauInterfbar = zeros(m,p+1);
    end
    if PostProcessingPressure
        stdpresbar = zeros(m,p+1);
    end
end

for t=0:p
    switch usePCA
        case 'no'
            if g<2^7
                Yt = Y(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
            end
            stdYt = std(Yt);
            % ut = Yt(:,1:dim,:);
            % pt = Yt(:,dim+1,:);
            % Ct = Yt(:,dim+2,:);
            % clear Yt
            % stdU(:,t+1) = std(ut(:,:))';
            % stdp(:,t+1) = std(pt(:,:))';
            % stdC(:,t+1) = std(Ct(:,:))';
            % clear ut pt Ct
            
        case 'single'
            %[Yt,mYt,Vt] = s.reconstructionAtStep(mY,V,sv,X,t+1);
            %CYt = cov(Yt');
            %CYt = s.cov(Vt,sv);
            %stdYt = sqrt(diag(CYt));
            % Yt = s.reconstructionAtStep(mY,V,sv,X,t+1);
            % stdYt = std(Yt,0,2);
            % stdYt = s.unscaling(stdYt,Ya(:,t+1),zeros(size(Yb(:,t+1))))';
            Vt = V(:,:,t+1);
            stdYt = s.std(Vt,sv);
            stdYt = s.unscaling(stdYt,Ya(:,t+1),zeros(size(Yb(:,t+1))))';
            
        case 'double'
            % if g<2^7
            %     [Yt,mYt,Vt,svt,mZt,Wt,Rt] = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
            % else
            %     [Yt,mYt,Vt,svt,mZt,Wt,Rt] = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
            % end
            %CYt = cov(Yt');
            %stdYt = sqrt(diag(CYt));
            % stdYt = std(Yt,0,2);
            % stdYt = sSpace.unscaling(stdYt,Ya(:,t+1),zeros(size(Yb(:,t+1))))';
            if g<2^7
                [Vt,svt,Wt,Rt] = sSpace.getPrincipalComponentsDoubleAtStep(V,sv,W,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
            else
                [Vt,svt,Wt,Rt] = sSpace.getPrincipalComponentsDoubleAtStep([],sv,W,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
            end
            stdYt = sSpace.stdDouble(Vt,svt,Wt,sw);
            stdYt = sSpace.unscaling(stdYt,Ya(:,t+1),zeros(size(Yb(:,t+1))))';
    end
    stdYt = reshape(stdYt,[n,m]);
    stdU(:,t+1) = reshape(stdYt(1:dim,:),[dim*m,1]);
    stdp(:,t+1) = reshape(stdYt(dim+1,:),[m,1]);
    stdC(:,t+1) = reshape(stdYt(dim+2,:),[m,1]);
    clear stdYt
    if PostProcessingTau
        if g<2^7
            Taut = Tau(:,:,:,t+1);
        else
            load(fullfile(gridpathname,[prefix 'data_tau_t' num2str(t) '.mat']),'Taut');
        end
        stdTaut = reshape(std(Taut),[ntau,m]);
        clear Taut
        stdtauTime(:,t+1) = reshape(stdTaut(1:dim,:),[dim*m,1]);
        stddivtauConv(:,t+1) = reshape(stdTaut(dim+(1:dim),:),[dim*m,1]);
        stddivtauDiff(:,t+1) = reshape(stdTaut(2*dim+(1:dim),:),[dim*m,1]);
        stdtauSurf(:,t+1) = reshape(stdTaut(3*dim+(1:dim),:),[dim*m,1]);
        stdtauInterf(:,t+1) = reshape(stdTaut(4*dim+1,:),[m,1]);
        clear stdTaut
    end
    if PostProcessingPressure
        if g<2^7
            prest = pres(:,:,:,t+1);
        else
            load(fullfile(gridpathname,[prefix 'data_pressure_t' num2str(t) '.mat']),'prest');
        end
        stdpres(:,t+1) = std(prest(:,:))';
        clear prest
    end
    
    if g~=gref && Filtering
        if g<2^7
            Ybart = Ybar(:,:,:,t+1);
        else
            load(fullfile(gridpathname,[prefix 'data_filtered_t' num2str(t) '.mat']),'Ybart');
        end
        stdYbart = reshape(std(Ybart),[n,m]);
        stdUbar(:,t+1) = reshape(stdYbart(1:dim,:),[dim*m,1]);
        stdpbar(:,t+1) = reshape(stdYbart(dim+1,:),[m,1]);
        stdCbar(:,t+1) = reshape(stdYbart(dim+2,:),[m,1]);
        clear stdYbart
        if PostProcessingTau
            if g<2^7
                Taubart = Taubar(:,:,:,t+1);
            else
                load(fullfile(gridpathname,[prefix 'data_tau_filtered_t' num2str(t) '.mat']),'Taubart');
            end
            stdTaubart = reshape(std(Taubart),[ntau,m]);
            stdtauTimebar(:,t+1) = reshape(stdTaubart(1:dim,:),[dim*m,1]);
            stddivtauConvbar(:,t+1) = reshape(stdTaubart(dim+(1:dim),:),[dim*m,1]);
            stddivtauDiffbar(:,t+1) = reshape(stdTaubart(2*dim+(1:dim),:),[dim*m,1]);
            stdtauSurfbar(:,t+1) = reshape(stdTaubart(3*dim+(1:dim),:),[dim*m,1]);
            stdtauInterfbar(:,t+1) = reshape(stdTaubart(4*dim+1,:),[m,1]);
            clear stdTaubart
        end
        if PostProcessingPressure
            if g<2^7
                Presbart = Presbar(:,:,:,t+1);
            else
                load(fullfile(gridpathname,[prefix 'data_pressure_filtered_t' num2str(t) '.mat']),'Presbart');
            end
            stdpresbar(:,t+1) = std(Presbart(:,:))';
            clear Presbart
        end
    end
end

if saveStd
fprintf('\nSaving standard deviation of random solution');
t_save = tic;
for t=0:p
    stdUt = stdU(:,t+1);
    stdpt = stdp(:,t+1);
    stdCt = stdC(:,t+1);
    fields = {stdUt,stdpt,stdCt};
    fieldnames = {'velocity','pressure','phase'};
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
        fieldnames = [fieldnames,'pressure postprocessed'];
    end
    
    if g~=gref && Filtering
        stdUbart = stdUbar(:,t+1);
        stdpbart = stdpbar(:,t+1);
        stdCbart = stdCbar(:,t+1);
        fields = [fields,stdUbart,stdpbart,stdCbart];
        fieldnames = [fieldnames,'velocity filtered','pressure filtered','phase filtered'];
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
            fieldnames = [fieldnames,'pressure postprocessed filtered'];
        end
    end
    write_vtk_mesh(M,fields,[],fieldnames,[],gridpathname,[prefix 'diphasic_fluids_grid' num2str(g) '_std'],1,t);
end
make_pvd_file(gridpathname,[prefix 'diphasic_fluids_grid' num2str(g) '_std'],1,p+1);
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
                Yt = sSpace.reconstructionDoubleAtStep(mY,V,sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
            else
                Yt = sSpace.reconstructionDoubleAtStep(mY,[],sv,mZ,W,sw,X,R,t+1,'pathname',gridpathname,'filename',[prefix 'space_data_t' num2str(t) '.mat']);
            end
            Yt = sSpace.unscaling(Yt,Ya(:,t+1),Yb(:,t+1))';
    end
    Yt = reshape(Yt,[N,n,m]);
    Ut = reshape(Yt(:,1:dim,:),[N,dim*m])';
    pt = reshape(Yt(:,dim+1,:),[N,m])';
    Ct = reshape(Yt(:,dim+2,:),[N,m])';
    clear Yt
    if PostProcessingTau
        Taut = Tau(:,:,:,t+1);
        tauTimet = reshape(Taut(:,1:dim,:),[N,dim*m])';
        divtauConvt = reshape(Taut(:,dim+(1:dim),:),[N,dim*m])';
        divtauDifft = reshape(Taut(:,2*dim+(1:dim),:),[N,dim*m])';
        tauSurft = reshape(Taut(:,3*dim+(1:dim),:),[N,dim*m])';
        tauInterft = reshape(Taut(:,4*dim+1,:),[N,m])';
        clear Taut
    end
    if PostProcessingPressure
        prest = reshape(pres(:,1,:,t+1),[N,m])';
    end
    
    if g~=gref && Filtering
        if g<2^7
            Ybart = Ybar(:,:,:,t+1);
        else
            load(fullfile(gridpathname,[prefix 'data_filtered_t' num2str(t) '.mat']),'Ybart');
        end
        Ubart = reshape(Ybart(:,1:dim,:),[N,dim*m])';
        pbart = reshape(Ybart(:,dim+1,:),[N,m])';
        Cbart = reshape(Ybart(:,dim+2,:),[N,m])';
        clear Ybart
        if PostProcessingTau
            Taubart = Taubar(:,:,:,t+1);
            tauTimebart = reshape(Taubart(:,1:dim,:),[N,dim*m])';
            divtauConvbart = reshape(Taubart(:,dim+(1:dim),:),[N,dim*m])';
            divtauDiffbart = reshape(Taubart(:,2*dim+(1:dim),:),[N,dim*m])';
            tauSurfbart = reshape(Taubart(:,3*dim+(1:dim),:),[N,dim*m])';
            tauInterfbart = reshape(Taubart(:,4*dim+1,:),[N,m])';
            clear Taubart
        end
        if PostProcessingPressure
            presbart = reshape(Presbar(:,1,:,t+1),[N,m])';
        end
    end
    
    for l=1:N
        Ult = Ut(:,l);
        plt = pt(:,l);
        Clt = Ct(:,l);
        fields = {Ult,plt,Clt};
        fieldnames = {'velocity','pressure','phase'};
        if PostProcessingTau
            tauTimelt = tauTimet(:,l);
            divtauConvlt = divtauConvt(:,l);
            divtauDifflt = divtauDifft(:,l);
            tauSurflt = tauSurft(:,l);
            tauInterflt = tauInterft(:,l);
            fields = [fields,tauTimelt,divtauConvlt,divtauDifflt,tauSurflt,tauInterflt];
            fieldnames = [fieldnames,'tauTime','div(tauConv)','div(tauDiff)','tauSurf','tauInterf'];
        end
        if PostProcessingPressure
            preslt = prest(:,l);
            fields = [fields,preslt];
            fieldnames = [fieldnames,'pressure postprocessed'];
        end
        
        if g~=gref && Filtering
            Ubarlt = Ubart(:,l);
            pbarlt = pbart(:,l);
            Cbarlt = Cbart(:,l);
            fields = [fields,Ubarlt,Cbarlt];
            fieldnames = [fieldnames,'velocity filtered','pressure filtered','phase filtered'];
            if PostProcessingTau
                tauTimebarlt = tauTimebart(:,l);
                divtauConvbarlt = divtauConvbart(:,l);
                divtauDiffbarlt = divtauDiffbart(:,l);
                tauSurfbarlt = tauSurfbart(:,l);
                tauInterfbarlt = tauInterfbart(:,l);
                fields = [fields,tauTimebarlt,divtauConvbarlt,divtauDiffbarlt,tauSurfbarlt,tauInterfbarlt];
                fieldnames = [fieldnames,'tauTime filtered','div(tauConv) filtered','div(tauDiff) filtered','tauSurf filtered','tauInterf filtered'];
            end
            if PostProcessingPressure
                presbarlt = presbart(:,l);
                fields = [fields,presbarlt];
                fieldnames = [fieldnames,'pressure postprocessed filtered'];
            end
        end
        write_vtk_mesh(M,fields,[],fieldnames,[],gridpathname,[prefix 'diphasic_fluids_grid' num2str(g) '_sample_' num2str(l)],1,t);
    end
    time_savet = toc(t_savet);
    fprintf('\nTime step %2d/%d : elapsed time = %f s',t,p,time_savet);
end
for l=1:N
    make_pvd_file(gridpathname,[prefix 'diphasic_fluids_grid' num2str(g) '_sample_' num2str(l)],1,p+1);
end
time_save = toc(t_save);
fprintf('\nelapsed time = %f s',time_save);
fprintf('\n');
end

time_Total = toc(t_Total);
fprintf('\nElapsed time = %f s\n',time_Total);

myparallel('stop');
