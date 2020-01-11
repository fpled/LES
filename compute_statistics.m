clc
clearvars
close all

solveProblem = true;
displaySolution = false;
displayEigenvales = true;
displayCovariance = false;

% index = 'time';
index = 'coord';

cmap = 'default';
% cmap = 'gray';
framerate = 5;
fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc'};
renderer = 'OpenGL';

pathname = fileparts(mfilename('fullpath'));

L = 1; % domain size (m)
sigma = 0.45; % surface tension (N/m)
mu = [0.1,0.1]; % dynamic viscosity of [water,oil] (Pa.s)
rho = [1000,900]; % mass density of [water,oil] (kg/m3)
gravity = 9.81; % gravity (m/s2)

% Time scheme
dt = 5e-3; % time step (s)
dt = 100*dt; % physical time step stored every 100 time iterations

order = [3 1 2]; % dimension order
perm = @(u) permute(u,[2,order+2,1]);
iperm = @(u) ipermute(u,[2,order+2,1]);

tolsvdYc = eps; % relative precision for truncated SVD of Yc
tolsvdZc = 1e-6; % relative precision for truncated SVD of Zc

% for g=2.^(4:8)
for g=2^4
    tic
    gridname = ['Grid' num2str(g)];
    disp(gridname)
    gridpathname = fullfile(pathname,gridname);
    load(fullfile(gridpathname,'data.mat'),'N','n','m','p');
    r = n*m;
    
    fprintf('\nn = %d variables',n);
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
    
    % Explicit backward Euler time scheme (first-order accurate, conditionally stable)
    % Dt = spdiags(repmat([1 -1],p+1,1),[0 -1],p+1,p+1)/dt;
    % Implicit forward Euler time scheme (first-order accurate, unconditionally stable)
    % Dt = spdiags(repmat([1 -1],p+1,1),[1 0],p+1,p+1)/dt;
    % Implicit central-difference time scheme (second-order accurate, unconditionally stable)
    Dt = spdiags(repmat([1 -1],p+1,1),[1 -1],p+1,p+1)/(2*dt);
    
    if solveProblem
        %% First reduction step
        fprintf('\nPCA in space');
        Rinit = min(r,N);
        if g<2^7
            load(fullfile(gridpathname,'data.mat'),'Y');
            mY = mean(Y,1);
            Yc = Y - repmat(mY,N,1,1,1); % Yc = Y - mY.*ones(N,1,1,1);
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
        Rmax = 1;
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
            Rmax = max(Rmax,Rt);
            fprintf('\nTime %2d, t = %4g s : rank R = %d, error = %.3e for Y',t,t*dt,Rt,errYct);
            
            norm2Yct = sum(var(Yct,0,2));
            % norm2Yct_approx = sum(var(Yct_approx,0,2));
            err2Yct = errYct^2*norm2Yct;
            % err2Yct = norm2Yct-norm2Yct_approx;
            % sigtf = svdtruncate(Yct,eps);
            % sigtf = sigtf/sqrt(N-1);
            % err2Yct = sum(sigtf.^2)-sum(sigt.^2);
            
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
            
            %if verLessThan('matlab','9.1') % compatibility (<R2016b)
            %    CYt_approx = Vt*Sigt.^2*Vt';
            %else
            %    CYt_approx = Vt*(sigt.^2.*Vt');
            %end 
            %CYt = cov(Yct'); % CYt = 1/(N-1)*Yct*Yct';
            %errCYt = norm(CYt_approx-CYt)/norm(CYt);
            %fprintf('\n                                                  error = %.3e for CY',errCYt);
            
            % mZt = mean(Zt,1)';
            % CZt = cov(Zt); % CZt = 1/(N-1)*Zt'*Zt;
            % norm(mZt)
            % norm(CZt-eye(Rt))
            % norm(Vt'*Vt-eye(Rt))
        end
        fprintf('\n');
        
        t = (0:p)*dt;
        errL2 = trapz(t,err2Yc,2)/trapz(t,norm2Yc,2);
        fprintf('\nerror = %.3e for Y',errL2);
        fprintf('\n');
        
        sig = sig(1:Rmax,:);
        Z = Z(:,1:Rmax,:);
        errsvdYc = errsvdYc(1:Rmax,:);
        save(fullfile(gridpathname,'PCA_space.mat'),'mY','sig','Z','errsvdYc','R','Rmax');
        if g<2^7
            V = V(:,1:Rmax,:);
            save(fullfile(gridpathname,'PCA_space.mat'),'V','-append');
        end
        
        %% Second reduction step for each coordinate
%         fprintf('\nPCA in time');
%         Q = min(p+1,N);
%         s = zeros(Q,Rmax);
%         W = zeros(p+1,Q,Rmax);
%         X = zeros(N,Q,Rmax);
%         errsvdZc = zeros(Q,Rmax);
%         err2Zc = zeros(Rmax,1);
%         Q = zeros(1,Rmax);
%         Qmax = 1;
%         for a=1:Rmax
%             Zca = Z(:,a,:);
%             Zca = Zca(:,:)';
%             [Wa,Sa,Xa,errsvdZca] = svdtruncate(Zca,tolsvdZc);
%             Sa = Sa/sqrt(N-1);
%             sa = diag(Sa);
%             Xa = Xa*sqrt(N-1);
%             if verLessThan('matlab','9.1') % compatibility (<R2016b)
%                 Zca_approx = Wa*Sa*Xa';
%             else
%                 Zca_approx = Wa*(sa.*Xa');                
%             end
%             % errZca = norm(Zca_approx-Zca)/norm(Zca);
%             errZca = errsvdZca(end);
%             Qa = length(sa);
%             Qmax = max(Qmax,Qa);
%             fprintf('\nCoordinate alpha = %2.f : rank Q = %d, error = %.3e for Z',a,Qa,errZca);
%             
%             % norm2Zca = sum(var(Zca,0,2));
%             % norm2Zca_approx = sum(var(Zca_approx,0,2));
%             % err2Zca = errZca^2*norm2Zca;
%             % err2Zca = norm2Zca-norm2Zca_approx;
%             % saf = svdtruncate(Zca,eps);
%             % saf = saf/sqrt(N-1);
%             % err2Zca = sum(saf.^2)-sum(sa.^2);
%             
%             s(1:Qa,a) = sa;
%             W(:,1:Qa,a) = Wa;
%             X(:,1:Qa,a) = Xa;
%             Q(a) = Qa;
%             errsvdZc(1:Qa,a) = errsvdZca;
%             
%             if verLessThan('matlab','9.1') % compatibility (<R2016b)
%                 CZa_approx = Wa*Sa.^2*Wa';
%             else
%                 CZa_approx = Wa*(sa.^2.*Wa');                
%             end
%             CZa = cov(Zca'); % CZa = 1/(N-1)*Zca*Zca';
%             errCZa = norm(CZa_approx-CZa)/norm(CZa);
%             fprintf('\n                                     error = %.3e for CZ',errCZa);
%             
%             if displayCovariance
%                 figure('Name','Covariance matrix')
%                 clf
%                 imagesc(CZa)
%                 colorbar
%                 axis image
%                 set(gca,'FontSize',fontsize)
%                 xlabel('$k$','Interpreter','latex')
%                 ylabel('$k''$','Interpreter','latex')
%                 title(['Covariance matrix $[C_{\zeta}(t^k,t^{k''})]_{\alpha,\alpha} = [C_{Z_{\alpha}}]_{k,k''}$ for $\alpha=$' num2str(a)],'Interpreter','latex')
%                 mysaveas(gridpathname,['covariance_CZ_a' num2str(a)],formats,renderer);
%                 mymatlab2tikz(gridpathname,['covariance_CZ_a' num2str(a) '.tex']);
%             end
%             
%             % mXa = mean(Xa,1)';
%             % CXa = cov(Xa); % CXa = 1/(N-1)*Xa'*Xa;
%             % norm(mXa)
%             % norm(CXa-eye(Qa))
%             % norm(Wa'*Wa-eye(Qa))
%         end
%         fprintf('\n');
%         s = s(1:Qmax,:);
%         W = W(:,1:Qmax,:);
%         X = X(:,1:Qmax,:);
%         errsvdZc = errsvdZc(1:Qmax,:);
%         Q = Qmax;
        
        %% Second reduction step
        fprintf('\nPCA in time');
        q = (p+1)*Rmax;
        switch index
            case 'time'
                Zc = permute(Z,[1,3,2]);
                Zc = Zc(:,:)';
            case 'coord'
                Zc = Z(:,:)';
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
        
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            CZ_approx = W*S.^2*W';
        else
            CZ_approx = W*(s.^2.*W');
        end
        CZ = cov(Zc'); % CZ = 1/(N-1)*Zc*Zc';
        errCZ = norm(CZ_approx-CZ)/norm(CZ);
        fprintf('\n                          error = %.3e for CZ',errCZ);
        fprintf('\n');
        
        save(fullfile(gridpathname,'PCA_time.mat'),'s','W','errsvdZc','Q');
        
        % mX = mean(X,1)';
        % CX = cov(X); % CX = 1/(N-1)*X'*X;
        % norm(mX)
        % norm(CX-eye(Q))
        % norm(W'*W-eye(Q))
        % norm(CZ_approx*W-(p+1)*W)
        % norm(CZ*W-(p+1)*W)
        % norm(abs(s.^2-(p+1)))
        
        if displayCovariance
            figure('Name','Covariance matrix')
            clf
            imagesc(CZ)
            colorbar
            axis image
            set(gca,'FontSize',fontsize)
            switch index
                case 'time'
                    xlabel('$K=\alpha p+k$','Interpreter','latex')
                    ylabel('$K''=\alpha'' p+k''$','Interpreter','latex')
                case 'coord'
                    xlabel('$K=(k-1)R+\alpha$','Interpreter','latex')
                    ylabel('$K''=(k''-1)R+\alpha''$','Interpreter','latex')
            end
            title('Covariance matrix $[C_{\zeta}(t^k,t^{k''})]_{\alpha,\alpha''} = [C_Z]_{K,K''}$','Interpreter','latex')
            mysaveas(gridpathname,'covariance_CZ',formats,renderer);
            mymatlab2tikz(gridpathname,'covariance_CZ.tex');
            
%             for t=0:p
%                 switch index
%                     case 'time'
%                         ind = (0:Rmax-1)*(p+1)+t+1;
%                     case 'coord'
%                         ind = t*Rmax+(1:Rmax);
%                 end
%                 CZt_approx = CZ_approx(ind,ind);
%                 CZt = CZ(ind,ind);
%                 
%                 figure('Name','Covariance matrix')
%                 clf
%                 imagesc(CZt)
%                 colorbar
%                 axis image
%                 set(gca,'FontSize',fontsize)
%                 xlabel('$\alpha$','Interpreter','latex')
%                 ylabel('$\alpha''$','Interpreter','latex')
%                 title(['Covariance matrix $[C_{\zeta}(t^k,t^k)]_{\alpha,\alpha''} = [C_{Z_k}]_{\alpha,\alpha''}$ for $k=$' num2str(t)],'Interpreter','latex')
%                 mysaveas(gridpathname,['covariance_CZ_t' num2str(t*100)],formats,renderer);
%                 mymatlab2tikz(gridpathname,['covariance_CZ_t' num2str(t*100) '.tex']);
%             end
        end
        
        %% Post-processing data
        fprintf('\nPost-processing data');
        nt = 13; % number of post-processed variables
        ne = 9; % number of energy variables
        nvar = n+nt+ne;
        fprintf('\nn = %d post-processed variables',nt);
        fprintf('\n  = %d energy variables',ne);
        fprintf('\n');
        
        if g<2^7
            Tau = zeros(N,nt,m,p+1);
            E = zeros(N,ne,m,p+1);
        end
        I = zeros(N,nt+ne,2,p+1);
        
        for t=0:p
            time = ['Time ' num2str(t)];
            disp(time)
            
            if t==0
                mYt_old = zeros(1,n,m);
                Rt_old = Rmax;
                sigt_old = zeros(Rmax,1);
                % Zt_old = zeros(N,Rmax);
                Vt_old = zeros(n*m,Rmax);
                Wt_old = zeros(Rmax,Q);
            else
                mYt_old = mY(1,:,:,t);
                Rt_old = R(t);
                sigt_old = sig(1:Rt_old,t);
                % Zt_old = Z(:,1:Rt_old,t);
                if g<2^7
                    Vt_old = V(:,1:Rt_old,t);
                else
                    load(fullfile(gridpathname,['PCAspace_t' num2str(t-1) '.mat']),'Vt');
                    Vt_old = Vt;
                end
                switch index
                    case 'time'
                        Wt_old = W(t:(p+1):end,:);
                        Wt_old = Wt_old(1:Rt_old,:);
                    case 'coord'
                        Wt_old = W(Rmax*(t-1)+(1:Rt_old),:);
                end
            end
            
            if t==p
                mYt_new = zeros(1,n,m);
                Rt_new = Rmax;
                sigt_new = zeros(Rmax,1);
                % Zt_new = zeros(N,Rmax);
                Vt_new = zeros(n*m,Rmax);
                Wt_new = zeros(Rmax,Q);
            else
                mYt_new = mY(1,:,:,t+2);
                Rt_new = R(t+2);
                sigt_new = sig(1:Rt_new,t+2);
                % Zt_new = Z(:,1:Rt_new,t+2);
                if g<2^7
                    Vt_new = V(:,1:Rt_new,t+2);
                else
                    load(fullfile(gridpathname,['PCAspace_t' num2str(t+1) '.mat']),'Vt');
                    Vt_new = Vt;
                end
                switch index
                    case 'time'
                        Wt_new = W((t+2):(p+1):end,:);
                        Wt_new = Wt_new(1:Rt_new,:);
                    case 'coord'
                        Wt_new = W(Rmax*(t+1)+(1:Rt_new),:);
                end
            end
            
            mYt = mY(1,:,:,t+1);
            Rt = R(t+1);
            sigt = sig(1:Rt,t+1);
            % Zt = Z(:,1:Rt,t+1);
            if g<2^7
                Vt = V(:,1:Rt,t+1);
            else
                load(fullfile(gridpathname,['PCAspace_t' num2str(t) '.mat']),'Vt');
            end
            switch index
                case 'time'
                    Wt = W((t+1):(p+1):end,:);
                    Wt = Wt(1:Rt,:);
                case 'coord'
                    Wt = W(Rmax*t+(1:Rt),:);
            end
            
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                Zct = Wt*diag(s)*X'; % Zct = Zt';
                Zct_old = Wt_old*diag(s)*X'; % Zct_old = Zt_old';
                Zct_new = Wt_new*diag(s)*X'; % Zct_new = Zt_new';
                Yct = Vt*diag(sigt)*Zct;
                Yct_old = Vt_old*diag(sigt_old)*Zct_old;
                Yct_new = Vt*diag(sigt_new)*Zct_new;
            else
                Zct = Wt*(s.*X'); % Zct = Zt';
                Zct_old = Wt_old*(s.*X'); % Zct_old = Zt_old';
                Zct_new = Wt_new*(s.*X'); % Zct_new = Zt_new';
                Yct = Vt*(sigt.*Zct);
                Yct_old = Vt_old*(sigt_old.*Zct_old);
                Yct_new = Vt*(sigt_new.*Zct_new);
            end
            
            mYt = perm(reshape(mYt,[1,n,sx]));
            mYt_old = perm(reshape(mYt_old,[1,n,sx]));
            mYt_new = perm(reshape(mYt_new,[1,n,sx]));
            
            Yt = repmat(mYt,[ones(1,4),N]) + perm(reshape(Yct',[N,n,sx]));
            Yt_old = repmat(mYt_old,[ones(1,4),N]) + perm(reshape(Yct_old',[N,n,sx]));
            Yt_new = repmat(mYt_new,[ones(1,4),N]) + perm(reshape(Yct_new',[N,n,sx]));
            
%             vt = perm(reshape(Vt',[Rt,n,sx]));
%             vt_old = perm(reshape(Vt_old',[Rt_old,n,sx]));
%             vt_new = perm(reshape(Vt_new',[Rt_new,n,sx]));
            
%             rep = @(Yt) [3,ones(1,ndims(Yt))];
%             get_ut = @(Yt) Yt(1:3,:,:,:,:);
%             get_uti = @(Yt,i) Yt(i,:,:,:,:);
%             get_Ct = @(Yt) Yt(4,:,:,:,:);
%             get_rhot = @(Yt) rho(1) + (rho(2)-rho(1))*get_Ct(Yt);
%             get_mut = @(Yt) mu(1) + (mu(2)-mu(1))*get_Ct(Yt);
%             get_rhout = @(Yt) repmat(get_rhot(Yt),[3,ones(1,ndims(Yt)-1)]).*get_ut(Yt);
%             get_rhouti = @(Yt,i) get_rhot(Yt).*get_uti(Yt,i);
%             get_tauTime = @(Yt_new,Yt_old) (get_rhout(Yt_new)-get_rhout(Yt_old))/(2*dt);
%             get_u2t = @(Yt) dot(get_ut(Yt),get_ut(Yt),1);
%             get_Ek = @(Yt) 1/2*get_rhot(Yt).*get_u2t(Yt);
%             get_energyKinTime = @(Yt_new,Yt_old) (get_Ek(Yt_new)-get_Ek(Yt_old))/(2*dt);
%             
%             get_divut = @(Yt) div(get_ut(Yt),Dx);
%             get_gradut = @(Yt) grad(get_ut(Yt),Dx);
%             get_gradrhout = @(Yt) grad(repmat(get_rhot(Yt),rep(Yt)).*get_ut(Yt),Dx);
%             get_St = @(Yt) (get_gradut(Yt)+permute(get_gradut(Yt),[2,1,3:ndims(Yt)+1]))/2;
%             get_muSt = @(Yt) repmat(shiftdim(get_mut(Yt),-1),[3,3,ones(1,ndims(Yt)-1)]).*get_St(Yt);
%             get_gradCt = @(Yt) grad(get_Ct(Yt),Dx);
%             get_rhou2t = @(Yt) get_rhot(Yt).*get_u2t(Yt);
%             
%             ut = get_ut(Yt); 
%             ut_old = get_ut(Yt_old);
%             ut_new = get_ut(Yt_new);
%             
%             Ct = get_Ct(Yt);
%             Ct_old = get_Ct(Yt_old);
%             Ct_new = get_Ct(Yt_new);
%             
%             % rhot = get_rhot(Yt);
%             rhot_old = get_rhot(Yt_old);
%             rhot_new = get_rhot(Yt_new);
%             % rhout = get_rhout(Yt);
%             rhout_old = get_rhout(Yt_old);
%             rhout_new = get_rhout(Yt_new);
%             tauTimet = get_tauTime(Yt_new,Yt_old) ;
%             % u2t = get_u2t(Yt);
%             u2t_old = get_u2t(Yt_old);
%             u2t_new = get_u2t(Yt_new);
%             % Ek = get_Ek(Yt);
%             Ek_old = get_Ek(Yt_old);
%             Ek_new = get_Ek(Yt_new);
%             energyKinTimet = get_energyKinTime(Yt_new,Yt_old);
%             
%             rhot = get_rhot(Yt);
%             % mut =  get_mut(Yt);
%             gradut = get_gradut(Yt);
%             rhout = get_rhout(Yt);
%             gradrhout = get_gradrhout(Yt);
%             % St = get_St(Yt);
%             muSt = get_muSt(Yt);
%             gradCt = get_gradCt(Yt);
%             ngradCt = normal(get_gradCt(Yt));
%             kappa = div(normal(get_gradCt(Yt)),Dx);
%             u2t = get_u2t(Yt);
%             rhou2t = get_rhou2t(Yt);
            
%             get_tauConv = @(Yt) squeeze(sum(get_gradrhout(Yt).*repmat(shiftdim(get_ut(Yt),-1),rep(Yt)),2));
%             % get_tauConv = @(Yt) div(permute(repmat(shiftdim(get_rhout(Yt),-1),rep(Yt)),[2,1,3:ndims(Yt)+1]).*repmat(shiftdim(get_ut(Yt),-1),rep(Yt)),Dx);
%             get_tauDiff = @(Yt) div(2*get_muSt(Yt),Dx);
%             get_tauSurf = @(Yt) sigma*repmat(shiftdim(kappa,-1),[3,ones(1,ndims(Yt)-1)]).*get_gradCt(Yt);
%             get_tauInterf = @(Yt) dot(get_ut(Yt),get_gradCt(Yt),1);
%             get_energyConv = @(Yt) shiftdim(div(repmat(get_rhou2t(Yt),[3,ones(1,ndims(Yt)-1)]).*get_ut(Yt),Dx),-1);
%             get_energyGrav = @(Yt) gravity.*get_rhouti(Yt,3);
%             get_energyKinSpace = @(Yt) dot(get_rhout(Yt),grad(get_u2t(Yt)/2,Dx),1);
%             get_energyDiff = @(Yt) shiftdim(div(squeeze(dot(2*get_muSt(Yt),repmat(shiftdim(get_ut(Yt),-1),rep(Yt)),2)),Dx),-1);
%             get_energyVisc = @(Yt) shiftdim(sum(sum(2*get_muSt(Yt).*get_gradut(Yt),1),2),1);
%             get_energySurf = @(Yt) dot(get_tauSurf(Yt),get_ut(Yt),1);
            
%             tauConvt = get_tauConv(Yt);
%             tauDifft = get_tauDiff(Yt);
%             tauSurft = get_tauSurf(Yt);
%             tauInterft = get_tauInterf(Yt);
%             energyConvt = get_energyConv(Yt);
%             energyGravt = get_energyGrav(Yt);
%             energyPrest = zeros(1,g+1,g+1,g+1,N);
%             energyPresDilt = zeros(1,g+1,g+1,g+1,N);
%             energyKinSpacet = get_energyKinSpace(Yt);
%             energyDifft = get_energyDiff(Yt);
%             energyVisct = get_energyVisc(Yt);
%             energySurft = get_energySurf(Yt);

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
            
            tauConvt = squeeze(sum(gradrhout.*repmat(shiftdim(ut,-1),[3,ones(1,5)]),2));
            % tauConvt = div(permute(repmat(shiftdim(rhout,-1),[3,ones(1,5)]),[2,1,3:6]).*repmat(shiftdim(ut,-1),[3,ones(1,5)]),Dx);
            tauDifft = div(2*muSt,Dx);
            tauSurft = sigma*repmat(shiftdim(kappa,-1),[3,ones(1,4)]).*gradCt;
            tauInterft = dot(ut,gradCt,1);
            energyConvt = shiftdim(div(repmat(rhou2t,[3,ones(1,4)]).*ut,Dx),-1);
            energyGravt = gravity.*rhout(3,:,:,:,:);
            energyPrest = zeros(1,g+1,g+1,g+1,N);
            energyPresDilt = zeros(1,g+1,g+1,g+1,N);
            energyKinSpacet = dot(rhout,grad(u2t/2,Dx),1);
            energyDifft = shiftdim(div(squeeze(dot(2*muSt,repmat(shiftdim(ut,-1),[3,ones(1,5)]),2)),Dx),-1);
            energyVisct = shiftdim(sum(sum(2*muSt.*gradut,1),2),1);
            energySurft = dot(tauSurft,ut,1);
            
%             mut = get_ut(mYt);
%             mut_old = get_ut(mYt_old);
%             mut_new = get_ut(mYt_new);
%             
%             mCt = get_Ct(mYt);
%             mCt_old = get_Ct(mYt_old);
%             mCt_new = get_Ct(mYt_new);
%             
%             % mrhot = get_rhot(mYt); % mrhot = mCt*rho(2) + (1-mCt)*rho(1);
%             mrhot_old = get_rhot(mYt_old); % rhot_old = mCt_old*rho(2) + (1-mCt_old)*rho(1);
%             mrhot_new = get_rhot(mYt_new); % rhot_new = mCt_new*rho(2) + (1-mCt_new)*rho(1);
%             % rhout = get_rhout(Yt); % rhout = repmat(rhot,[3,ones(1,4)]).*ut;
%             rhout_old = get_rhout(Yt_old); % rhout_old = repmat(rhot_old,[3,ones(1,4)]).*ut_old;
%             rhout_new = get_rhout(Yt_new); % rhout_new = repmat(rhot_new,[3,ones(1,4)]).*ut_new;
%             tauTime = get_tauTime(Yt_new,Yt_old) ; % tauTime = (rhout_new-rhout_old)/(2*dt);
%             % u2t = get_u2t(Yt); % u2t = dot(ut,ut,1);
%             u2t_old = get_u2t(Yt_old); % u2t_old = dot(ut_old,ut_old,1);
%             u2t_new = get_u2t(Yt_new); % u2t_new = dot(ut_new,ut_new,1);
%             % Ek = get_Ek(Yt); % Ek = 1/2*rhot.*u2t;
%             Ek_old = get_Ek(Yt_old); % Ek_old = 1/2*rhot_old.*u2t_old;
%             Ek_new = get_Ek(Yt_new); % Ek_new = 1/2*rhot_new.*u2t_new;
%             energyKinTimet = get_energyKinTime(Yt_new,Yt_old); % energyKinTime = (Ek_new-Ek_old)/(2*dt);
%             
%             mrhot = get_rhot(mYt); % mrhot = mCt*rho(2) + (1-mCt)*rho(1);
%             % mmut =  get_mut(mYt); % mmut = mCt*mu(2) + (1-mCt)*mu(1);
%             gradmut = get_gradut(mYt); % gradut = grad(ut,Dx);
%             mrhout = get_rhout(mYt); % rhout = repmat(rhot,[3,ones(1,4)]).*ut;
%             gradmrhout = get_gradrhout(mYt); % gradrhout = grad(rhout,Dx);
%             % mSt = get_St(mYt); % St = (gradmut+permute(gradmut,[2,1,3:6]))/2;
%             mmuSt = get_muSt(mYt); % mmuSt = repmat(shiftdim(mmut,-1),[3,3,ones(1,4)]).*St;
%             gradmCt = get_gradCt(mYt); % gradmCt = grad(mCt,Dx);
%             ngradmCt = normal(get_gradCt(mYt)); % ngradmCt = normal(gradmCt);
%             mkappa = div(normal(get_gradCt(mYt)),Dx); % mkappa = div(ngradmCt,Dx);
%             mu2t = get_u2t(mYt); % mu2t = dot(mut,mut,1);
%             mrhou2t = get_rhou2t(mYt); % mrhou2t = mrhot.*mu2t;
%             
%             mtauConv = div(permute(repmat(shiftdim(mrhout,-1),rep(mYt)),[2,1,3:ndims(mYt)+1]).*repmat(shiftdim(mut,-1),rep(mYt)),Dx);
            
            Taut = cat(1,tauTimet,tauConvt,tauDifft,tauSurft,tauInterft);
            Et = cat(1,energyKinTimet,energyConvt,energyGravt,energyPrest,energyPresDilt,energyKinSpacet,energyDifft,energyVisct,energySurft);
            
            It = cat(1,Taut,Et);
            It = cat(2,reshape(trapz(x,trapz(x,trapz(x,repmat(1-Ct,[nt+ne,ones(1,4)]).*It,2),3),4),[nt+ne,1,N]),...
                reshape(trapz(x,trapz(x,trapz(x,repmat(Ct,[nt+ne,ones(1,4)]).*It,2),3),4),[nt+ne,1,N]));
            It = shiftdim(It,2);
            
            Taut = iperm(Taut);
            Et = iperm(Et);
            Taut = Taut(:,:,:);
            Et = Et(:,:,:);
            
            if g<2^7
                Tau(:,:,:,t+1) = Taut;
                E(:,:,:,t+1) = Et;
            else
                save(fullfile(gridpathname,['PP_t' num2str(t) '.mat']),'Taut','Et');
            end
            I(:,:,:,t+1) = It;
        end
        fprintf('\n');
        if g<2^7
            save(fullfile(gridpathname,'PP.mat'),'Tau','E');
        end
        save(fullfile(gridpathname,'PP.mat'),'I','nt','ne','-append');
        
        %% Post processing PCA
%         for t=0:p
%             time = ['Time ' num2str(t)];
%             disp(time)
%             
%             mYt = mY(1,:,:,t+1);
%             Rt = R(t+1);
%             sigt = sig(1:Rt,t+1);
%             Zt = Z(:,1:Rt,t+1);
%             if g<2^7
%                 Vt = V(:,1:Rt,t+1);
%             else
%                 load(fullfile(gridpathname,['PCAspace_t' num2str(t) '.mat']),'Vt');
%             end
%             if verLessThan('matlab','9.1') % compatibility (<R2016b)
%                 Yct = Vt*(sigt.*Zt');
%             else
%                 Yct = Vt*diag(sigt)*Zt';
%             end
%             % Yct = permute(reshape(Yct,[n,sx,N]),[1,order+1,1]);
%             
%             % Vt = permute(reshape(Vt,[n,sx,Rt]),[1,order+1,1]);
%             
%             switch index
%                 case 'time'
%                     Wt = W((t+1):(p+1):end,:);
%                     Wt = Wt(1:Rt,:);
%                 case 'coord'
%                     Wt = W(Rmax*t+(1:Rt),:);
%             end
%             
%             if verLessThan('matlab','9.1') % compatibility (<R2016b)
%                 % Zct_approx = Wt*diag(s)*X'; % Zct_approx = Z_approx(:,1:Rt,t+1)';
%                 % Yct_approx = Vt*diag(sigt)*Zct_approx;
%                 Yct_approx = Vt*diag(sigt)*Wt*diag(s)*X';
%             else
%                 % Zct_approx = Wt*(s.*X'); % Zct_approx = Z_approx(:,1:Rt,t+1)';
%                 % Yct_approx = Vt*(sigt.*Zct_approx);
%                 Yct_approx = Vt*(sigt.*Wt*(s.*X'));
%             end
%         end
    
    else
        if g<2^7
            load(fullfile(gridpathname,'PCA_space.mat'),'V');
            load(fullfile(gridpathname,'PP.mat'),'Tau','E'); 
        end
        load(fullfile(gridpathname,'PCA_space.mat'),'mY','sig','Z','errsvdYc','R','Rmax');
        load(fullfile(gridpathname,'PCA_time.mat'),'s','W','Q','errsvdZc');
        load(fullfile(gridpathname,'PP.mat'),'I','nt','ne');
    end
    
    %% Outputs
    % Display eigenvalues
    if displayEigenvales
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
        semilogy(1:Q,s(1:Q).^2,'LineStyle','-','Color','b','LineWidth',1);
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
    
    % Mesh
    L = 1;
    D = DOMAIN(3,[0.0,0.0,0.0],[L,L,L]);
    elemtype = 'CUB8';
    nbelem = repmat(g,1,3);
    M = build_model(D,'nbelem',nbelem,'elemtype',elemtype);
    coord = getcoord(getnode(M));
    M = setnode(M,NODE(fliplr(coord)));
    M = final(M,DDL(DDLVECT('U',M.syscoord)));
    tf = p*dt;
    T = TIMEMODEL(0,tf,p);
    
    % Solution
    switch index
        case 'time'
            W = permute(reshape(W',[Q,p+1,Rmax]),[3,1,2]);
            % Z_approx = permute(reshape(Zc_approx',[N,p+1,Rmax]),[1,3,2]);
        case 'coord'
            W = permute(reshape(W',[Q,Rmax,p+1]),[2,1,3]);
            % Z_approx = reshape(Zc_approx',[N,Rmax,p+1]);
    end
    if g<2^7
        % Yc_approx = zeros(N,n,m,p+1);
        Uc_approx = zeros(Q,n,m,p+1);
    end
    
    if g<2^7
        load(fullfile(gridpathname,'data.mat'),'Y');
        load(fullfile(gridpathname,'PP.mat'),'Tau','E');
        Y = cat(2,Y,Tau,E);
        clear Tau E
        mY = mean(Y,1);
        Yc = Y - repmat(mY,N,1,1,1); % Yc = Y - mY.*ones(N,1,1,1);
        clear Y
    end
    
    % Compute mean
    disp('Mean')
    mu = reshape(mY(1,1:3,:,:),[3*m,p+1]);
    mC = reshape(mY(1,4,:,:),[m,p+1]);
    mtauTime = reshape(mY(1,5:7,:,:),[3*m,p+1]);
    mtauConv = reshape(mY(1,8:10,:,:),[3*m,p+1]);
    mtauDiff = reshape(mY(1,11:13,:,:),[3*m,p+1]);
    mtauSurf = reshape(mY(1,14:16,:,:),[3*m,p+1]);
    mtauInterf = reshape(mY(1,17,:,:),[m,p+1]);
    menergyKinTime = reshape(mY(1,18,:,:),[m,p+1]);
    menergyConv = reshape(mY(1,19,:,:),[m,p+1]);
    menergyGrav = reshape(mY(1,20,:,:),[m,p+1]);
    menergyPres = reshape(mY(1,21,:,:),[m,p+1]);
    menergyPresDil = reshape(mY(1,22,:,:),[m,p+1]);
    menergyKinSpace = reshape(mY(1,23,:,:),[m,p+1]);
    menergyDiff = reshape(mY(1,24,:,:),[m,p+1]);
    menergyVisc = reshape(mY(1,25,:,:),[m,p+1]);
    menergySurf = reshape(mY(1,26,:,:),[m,p+1]);
    
    % Compute variance
    disp('Variance')
    vu = zeros(3*m,p+1);
    vC = zeros(m,p+1);
    vtauTime = zeros(3*m,p+1);
    vtauConv = zeros(3*m,p+1);
    vtauDiff = zeros(3*m,p+1);
    vtauSurf = zeros(3*m,p+1);
    vtauInterf = zeros(m,p+1);
    venergyKinTime = zeros(m,p+1);
    venergyConv = zeros(m,p+1);
    venergyGrav = zeros(m,p+1);
    venergyPres = zeros(m,p+1);
    venergyPresDil = zeros(m,p+1);
    venergyKinSpace = zeros(m,p+1);
    venergyDiff = zeros(m,p+1);
    venergyVisc = zeros(m,p+1);
    venergySurf = zeros(m,p+1);
    for t=0:p
        if g<2^7
            Yct = Yc(:,:,:,t+1);
            Vt = V(:,1:Rt,t+1);
        else
            load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
            load(fullfile(gridpathname,['PP_t' num2str(t) '.mat']),'Taut','Et');
            Yt = cat(2,Yt,Taut,Et);
            clear Taut Et
            mYt = mean(Yt,1);
            Yct = Yt - repmat(mYt,[N,1,1]); % Yct = Yt - mYt.*ones(N,1,1);
            clear Yt
            load(fullfile(gridpathname,['PCA_space_t' num2str(t) '.mat']),'Vt');
        end
        Rt = R(t+1);
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
        if g<2^7
            % Yc_approx(:,:,:,t+1) = reshape(Yct_approx',[N,n,m]);
            Uc_approx(:,:,:,t+1) = reshape(Uct_approx',[Q,n,m]);
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
        
        indut = (0:m-1)*nvar+(1:3)';
        indCt = (0:m-1)*nvar+4;
        indtauTimet = (0:m-1)*nvar+(5:7)';
        indtauConvt = (0:m-1)*nvar+(8:10)';
        indtauDifft = (0:m-1)*nvar+(11:13)';
        indtauSurft = (0:m-1)*nvar+(14:16)';
        indtauInterft = (0:m-1)*nvar+17;
        indenergyKinTimet = (0:m-1)*nvar+18;
        indenergyConvt = (0:m-1)*nvar+19;
        indenergyGravt = (0:m-1)*nvar+20;
        indenergyPrest = (0:m-1)*nvar+21;
        indenergyPresDilt = (0:m-1)*nvar+22;
        indenergyKinSpacet = (0:m-1)*nvar+23;
        indenergyDifft = (0:m-1)*nvar+24;
        indenergyVisct = (0:m-1)*nvar+25;
        indenergySurft = (0:m-1)*nvar+26;
        
        tauTimet = Yct(:,5:7,:);
        tauConvt = Yct(:,8:10,:);
        tauDifft = Yct(:,11:13,:);
        tauSurft = Yct(:,14:16,:);
        tauInterft = Yct(:,17,:);
        energyKinTimet = Yct(:,18,:);
        energyConvt = Yct(:,19,:);
        energyGravt = Yct(:,20,:);
        energyPrest = Yct(:,21,:);
        energyPresDilt = Yct(:,22,:);
        energyKinSpacet = Yct(:,23,:);
        energyDifft = Yct(:,24,:);
        energyVisct = Yct(:,25,:);
        energySurft = Yct(:,26,:);
        
        indut_ = (0:m-1)*n+(1:3)';
        indCt_ = (0:m-1)*n+4;
        
        vut = vYt_approx(indut_(:));
        vCt = vYt_approx(indCt_(:));
        vtauTimet = 1/(N-1)*sum(tauTimet(:,:).^2)';
        vtauConvt = 1/(N-1)*sum(tauConvt(:,:).^2)';
        vtauDifft = 1/(N-1)*sum(tauDifft(:,:).^2)';
        vtauSurft = 1/(N-1)*sum(tauSurft(:,:).^2)';
        vtauInterft = 1/(N-1)*sum(tauInterft(:,:).^2)';
        venergyKinTimet = 1/(N-1)*sum(energyKinTimet(:,:).^2)';
        venergyConvt = 1/(N-1)*sum(energyConvt(:,:).^2)';
        venergyGravt = 1/(N-1)*sum(energyGravt(:,:).^2)';
        venergyPrest = 1/(N-1)*sum(energyPrest(:,:).^2)';
        venergyPresDilt = 1/(N-1)*sum(energyPresDilt(:,:).^2)';
        venergyKinSpacet = 1/(N-1)*sum(energyKinSpacet(:,:).^2)';
        venergyDifft = 1/(N-1)*sum(energyDifft(:,:).^2)';
        venergyVisct = 1/(N-1)*sum(energyVisct(:,:).^2)';
        venergySurft = 1/(N-1)*sum(energySurft(:,:).^2)';
        
        vYt_approx = zeros(nvar*m,1);
        vYt_approx(indut(:)) = vut;
        vYt_approx(indCt(:)) = vCt;
        vYt_approx(indtauTimet(:)) = vtauTimet;
        vYt_approx(indtauConvt(:)) = vtauConvt;
        vYt_approx(indtauDifft(:)) = vtauDifft;
        vYt_approx(indtauSurft(:)) = vtauSurft;
        vYt_approx(indtauInterft(:)) = vtauInterft;
        vYt_approx(indenergyKinTimet(:)) = venergyKinTimet;
        vYt_approx(indenergyConvt(:)) = venergyConvt;
        vYt_approx(indenergyGravt(:)) = venergyGravt;
        vYt_approx(indenergyPrest(:)) = venergyPrest;
        vYt_approx(indenergyPresDilt(:)) = venergyPresDilt;
        vYt_approx(indenergyKinSpacet(:)) = venergyKinSpacet;
        vYt_approx(indenergyDifft(:)) = venergyDifft;
        vYt_approx(indenergyVisct(:)) = venergyVisct;
        vYt_approx(indenergySurft(:)) = venergySurft;
        
        % CYt = cov(Yct(:,:)); % CYt = 1/(N-1)*Yct(:,:)'*Yct(:,:);
        % vYt = diag(CYt);
        vYt = 1/(N-1)*sum(Yct(:,:).^2)';
        errvYt = norm(vYt_approx-vYt)/norm(vYt);
        fprintf('\nTime %2d, t = %4g s : error = %.3e for VY',t,t*dt,errvYt);
        
        vu(:,t+1) = vYt_approx(indut(:));
        vC(:,t+1) = vYt_approx(indCt(:));
        vtauTime(:,t+1) = vYt_approx(indtauTimet(:));
        vtauConv(:,t+1) = vYt_approx(indtauConvt(:));
        vtauDiff(:,t+1) = vYt_approx(indtauDifft(:));
        vtauSurf(:,t+1) = vYt_approx(indtauSurft(:));
        vtauInterf(:,t+1) = vYt_approx(indtauInterft(:));
        venergyKinTimet(:,t+1) = vYt_approx(indenergyKinTimet(:));
        venergyConvt(:,t+1) = vYt_approx(indenergyConvt(:));
        venergyGravt(:,t+1) = vYt_approx(indenergyGravt(:));
        venergyPrest(:,t+1) = vYt_approx(indenergyPrest(:));
        venergyPresDilt(:,t+1) = vYt_approx(indenergyPresDilt(:));
        venergyKinSpacet(:,t+1) = vYt_approx(indenergyKinSpacet(:));
        venergyDifft(:,t+1) = vYt_approx(indenergyDifft(:));
        venergyVisct(:,t+1) = vYt_approx(indenergyVisct(:));
        venergySurft(:,t+1) = vYt_approx(indenergySurft(:));
    end
    
    fprintf('\n');
    mu = TIMEMATRIX(mu,T);
    mC = TIMEMATRIX(mC,T);
    mtauTime = TIMEMATRIX(mtauTime,T);
    mtauConv = TIMEMATRIX(mtauConv,T);
    mtauDiff = TIMEMATRIX(mtauDiff,T);
    mtauSurf = TIMEMATRIX(mtauSurf,T);
    mtauInterf = TIMEMATRIX(mtauInterf,T);
    menergyKinTime = TIMEMATRIX(menergyKinTime,T);
    menergyConv = TIMEMATRIX(menergyConv,T);
    menergyGrav = TIMEMATRIX(menergyGrav,T);
    menergyPres = TIMEMATRIX(menergyPres,T);
    menergyPresDil = TIMEMATRIX(menergyPresDil,T);
    menergyKinSpace = TIMEMATRIX(menergyKinSpace,T);
    menergyDiff = TIMEMATRIX(menergyDiff,T);
    menergyVisc = TIMEMATRIX(menergyVisc,T);
    menergySurf = TIMEMATRIX(menergySurf,T);
    
    vu = TIMEMATRIX(vu,T);
    vC = TIMEMATRIX(vC,T);
    vtauTime = TIMEMATRIX(vtauTime,T);
    vtauConv = TIMEMATRIX(vtauConv,T);
    vtauDiff = TIMEMATRIX(vtauDiff,T);
    vtauSurf = TIMEMATRIX(vtauSurf,T);
    vtauInterf = TIMEMATRIX(vtauInterf,T);
    venergyKinTime = TIMEMATRIX(venergyKinTime,T);
    venergyConv = TIMEMATRIX(venergyConv,T);
    venergyGrav = TIMEMATRIX(venergyGrav,T);
    venergyPres = TIMEMATRIX(venergyPres,T);
    venergyPresDil = TIMEMATRIX(venergyPresDil,T);
    venergyKinSpace = TIMEMATRIX(venergyKinSpace,T);
    venergyDiff = TIMEMATRIX(venergyDiff,T);
    venergyVisc = TIMEMATRIX(venergyVisc,T);
    venergySurf = TIMEMATRIX(venergySurf,T);
    
    if g<2^7
        % CY_approx = cov(Yc_approx(:,:)); % CY_approx = 1/(N-1)*Yc_approx(:,:)'*Yc_approx(:,:);
        % if verLessThan('matlab','9.1') % compatibility (<R2016b)
        %     CY_approx = Uc_approx(:,:)'*diag(s).^2*Uc_approx(:,:);
        % else
        %     CY_approx = Uc_approx(:,:)'*(s.^2.*Uc_approx(:,:));
        % end
        % vY_approx = diag(CY_approx);
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            vY_approx = sum((diag(s)*Uc_approx(:,:)).^2)';
        else
            vY_approx = sum((s.*Uc_approx(:,:)).^2)';
        end
        tauTime = permute(Yc(:,5:7,:,:),[1,4,2,3]);
        tauConv = permute(Yc(:,8:10,:,:),[1,4,2,3]);
        tauDiff = permute(Yc(:,11:13,:,:),[1,4,2,3]);
        tauSurf = permute(Yc(:,14:16,:,:),[1,4,2,3]);
        tauInterf = permute(Yc(:,17,:,:),[1,4,2,3]);
        energyKinTime = permute(Yc(:,18,:,:),[1,4,2,3]);
        energyConv = permute(Yc(:,19,:,:),[1,4,2,3]);
        energyGrav = permute(Yc(:,20,:,:),[1,4,2,3]);
        energyPres = permute(Yc(:,21,:,:),[1,4,2,3]);
        energyPresDil = permute(Yc(:,22,:,:),[1,4,2,3]);
        energyKinSpace = permute(Yc(:,23,:,:),[1,4,2,3]);
        energyDiff = permute(Yc(:,24,:,:),[1,4,2,3]);
        energyVisc = permute(Yc(:,25,:,:),[1,4,2,3]);
        energySurf = permute(Yc(:,26,:,:),[1,4,2,3]);
        indu_ = (0:m-1)*n+(1:3)';
        indC_ = (0:m-1)*n+4;
        indu_ = indu_(:)'+n*m*(0:p)';
        indC_ = indC_(:)'+n*m*(0:p)';
        varu = vY_approx(indu_(:));
        varC = vY_approx(indC_(:));
        vartauTime = 1/(N-1)*sum(tauTime(:,:).^2)';
        vartauConv = 1/(N-1)*sum(tauConv(:,:).^2)';
        vartauDiff = 1/(N-1)*sum(tauDiff(:,:).^2)';
        vartauSurf = 1/(N-1)*sum(tauSurf(:,:).^2)';
        vartauInterf = 1/(N-1)*sum(tauInterf(:,:).^2)';
        varenergyKinTime = 1/(N-1)*sum(energyKinTime(:,:).^2)';
        varenergyConv = 1/(N-1)*sum(energyConv(:,:).^2)';
        varenergyGrav = 1/(N-1)*sum(energyGrav(:,:).^2)';
        varenergyPres = 1/(N-1)*sum(energyPres(:,:).^2)';
        varenergyPresDil = 1/(N-1)*sum(energyPresDil(:,:).^2)';
        varenergyKinSpace = 1/(N-1)*sum(energyKinSpace(:,:).^2)';
        varenergyDiff = 1/(N-1)*sum(energyDiff(:,:).^2)';
        varenergyVisc = 1/(N-1)*sum(energyVisc(:,:).^2)';
        varenergySurf = 1/(N-1)*sum(energySurf(:,:).^2)';
        vY_approx = zeros(nvar*m*(p+1),1);
        indu = (0:m-1)*nvar+(1:3)';
        indC = (0:m-1)*nvar+4;
        indtauTime = (0:m-1)*nvar+(5:7)';
        indtauConv = (0:m-1)*nvar+(8:10)';
        indtauDiff = (0:m-1)*nvar+(11:13)';
        indtauSurf = (0:m-1)*nvar+(14:16)';
        indtauInterf = (0:m-1)*nvar+17;
        indenergyKinTime = (0:m-1)*nvar+18;
        indenergyConv = (0:m-1)*nvar+19;
        indenergyGrav = (0:m-1)*nvar+20;
        indenergyPres = (0:m-1)*nvar+21;
        indenergyPresDil = (0:m-1)*nvar+22;
        indenergyKinSpace = (0:m-1)*nvar+23;
        indenergyDiff = (0:m-1)*nvar+24;
        indenergyVisc = (0:m-1)*nvar+25;
        indenergySurf = (0:m-1)*nvar+26;
        
        indu = indu(:)'+nvar*m*(0:p)';
        indC = indC(:)'+nvar*m*(0:p)';
        indtauTime = indtauTime(:)'+nvar*m*(0:p)';
        indtauConv = indtauConv(:)'+nvar*m*(0:p)';
        indtauDiff = indtauDiff(:)'+nvar*m*(0:p)';
        indtauSurf = indtauSurf(:)'+nvar*m*(0:p)';
        indtauInterf = indtauInterf(:)'+nvar*m*(0:p)';
        indenergyKinTime = indenergyKinTime(:)'+nvar*m*(0:p)';
        indenergyConv = indenergyConv(:)'+nvar*m*(0:p)';
        indenergyGrav = indenergyGrav(:)'+nvar*m*(0:p)';
        indenergyPres = indenergyPres(:)'+nvar*m*(0:p)';
        indenergyPresDil = indenergyPresDil(:)'+nvar*m*(0:p)';
        indenergyKinSpace = indenergyKinSpace(:)'+nvar*m*(0:p)';
        indenergyDiff = indenergyDiff(:)'+nvar*m*(0:p)';
        indenergyVisc = indenergyVisc(:)'+nvar*m*(0:p)';
        indenergySurf = indenergySurf(:)'+nvar*m*(0:p)';
        
        vY_approx(indu(:)) = varu;
        vY_approx(indC(:)) = varC;
        vY_approx(indtauTime(:)) = vartauTime;
        vY_approx(indtauConv(:)) = vartauConv;
        vY_approx(indtauDiff(:)) = vartauDiff;
        vY_approx(indtauSurf(:)) = vartauSurf;
        vY_approx(indtauInterf(:)) = vartauInterf;
        vY_approx(indenergyKinTime(:)) = varenergyKinTime;
        vY_approx(indenergyConv(:)) = varenergyConv;
        vY_approx(indenergyGrav(:)) = varenergyGrav;
        vY_approx(indenergyPres(:)) = varenergyPres;
        vY_approx(indenergyPresDil(:)) = varenergyPresDil;
        vY_approx(indenergyKinSpace(:)) = varenergyKinSpace;
        vY_approx(indenergyDiff(:)) = varenergyDiff;
        vY_approx(indenergyVisc(:)) = varenergyVisc;
        vY_approx(indenergySurf(:)) = varenergySurf;
        
        % CY = cov(Yc(:,:)); % CY = 1/(N-1)*Yc(:,:)'*Yc(:,:);
        % vY = diag(CY);
        vY = 1/(N-1)*sum(Yc(:,:).^2)';
        errvY = norm(vY_approx-vY)/norm(vY);
        fprintf('\nerror = %.3e for VY',errvY);
        fprintf('\n');
    end
    
    mI = shiftdim(mean(I,1),1);
    mItauTime = mI(1:3,:,:);
    mItauConv = mI(4:6,:,:);
    mItauDiff = mI(7:9,:,:);
    mItauSurf = mI(10:12,:,:);
    mItauInterf = mI(13,:,:);
    mIenergyKinTime = mI(14,:,:);
    mIenergyConv = mI(15,:,:);
    mIenergyGrav = mI(16,:,:);
    mIenergyPres = mI(17,:,:);
    mIenergyPresDil = mI(18,:,:);
    mIenergyKinSpace = mI(19,:,:);
    mIenergyDiff = mI(20,:,:);
    mIenergyVisc = mI(21,:,:);
    mIenergySurf = mI(22,:,:);
    
    sz = size(I);
    sI = prod(sz(2:end-1));
    CI = zeros(sI,sI,p+1,p+1);
    for t=0:p
        It = reshape(I(:,:,:,t+1),[sz(1) sI]);
        for s=0:p
            Is = reshape(I(:,:,:,s+1),[sz(1) sI]);
            CI(:,:,t+1,s+1) = 1/(N-1)*It'*Is;
        end
    end
    
    t = (0:p)*dt;
    
    figure('Name','Mean of spatial average of tau time')
    clf
    hdl(1) = plot(t,squeeze(mItauTime(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
    hold on
    hdl(2) = plot(t,squeeze(mItauTime(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
    hdl(3) = plot(t,squeeze(mItauTime(2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
    hdl(4) = plot(t,squeeze(mItauTime(2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
    hdl(5) = plot(t,squeeze(mItauTime(3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
    hdl(6) = plot(t,squeeze(mItauTime(3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$t$ [s]','Interpreter','latex')
    ylabel('$\tau_{\mathrm{time}}$','Interpreter','latex')
    leg = {'component 1 in phase 1','component 1 in phase 2','component 2 in phase 1','component 2 in phase 2','component 3 in phase 1','component 3 in phase 2'};
    legend(leg{:},'Location','NorthEast')
    mysaveas(gridpathname,'mean_tau_time',formats,renderer);
    mymatlab2tikz(gridpathname,'mean_tau_time.tex');
    
    figure('Name','Mean of spatial average of tau conv')
    clf
    hdl(1) = plot(t,squeeze(mItauConv(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
    hold on
    hdl(2) = plot(t,squeeze(mItauConv(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
    hdl(3) = plot(t,squeeze(mItauConv(2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
    hdl(4) = plot(t,squeeze(mItauConv(2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
    hdl(5) = plot(t,squeeze(mItauConv(3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
    hdl(6) = plot(t,squeeze(mItauConv(3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$t$ [s]','Interpreter','latex')
    ylabel('$\tau_{\mathrm{conv}}$','Interpreter','latex')
    leg = {'component 1 in phase 1','component 1 in phase 2','component 2 in phase 1','component 2 in phase 2','component 3 in phase 1','component 3 in phase 2'};
    legend(leg{:},'Location','NorthEast')
    mysaveas(gridpathname,'mean_tau_conv',formats,renderer);
    mymatlab2tikz(gridpathname,'mean_tau_conv.tex');
    
    figure('Name','Mean of spatial average of tau diff')
    clf
    hdl(1) = plot(t,squeeze(mItauDiff(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
    hold on
    hdl(2) = plot(t,squeeze(mItauDiff(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
    hdl(3) = plot(t,squeeze(mItauDiff(2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
    hdl(4) = plot(t,squeeze(mItauDiff(2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
    hdl(5) = plot(t,squeeze(mItauDiff(3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
    hdl(6) = plot(t,squeeze(mItauDiff(3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$t$ [s]','Interpreter','latex')
    ylabel('$\tau_{\mathrm{diff}}$','Interpreter','latex')
    leg = {'component 1 in phase 1','component 1 in phase 2','component 2 in phase 1','component 2 in phase 2','component 3 in phase 1','component 3 in phase 2'};
    legend(leg{:},'Location','NorthEast')
    mysaveas(gridpathname,'mean_tau_diff',formats,renderer);
    mymatlab2tikz(gridpathname,'mean_tau_diff.tex');
    
    figure('Name','Mean of spatial average of tau surf')
    clf
    hdl(1) = plot(t,squeeze(mItauSurf(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
    hold on
    hdl(2) = plot(t,squeeze(mItauSurf(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
    hdl(3) = plot(t,squeeze(mItauSurf(2,1,:)),'LineStyle','-','Color','r','LineWidth',1);
    hdl(4) = plot(t,squeeze(mItauSurf(2,2,:)),'LineStyle','--','Color','r','LineWidth',1);
    hdl(5) = plot(t,squeeze(mItauSurf(3,1,:)),'LineStyle','-','Color','g','LineWidth',1);
    hdl(6) = plot(t,squeeze(mItauSurf(3,2,:)),'LineStyle','--','Color','g','LineWidth',1);
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$t$ [s]','Interpreter','latex')
    ylabel('$\tau_{\mathrm{surf}}$','Interpreter','latex')
    leg = {'component 1 in phase 1','component 1 in phase 2','component 2 in phase 1','component 2 in phase 2','component 3 in phase 1','component 3 in phase 2'};
    legend(leg{:},'Location','NorthEast')
    mysaveas(gridpathname,'mean_tau_surf',formats,renderer);
    mymatlab2tikz(gridpathname,'mean_tau_surf.tex');
    
    figure('Name','Mean of spatial average of tau interf')
    clf
    hdl(1) = plot(t,squeeze(mItauInterf(1,1,:)),'LineStyle','-','Color','b','LineWidth',1);
    hold on
    hdl(2) = plot(t,squeeze(mItauInterf(1,2,:)),'LineStyle','--','Color','b','LineWidth',1);
    hold off
    grid on
    box on
    set(gca,'FontSize',10)
    xlabel('$t$ [s]','Interpreter','latex')
    ylabel('$\tau_{\mathrm{interf}}$','Interpreter','latex')
    leg = {'phase 1','phase 2'};
    legend(leg{:},'Location','NorthEast')
    mysaveas(gridpathname,'mean_tau_interf',formats,renderer);
    mymatlab2tikz(gridpathname,'mean_tau_interf.tex');
    
    % Display solution
    if displaySolution
        t = 11;
        ampl = 0;
        % ampl = getsize(M)/max(abs(mut))/5;
        for i=1:3
            evolSolution(M,mu,'displ',i,'colormap',cmap,'filename',['evol_mean_u' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
                'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            evolSolution(M,mtauTime,'displ',i,'colormap',cmap,'filename',['evol_mean_tauTime' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
                'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            evolSolution(M,mtauConv,'displ',i,'colormap',cmap,'filename',['evol_mean_tauConv' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
                'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            evolSolution(M,mtauDiff,'displ',i,'colormap',cmap,'filename',['evol_mean_tauDiff' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
                'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            evolSolution(M,mtauSurf,'displ',i,'colormap',cmap,'filename',['evol_mean_tauSurf' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
                'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            
%             evolSolution(M,vu,'displ',i,'colormap',cmap,'filename',['evol_var_u' num2str(i)],'pathname',gridpathname,'FrameRate',framerate);,...
%                 'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%             evolSolution(M,vtauTime,'displ',i,'colormap',cmap,'filename',['evol_var_tauTime' num2str(i)],'pathname',gridpathname,'FrameRate',framerate);,...
%                 'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%             evolSolution(M,vtauConv,'displ',i,'colormap',cmap,'filename',['evol_var_tauConv' num2str(i)],'pathname',gridpathname,'FrameRate',framerate);,...
%                 'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%             evolSolution(M,vtauDiff,'displ',i,'colormap',cmap,'filename',['evol_var_tauDiff' num2str(i)],'pathname',gridpathname,'FrameRate',framerate);,...
%                 'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%             evolSolution(M,vtauSurf,'displ',i,'colormap',cmap,'filename',['evol_var_tauSurf' num2str(i)],'pathname',gridpathname,'FrameRate',framerate);,...
%                 'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            
            figure('Name',['Mean of velocity u' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mu,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(gridpathname,['mean_u' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau time ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauTime,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(gridpathname,['mean_tauTime' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau conv ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauConv,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(gridpathname,['mean_tauConv' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau diff ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauDiff,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(gridpathname,['mean_tauDiff' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau surf ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauSurf,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(gridpathname,['mean_tauSurf' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
%             figure('Name',['Variance of velocity u' num2str(i)])
%             clf
%             plot_sol(M,getmatrixatstep(vu,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
%             colormap(cmap)
%             colorbar
%             axis on
%             box on
%             set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%             mysaveas(gridpathname,['var_u' num2str(i) '_t' num2str(t*100)],formats,renderer);
%             
%             figure('Name',['Variance of tau time ' num2str(i)])
%             clf
%             plot_sol(M,getmatrixatstep(vtauTime,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
%             colormap(cmap)
%             colorbar
%             axis on
%             box on
%             set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%             mysaveas(gridpathname,['var_tauTime' num2str(i) '_t' num2str(t*100)],formats,renderer);
%             
%             figure('Name',['Variance of tau conv ' num2str(i)])
%             clf
%             plot_sol(M,getmatrixatstep(vtauConv,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
%             colormap(cmap)
%             colorbar
%             axis on
%             box on
%             set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%             mysaveas(gridpathname,['var_tauConv' num2str(i) '_t' num2str(t*100)],formats,renderer);
%             
%             figure('Name',['Variance of tau diff ' num2str(i)])
%             clf
%             plot_sol(M,getmatrixatstep(vtauDiff,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
%             colormap(cmap)
%             colorbar
%             axis on
%             box on
%             set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%             mysaveas(gridpathname,['var_tauDiff' num2str(i) '_t' num2str(t*100)],formats,renderer);
%             
%             figure('Name',['Variance of tau surf ' num2str(i)])
%             clf
%             plot_sol(M,getmatrixatstep(vtauSurf,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
%             colormap(cmap)
%             colorbar
%             axis on
%             box on
%             set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%             mysaveas(gridpathname,['var_tauSurf' num2str(i) '_t' num2str(t*100)],formats,renderer);
        end
        
        Mscal = final(M,DDL('C'));
        
        evolSolution(Mscal,mC,'colormap',cmap,'filename','evol_mean_C','pathname',gridpathname,'FrameRate',framerate,...
            'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        evolSolution(Mscal,mtauInterf,'colormap',cmap,'filename','evol_mean_tauInterf','pathname',gridpathname,'FrameRate',framerate,...
            'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        
%         evolSolution(Mscal,vC,'colormap',cmap,'filename','evol_var_C','pathname',gridpathname,'FrameRate',framerate,...
%             'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%         evolSolution(Mscal,vtauInterf,'colormap',cmap,'filename','evol_var_tauInterf','pathname',gridpathname,'FrameRate',framerate,...
%             'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        
        figure('Name','Mean of indicator function C')
        clf
        plot_sol(Mscal,getmatrixatstep(mC,t+1));
        title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(gridpathname,['mean_C_t' num2str(t*100)],formats,renderer);
        
        figure('Name','Mean of tau interf')
        clf
        plot_sol(Mscal,getmatrixatstep(mtauInterf,t+1));
        title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(gridpathname,['mean_tauInterf_t' num2str(t*100)],formats,renderer);
        
%         figure('Name','Variance of indicator function C')
%         clf
%         plot_sol(Mscal,getmatrixatstep(vC,t+1));
%         title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
%         colormap(cmap)
%         colorbar
%         box on
%         set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%         mysaveas(gridname,['var_C_t' num2str(t*100)],formats,renderer);
%         
%         figure('Name','Variance of tau interf')
%         clf
%         plot_sol(Mscal,getmatrixatstep(vtauInterf,t+1));
%         title(['time ' num2str(t*dt,'%f') ' s'],'FontSize',fontsize)
%         colormap(cmap)
%         colorbar
%         box on
%         set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%         mysaveas(gridname,['var_tauInterf_t' num2str(t*100)],formats,renderer);
    end
    
    for t=0:p
        mut = getmatrixatstep(mu,t+1);
        mCt = getmatrixatstep(mC,t+1);
        mtauTimet = getmatrixatstep(mtauTime,t+1);
        mtauConvt = getmatrixatstep(mtauConv,t+1);
        mtauDifft = getmatrixatstep(mtauDiff,t+1);
        mtauSurft = getmatrixatstep(mtauSurf,t+1);
        mtauInterft = getmatrixatstep(mtauInterf,t+1);
        menergyKinTimet = getmatrixatstep(menergyKinTime,t+1);
        menergyConvt = getmatrixatstep(menergyConv,t+1);
        menergyGravt = getmatrixatstep(menergyGrav,t+1);
        menergyPrest = getmatrixatstep(menergyPres,t+1);
        menergyPresDilt = getmatrixatstep(menergyPresDil,t+1);
        menergyKinSpacet = getmatrixatstep(menergyKinSpace,t+1);
        menergyDifft = getmatrixatstep(menergyDiff,t+1);
        menergyVisct = getmatrixatstep(menergyVisc,t+1);
        menergySurft = getmatrixatstep(menergySurf,t+1);
        
        vut = getmatrixatstep(vu,t+1);
        vCt = getmatrixatstep(vC,t+1);
        vtauTimet = getmatrixatstep(vtauTime,t+1);
        vtauConvt = getmatrixatstep(vtauConv,t+1);
        vtauDifft = getmatrixatstep(vtauDiff,t+1);
        vtauSurft = getmatrixatstep(vtauSurf,t+1);
        vtauInterft = getmatrixatstep(vtauInterf,t+1);
        venergyKinTimet = getmatrixatstep(venergyKinTime,t+1);
        venergyConvt = getmatrixatstep(venergyConv,t+1);
        venergyGravt = getmatrixatstep(venergyGrav,t+1);
        venergyPrest = getmatrixatstep(venergyPres,t+1);
        venergyPresDilt = getmatrixatstep(venergyPresDil,t+1);
        venergyKinSpacet = getmatrixatstep(venergyKinSpace,t+1);
        venergyDifft = getmatrixatstep(venergyDiff,t+1);
        venergyVisct = getmatrixatstep(venergyVisc,t+1);
        venergySurft = getmatrixatstep(venergySurf,t+1);
        
        write_vtk_mesh(M,{mut,mCt,mtauTimet,mtauConvt,mtauDifft,mtauSurft,mtauInterft,...
            menergyKinTimet,menergyConvt,menergyGravt,menergyPrest,menergyPresDilt,...
            menergyKinSpacet,menergyDifft,menergyVisct,menergySurft},[],...
            {'velocity','phase','tau time','tau conv','tau diff','tau surf','tau interf',...
            'kinetic energy','convection energy','gravity energy','power of external pressure forces','pressure-dilatation energy transfer',...
            'transport of gradient of kinetic energy','energy exchange with kinetic energy','power of external viscous stresses','capillary kinetic energy'},[],...
            gridpathname,'diphasic_fluids_mean',1,t);
        write_vtk_mesh(M,{vut,vCt,vtauTimet,vtauConvt,vtauDifft,vtauSurft,vtauInterft,...
            venergyKinTimet,venergyConvt,venergyGravt,venergyPrest,venergyPresDilt,...
            venergyKinSpacet,venergyDifft,venergyVisct,venergySurft},[],...
            {'velocity','phase','tau time','tau conv','tau diff','tau surf','tau interf',...
            'kinetic energy','convection energy','gravity energy','power of external pressure forces','pressure-dilatation energy transfer',...
            'transport of gradient of kinetic energy','energy exchange with kinetic energy','power of external viscous stresses','capillary kinetic energy'},[],...
            gridpathname,'diphasic_fluids_variance',1,t);
    end
    make_pvd_file(gridpathname,'diphasic_fluids_mean',1,p+1);
    make_pvd_file(gridpathname,'diphasic_fluids_variance',1,p+1);
    
    toc

end
