clc
clearvars
close all

solveProblem = true;
postProcess = true;
displaySolution = false;
displayEigenvales = true;
displayCovariance  = false;

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

% Time scheme
dt = 5e-3; % time step (s)
dt = 100*dt; % physical time step stored every 100 time iterations

order = [3 1 2]; % dimension order

tolsvdYc = eps; % relative precision for truncated SVD of Yc
tolsvdZc = 1e-6; % relative precision for truncated SVD of Zc

% for g=2.^(4:8)
for g=2^4
    tic
    gridname = ['Grid' num2str(g)];
    disp(gridname)
    gridpathname = fullfile(pathname,gridname);
    load(fullfile(gridpathname,'data.mat'),'N','n','m','p');
    if ~postProcess
        load(fullfile(gridpathname,'data_post.mat'),'nn');
        n = n+nn;
    end
    r = n*m;
    
    fprintf('\nn = %d variables',n);
    fprintf('\nN = %d samples',N);
    fprintf('\nm = %d spatial points',m);
    fprintf('\np+1 = %d time steps',p+1);
    fprintf('\n');
    
    s = [g+1,g+1,g+1]; % spatial dimensions
    
    % Spatial scheme
    dx = L/g; % spatial step (m)
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
        fprintf('\nFirst PCA');
        Rinit = min(r,N);
        if g<2^7
            load(fullfile(gridpathname,'data.mat'),'Y');
            if ~postProcess
                load(fullfile(gridpathname,'data_post.mat'),'YY');
                Y = cat(2,Y,YY);
                clear YY
            end
            mY = mean(Y,1);
            Yc = Y - repmat(mY,N,1,1,1); % Yc = Y - mY.*ones(N,1,1,1);
            clear Y
            Sig = zeros(Rinit,p+1);
            V = zeros(r,Rinit,p+1);
            R = zeros(p+1,1);
            errsvdYc = zeros(Rinit,p+1);
        else
            mY = zeros(1,n,m,p+1);
        end
        Z = zeros(N,Rinit,p+1);
        
        Rmax = 1;
        for t=0:p
            if g<2^7
                Yct = Yc(:,:,:,t+1);
            else
                load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
                if ~postProcess
                    load(fullfile(gridpathname,['data_post_t' num2str(t) '.mat']),'YYt');
                    Yt = cat(2,Yt,YYt);
                    clear YYt
                end
                mYt = mean(Yt,1);
                mY(1,:,:,t+1) = mYt;
                Yct = Yt - repmat(mYt,N,1,1); % Yct = Yt - mYt.*ones(N,1,1);
                clear Yt
            end
            Yct = Yct(:,:)';
            % if t==0
            %     [Vt,Sigt,Zt,errt] = svdtruncate(Yct,Rinit-1);
            % else
                [Vt,Sigt,Zt,errsvdYct] = svdtruncate(Yct,tolsvdYc);
            % end
            Sigt = Sigt/sqrt(N-1);
            Zt = Zt*sqrt(N-1);
            Yct_approx = Vt*(Sigt.*Zt'); % Yct_approx = Vt*diag(Sigt)*Zt';
            errYct = norm(Yct_approx-Yct)/norm(Yct);
            Rt = length(Sigt);
            Rmax = max(Rmax,Rt);
            fprintf('\nTime %2d, t = %4g s : rank R = %d, error = %.3e for Y',t,t*100*dt,Rt,errYct);
            %CYt_approx = Vt*(Sigt.^2.*Vt'); %CYt_approx = Vt*diag(Sigt).^2*Vt';
            %CYt = cov(Yct'); % CYt = 1/(N-1)*Yct*Yct';
            %errCYt = norm(CYt_approx-CYt)/norm(CYt);
            %fprintf('\n                                                  error = %.3e for CY',errCYt);
            
            if g<2^7
                Sig(1:Rt,t+1) = Sigt;
                V(:,1:Rt,t+1) = Vt;
                R(t+1) = Rt;
                errsvdYc(1:Rt,t+1) = errsvdYct;
            else
                save(fullfile(gridpathname,['PCA_t' num2str(t) '.mat']),'Sigt','Vt','Rt','errsvdYct');
            end
            Z(:,1:Rt,t+1) = Zt;
            
            % mZt = mean(Zt,1)';
            % CZt = cov(Zt); % CZt = 1/(N-1)*Zt'*Zt;
            % norm(mZt)
            % norm(CZt-eye(Rt))
            % norm(Vt'*Vt-eye(Rt))
        end
        fprintf('\n');
        if g<2^7
            Sig = Sig(1:Rmax,:);
            V = V(:,1:Rmax,:);
            errsvdYc = errsvdYc(1:Rmax,:);
        end
        Z = Z(:,1:Rmax,:);
        
        %% Second reduction step for each coordinate
%         fprintf('\nSecond PCA');
%         Q = min(p+1,N);
%         S = zeros(Q,Rmax);
%         W = zeros(p+1,Q,Rmax);
%         X = zeros(N,Q,Rmax);
%         errsvdZc = zeros(Q,Rmax);
%         Q = zeros(1,Rmax);
%         Qmax = 1;
%         for a=1:Rmax
%             Zca = Z(:,a,:);
%             Zca = Zca(:,:)';
%             [Wa,Sa,Xa,errsvdZca] = svdtruncate(Zca,tolsvdZc);
%             Sa = Sa/sqrt(N-1);
%             Xa = Xa*sqrt(N-1);
%             Zca_approx = Wa*(Sa.*Xa'); % Zca_approx = Wa*diag(Sa)*Xa';
%             errZa = norm(Zca_approx-Zca)/norm(Zca);
%             Qa = length(Sa);
%             Qmax = max(Qmax,Qa);
%             fprintf('\nCoordinate alpha = %2.f : rank Q = %d, error = %.3e for Z',a,Qa,errZa);
%             S(1:Qa,a) = Sa;
%             W(:,1:Qa,a) = Wa;
%             X(:,1:Qa,a) = Xa;
%             Q(a) = Qa;
%             errsvdZc(1:Qa,a) = errsvdZca;
%             
%             CZa_approx = Wa*(Sa.^2.*Wa'); % CZa_approx = Wa*diag(Sa).^2*Wa';
%             CZa = cov(Zca'); % CZa = 1/(N-1)*Zca*Zca';
%             errCZa = norm(CZa_approx-CZa)/norm(CZa);
%             fprintf('\n                                     error = %.3e for CZ',errCZa);
%             
%             % mXa = mean(Xa,1)';
%             % CXa = cov(Xa); % CXa = 1/(N-1)*Xa'*Xa;
%             % norm(mXa)
%             % norm(CXa-eye(Qa))
%             % norm(Wa'*Wa-eye(Qa))
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
%         end
%         fprintf('\n');
%         S = S(1:Qmax,:);
%         W = W(:,1:Qmax,:);
%         X = X(:,1:Qmax,:);
%         errsvdZc = errsvdZc(1:Qmax,:);
%         Q = Qmax;
        
        %% Second reduction step
        fprintf('\nSecond PCA');
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
        X = X*sqrt(N-1);
        Zc_approx = W*(S.*X'); % Zc_approx = W*diag(S)*X';
        errZ = norm(Zc_approx-Zc)/norm(Zc);
        Q = length(S);
        fprintf('\nrank R = %d, rank Q = %d, error = %.3e for Z',Rmax,Q,errZ);
        
        CZ_approx = W*(S.^2.*W'); % CZ_approx = W*diag(S).^2*W';
        CZ = cov(Zc'); % CZ = 1/(N-1)*Zc*Zc';
        errCZ = norm(CZ_approx-CZ)/norm(CZ);
        fprintf('\n                          error = %.3e for CZ',errCZ);
        fprintf('\n');
        
        % mX = mean(X,1)';
        % CX = cov(X); % CX = 1/(N-1)*X'*X;
        % norm(mX)
        % norm(CX-eye(Q))
        % norm(W'*W-eye(Q))
        % norm(CZ_approx*W-(p+1)*W)
        % norm(CZ*W-(p+1)*W)
        % norm(abs(S.^2-(p+1)))
        
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
        
        if g<2^7
            save(fullfile(gridpathname,'solution.mat'),'mY','Sig','S','V','W','Z','Rmax','R','Q','errsvdYc','errsvdZc');
        else
            save(fullfile(gridpathname,'solution.mat'),'mY','S','W','Z','Rmax','Q','errsvdYc','errsvdZc');
        end
    else
        if g<2^7
            load(fullfile(gridpathname,'solution.mat'),'mY','Sig','S','V','W','Z','Rmax','R','Q','errsvdYc','errsvdZc');
        else
            load(fullfile(gridpathname,'solution.mat'),'mY','S','W','Z','Rmax','Q','errsvdYc','errsvdZc');
        end
    end
    
    %% Post-processing
    if postProcess
        %% Post-processing data
        fprintf('\nPost-processing data');
        nn = 13; % number of post-processed variables
        fprintf('\nn = %d post-processed variables',nn);
        fprintf('\n');
        
        if g<2^7
            YY = zeros(N,nn,m,p+1);
        end
        
        for t=0:p
            time = ['Time ' num2str(t)];
            disp(time)
            
            mYt = mY(1,:,:,t+1);
            mYt_old = zeros(1,n,m);
            mYt_new = zeros(1,n,m);
            Rt_old = Rmax;
            Rt_new = Rmax;
            Sigt_old = zeros(Rmax,1);
            Sigt_new = zeros(Rmax,1);
            Vt_old = zeros(n*m,Rmax);
            Vt_new = zeros(n*m,Rmax);
            Zt_old = zeros(N,Rmax,1);
            Zt_new = zeros(N,Rmax,1);
            if t>0
                mYt_old = mY(1,:,:,t);
            end
            if t<p
                mYt_new = mY(1,:,:,t+2);
            end
            if g<2^7
                Rt = R(t+1);
                Sigt = Sig(1:Rt,t+1);
                Vt = V(:,1:Rt,t+1);
                if t>0
                    Rt_old = R(t);
                    Sigt_old = Sig(1:Rt_old,t);
                    Vt_old = V(:,1:Rt_old,t);
                end
                if t<p
                    Rt_new = R(t+2);
                    Sigt_new = Sig(1:Rt_new,t+2);
                    Vt_new = V(:,1:Rt_new,t+2);
                end
            else
                if t>0
                    load(fullfile(gridpathname,['PCA_t' num2str(t-1) '.mat']),'Sigt','Vt','Rt');
                    Rt_old = Rt;
                    Sigt_old = Sigt;
                    Vt_old = Vt;
                end
                if t<p
                    load(fullfile(gridpathname,['PCA_t' num2str(t+1) '.mat']),'Sigt','Vt','Rt');
                    Rt_new = Rt;
                    Sigt_new = Sigt;
                    Vt_new = Vt;
                end
                load(fullfile(gridpathname,['PCA_t' num2str(t) '.mat']),'Sigt','Vt','Rt');
            end
            Zt = Z(:,1:Rt,t+1);
            if t>0
                Zt_old = Z(:,1:Rt_old,t);
            end
            if t<p
                Zt_new = Z(:,1:Rt_new,t+2);
            end
            
            Yct = Vt*(Sigt.*Zt');
            Yct_old = Vt_old*(Sigt_old.*Zt_old');
            Yct_new = Vt*(Sigt_new.*Zt_new');
            
            Yt = reshape(mYt(:)' + Yct',[N,n,m]);
            Yt_old = reshape(mYt_old(:)' + Yct_old',[N,n,m]);
            Yt_new = reshape(mYt_new(:)' + Yct_new',[N,n,m]);
            
            ut = Yt(:,1:3,:);
            ut_old = Yt_old(:,1:3,:);
            ut_new = Yt_new(:,1:3,:);
            Ct = Yt(:,4,:);
            Ct_old = Yt_old(:,4,:);
            Ct_new = Yt_new(:,4,:);
            rhot = Ct*rho(2) + (1-Ct)*rho(1);
            rhot_old = Ct_old*rho(2) + (1-Ct_old)*rho(1);
            rhot_new = Ct_new*rho(2) + (1-Ct_new)*rho(1);
            rhout = rhot.*ut;
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
            
%             mYt = mean(YYt,1);
%             
%             mut = mYt(1,1:3,:);
%             mut_old = mYt_old(1,1:3,:);
%             mut_new = mYt_new(1,1:3,:);
%             mCt = mYt(1,4,:);
%             mCt_old = mYt_old(1,4,:);
%             mCt_new = mYt_new(1,4,:);
%             
%             Vt = reshape(Vt',[Rt,n,m]);
%             Vt_old = reshape(Vt_old',[Rt_old,n,m]);
%             Vt_new = reshape(Vt_new',[Rt_new,n,m]);
%             Vut = Vt(:,1:3,:);
%             Vut_old = Vt_old(:,1:3,:);
%             Vut_new = Vt_new(:,1:3,:);
%             VCt = Vt(:,4,:);
%             VCt_old = Vt_old(:,4,:);
%             VCt_new = Vt_new(:,4,:);
%             mrhot = (mCt*rho(2) + (1-mCt)*rho(1));
%             mrhot_old = (mCt_old*rho(2) + (1-mCt_old)*rho(1));
%             mrhot_new = (mCt_new*rho(2) + (1-mCt_new)*rho(1));
%             mrhout  = mrhot.*mut + (rho(2)-rho(1))*Vut.*reshape(Sigt.^2.*VCt(:,:),[Rt,1,m]);
%             mrhout_old  = mrhot_old.*mut_old + (rho(2)-rho(1))*Vut_old.*reshape(Sigt_old.^2.*VCt_old(:,:),[Rt,1,m]);
%             mrhout_new  = mrhot_new.*mut_new + (rho(2)-rho(1))*Vut_new.*reshape(Sigt_new.^2.*VCt_new(:,:),[Rt,1,m]);
%             mtauTime = (mrhout_new-mrhout_old)/(2*dt);
%             
%             mut = permute(reshape(mut,[1,3,s]),[2,order+2,1]);
%             mCt = permute(reshape(mCt,[1,s]),[order+1,1]);
%             gradmut = grad(mut,Dx);
%             mrhot = (mCt*rho(2) + (1-mCt)*rho(1));
%             mmut = (mCt*mu(2) + (1-mCt)*mu(1));
%             
%             Vut = permute(reshape(Vut,[Rt,3,s]),[2,order+2,1]);
%             VCt = permute(reshape(VCt,[Rt,1,s]),[2,order+2,1]);
%             gradVCVut = grad(VCt.*Vut,Dx);
%             gradVCmut = grad(VCt.*mut,Dx);
%             gradmrhoVut = grad(shiftdim(mrhot,-1).*Vut,Dx);
%             gradVCVumut = reshape(compute_tauConv(mut,gradVCVut),[3*m,Rt])';
%             gradVCmuVut = reshape(compute_tauConv(Vut,gradVCmut),[3*m,Rt])';
%             gradmrhoVuVut = reshape(compute_tauConv(Vut,gradmrhoVut),[3*m,Rt])';
%             mtauConv = shiftdim(mrhot,-1).*compute_tauConv(mut,gradmut) + (rho(2)-rho(1))*reshape(sum(Sigt.*(gradVCVumut+gradVCVumut),1),[3,s])...
%                 + reshape(sum(Sigt.*gradmrhoVuVut,1),[3,s]);
            
%             tTime = yy(:,1:3,:,t+1);
%             tConv = yy(:,4:6,:,t+1);
%             tDiff = yy(:,7:9,:,t+1);
%             tSurf = yy(:,10:12,:,t+1);
%             tInterf = yy(:,13,:,t+1);
            
%             yyt = yy(:,:,:,t+1);
%             dY = YYt(:,:)-yyt(:,:);
%             dTime = tTime(:,:)-tauTime(:,:);
%             dConv = tConv(:,:)-tauConv(:,:);
%             dDiff = tDiff(:,:)-tauDiff(:,:);
%             dSurf = tSurf(:,:)-tauSurf(:,:);
%             dInterf = tInterf(:,:)-tauInterf(:,:);
            
%             norm(dY)/norm(yyt(:,:))
%             norm(dTime)/norm(tTime(:,:))
%             norm(dConv)/norm(tConv(:,:))
%             norm(dDiff)/norm(tDiff(:,:))
%             norm(dSurf)/norm(tSurf(:,:))
%             norm(dInterf)/norm(tInterf(:,:))
            
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
    end
    
    %% Outputs
    % Display eigenvalues
    if displayEigenvales
        times = [1,10:10:p];
        plot_eigenvalues_CY(times,dt,5,g,gridpathname,9)
        mysaveas(gridpathname,'eigenvalues_CY',formats,renderer);
        mymatlab2tikz(gridpathname,'eigenvalues_CY.tex');
        
        plot_error_svdYc(times,dt,5,g,gridpathname,9)
        mysaveas(gridpathname,'error_svdYc',formats,renderer);
        mymatlab2tikz(gridpathname,'error_svdYc.tex');
        
        figure('Name','Evolution of eigenvalues')
        clf
        semilogy(1:Q,S(1:Q).^2,'LineStyle','-','Color','b','LineWidth',1);
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
    tf = p*100*dt;
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
        load(fullfile(gridpathname,'data_post.mat'),'YY');
        Y = cat(2,Y,YY);
        clear YY
        mY = mean(Y,1);
        Yc = Y - repmat(mY,N,1,1,1); % Yc = Y - mY.*ones(N,1,1,1);
        clear Y
    end
    
    mu = zeros(3*m,p+1);
    vu = zeros(3*m,p+1);
    mC = zeros(m,p+1);
    vC = zeros(m,p+1);
    mtauTime = zeros(3*m,p+1);
    vtauTime = zeros(3*m,p+1);
    mtauConv = zeros(3*m,p+1);
    vtauConv = zeros(3*m,p+1);
    mtauDiff = zeros(3*m,p+1);
    vtauDiff = zeros(3*m,p+1);
    mtauSurf = zeros(3*m,p+1);
    vtauSurf = zeros(3*m,p+1);
    mtauInterf = zeros(m,p+1);
    vtauInterf = zeros(m,p+1);
    for t=0:p
        mYt = mY(1,:,:,t+1);
        if ~postProcess
            mYt = reshape(mYt,[n,m]);
        else
            mYt = reshape(mYt,[n+nn,m]);
        end
        mut = mYt(1:3,:);
        mCt = mYt(4,:);
        mtauTimet = mYt(5:7,:);
        mtauConvt = mYt(8:10,:);
        mtauDifft = mYt(11:13,:);
        mtauSurft = mYt(14:16,:);
        mtauInterft = mYt(17,:);
        mu(:,t+1) = mut(:);
        mC(:,t+1) = mCt(:);
        mtauTime(:,t+1) = mtauTimet(:);
        mtauConv(:,t+1) = mtauConvt(:);
        mtauDiff(:,t+1) = mtauDifft(:);
        mtauSurf(:,t+1) = mtauSurft(:);
        mtauInterf(:,t+1) = mtauInterft(:);
        
        if g<2^7
            Yct = Yc(:,:,:,t+1);
            Rt = R(t+1);
            Sigt = Sig(1:Rt,t+1);
            Vt = V(:,1:Rt,t+1);
        else
            load(fullfile(gridpathname,['data_t' num2str(t) '.mat']),'Yt');
            load(fullfile(gridpathname,['data_post_t' num2str(t) '.mat']),'YYt');
            Yt = cat(2,Yt,YYt);
            clear YYt
            mYt = mean(Yt,1);
            Yct = Yt - repmat(mYt,N,1,1); % Yct = Yt - mYt.*ones(N,1,1);
            clear Yt
            load(fullfile(gridpathname,['PCA_t' num2str(t) '.mat']),'Sigt','Vt','Rt');
        end
        Wt = W(1:Rt,:,t+1);
        % Zct_approx = Wt*(S.*X'); % Zct_approx = Wt*diag(S)*X'; % Zct_approx = Z_approx(:,1:Rt,t+1)';
        % Yct_approx = Vt*(Sigt.*Zct_approx'); % Yct_approx = Vt*diag(Sigt)*Zct_approx';
        Uct_approx = Vt*(Sigt.*Wt); % Uct_approx = Vt*diag(Sigt)*Wt;
        if g<2^7
            % Yc_approx(:,:,:,t+1) = reshape(Yct_approx',[N,n,m]);
            Uc_approx(:,:,:,t+1) = reshape(Uct_approx',[Q,n,m]);
        end
        
        % CYt_approx = cov(Yct_approx'); % CYt_approx = 1/(N-1)*Yct_approx*Yct_approx';
        % CYt_approx = Uct_approx*(S.^2.*Uct_approx'); % CYt_approx = Uct_approx*diag(S).^2*Uct_approx';
        % vYt_approx = diag(CYt_approx);
        vYt_approx = sum((S.*Uct_approx').^2)'; % vYt_approx = sum((diag(S)*Uct_approx').^2)';
        
        s = n;
        if postProcess
            s = s+nn;
        end
        indut = (0:m-1)*s+(1:3)';
        indCt = 4:s:(s*m); % indCt = (0:m-1)*s+4;
        indtauTimet = (0:m-1)*s+(5:7)';
        indtauConvt = (0:m-1)*s+(8:10)';
        indtauDifft = (0:m-1)*s+(11:13)';
        indtauSurft = (0:m-1)*s+(14:16)';
        indtauInterft = 17:s:(s*m); % indtauInterft = (0:m-1)*s+17;
        
        if postProcess
            tauTimet = Yct(:,5:7,:);
            tauConvt = Yct(:,8:10,:);
            tauDifft = Yct(:,11:13,:);
            tauSurft = Yct(:,14:16,:);
            tauInterft = Yct(:,17,:);
            indut_ = (0:m-1)*n+(1:3)';
            indCt_ = 4:n:n*m;
            vut = vYt_approx(indut_(:));
            vCt = vYt_approx(indCt_(:));
            vtauTimet = 1/(N-1)*sum(tauTimet(:,:).^2)';
            vtauConvt = 1/(N-1)*sum(tauConvt(:,:).^2)';
            vtauDifft = 1/(N-1)*sum(tauDifft(:,:).^2)';
            vtauSurft = 1/(N-1)*sum(tauSurft(:,:).^2)';
            vtauInterft = 1/(N-1)*sum(tauInterft(:,:).^2)';
            vYt_approx = zeros((n+nn)*m,1);
            vYt_approx(indut(:)) = vut;
            vYt_approx(indCt(:)) = vCt;
            vYt_approx(indtauTimet(:)) = vtauTimet;
            vYt_approx(indtauConvt(:)) = vtauConvt;
            vYt_approx(indtauDifft(:)) = vtauDifft;
            vYt_approx(indtauSurft(:)) = vtauSurft;
            vYt_approx(indtauInterft(:))  = vtauInterft;
        end
        
        % CYt = cov(Yct(:,:)); % CYt = 1/(N-1)*Yct(:,:)'*Yct(:,:);
        % vYt = diag(CYt);
        vYt = 1/(N-1)*sum(Yct(:,:).^2)';
        errvYt = norm(vYt_approx-vYt)/norm(vYt);
        fprintf('\nTime %2d, t = %4g s : error = %.3e for VY',t,t*100*dt,errvYt);
        
        vu(:,t+1) = vYt_approx(indut(:));
        vC(:,t+1) = vYt_approx(indCt(:));
        vtauTime(:,t+1) = vYt_approx(indtauTimet(:));
        vtauConv(:,t+1) = vYt_approx(indtauConvt(:));
        vtauDiff(:,t+1) = vYt_approx(indtauDifft(:));
        vtauSurf(:,t+1) = vYt_approx(indtauSurft(:));
        vtauInterf(:,t+1) = vYt_approx(indtauInterft(:));
    end
    fprintf('\n');
    mu = TIMEMATRIX(mu,T);
    mC = TIMEMATRIX(mC,T);
    mtauTime = TIMEMATRIX(mtauTime,T);
    mtauConv = TIMEMATRIX(mtauConv,T);
    mtauDiff = TIMEMATRIX(mtauDiff,T);
    mtauSurf = TIMEMATRIX(mtauSurf,T);
    mtauInterf = TIMEMATRIX(mtauInterf,T);
    vu = TIMEMATRIX(vu,T);
    vC = TIMEMATRIX(vC,T);
    vtauTime = TIMEMATRIX(vtauTime,T);
    vtauConv = TIMEMATRIX(vtauConv,T);
    vtauDiff = TIMEMATRIX(vtauDiff,T);
    vtauSurf = TIMEMATRIX(vtauSurf,T);
    vtauInterf = TIMEMATRIX(vtauInterf,T);
    
    if g<2^7
        % CY_approx = cov(Yc_approx(:,:)); % CY_approx = 1/(N-1)*Yc_approx(:,:)'*Yc_approx(:,:);
        % CY_approx = Uc_approx(:,:)'*(S.^2.*Uc_approx(:,:)); % CY_approx = Uc_approx(:,:)'*diag(S).^2*Uc_approx(:,:);
        % vY_approx = diag(CY_approx);
        vY_approx = sum((S.*Uc_approx(:,:)).^2)'; % vY_approx = sum((diag(S)*Uc_approx(:,:)).^2)';
        if postProcess
            tauTime = permute(Yc(:,5:7,:,:),[1,4,2,3]);
            tauConv = permute(Yc(:,8:10,:,:),[1,4,2,3]);
            tauDiff = permute(Yc(:,11:13,:,:),[1,4,2,3]);
            tauSurf = permute(Yc(:,14:16,:,:),[1,4,2,3]);
            tauInterf = permute(Yc(:,17,:,:),[1,4,2,3]);
            indu_ = (0:m-1)*n+(1:3)';
            indC_ = 4:n:n*m;
            indu_ = indu_(:)'+n*m*(0:p)';
            indC_ = indC_(:)'+n*m*(0:p)';
            varu = vY_approx(indu_(:));
            varC = vY_approx(indC_(:));
            vartauTime = 1/(N-1)*sum(tauTime(:,:).^2)';
            vartauConv = 1/(N-1)*sum(tauConv(:,:).^2)';
            vartauDiff = 1/(N-1)*sum(tauDiff(:,:).^2)';
            vartauSurf = 1/(N-1)*sum(tauSurf(:,:).^2)';
            vartauInterf = 1/(N-1)*sum(tauInterf(:,:).^2)';
            vY_approx = zeros((n+nn)*m*(p+1),1);
            s = n+nn;
            indu = (0:m-1)*s+(1:3)';
            indC = 4:s:(s*m); % indC = (0:m-1)*s+4;
            indtauTime = (0:m-1)*s+(5:7)';
            indtauConv = (0:m-1)*s+(8:10)';
            indtauDiff = (0:m-1)*s+(11:13)';
            indtauSurf = (0:m-1)*s+(14:16)';
            indtauInterf = 17:s:(s*m); % indtauInterf = (0:m-1)*s+17;
            indu = indu(:)'+s*m*(0:p)';
            indC = indC(:)'+s*m*(0:p)';
            indtauTime = indtauTime(:)'+s*m*(0:p)';
            indtauConv = indtauConv(:)'+s*m*(0:p)';
            indtauDiff = indtauDiff(:)'+s*m*(0:p)';
            indtauSurf = indtauSurf(:)'+s*m*(0:p)';
            indtauInterf = indtauInterf(:)'+s*m*(0:p)';
            
            vY_approx(indu(:)) = varu;
            vY_approx(indC(:)) = varC;
            vY_approx(indtauTime(:)) = vartauTime;
            vY_approx(indtauConv(:)) = vartauConv;
            vY_approx(indtauDiff(:)) = vartauDiff;
            vY_approx(indtauSurf(:)) = vartauSurf;
            vY_approx(indtauInterf(:))  = vartauInterf;
        end
        
        % CY = cov(Yc(:,:)); % CY = 1/(N-1)*Yc(:,:)'*Yc(:,:);
        % vY = diag(CY);
        vY = 1/(N-1)*sum(Yc(:,:).^2)';
        errvY = norm(vY_approx-vY)/norm(vY);
        fprintf('\nerror = %.3e for VY',errvY);
        fprintf('\n');
    end
    
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
            title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(gridpathname,['mean_u' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau time ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauTime,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(gridpathname,['mean_tauTime' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau conv ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauConv,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(gridpathname,['mean_tauConv' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau diff ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauDiff,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(gridpathname,['mean_tauDiff' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau surf ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauSurf,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(gridpathname,['mean_tauSurf' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
%             figure('Name',['Variance of velocity u' num2str(i)])
%             clf
%             plot_sol(M,getmatrixatstep(vu,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
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
%             title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
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
%             title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
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
%             title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
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
%             title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
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
        title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(gridpathname,['mean_C_t' num2str(t*100)],formats,renderer);
        
        figure('Name','Mean of tau interf')
        clf
        plot_sol(Mscal,getmatrixatstep(mtauInterf,t+1));
        title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(gridpathname,['mean_tauInterf_t' num2str(t*100)],formats,renderer);
        
%         figure('Name','Variance of indicator function C')
%         clf
%         plot_sol(Mscal,getmatrixatstep(vC,t+1));
%         title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
%         colormap(cmap)
%         colorbar
%         box on
%         set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%         mysaveas(gridname,['var_C_t' num2str(t*100)],formats,renderer);
%         
%         figure('Name','Variance of tau interf')
%         clf
%         plot_sol(Mscal,getmatrixatstep(vtauInterf,t+1));
%         title(['time ' num2str(t*100*dt,'%f') ' s'],'FontSize',fontsize)
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
        vut = getmatrixatstep(vu,t+1);
        vCt = getmatrixatstep(vC,t+1);
        vtauTimet = getmatrixatstep(vtauTime,t+1);
        vtauConvt = getmatrixatstep(vtauConv,t+1);
        vtauDifft = getmatrixatstep(vtauDiff,t+1);
        vtauSurft = getmatrixatstep(vtauSurf,t+1);
        vtauInterft = getmatrixatstep(vtauInterf,t+1);
        write_vtk_mesh(M,{mut,mCt,mtauTimet,mtauConvt,mtauDifft,mtauSurft,mtauInterft},[],...
            {'velocity','phase','tau time','tau conv','tau diff','tau surf','tau interf'},[],...
            gridpathname,'diphasic_fluids_mean',1,t);
        write_vtk_mesh(M,{vut,vCt,vtauTimet,vtauConvt,vtauDifft,vtauSurft,vtauInterft},[],...
            {'velocity','phase','tau time','tau conv','tau diff','tau surf','tau interf'},[],...
            gridpathname,'diphasic_fluids_variance',1,t);
    end
    make_pvd_file(gridpathname,'diphasic_fluids_mean',1,p+1);
    make_pvd_file(gridpathname,'diphasic_fluids_variance',1,p+1);
    
    toc

end
