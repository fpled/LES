clc
clearvars
close all

solveProblem = true;
displaySolution = false;
displayEigenvales = false;
displayCovariance  = false;

% index = 'time';
index = 'coord';

cmap = 'default';
% cmap = 'gray';
framerate = 5;
fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc2'};
renderer = 'OpenGL';

pathname = fileparts(mfilename('fullpath'));

for g=2.^(4:8)
% for g=2^4
    gridname = ['Grid' num2str(g)];
    disp(gridname)
    pathnamegrid = fullfile(pathname,gridname);
    load(fullfile(pathnamegrid,'data.mat'),'N','n','m','p');
    load(fullfile(pathnamegrid,'data_post.mat'),'nn');
    n = n+nn;
    r = n*m;
    
    fprintf('\nn = %d variables',n);
    fprintf('\nN = %d samples',N);
    fprintf('\nm = %d spatial points',m);
    fprintf('\np+1 = %d time steps',p+1);
    fprintf('\n');
    
    tol = eps;
    
    if solveProblem
        %% First reduction step
        if g<2^7
            load(fullfile(pathnamegrid,'data.mat'),'Y');
            load(fullfile(pathnamegrid,'data_post.mat'),'YY');
            Y = cat(2,Y,YY);
            clear YY
            mY = mean(Y,1);
            Yc = Y - repmat(mY,N,1,1,1); % Yc = Y - mY.*ones(N,1,1,1);
            clear Y
        else
            mY = zeros(1,n,m,p+1);
        end
        Rinit = min(r,N);
        if g<2^7
            Sig = zeros(Rinit,p+1);
            V = zeros(r,Rinit,p+1);
            R = zeros(p+1,1);
        end
        Z = zeros(N,Rinit,p+1);
        
        Rmax = 1;
        for t=0:p
            if g<2^7
                Yct = Yc(:,:,:,t+1);
            else
                load(fullfile(pathnamegrid,['data_t' num2str(t) '.mat']),'Yt');
                load(fullfile(pathnamegrid,['data_post_t' num2str(t) '.mat']),'YYt');
                Yt = cat(2,Yt,YYt);
                clear YYt
                mYt = mean(Yt,1);
                mY(1,:,:,t+1) = mYt;
                Yct = Yt - repmat(mYt,N,1,1); % Yct = Yt - mYt.*ones(N,1,1);
                clear Yt
            end
            Yct = Yct(:,:)';
            if t==0
                [Vt,Sigt,Zt] = svdtruncate(Yct,Rinit-1);
            else
                [Vt,Sigt,Zt] = svdtruncate(Yct,tol);
            end
            Sigt = Sigt/sqrt(N-1);
            Zt = Zt*sqrt(N-1);
            Yct_approx = Vt*diag(Sigt)*Zt';
            errYct = norm(Yct_approx-Yct)/norm(Yct);
            Rt = length(Sigt);
            Rmax = max(Rmax,Rt);
            fprintf('\nTime t = %4.f s : rank R = %d, error = %.3e for Y',t*100,Rt,errYct);
            if g<2^7
                Sig(1:Rt,t+1) = Sigt;
                V(:,1:Rt,t+1) = Vt;
                R(t+1) = Rt;
            else
                save(fullfile(pathnamegrid,['PCA_t' num2str(t) '.mat']),'Sigt','Vt','Rt');
            end
            Z(:,1:Rt,t+1) = Zt;
            
%             mZt = mean(Zt,1)';
%             CZt = cov(Zt); % CZt = 1/(N-1)*Zt(:,:)'*Zt(:,:);
%             norm(mZt)
%             norm(CZt-eye(Rt))
%             norm(Vt'*Vt-eye(Rt))
        end
        fprintf('\n');
        if g<2^7
            Sig = Sig(1:Rmax,:);
            V = V(:,1:Rmax,:);
        end
        Z = Z(:,1:Rmax,:);
        
        %% Second reduction step for each coordinate
%         Q = min(p+1,N);
%         S = zeros(Q,Rmax);
%         W = zeros(p+1,Q,Rmax);
%         X = zeros(N,Q,Rmax);
%         Q = zeros(1,Rmax);
%         Qmax = 1;
%         for a=1:Rmax
%             Zca = Z(:,a,:);
%             Zca = Zca(:,:)';
%             [Wa,Sa,Xa] = svdtruncate(Zca,tol);
%             Sa = Sa/sqrt(N-1);
%             Xa = Xa*sqrt(N-1);
%             Zca_approx = Wa*diag(Sa)*Xa';
%             errZa = norm(full(Zca_approx-Zca))/norm(full(Zca));
%             Qa = length(Sa);
%             Qmax = max(Qmax,Qa);
%             fprintf('\nCoordinate alpha = %2.f : rank Q = %d, error = %.3e for Z',a,Qa,errZa);
%             S(1:Qa,a) = Sa;
%             W(:,1:Qa,a) = Wa;
%             X(:,1:Qa,a) = Xa;
%             Q(a) = Qa;
%             
%             % mXa = mean(Xa,1)';
%             % CXa = cov(Xa); % CXa = 1/(N-1)*Xa(:,:)'*Xa(:,:);
%             % norm(mXa)
%             % norm(CXa-eye(Qa))
%             % norm(Wa'*Wa-eye(Qa))
%             
%             CZa_approx = Wa*diag(Sa).^2*Wa';
%             CZa = cov(Zca'); % CZa = 1/(N-1)*Zca(:,:)*Zca(:,:)';
%             errCZa = norm(full(CZa_approx-CZa))/norm(full(CZa));
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
%                 mysaveas(pathnamegrid,['covariance_CZ_a' num2str(a)],formats,renderer);
%                 mymatlab2tikz(pathnamegrid,['covariance_CZ_a' num2str(a) '.tex']);
%             end
%         end
%         fprintf('\n');
%         S = S(1:Qmax,:);
%         W = W(:,1:Qmax,:);
%         X = X(:,1:Qmax,:);
        
        %% Second reduction step
        % q = (p+1)*Rmax;
        switch index
            case 'time'
                Zc = permute(Z,[1,3,2]);
                Zc = Zc(:,:)';
            case 'coord'
                Zc = Z(:,:)';
        end
        [W,S,X] = svdtruncate(Zc,tol);
        S = S/sqrt(N-1);
        X = X*sqrt(N-1);
        Zc_approx = W*diag(S)*X';
        errZ = norm(full(Zc_approx-Zc))/norm(full(Zc));
        Q = length(S);
        fprintf('\nrank R = %d, rank Q = %d, error = %.3e for Z',Rmax,Q,errZ);
        
        % mX = mean(X,1)';
        % CX = cov(X); % CX = 1/(N-1)*X(:,:)'*X(:,:);
        % norm(mX)
        % norm(CX-eye(Q))
        % norm(W'*W-eye(Q))
        
        CZ_approx = W*diag(S).^2*W';
        CZ = cov(Zc'); % CZ = 1/(N-1)*Zc(:,:)*Zc(:,:)';
        errCZ = norm(full(CZ_approx-CZ))/norm(full(CZ));
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
                case 'time'
                    xlabel('$K=\alpha p+k$','Interpreter','latex')
                    ylabel('$K''=\alpha'' p+k''$','Interpreter','latex')
                case 'coord'
                    xlabel('$K=(k-1)R+\alpha$','Interpreter','latex')
                    ylabel('$K''=(k''-1)R+\alpha''$','Interpreter','latex')
            end
            title('Covariance matrix $[C_{\zeta}(t^k,t^{k''})]_{\alpha,\alpha''} = [C_Z]_{K,K''}$','Interpreter','latex')
            mysaveas(pathnamegrid,'covariance_CZ',formats,renderer);
            mymatlab2tikz(pathnamegrid,'covariance_CZ.tex');
            
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
%                 mysaveas(pathnamegrid,['covariance_CZ_t' num2str(t*100)],formats,renderer);
%                 mymatlab2tikz(pathnamegrid,['covariance_CZ_t' num2str(t*100) '.tex']);
%             end
        end
        
        if g<2^7
            save(fullfile(pathnamegrid,'solution.mat'),'mY','Sig','S','V','W','Rmax','R','Q');
        else
            save(fullfile(pathnamegrid,'solution.mat'),'mY','S','W','Rmax','Q');
        end
    else
        if g<2^7
            load(fullfile(pathnamegrid,'solution.mat'),'mY','Sig','S','V','W','Rmax','R','Q');
        else
            load(fullfile(pathnamegrid,'solution.mat'),'mY','S','W','Rmax','Q');
        end
    end
    
    %% Outputs
    % Display eigenvalues
    if displayEigenvales
        time = 0:10:p;
        for i=1:length(time)-1
            figure('Name','Evolution of eigenvalues')
            clf
            leg = cell(1,11);
            c = 0;
            for t=time(i):time(i+1)
                c = c+1;
                if g<2^7
                    Rt = R(t+1);
                    semilogy(1:Rt,Sig(1:Rt,t+1).^2,'LineStyle','-','Color',getfacecolor(c),'LineWidth',1);
                else
                    load(fullfile(pathnamegrid,['PCA_t' num2str(t) '.mat']),'Sigt','Rt');
                    semilogy(1:Rt,Sigt(:).^2,'LineStyle','-','Color',getfacecolor(c),'LineWidth',1);
                end
                leg{c} = ['t = ' num2str(t*100) ' s'];
                hold on
            end
            hold off
            grid on
            box on
            set(gca,'FontSize',fontsize)
            xlabel('$\alpha$','Interpreter','latex')
            ylabel('$\lambda_{\alpha}(t)$','Interpreter','latex')
            legend(leg{:},'Location','NorthEastOutside')
            mysaveas(pathnamegrid,['eigenvalues_CY_t' num2str(time(i)*100) '_t' num2str(time(i+1)*100)],formats,renderer);
            mymatlab2tikz(pathnamegrid,['eigenvalues_CY_t' num2str(time(i)*100) '_t' num2str(time(i+1)*100) '.tex']);
        end
        
%         figure('Name','Evolution of eigenvalues')
%         clf
%         leg = cell(1,Rmax);
%         for a=1:Rmax
%             semilogy(1:Q(a),S(1:Q(a),a).^2,'LineStyle','-','Color',getfacecolor(a),'LineWidth',1);
%             leg{a} = ['$\alpha$ = ' num2str(a)];
%             hold on
%         end
%         hold off
%         grid on
%         box on
%         set(gca,'FontSize',fontsize)
%         xlabel('$\beta$','Interpreter','latex')
%         ylabel('$\Lambda_{\beta}$','Interpreter','latex')
%         l = legend(leg{:},'Location','NorthEastOutside');
%         set(l,'Interpreter','latex');
%         mysaveas(pathnamegrid,'eigenvalues_CZa',formats,renderer);
%         mymatlab2tikz(pathnamegrid,'eigenvalues_CZa.tex');
        
        figure('Name','Evolution of eigenvalues')
        clf
        semilogy(1:Q,S(:).^2,'LineStyle','-','Color','b','LineWidth',1);
        grid on
        box on
        set(gca,'FontSize',fontsize)
        xlabel('$\beta$','Interpreter','latex')
        ylabel('$\Lambda_{\beta}$','Interpreter','latex')
        mysaveas(pathnamegrid,'eigenvalues_CZ',formats,renderer);
        mymatlab2tikz(pathnamegrid,'eigenvalues_CZ.tex');
        
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
    T = TIMEMODEL(0,p*100,p);
    
    % Solution
    switch index
        case 'time'
            W = permute(reshape(W',[Q,p+1,Rmax]),[1,3,2]);
            % Z_approx = permute(reshape(Zc_approx',[N,p+1,Rmax]),[1,3,2]);
        case 'coord'
            W = reshape(W',[Q,Rmax,p+1]);
            % Z_approx = reshape(Zc_approx',[N,Rmax,p+1]);
    end
    if g<2^7
        % Yc_approx = zeros(N,n,m,p+1);
        Uc_approx = zeros(Q,n,m,p+1);
    end
    
    if g<2^7
        load(fullfile(pathnamegrid,'data.mat'),'Y');
        load(fullfile(pathnamegrid,'data_post.mat'),'YY');
        Y = cat(2,Y,YY);
        clear YY
        Yc = Y - repmat(mY,N,1,1,1); % Yc = Y - mY.*ones(N,1,1,1);
        clear Y
    end
    
    mu = zeros(3*m,p+1);
    Vu = zeros(3*m,p+1);
    mC = zeros(m,p+1);
    VC = zeros(m,p+1);
    mtauTime = zeros(3*m,p+1);
    VtauTime = zeros(3*m,p+1);
    mtauConv = zeros(3*m,p+1);
    VtauConv = zeros(3*m,p+1);
    mtauDiff = zeros(3*m,p+1);
    VtauDiff = zeros(3*m,p+1);
    mtauSurf = zeros(3*m,p+1);
    VtauSurf = zeros(3*m,p+1);
    mtauInterf = zeros(m,p+1);
    VtauInterf = zeros(m,p+1);
    for t=0:p
        mYt = reshape(mY(1,:,:,t+1),[n,m]);
        mut = mYt(1:3,:);
        mCt = mYt(4,:);
        mtauTimet = mYt(4+(1:3),:);
        mtauConvt = mYt(4+(4:6),:);
        mtauDifft = mYt(4+(7:9),:);
        mtauSurft = mYt(4+(10:12),:);
        mtauInterft = mYt(4+13,:);
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
            Vt = reshape(V(:,:,t+1),[r,Rt]);
        else
            load(fullfile(pathnamegrid,['data_t' num2str(t) '.mat']),'Yt');
            load(fullfile(pathnamegrid,['data_post_t' num2str(t) '.mat']),'YYt');
            Yt = cat(2,Yt,YYt);
            clear YYt
            % mYt = mean(Yt,1);
            mYt = mY(1,:,:,t+1);
            Yct = Yt - repmat(mYt,N,1,1); % Yct = Yt - mYt.*ones(N,1,1);
            clear Yt
            load(fullfile(pathnamegrid,['PCA_t' num2str(t) '.mat']),'Sigt','Vt','Rt');
        end
        Wt = reshape(W(:,:,t+1),[Q,Rt]);
        % Zt_approx = reshape(Z_approx(:,:,t+1),[N,Rt]);
        % Zt_approx = X*diag(S)*Wt;
        % Yct_approx = Vt*diag(Sigt)*Zt_approx';
        % Yc_approx(:,:,:,t+1) = reshape(Yct_approx',[N,n,m]);
        Uct_approx = Wt*diag(Sigt)*Vt';
        if g<2^7
            Uc_approx(:,:,:,t+1) = reshape(Uct_approx,[Q,n,m]);
        end
        
        % CYt_approx = cov(Yct_approx(:,:)'); % CYt_approx = 1/(N-1)*Yct_approx(:,:)*Yct_approx(:,:)';
        % CYt_approx = Uct_approx(:,:)'*diag(S).^2*Uct_approx(:,:);
        % VYt_approx = diag(CYt_approx);
        VYt_approx = sum((diag(S)*Uct_approx(:,:)).^2)';
        
        % CYt = cov(Yct(:,:)); % CYt = 1/(N-1)*Yct(:,:)'*Yct(:,:);
        % VYt = diag(CYt);
        VYt = 1/(N-1)*sum(Yct(:,:).^2)';
        errVYt = norm(full(VYt_approx-VYt))/norm(full(VYt));
        fprintf('\nTime t = %4.f s : error = %.3e for VY',t*100,errVYt);
        
        indu = zeros(3*m,1);
        indtauTime = zeros(3*m,1);
        indtauConv = zeros(3*m,1);
        indtauDiff = zeros(3*m,1);
        indtauSurf = zeros(3*m,1);
        for i=1:m
            indu(3*(i-1)+(1:3)) = n*(i-1)+(1:3);
            indtauTime(3*(i-1)+(1:3)) = n*(i-1)+(5:7);
            indtauConv(3*(i-1)+(1:3)) = n*(i-1)+(8:10);
            indtauDiff(3*(i-1)+(1:3)) = n*(i-1)+(11:13);
            indtauSurf(3*(i-1)+(1:3)) = n*(i-1)+(14:16);
        end
        Vu(:,t+1) = VYt_approx(indu);
        VC(:,t+1) = VYt_approx(4:n:end);
        VtauTime(:,t+1) = VYt_approx(indtauTime);
        VtauConv(:,t+1) = VYt_approx(indtauConv);
        VtauDiff(:,t+1) = VYt_approx(indtauDiff);
        VtauSurf(:,t+1) = VYt_approx(indtauSurf);
        VtauInterf(:,t+1) = VYt_approx(n:n:end);
    end
    fprintf('\n');
    mu = TIMEMATRIX(mu,T);
    mC = TIMEMATRIX(mC,T);
    mtauTime = TIMEMATRIX(mtauTime,T);
    mtauConv = TIMEMATRIX(mtauConv,T);
    mtauDiff = TIMEMATRIX(mtauDiff,T);
    mtauSurf = TIMEMATRIX(mtauSurf,T);
    mtauInterf = TIMEMATRIX(mtauInterf,T);
    Vu = TIMEMATRIX(Vu,T);
    VC = TIMEMATRIX(VC,T);
    VtauTime = TIMEMATRIX(VtauTime,T);
    VtauConv = TIMEMATRIX(VtauConv,T);
    VtauDiff = TIMEMATRIX(VtauDiff,T);
    VtauSurf = TIMEMATRIX(VtauSurf,T);
    VtauInterf = TIMEMATRIX(VtauInterf,T);
    
    if g<2^7
        % CY_approx = cov(Yc_approx(:,:)); % CY_approx = 1/(N-1)*Yc_approx(:,:)'*Yc_approx(:,:);
        % CY_approx = Uc_approx(:,:)'*diag(S).^2*Uc_approx(:,:);
        % VY_approx = diag(CY_approx);
        VY_approx = sum((diag(S)*Uc_approx(:,:)).^2)';
        
        % CY = cov(Yc(:,:)); % CY = 1/(N-1)*Yc(:,:)'*Yc(:,:);
        % VY = diag(CY);
        VY = 1/(N-1)*sum(Yc(:,:).^2)';
        errVY = norm(full(VY_approx-VY))/norm(full(VY));
        fprintf('\nerror = %.3e for VY',errVY);
        fprintf('\n');
    end
    
    % Display solution
    if displaySolution
        t = 11;
        ampl = 0;
        % ampl = getsize(M)/max(abs(mut))/5;
        for i=1:3
            evolSolution(M,mu,'displ',i,'colormap',cmap,'filename',['evol_mean_u' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate,...
                'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            evolSolution(M,mtauTime,'displ',i,'colormap',cmap,'filename',['evol_mean_tauTime' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate,...
                'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            evolSolution(M,mtauConv,'displ',i,'colormap',cmap,'filename',['evol_mean_tauConv' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate,...
                'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            evolSolution(M,mtauDiff,'displ',i,'colormap',cmap,'filename',['evol_mean_tauDiff' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate,...
                'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            evolSolution(M,mtauSurf,'displ',i,'colormap',cmap,'filename',['evol_mean_tauSurf' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate,...
                'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            
%             evolSolution(M,Vu,'displ',i,'colormap',cmap,'filename',['evol_var_u' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate);,...
%                 'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%             evolSolution(M,VtauTime,'displ',i,'colormap',cmap,'filename',['evol_var_tauTime' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate);,...
%                 'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%             evolSolution(M,VtauConv,'displ',i,'colormap',cmap,'filename',['evol_var_tauConv' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate);,...
%                 'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%             evolSolution(M,VtauDiff,'displ',i,'colormap',cmap,'filename',['evol_var_tauDiff' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate);,...
%                 'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%             evolSolution(M,VtauSurf,'displ',i,'colormap',cmap,'filename',['evol_var_tauSurf' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate);,...
%                 'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
            
            figure('Name',['Mean of velocity u' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mu,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(pathnamegrid,['mean_u' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau time ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauTime,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(pathnamegrid,['mean_tauTime' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau conv ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauConv,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(pathnamegrid,['mean_tauConv' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau diff ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauDiff,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(pathnamegrid,['mean_tauDiff' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
            figure('Name',['Mean of tau surf ' num2str(i)])
            clf
            plot_sol(M,getmatrixatstep(mtauSurf,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            axis on
            box on
            set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
            mysaveas(pathnamegrid,['mean_tauSurf' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
%             figure('Name',['Variance of velocity u' num2str(i)])
%             clf
%             plot_sol(M,getmatrixatstep(Vu,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
%             colormap(cmap)
%             colorbar
%             axis on
%             box on
%             set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%             mysaveas(pathnamegrid,['var_u' num2str(i) '_t' num2str(t*100)],formats,renderer);
%             
%             figure('Name',['Variance of tau time ' num2str(i)])
%             clf
%             plot_sol(M,getmatrixatstep(VtauTime,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
%             colormap(cmap)
%             colorbar
%             axis on
%             box on
%             set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%             mysaveas(pathnamegrid,['var_tauTime' num2str(i) '_t' num2str(t*100)],formats,renderer);
%             
%             figure('Name',['Variance of tau conv ' num2str(i)])
%             clf
%             plot_sol(M,getmatrixatstep(VtauConv,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
%             colormap(cmap)
%             colorbar
%             axis on
%             box on
%             set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%             mysaveas(pathnamegrid,['var_tauConv' num2str(i) '_t' num2str(t*100)],formats,renderer);
%             
%             figure('Name',['Variance of tau diff ' num2str(i)])
%             clf
%             plot_sol(M,getmatrixatstep(VtauDiff,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
%             colormap(cmap)
%             colorbar
%             axis on
%             box on
%             set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%             mysaveas(pathnamegrid,['var_tauDiff' num2str(i) '_t' num2str(t*100)],formats,renderer);
%             
%             figure('Name',['Variance of tau surf ' num2str(i)])
%             clf
%             plot_sol(M,getmatrixatstep(VtauSurf,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
%             colormap(cmap)
%             colorbar
%             axis on
%             box on
%             set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%             mysaveas(pathnamegrid,['var_tauSurf' num2str(i) '_t' num2str(t*100)],formats,renderer);
        end
        
        Mscal = final(M,DDL('C'));
        
        evolSolution(Mscal,mC,'colormap',cmap,'filename','evol_mean_C','pathname',pathnamegrid,'FrameRate',framerate,...
            'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        evolSolution(Mscal,mtauInterf,'colormap',cmap,'filename','evol_mean_tauInterf','pathname',pathnamegrid,'FrameRate',framerate,...
            'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        
%         evolSolution(Mscal,VC,'colormap',cmap,'filename','evol_var_C','pathname',pathnamegrid,'FrameRate',framerate,...
%             'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%         evolSolution(Mscal,VtauInterf,'colormap',cmap,'filename','evol_var_tauInterf','pathname',pathnamegrid,'FrameRate',framerate,...
%             'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        
        figure('Name','Mean of indicator function C')
        clf
        plot_sol(Mscal,getmatrixatstep(mC,t+1));
        title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(pathnamegrid,['mean_C_t' num2str(t*100)],formats,renderer);
        
        figure('Name','Mean of tau interf')
        clf
        plot_sol(Mscal,getmatrixatstep(mtauInterf,t+1));
        title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(pathnamegrid,['mean_tauInterf_t' num2str(t*100)],formats,renderer);
        
%         figure('Name','Variance of indicator function C')
%         clf
%         plot_sol(Mscal,getmatrixatstep(VC,t+1));
%         title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
%         colormap(cmap)
%         colorbar
%         box on
%         set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%         mysaveas(gridname,['var_C_t' num2str(t*100)],formats,renderer);
%         
%         figure('Name','Variance of tau interf')
%         clf
%         plot_sol(Mscal,getmatrixatstep(VtauInterf,t+1));
%         title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
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
        Vut = getmatrixatstep(Vu,t+1);
        VCt = getmatrixatstep(VC,t+1);
        VtauTimet = getmatrixatstep(VtauTime,t+1);
        VtauConvt = getmatrixatstep(VtauConv,t+1);
        VtauDifft = getmatrixatstep(VtauDiff,t+1);
        VtauSurft = getmatrixatstep(VtauSurf,t+1);
        VtauInterft = getmatrixatstep(VtauInterf,t+1);
        write_vtk_mesh(M,mut,mCt,mtauTimet,mtauConvt,mtauDifft,mtauSurft,mtauInterft,pathnamegrid,'diphasic_fluids_mean',1,t);
        write_vtk_mesh(M,Vut,VCt,VtauTimet,VtauConvt,VtauDifft,VtauSurft,VtauInterft,pathnamegrid,'diphasic_fluids_variance',1,t);
    end
    make_pvd_file(pathnamegrid,'diphasic_fluids_mean',1,p+1);
    make_pvd_file(pathnamegrid,'diphasic_fluids_variance',1,p+1);

end
