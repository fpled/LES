clc
clearvars
close all

displaySolution = true;

cmap = 'default';
% cmap = 'gray';
framerate = 5;
fontsize = 16;
interpreter = 'latex';
formats = {'fig','epsc2'};
renderer = 'OpenGL';

% for g=2.^(4:8)
% for g=2.^(4:6)
for g=2.^4
    gridname = ['Grid' num2str(g)];
    disp(gridname)
    load(fullfile(gridname,'data.mat'),'Y');
    
    N = size(Y,1);
    n = size(Y,2);
    m = size(Y,3);
    p = size(Y,4)-1;
    r = n*m;
    R = min(r,N);
    fprintf('\nn = %d variables',n);
    fprintf('\nN = %d samples',N);
    fprintf('\nm = %d spatial points',m);
    fprintf('\np+1 = %d time steps',p+1);
    fprintf('\n');
    
    tol = eps;
    
    % First reduction step
    mY = mean(Y,1);
    Yc = Y - mY.*ones(N,1,1,1); % Yc = Y - repmat(mY,N,1,1,1);
    Sig = zeros(R,p+1);
    V = zeros(r,R,p+1);
    Z = zeros(N,R,p+1);
    Rmax = 1;
    for t=0:p
        Yct = Yc(:,:,:,t+1);
        Yct = Yct(:,:)';
        [Vt,Sigt,Zt] = svdtruncate(Yct,tol);
        Sigt = Sigt/sqrt(N-1);
        Zt = Zt*sqrt(N-1);
        Yct_approx = Vt*diag(Sigt)*Zt';
        errYct = norm(Yct_approx-Yct)/norm(Yct);
        Rt = length(Sigt);
        Rmax = max(Rmax,Rt);
        fprintf('\nTime t = %4.f s : rank R = %d, error = %.3e for Y',t*100,Rt,errYct);
        Sig(1:Rt,t+1) = Sigt;
        V(:,1:Rt,t+1) = Vt;
        Z(:,1:Rt,t+1) = Zt;
        
        mZt = mean(Zt,1)';
        CZt = cov(Zt); % CZt = 1/(N-1)*Zt(:,:)'*Zt(:,:);
%         norm(mZt)
%         norm(CZt-eye(Rt))
%         norm(Vt'*Vt-eye(Rt))
    end
    fprintf('\n');
    Sig = Sig(1:Rmax,:);
    V = V(:,1:Rmax,:);
    Z = Z(:,1:Rmax,:);
    
    time = 0:10:p;
    for i=1:length(time)-1
        figure('Name','Evolution of eigenvalues')
        clf
        leg = cell(1,11);
        c = 0;
        for t=time(i):time(i+1)
            c = c+1;
            semilogy(1:Rmax,Sig(:,t+1).^2,'LineStyle','-','Color',getfacecolor(c),'LineWidth',1);
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
        mysaveas(gridname,['eigenvalues_CY_t' num2str(time(i)*100) '_t' num2str(time(i+1)*100)],formats,renderer);
        mymatlab2tikz(gridname,['eigenvalues_CY_t' num2str(time(i)*100) '_t' num2str(time(i+1)*100) '.tex']);
    end
    
    % Second reduction step
    Q = min(p+1,N);
    S = zeros(Q,Rmax);
    W = zeros(p+1,Q,Rmax);
    X = zeros(N,Q,Rmax);
    Q = zeros(1,Rmax);
    Qmax = 1;
    for a=1:Rmax
        Zca = Z(:,a,:);
        Zca = Zca(:,:)';
        [Wa,Sa,Xa] = svdtruncate(Zca,tol);
        Sa = Sa/sqrt(N-1);
        Xa = Xa*sqrt(N-1);
        Zca_approx = Wa*diag(Sa)*Xa';
        errZa = norm(full(Zca_approx-Zca))/norm(full(Zca));
        Qa = length(Sa);
        Qmax = max(Qmax,Qa);
        fprintf('\nCoordinate alpha = %2.f : rank Q = %d, error = %.3e for Z',a,Qa,errZa);
        S(1:Qa,a) = Sa;
        W(:,1:Qa,a) = Wa;
        X(:,1:Qa,a) = Xa;
        Q(a) = Qa;
        
        mXa = mean(Xa,1)';
        CXa = cov(Xa); % CXa = 1/(N-1)*Xa(:,:)'*Xa(:,:);
%         norm(mXa)
%         norm(CXa-eye(Qa))
%         norm(Wa'*Wa-eye(Qa))
        
        CZa_approx = Wa*diag(Sa).^2*Wa';
        CZa = cov(Zca'); % CZa = 1/(N-1)*Zca(:,:)*Zca(:,:)';
        errCZa = norm(full(CZa_approx-CZa))/norm(full(CZa));
        fprintf('\n                                     error = %.3e for CZ',errCZa);
        
        figure('Name','Covariance matrix')
        clf
        % surf(CZa)
        imagesc(CZa)
        colorbar
        set(gca,'FontSize',fontsize)
        xlabel('$k$','Interpreter','latex')
        ylabel('$k''$','Interpreter','latex')
        title(['Covariance matrix $[C_{Z_{\alpha}}]_{k,k''}$ for $\alpha=$' num2str(a)],'Interpreter','latex')
        mysaveas(gridname,['covariance_CZ_a' num2str(a)],formats,renderer);
        mymatlab2tikz(gridname,['covariance_CZ_a' num2str(a) '.tex']);
    end
    fprintf('\n');
    S = S(1:Qmax,:);
    W = W(:,1:Qmax,:);
    X = X(:,1:Qmax,:);
    
    figure('Name','Evolution of eigenvalues')
    clf
    leg = cell(1,Rmax);
    for a=1:Rmax
        semilogy(1:Q(a),S(1:Q(a),a).^2,'LineStyle','-','Color',getfacecolor(a),'LineWidth',1);
        leg{a} = ['$\alpha$ = ' num2str(a)];
        hold on
    end
    hold off
    grid on
    box on
    set(gca,'FontSize',fontsize)
    xlabel('$\beta$','Interpreter','latex')
    ylabel('$\Lambda_{\beta}$','Interpreter','latex')
    l = legend(leg{:},'Location','NorthEastOutside');
    set(l,'Interpreter','latex');
    mysaveas(gridname,'eigenvalues_CZ',formats,renderer);
    mymatlab2tikz(gridname,'eigenvalues_CZ');
    
%     q = (p+1)*Rmax;
%     Zc = Z(:,:)';
%     [W,S,X] = svdtruncate(Zc,tol);
%     S = S/sqrt(N-1);
%     X = X*sqrt(N-1);
%     Zc_approx = W*diag(S)*X';
%     errZ = norm(full(Zc_approx-Zc))/norm(full(Zc));
%     Q = length(S);
%     fprintf('\nrank R = %d, rank Q = %d, error = %.3e for Z',Rmax,Q,errZ);
%     
%     mX = mean(X,1)';
%     CX = cov(X); % CX = 1/(N-1)*X'*X;
%     norm(mX)
%     norm(CX-eye(Q))
%     norm(W'*W-eye(Q))
%     
%     CZ_approx = W*diag(S).^2*W';
%     CZ = cov(Zc'); % CZ = 1/(N-1)*Zc(:,:)*Zc(:,:)';
%     errCZ = norm(full(CZ_approx-CZ))/norm(full(CZ));
%     fprintf('\nerror = %.3e for CZ',errCZ);
%     fprintf('\n');
%     for t=0:p
%         CZt_approx = CZ_approx(t*Rmax+(1:Rmax),t*Rmax+(1:Rmax));
%         CZt = CZ(t*Rmax+(1:Rmax),t*Rmax+(1:Rmax));
%         norm(CZt_approx-eye(Q,Q))
%         norm(CZt-eye(Q,Q))
%     end
%     
%     figure('Name','Evolution of eigenvalues')
%     clf
%     semilogy(1:Q,S(:).^2,'LineStyle','-','Color','b','LineWidth',1);
%     grid on
%     box on
%     set(gca,'FontSize',fontsize)
%     xlabel('$\beta$','Interpreter','latex')
%     ylabel('$\Lambda_{\beta}$','Interpreter','latex')
%     mysaveas(gridname,'eigenvalues_CZ',formats,renderer);
%     mymatlab2tikz(gridname,'eigenvalues_CZ');
    
    
    % save(fullfile(gridname,'res.mat'),'mY','V','W','Sig','S');
    
    if displaySolution
        L = 1;
        D = DOMAIN(3,[0.0,0.0,0.0],[L,L,L]);
        elemtype = 'CUB8';
        nbelem = repmat(g,1,3);
        Mesh = build_model(D,'nbelem',nbelem,'elemtype',elemtype);
        coord = getcoord(getnode(Mesh));
        Mesh = setnode(Mesh,NODE(fliplr(coord)));
        Mesh = final(Mesh,DDL(DDLVECT('U',Mesh.syscoord)));
        T = TIMEMODEL(0,p*100,p);
        
        W = reshape(W',[Q,Rmax,p+1]);
        % Z_approx = reshape(Zc_approx',[N,Rmax,p+1]);
        % Yc_approx = zeros(N,n,m,p+1);
        Uc_approx = zeros(Q,n,m,p+1);
        mu = zeros(3*m,p+1);
        Vu = zeros(3*m,p+1);
        mphase = zeros(m,p+1);
        Vphase = zeros(m,p+1);
        for t=0:p
            mYt = reshape(mY(1,:,:,t+1),[n,m]);
            mut = mYt(1:3,:);
            mphaset = mYt(4,:);
            mu(:,t+1) = mut(:);
            mphase(:,t+1) = mphaset(:);
            
            Yct = Yc(:,:,:,t+1);
            Vt = reshape(V(:,:,t+1),[r,Rmax]);
            Sigt = Sig(:,t+1);
            Wt = reshape(W(:,:,t+1),[Q,Rmax]);
            % Zt_approx = reshape(Z_approx(:,:,t+1),[N,Rmax]);
            % Zt_approx = X*diag(S)*Wt;
            % Yct_approx = Vt*diag(Sigt)*Zt_approx';
            % Yc_approx(:,:,:,t+1) = reshape(Yct_approx',[N,n,m]);
            % Yc_approx(:,:,:,t+1) = reshape(Zt_approx*diag(Sigt)*Vt',[N,n,m]);
            % Yc_approx(:,:,:,t+1) = reshape(X*diag(S)*Wt*diag(Sigt)*Vt',[N,n,m]);
            Uct_approx = reshape(Wt*diag(Sigt)*Vt',[Q,n,m]);
            Uc_approx(:,:,:,t+1) = Uct_approx;
            
%             % CYt_approx = 1/(N-1)*Yct_approx(:,:)'*Yct_approx(:,:);
%             CYt_approx = 1/(N-1)*Uct_approx(:,:)'*diag(S).^2*Uct_approx(:,:);
%             VYt_approx = diag(CYt_approx);
            VYt_approx = 1/(N-1)*sum((diag(S)*Uct_approx(:,:)).^2)';
%             CYt = 1/(N-1)*Yct(:,:)'*Yct(:,:);
%             VYt = diag(CYt);
            VYt = 1/(N-1)*sum(Yct(:,:).^2)';
            errVYt = norm(full(VYt_approx-VYt))/norm(full(VYt));
            fprintf('\nTime t = %4.f s : error = %.3e for VY',t*100,errVYt);
            Vu(:,t+1) = VYt(setdiff(1:end,n:n:end));
            Vphase(:,t+1) = VYt(n:n:end);
        end
        fprintf('\n');
        mu = TIMEMATRIX(mu,T);
        mphase = TIMEMATRIX(mphase,T);
        Vu = TIMEMATRIX(Vu,T);
        Vphase = TIMEMATRIX(Vphase,T);
        
%         % CY_approx = 1/(N-1)*Yc_approx(:,:)'*Yc_approx(:,:); % CY_approx = cov(Yc_approx(:,:));
%         CY_approx = 1/(N-1)*Uc_approx(:,:)'*diag(S).^2*Uc_approx(:,:); % CY_approx = cov(diag(S)*Uc_approx(:,:));
%         CY = 1/(N-1)*Yc(:,:)'*Yc(:,:); % CY = cov(Yc(:,:));
%         VY_approx = diag(CY_approx);
%         VY = diag(CY);
%         errVY = norm(full(VY_approx-VY))/norm(full(VY));
%         fprintf('\nerror = %.3e for VY',errVY);
        
        t = 11;
        ampl = 0;
        % ampl = getsize(Mesh)/max(abs(mut))/5;
        for i=1:3
            evolSolution(Mesh,mu,'displ',i,'colormap',cmap,'filename',['evol_mean_u' num2str(i)],'pathname',gridname,'FrameRate',framerate);
            %         evolSolution(Mesh,Vu,'displ',i,'colormap',cmap,'filename',['evol_var_u' num2str(i)],'pathname',gridname,'FrameRate',framerate);
            
            figure('Name',['Mean of velocity u' num2str(i)])
            clf
            plot_sol(Mesh,getmatrixatstep(mu,t+1),'displ',i,'ampl',ampl);
            title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
            colormap(cmap)
            colorbar
            set(gca,'FontSize',fontsize)
            mysaveas(gridname,['mean_u' num2str(i) '_t' num2str(t*100)],formats,renderer);
            
%             figure('Name',['Variance of velocity u' num2str(i)])
%             clf
%             plot_sol(Mesh,getmatrixatstep(Vu,t+1),'displ',i,'ampl',ampl);
%             title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
%             colormap(cmap)
%             colorbar
%             set(gca,'FontSize',fontsize)
%             mysaveas(gridname,['var_u' num2str(i) '_t' num2str(t*100)],formats,renderer);
        end
        
        Mesh = final(Mesh,DDL('Phase'));
        evolSolution(Mesh,mphase,'colormap',cmap,'filename','evol_mean_phase','pathname',gridname,'FrameRate',framerate);
%         evolSolution(Mesh,Vphase,'colormap',cmap,'filename','evol_var_phase','pathname',gridname,'FrameRate',framerate);
        
        figure('Name','Mean of phase')
        clf
        plot_sol(Mesh,getmatrixatstep(mphase,t+1));
        title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        set(gca,'FontSize',fontsize)
        mysaveas(gridname,['mean_phase_t' num2str(t*100)],formats,renderer);
        
%         figure('Name','Variance of phase')
%         clf
%         plot_sol(Mesh,getmatrixatstep(Vphase,t+1));
%         title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
%         colormap(cmap)
%         colorbar
%         set(gca,'FontSize',fontsize)
%         mysaveas(gridname,['var_phase_t' num2str(t*100)],formats,renderer);
    end
end
