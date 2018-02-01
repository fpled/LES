clc
clearvars
% close all

solveProblem = true;
displaySolution = false;
displayEigenvales = false;
displayCovariance  = true;

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

% for g=2.^(4:8)
% for g=2.^(4:7)
for g=2.^4
    gridname = ['Grid' num2str(g)];
    disp(gridname)
    pathnamegrid = fullfile(pathname,gridname);
    load(fullfile(pathnamegrid,'data.mat'),'Y');
    
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
    
    if solveProblem
        %% First reduction step
        mY = mean(Y,1);
        Yc = Y - repmat(mY,N,1,1,1); % Yc = Y - mY.*ones(N,1,1,1);
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
            
            % mZt = mean(Zt,1)';
            % CZt = cov(Zt); % CZt = 1/(N-1)*Zt(:,:)'*Zt(:,:);
            % norm(mZt)
            % norm(CZt-eye(Rt))
            % norm(Vt'*Vt-eye(Rt))
        end
        fprintf('\n');
        Sig = Sig(1:Rmax,:);
        V = V(:,1:Rmax,:);
        Z = Z(:,1:Rmax,:);
        
        if displayEigenvales
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
                mysaveas(pathnamegrid,['eigenvalues_CY_t' num2str(time(i)*100) '_t' num2str(time(i+1)*100)],formats,renderer);
                mymatlab2tikz(pathnamegrid,['eigenvalues_CY_t' num2str(time(i)*100) '_t' num2str(time(i+1)*100) '.tex']);
            end
        end
        
        %% Second reduction step for each coordinate
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
            
            % mXa = mean(Xa,1)';
            % CXa = cov(Xa); % CXa = 1/(N-1)*Xa(:,:)'*Xa(:,:);
            % norm(mXa)
            % norm(CXa-eye(Qa))
            % norm(Wa'*Wa-eye(Qa))
            
            CZa_approx = Wa*diag(Sa).^2*Wa';
            CZa = cov(Zca'); % CZa = 1/(N-1)*Zca(:,:)*Zca(:,:)';
            errCZa = norm(full(CZa_approx-CZa))/norm(full(CZa));
            fprintf('\n                                     error = %.3e for CZ',errCZa);
            
            if displayCovariance
                figure('Name','Covariance matrix')
                clf
                imagesc(CZa)
                colorbar
                axis image
                set(gca,'FontSize',fontsize)
                xlabel('$k$','Interpreter','latex')
                ylabel('$k''$','Interpreter','latex')
                title(['Covariance matrix $[C_{\zeta}(t^k,t^{k''})]_{\alpha,\alpha} = [C_{Z_{\alpha}}]_{k,k''}$ for $\alpha=$' num2str(a)],'Interpreter','latex')
                mysaveas(pathnamegrid,['covariance_CZ_a' num2str(a)],formats,renderer);
                mymatlab2tikz(pathnamegrid,['covariance_CZ_a' num2str(a) '.tex']);
            end
        end
        fprintf('\n');
        S = S(1:Qmax,:);
        W = W(:,1:Qmax,:);
        X = X(:,1:Qmax,:);
        
        if displayEigenvales
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
            mysaveas(pathnamegrid,'eigenvalues_CZa',formats,renderer);
            mymatlab2tikz(pathnamegrid,'eigenvalues_CZa.tex');
        end
        
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
        
        if displayEigenvales
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
        
        save(fullfile(pathnamegrid,'solution.mat'),'mY','V','W','Sig','S','Rmax','Q');
    else
        load(fullfile(pathnamegrid,'solution.mat'),'mY','V','W','Sig','S','Rmax','Q');
    end
    
    %% Outputs
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
    % Yc_approx = zeros(N,n,m,p+1);
    Uc_approx = zeros(Q,n,m,p+1);
    
%     u = repmat({zeros(3*m,p+1)},1,N);
%     phase = repmat({zeros(m,p+1)},1,N);
    mu = zeros(3*m,p+1);
    Vu = zeros(3*m,p+1);
    mphase = zeros(m,p+1);
    Vphase = zeros(m,p+1);
%     v = repmat({zeros(3*m,p+1)},1,N);
%     vphase = repmat({zeros(m,p+1)},1,N);
    for t=0:p
%         for l=1:N
%             Ylt = reshape(Y(l,:,:,t+1),[n,m]);
%             ult = Ylt(1:3,:);
%             phaselt = Ylt(4,:);
%             u{l}(:,t+1) = ult(:);
%             phase{l}(:,t+1) = phaselt(:);
%         end
%         for a=1:Rmax
%             Vat = reshape(V(:,a,t+1),[n,m]);
%             vat = Vat(1:3,:);
%             vphaseat = Vat(4,:);
%             v{a}(:,t+1) = vat(:);
%             vphase{a}(:,t+1) = vphaseat(:);
%         end
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
        Uct_approx = Wt*diag(Sigt)*Vt';
        Uc_approx(:,:,:,t+1) = reshape(Uct_approx,[Q,n,m]);
        
        % CYt_approx = cov(Yct_approx(:,:)'); % CYt_approx = 1/(N-1)*Yct_approx(:,:)*Yct_approx(:,:)';
        % CYt_approx = Uct_approx(:,:)'*diag(S).^2*Uct_approx(:,:);
        % VYt_approx = diag(CYt_approx);
        VYt_approx = sum((diag(S)*Uct_approx(:,:)).^2)';
        % CYt = cov(Yct(:,:)); % CYt = 1/(N-1)*Yct(:,:)'*Yct(:,:);
        % VYt = diag(CYt);
        VYt = 1/(N-1)*sum(Yct(:,:).^2)';
        errVYt = norm(full(VYt_approx-VYt))/norm(full(VYt));
        fprintf('\nTime t = %4.f s : error = %.3e for VY',t*100,errVYt);
        
        Vu(:,t+1) = VYt(setdiff(1:end,n:n:end));
        Vphase(:,t+1) = VYt(n:n:end);
    end
    fprintf('\n');
%     for l=1:N
%         u{l} = TIMEMATRIX(u{l},T);
%         phase{l} = TIMEMATRIX(phase{l},T);
%     end
%     for a=1:Rmax
%         v{a} = TIMEMATRIX(v{a},T);
%         vphase{a} = TIMEMATRIX(vphase{a},T);
%     end
    mu = TIMEMATRIX(mu,T);
    mphase = TIMEMATRIX(mphase,T);
    Vu = TIMEMATRIX(Vu,T);
    Vphase = TIMEMATRIX(Vphase,T);
    
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
    
    % Display solution
    if displaySolution
        t = 11;
        ampl = 0;
        % ampl = getsize(M)/max(abs(mut))/5;
        for i=1:3
            evolSolution(M,mu,'displ',i,'colormap',cmap,'filename',['evol_mean_u' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate,...
                'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%             evolSolution(M,Vu,'displ',i,'colormap',cmap,'filename',['evol_var_u' num2str(i)],'pathname',pathnamegrid,'FrameRate',framerate);,...
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
        end
        
        Mphase = final(M,DDL('Phase'));
        
        evolSolution(Mphase,mphase,'colormap',cmap,'filename','evol_mean_phase','pathname',pathnamegrid,'FrameRate',framerate,...
            'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
%         evolSolution(Mphase,Vphase,'colormap',cmap,'filename','evol_var_phase','pathname',pathnamegrid,'FrameRate',framerate,...
%             'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
        
        figure('Name','Mean of phase')
        clf
        plot_sol(Mphase,getmatrixatstep(mphase,t+1));
        title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
        colormap(cmap)
        colorbar
        axis on
        box on
        set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
        mysaveas(pathnamegrid,['mean_phase_t' num2str(t*100)],formats,renderer);
        
%         figure('Name','Variance of phase')
%         clf
%         plot_sol(Mphase,getmatrixatstep(Vphase,t+1));
%         title(['time ' num2str(t*100,'%.2f') ' s'],'FontSize',fontsize)
%         colormap(cmap)
%         colorbar
%         box on
%         set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
%         mysaveas(gridname,['var_phase_t' num2str(t*100)],formats,renderer);
    
    end
    
%     for l=1:N
%         for t=0:p
%             ult = getmatrixatstep(u{l},t+1);
%             phaselt = getmatrixatstep(phase{l},t+1);
%             write_vtk_mesh(M,ult,phaselt,pathnamegrid,['diphasic_fluids_sample' num2str(l)],1,t);
%         end
%         make_pvd_file(pathnamegrid,['diphasic_fluids_sample' num2str(l)],1,p+1);
%     end
%     for a=1:Rmax
%         for t=0:p
%             vat = getmatrixatstep(v{a},t+1);
%             vphaseat = getmatrixatstep(vphase{a},t+1);
%             write_vtk_mesh(M,vat,vphaseat,pathnamegrid,['diphasic_fluids_eigenvector' num2str(a)],1,t);
%         end
%         make_pvd_file(pathnamegrid,['diphasic_fluids_eigenvector' num2str(a)],1,p+1);
%     end
    
    for t=0:p
        mut = getmatrixatstep(mu,t+1);
        mphaset = getmatrixatstep(mphase,t+1);
        Vut = getmatrixatstep(mu,t+1);
        Vphaset = getmatrixatstep(mphase,t+1);
        write_vtk_mesh(M,mut,mphaset,pathnamegrid,'diphasic_fluids_mean',1,t);
        write_vtk_mesh(M,Vut,Vphaset,pathnamegrid,'diphasic_fluids_variance',1,t);
    end
    make_pvd_file(pathnamegrid,'diphasic_fluids_mean',1,p+1);
    make_pvd_file(pathnamegrid,'diphasic_fluids_variance',1,p+1);

end
