% Class PrincipalComponentAnalysis

% Copyright (c) 2020, Florent Pled

classdef PrincipalComponentAnalysis
    
    properties
        tol = eps; % relative precision
        maxRank = Inf; % maximum rank
        thresholdingParameter
        thresholdingType = 'hard';
        checkOrthonormality = false;
        display = false;
    end
    
    methods
        
        function s = PrincipalComponentAnalysis(varargin)
            % s = PrincipalComponentAnalysis(varargin)
            % Principal component analysis (PCA) of a vector-valued random variable
            % s.tol: prescribed relative precision, positive number <1 (eps by default)
            % s.maxRank: maximal rank (Inf by default)
            % s.checkOrthonormality (true or false): check orthonormality for the coefficient matrix and the principal components (false by default)
            % s.thresholdingParameter: thresholding parameter for determining the minimal singular values ([] by default)
            % s.thresholdingType ('hard' or 'soft'): thresholding parameter for determining the minimal singular values ('hard' by default)
            % s.display (true or false): display (false by default)
            
            expectedThresholdingTypes = {'hard','soft'};
            
            p = ImprovedInputParser();
            addParamValue(p,'tol',eps,@isscalar);
            addParamValue(p,'maxRank',inf,@isscalar);
            addParamValue(p,'thresholdingParameter',[]);
            addParamValue(p,'thresholdingType','hard',...
                @(x) any(validatestring(x,expectedThresholdingTypes)));
            addParamValue(p,'checkOrthonormality',false,@islogical);
            addParamValue(p,'display',false,@islogical);
            
            parse(p,varargin{:});
            s = passMatchedArgsToProperties(p,s);
        end
        
        function [V,sv,Z,err,mu] = principalComponents(s,X)
            % [V,sv,Z,err,mu] = principalComponents(s,X)
            % computes the principal component coefficients V associated to the largest singular values sv 
            % and the principal components Z of the n-by-N data matrix X such that Y = mu + V*diag(sv)*Z'
            % approximates X with relative precision s.tol and maximal rank s.maxRank
            %
            % X: n-by-N matrix whose rows correspond to variables and columns correspond to observations
            % V: n-by-r matrix whose columns are the principal components such that V'*V = eye(r)
            % sv: column vector containing the corresponding singular values 
            % Z: N-by-r matrix whose rows correspond to observations and columns correspond to variables
            % err: column vector containing the corresponding L2-error
            % mu: mean vector of X returned as a column vector
            
            N = size(X,2);
            [Xc,mu] = s.center(X);
            clear X
            
            [V,S,Z,err] = s.svd(Xc);
            S = S/sqrt(N-1);
            sv = diag(S);
            Z = Z*sqrt(N-1);
            % if verLessThan('matlab','9.1') % compatibility (<R2016b)
            %     Yc = V*S*Z';
            % else
            %     Yc = V*(sv.*Z');
            % end
                
            % [coeff,score,latent] = pca(X');
            % V = coeff;
            % sv = sqrt(latent);
            % S = diag(sv);
            % if verLessThan('matlab','9.1') % compatibility (<R2016b)
            %     Z = score*diag(1./sv);
            % else
            %     Z = score./sv';
            % end
            % % Yc = coeff*score'; % Yc = V*(Z*S)' = V*S*Z'; % score = Z*S;
            
            if s.checkOrthonormality
                r = length(sv);
                mZ = mean(Z,1)';
                CZ = cov(Z); % CZ = 1/(N-1)*Z'*Z;
                fprintf('norm(mZ) = %f\n',norm(mZ));
                fprintf('norm(CZ-I) = %f\n',norm(CZ-eye(r)));
                fprintf('norm(V''*V-I) = %f\n',norm(V'*V-eye(r)));
            end
        end
        
        function sv = singularValues(s,X)
            % function sv = singularValues(X)
            % computes the singular values of the centered matrix Xc = X - mu of X with mean vector mu 
            % corresponding to the square root of the eigenvalues of the covariance matrix of X
            
            N = size(X,2);
            Xc = s.center(X);
            
            sv = svd(full(Xc),'econ');
            sv = sv/sqrt(N-1);
        end
        
        function la = eigenValues(s,X)
            % function la = eigenValues(X)
            % computes the eigenvalues of the covariance matrix of X
            
            sv = s.singularValues(X);
            la = sv.^2;
        end
        
        function [U,S,V,err] = svd(s,X,tol,p)
            % [U,S,V,err] = svd(s,X)
            % performs a truncated economy-size singular value decomposition (SVD) of the matrix X
            % such that Y = U*S*V' approximates X with relative precision s.tol and maximal rank s.maxRank
            %
            % [U,S,V,err] = svd(s,X,tol)
            % uses the argument tol for tolerance, instead of s.tol
            %
            % [U,S,V,err] = svd(x,tol,p)
            % performs a truncated economy-size SVD of the matrix X with relative precision tol
            % in schatten-p norm (1<=p<=Inf, p=2 for Frobenius), p=2 by default
            %
            % function [sv,err] = svd(s,X,tol,p)
            % function sv = svd(s,X,tol,p)
            % returns the first singular values of matrix X in descending order as a column vector
            
            if nargin==2 || isempty(tol)
                tol = s.tol;
            end
            if nargin<=3
                p = 2;
            end
            [U,S,V] = svd(full(X),'econ');
            sv = diag(S);
            
            %if p==Inf
            %    err = sv/max(sv);
            %else
            %    if verLessThan('matlab','8.2') % compatibility (<R2013b)
            %        err = (flipdim(cumsum(flipdim(sv.^p,1)),1)/sum(sv.^p)).^(1/p);
            %    else
            %        err = (flip(cumsum(flip(sv.^p)))/sum(sv.^p)).^(1/p);
            %    end
            %end
            %err = [err(2:end);0];
            if p==Inf
                err = sv/max(sv);
                err = [err(2:end);0];
            else
                err = (1-cumsum(sv.^p)/sum(sv.^p)).^(1/p);
            end
            m = find(err<tol);
            if isempty(m)
                m = min(size(X));
            else
                m = min(m);
            end
            %if m>s.maxRank
            %    warning('maxRank reached, tolerance not achieved')
            %end
            m = min(m,s.maxRank);
            if ~isempty(s.thresholdingParameter) && s.thresholdingParameter~=0
                switch lower(s.thresholdingType)
                    case 'soft'
                        sv = sv - s.thresholdingParameter;
                        m = min(m,find(sv>=0,1,'last'));
                    case 'hard'
                        m = min(m,find(sv>=s.thresholdingParameter,1,'last'));
                end
            end
            
            U = U(:,1:m);
            S = S(1:m,1:m);
            V = V(:,1:m);
            err = err(1:m);
            
            if nargout<=2
                U = sv(1:m);
                S = err(1:m);
            end
        end
        
        function [Yt,mYt,Vt,svt,mZt,Wt,Rt] = reconstructionDoubleAtStep(s,mY,V,sv,mZ,W,sw,X,R,t,varargin)
            % function [Yt,mYt,Vt,svt,mZt,Wt,Rt] = reconstructionDoubleAtStep(s,mY,V,sv,mZ,W,sw,X,R,t,varargin)
            % reconstructs the double PCA representation Y(t) = mY(t) + V(t)*diag(sv(t))*Z(t)' 
            % with Z(t) = mZ(t) + W(t)*diag(sw)*X' at step t
            
            % N = size(X,1);
            [Vt,svt,Wt,Rt] = s.getPrincipalComponentsDoubleAtStep(V,sv,W,R,t,varargin{:});
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                Zct = Wt*diag(sw)*X';
            else
                Zct = Wt*(sw.*X');
            end
            if isempty(mZ)
                mZt  = [];
                Zt = Zct;
            else
                mZt = mZ(:,1:Rt,t)';
                if verLessThan('matlab','9.1') % compatibility (<R2016b)
                    Zt = bsxfun(@plus,mZt,Zct);
                else
                    Zt = mZt + Zct;
                end
                % Zt = repmat(mZt,[1,N]) + Zct;
            end
            
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                Yct = Vt*diag(svt)*Zt;
            else
                Yct = Vt*(svt.*Zt);
            end
            
            if isempty(mY)
                mYt = [];
                Yt = Yct;
            else
                mYt = mY(:,:,t);
                if verLessThan('matlab','9.1') % compatibility (<R2016b)
                    Yt = bsxfun(@plus,mYt,Yct);
                else
                    Yt = mYt + Yct;
                end
                % Yt = repmat(mYt,[1,N]) + Yct;
            end
        end
        
        function [Xc,mu] = center(s,X)
            % function [Xc,mu] = center(X)
            % computes the centered matrix Xc = X - mu of X 
            % with mean vector mu
            
            mu = s.mean(X);
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                Xc = bsxfun(@minus,X,mu);
            else
                Xc = X - mu;
            end
        end
        
        function sig = std(s,V,sv)
            % function sig = std(V,sv)
            % returns the estimated standard deviation of the PCA representation
            
            v = s.var(V,sv);
            sig = sqrt(v);
        end
        
        function sig = stdDouble(s,V,sv,W,sw)
            % function sig = stdDouble(s,V,sv,W,sw)
            % returns the estimated standard deviation of a double PCA representation
            
            v = s.varDouble(V,sv,W,sw);
            sig = sqrt(v);
        end
        
    end
    
    methods (Static)
        function mu = mean(X)
            % function mu = mean(X)
            % computes the mean vector of X 
            
            mu = mean(X,2);
        end
        
        function C = cov(V,sv)
            % function C = cov(V,sv)
            % returns the estimated covariance of the PCA representation
            
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
               C = V*diag(sv).^2*V';
            else
               C = V*(sv.^2.*V');
            end
        end
        
        function v = var(V,sv)
            % function v = var(V,sv)
            % returns the estimated variance of the PCA representation
            
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                v = sum((diag(sv)*V').^2)';
            else
                v = sum((sv.*V').^2)';
            end
        end
        
        function v = varDouble(V,sv,W,sw)
            % function v = varDouble(V,sv,W,sw)
            % returns the estimated variance of a double PCA representation
            
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
               U = V*diag(sv)*W';
            else
               U = V*(sv.*W');
            end
            v = sum((sw.*U').^2)';
        end
        
        function Y = reconstruction(mY,V,sv,Z,r)
            % function Y = reconstruction(mY,V,sv,Z,r)
            % reconstructs the PCA representation Y = mY + V*diag(sv)*Z'
            
            if nargin==5
                V = V(:,1:r);
                sv = sv(1:r);
                Z = Z(:,1:r);
            end
            
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                Yc = V*diag(sv)*Z';
            else
                Yc = V*(sv.*Z');
            end
            
            if isempty(mY)
                Y = Yc;
            else
                if verLessThan('matlab','9.1') % compatibility (<R2016b)
                    Y = bsxfun(@plus,mY,Yc);
                else
                    Y = mY + Yc;
                end
                % Y = repmat(mY,[1,N]) + Yc;
            end
        end
        
        function [Yt,mYt,Vt] = reconstructionAtStep(mY,V,sv,Z,t,varargin)
            % function [Yt,mYt,Vt] = reconstructionAtStep(mY,V,sv,Z,t,varargin)
            % reconstructs the PCA representation Y(t) = mY(t) + V(t)*diag(sv)*Z' 
            % at step t
            
            % R = length(sv);
            % N = size(Z,1);
            Vt = V(:,:,t);
            
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                Yct = Vt*diag(sv)*Z';
            else
                Yct = Vt*(sv.*Z');
            end
            
            if isempty(mY)
                mYt = [];
                Yt = Yct;
            else
                mYt = mY(:,:,t);
                if verLessThan('matlab','9.1') % compatibility (<R2016b)
                    Yt = bsxfun(@plus,mYt,Yct);
                else
                    Yt = mYt + Yct;
                end
                % Yt = repmat(mYt,[1,N]) + Yct;
            end
        end
        
        function [Vt,svt,Wt,Rt] = getPrincipalComponentsDoubleAtStep(V,sv,W,R,t,varargin)
            % function [Vt,svt,Wt,Rt] = getPrincipalComponentsDoubleAtStep(V,sv,W,R,t,varargin)
            % gets the principal component coefficients V(t), corresponding singular values sv(t) and rank R(t) of the first PCA, 
            % and the principal component coefficients W(t) of the second PCA at step t
            
            Rt = R(t);
            pathname = getcharin('pathname',varargin,'.');
            filename = getcharin('filename',varargin,['double_PCA_space_data_t' num2str(t-1) '.mat']);
            if isempty(V)
                load(fullfile(pathname,filename),'Vt');
            else
                Vt = V(:,1:Rt,t);
            end
            svt = sv(1:Rt,t);
            Wt = W(:,1:Rt,t)';
        end
        
        function [Xs,Xa,Xb] = scaling(X)
            % function [Xs,Xa,Xb] = scaling(X)
            % returns the scaled normalized data Xs = (X-Xb)/Xa from the 
            % initial data X with Xa = Xmax-Xmin and Xb = Xmin
            
            Xmin = min(X,[],2);
            Xmax = max(X,[],2);
            Xa = Xmax - Xmin;
            clear Xmax
            Xb = Xmin;
            clear Xmin
            
            %% Implementation with Nans
            % if verLessThan('matlab','9.1') % compatibility (<R2016b)
            %     Xs = bsxfun(@rdvide,bsxfun(@minus,X,Xb),Xa);
            % else
            %     Xs = (X-Xb)./Xa;
            % end
            % isunscaled = isnan(Xs);
            % Xs(isunscaled) = X(isunscaled);
            
            %% Implementation without Nans
            ratio = 1./Xa;
            isunscaled = ~isfinite(ratio);
            ratio(isunscaled) = 1;
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                Xs = bsxfun(@times,bsxfun(@minus,X,Xb),ratio);
            else
                Xs = (X-Xb).*ratio;
            end
            Xs(isunscaled,:) = X(isunscaled,:);
        end
        
        function X = unscaling(Xs,Xa,Xb)
            % function X = unscaling(Xs,Xa,Xb)
            % returns back to the initial data X = Xs*Xa + Xb from the 
            % scaled normalized data Xs with Xa = Xmax-Xmin and Xb = Xmin
            
            if verLessThan('matlab','9.1') % compatibility (<R2016b)
                X = bsxfun(@plus,bsxfun(@times,Xs,Xa),Xb);
            else
                X = Xs.*Xa + Xb;
            end
            
        end
	end
end
