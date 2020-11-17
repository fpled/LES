function Yt = get_single_PCA_at_step(mY,Phi,sig,X,p,t,varargin)
% function Yt = get_single_PCA_at_step(mY,Phi,sig,X,p,t,varargin)

% R = length(sig);
N = size(X,1);
r = size(Phi,1)/(p+1);
% Phir = reshape(Phi',[R,r,p+1]);
% Phit = Phir(:,:,t)';
Phit = Phi(r*(t-1)+(1:r),:);

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    Yct = Phit*diag(sig)*X';
else
    Yct = Phit*(sig.*X');
end

if isempty(mY)
    Yt = Yct;
else
    mYt = mY(1,:,:,t);
    n = size(mY,2);
    m = size(mY,3);
    Yct = reshape(Yct',[N,n,m]);
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        Yt = bsxfun(@plus,mYt,Yct);
    else
        Yt = mYt + Yct;
    end
    % Yt = repmat(mYt,[N,1,1]) + Yct;
end

end
