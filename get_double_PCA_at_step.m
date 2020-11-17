function Yt = get_double_PCA_at_step(mY,V,sig,mZ,W,s,X,R,t,varargin)
% function Yt = get_double_PCA_at_step(mY,V,sig,mZ,W,s,X,R,t,varargin)

Rmax = max(R);
N = size(X,1);
Rt = R(t);
gridpathname = getcharin('gridpathname',varargin,'.');
if isempty(V)
    load(fullfile(gridpathname,['PCAspace_t' num2str(t-1) '.mat']),'Vt');
else
    Vt = V(:,1:Rt,t);
end
sigt = sig(1:Rt,t);
% p = length(R)-1;
if ~isempty(mZ)
    % mZr = reshape(mZ,[1,Rmax,p+1]);
    % mZt = mZr(:,1:Rt,t)';
    mZt = mZ(Rmax*(t-1)+(1:Rt),:);
end
% Q = length(s);
% Wr = reshape(W',[Q,Rmax,p+1]);
% Wt = Wr(:,1:Rt,t)';
Wt = W(Rmax*(t-1)+(1:Rt),:);

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    Zct = Wt*diag(s)*X';
else
    Zct = Wt*(s.*X');
end

if isempty(mZ)
    Zt = Zct;
else
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        Zt = bsxfun(@plus,mZt,Zct);
    else
        Zt = mZt + Zct;
    end
    % Zt = repmat(mZt,[1,N]) + Zct;
end

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    Yct = Vt*diag(sigt)*Zt;
else
    Yct = Vt*(sigt.*Zt);
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
