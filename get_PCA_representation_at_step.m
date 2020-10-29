function Yt = get_PCA_representation_at_step(mY,sig,V,s,W,X,R,t,varargin)
% function Yt = get_PCA_representation_at_step(mY,sig,V,s,W,X,R,t,varargin)

Rmax = max(R);
Rt = R(t);
sigt = sig(1:Rt,t);
% Zt = Z(:,1:Rt,t);
gridpathname = getcharin('gridpathname',varargin,'coord');
if isempty(V)
    load(fullfile(gridpathname,['PCAspace_t' num2str(t-1) '.mat']),'Vt');
else
    Vt = V(:,1:Rt,t);
end
index = getcharin('index',varargin,'coord');
switch index
    case 'coord'
        Wt = W(Rmax*(t-1)+(1:Rt),:);
    case 'time'
        p = length(R)-1;
        Wt = W(t:(p+1):end,:);
        Wt = Wt(1:Rt,:);
end

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    Zct = Wt*diag(s)*X'; % Zct = Zt';
    Yct = Vt*diag(sigt)*Zct;
else
    Zct = Wt*(s.*X'); % Zct = Zt';
    Yct = Vt*(sigt.*Zct);
end

if isempty(mY)
    Yt = Yct;
else
    mYt = mY(1,:,:,t);
    N = size(Yct,2);
    n = size(mY,2);
    m = size(mY,3);
    Yt = repmat(mYt,[N,1,1]) + reshape(Yct',[N,n,m]);
end

end
