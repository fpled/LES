function [Xs,Xa,Xb] = scaling(X)
% function [Xs,Xa,Xb] = scaling(X)

Xmin = min(X);
Xmax = max(X);
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
Xs(:,isunscaled) = X(:,isunscaled);

end

