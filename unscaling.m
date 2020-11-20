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

