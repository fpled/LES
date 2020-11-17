function X = unscaling(Xs,Xa,Xb)
% function X = unscaling(Xs,Xa,Xb)

if verLessThan('matlab','9.1') % compatibility (<R2016b)
    X = bsxfun(@plus,bsxfun(@times,Xs,Xa),Xb);
else
    X = Xs.*Xa + Xb;
end

end

