function n = normal(gradC)

% s = size(gradC);
% n = zeros(s);
% ngradC = repmat(sqrt(sum(gradC.^2,1)),[s(1),ones(1,ndims(gradC)-1)]);
% ind = find(ngradC);
% n(ind) = gradC(ind)./ngradC(ind);

ngradC = sqrt(sum(gradC.^2,1));
if verLessThan('matlab','9.1') % compatibility (<R2016b)
    n = bsxfun(@rdivide,gradC,ngradC);
else
    n = gradC./ngradC;
end
n(isnan(n)) = 0;

end
