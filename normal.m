function n = normal(gradC)

s = size(gradC);
n = zeros(s);
ngradC = repmat(sqrt(sum(gradC.^2,1)),[s(1),ones(1,ndims(gradC)-1)]);
ind = find(ngradC);
n(ind) = gradC(ind)./ngradC(ind);

end
