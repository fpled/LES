function n = normal(gradC)

s = size(gradC);
n = zeros(s);
normgradC = repmat(sqrt(sum(gradC.^2,1)),[s(1),ones(1,ndims(gradC)-1)]);
ind = find(normgradC);
n(ind) = gradC(ind)./normgradC(ind);

end
