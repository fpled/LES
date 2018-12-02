function tauDiff = compute_tauDiff(gradu,Dx)

S = (gradu+permute(gradu,[2,1,3:ndims(gradu)]))/2;
tauDiff = div(S,Dx);

end
