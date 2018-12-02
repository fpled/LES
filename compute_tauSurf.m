function tauSurf = compute_tauSurf(gradC,kappa)

tauSurf = shiftdim(kappa,-1).*gradC;

end
