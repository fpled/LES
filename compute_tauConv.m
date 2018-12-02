function tauConv = compute_tauConv(u,gradu)

tauConv = sum(gradu.*shiftdim(u,-1),2);
s = size(tauConv);
tauConv = reshape(tauConv,[s(1),s(3:end)]);

end
