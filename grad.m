function gradu = grad(u,Dx)

s = size(u);
if s(1)==1
    gradu = zeros([3,s(2:end)]);
else
    gradu = zeros([3,s]);
end
for i=1:3
    order = [i+1,setdiff(1:ndims(u),i+1)];
    ui = permute(u,order);
    si = size(ui);
    gradui = reshape(Dx{i}*ui(:,:),si);
    gradui = ipermute(gradui,order);
    gradu(i,:) = gradui(:);
end

end
