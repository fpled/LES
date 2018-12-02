function gradu = grad(u,Dx)

s = size(u);
gradu = zeros([3,s]);
for i=1:3
    dims = circshift([1 2 3],1-i);
    if ndims(u)==4
        if s(1)==3
            dims = [dims+1,1];
        else
            dims = [dims,4];
        end
    elseif ndims(u)==5
        dims = [dims+1,1,5];
    end
    ui = permute(u,dims);
    si = size(ui);
    gradui = ipermute(reshape(Dx*ui(:,:),si),dims);
    gradu(i,:) = gradui(:);
end

end
