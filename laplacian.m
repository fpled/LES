function laplacianu = laplacian(u,Dxx)

s = size(u);
laplacianu = zeros(s);
for i=1:3
    order = [i+1,setdiff(1:ndims(u),i+1)];
    ui = permute(u,order);
    si = size(ui);
    laplacianui = reshape(Dxx*ui(:,:),si);
    laplacianui = ipermute(laplacianui,order);
    laplacianu(:,:) = laplacianu(:,:) + laplacianui(:,:);
end

end
