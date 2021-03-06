function divu = div(u,Dx)

s = size(u);
s = s(2:end);
divu = zeros(s);
for i=1:3
    if s(1)==3
        ui = reshape(u(:,i,:),s);
        order = [i+1,setdiff(1:ndims(ui),i+1)];
    else
        ui = reshape(u(i,:),s);
        order = [i,setdiff(1:ndims(ui),i)];
    end
    ui = permute(ui,order);
    si = size(ui);
    divui = reshape(Dx{i}*ui(:,:),si);
    divui = ipermute(divui,order);
    if s(1)==3
        divui = sum(divui);
        divu(i,:) = divui(:,:);
    else
        divu(:,:) = divu(:,:) + divui(:,:);
    end
end
if s(1)~=3
    divu = shiftdim(divu,-1);
end

end
