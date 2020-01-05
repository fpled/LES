function divS = div(S,Dx)

s = size(S);
s = s(2:end);
divS = zeros(s);
for i=1:3
    dims = circshift([1 2 3],1-i,2);
    if ndims(S)==4 || (ndims(S)==5 && s(1)~=3)
        Si = reshape(S(i,:),s);
        if ndims(S)==5
            dims = [dims,4];
        end
        Si = permute(Si,dims);
        divSi = ipermute(reshape(Dx*Si(:,:),s),dims);
        divS = divS + divSi;
    else
        si = s(2:end);
        Si = reshape(sum(S(i,:,:),2),si);
        if ndims(S)==6
            dims = [dims,4];
        end
        Si = permute(Si,dims);
        divSi = ipermute(reshape(Dx*Si(:,:),si),dims);
        divS(i,:) = divSi(:);
    end
end

end
