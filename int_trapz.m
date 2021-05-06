function Qu = int_trapz(x,C,u)

dim = numel(x);
% s = size(u);
% N = s(1);
% n = s(2);
% sx = s(2+(1:dim));
switch dim
    case 1
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            Qu = cat(3,trapz(x{1},bsxfun(@times,1-C,u),3),trapz(x{1},bsxfun(@times,C,u),3));
        else
            Qu = cat(3,trapz(x{1},(1-C).*u,3),trapz(x{1},C.*u,3));
        end
        % Qu = cat(3,trapz(x{1},repmat(1-C,[1,n,ones(1,dim)]).*u,3),trapz(x{1},repmat(C,[1,n,ones(1,dim)]).*u,3));
    case 2
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            Qu = cat(3,trapz(x{2},trapz(x{1},bsxfun(@times,1-C,u),3),4),...
                trapz(x{2},trapz(x{1},bsxfun(@times,C,u),3),4));
        else
            Qu = cat(3,trapz(x{2},trapz(x{1},(1-C).*u,3),4),...
                trapz(x{2},trapz(x{1},C.*u,3),4));
        end
        % Qu = cat(3,trapz(x{2},trapz(x{1},repmat(1-C,[1,n,ones(1,dim)]).*u,3),4),...
        %     trapz(x{2},trapz(x{1},repmat(C,[1,n,ones(1,dim)]).*u,3),4));
    case 3
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            Qu = cat(3,trapz(x{3},trapz(x{2},trapz(x{1},bsxfun(@times,1-C,u),3),4),5),...
                trapz(x{3},trapz(x{2},trapz(x{1},bsxfun(@times,C,u),3),4),5));
        else
            Qu = cat(3,trapz(x{3},trapz(x{2},trapz(x{1},(1-C).*u,3),4),5),...
                trapz(x{3},trapz(x{2},trapz(x{1},C.*u,3),4),5));
        end
        % Qu = cat(3,trapz(x{3},trapz(x{2},trapz(x{1},repmat(1-C,[1,n,ones(1,dim)]).*u,3),4),5),...
        %     trapz(x{3},trapz(x{2},trapz(x{1},repmat(C,[1,n,ones(1,dim)]).*u,3),4),5));
    otherwise
        error(['Wrong spatial dimension ' num2str(dim)])
end

end
