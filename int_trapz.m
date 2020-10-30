function Qu = int_trapz(x,C,u,dim)

% s = size(u);
% N = s(1);
% n = s(2);
% sx = s(3+(0:dim));
switch dim
    case 1
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            Qu = cat(3,trapz(x,bsxfun(@times,1-C,u),3),trapz(x,bsxfun(@times,C,u),3));
        else
            Qu = cat(3,trapz(x,(1-C).*u,3),trapz(x,C.*u,3));
        end
        % Qu = cat(3,trapz(x,repmat(1-C,[1,n,ones(1,dim)]).*u,3),trapz(x,repmat(C,[1,n,ones(1,dim)]).*u,3));
    case 2
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            Qu = cat(3,trapz(x,trapz(x,bsxfun(@times,1-C,u),3),4),...
                trapz(x,trapz(x,bsxfun(@times,C,u),3),4));
        else
            Qu = cat(3,trapz(x,trapz(x,(1-C).*u,3),4),...
                trapz(x,trapz(x,C.*u,3),4));
        end
        % Qu = cat(3,trapz(x,trapz(x,repmat(1-C,[1,n,ones(1,dim)]).*u,3),4),...
        %     trapz(x,trapz(x,repmat(C,[1,n,ones(1,dim)]).*u,3),4));
    case 3
        if verLessThan('matlab','9.1') % compatibility (<R2016b)
            Qu = cat(3,trapz(x,trapz(x,trapz(x,bsxfun(@times,1-C,u),3),4),5),...
                trapz(x,trapz(x,trapz(x,bsxfun(@times,C,u),3),4),5));
        else
            Qu = cat(3,trapz(x,trapz(x,trapz(x,(1-C).*u,3),4),5),...
                trapz(x,trapz(x,trapz(x,C.*u,3),4),5));
        end
        % Qu = cat(3,trapz(x,trapz(x,trapz(x,repmat(1-C,[1,n,ones(1,dim)]).*u,3),4),5),...
        %     trapz(x,trapz(x,trapz(x,repmat(C,[1,n,ones(1,dim)]).*u,3),4),5));
    otherwise
        error(['Wrong spatial dimension ' num2str(dim)])
end

end
