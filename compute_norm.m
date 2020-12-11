function normu = compute_norm(x,u,dim)

% s = size(u);
% N = s(1);
% sx = s(1+(1:dim));
switch dim
    case 1
        normu = sqrt(mean(trapz(x,u,2),1));
    case 2
        normu = sqrt(mean(trapz(trapz(x,u,2),3),1));
    case 3
        normu = sqrt(mean(trapz(trapz(trapz(x,u,2),3),4),1));
    otherwise
        error(['Wrong spatial dimension ' num2str(dim)])
end

end
