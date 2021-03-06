function Yt = apply_filter(Yt,filterType,h,varargin)

% N = size(Yt,1);
% n = size(Yt,2);
% parfor l=1:N
%     for i=1:n
%         Ytl = squeeze(Yt(l,i,:,:,:));
%         switch filterType
%             case {'box','mean','average'}
%                 % filterSize = h;
%                 % padding = getcharin('padding',varargin,'replicate');
%                 % normalizationfactor = getcharin('normalizationfactor',varargin,1/prod(filterSize));
%                 % Ytl = imboxfilt3(Ytl,filterSize,'Padding',padding,'NormalizationFactor',normalizationfactor);
%                 Ytl = imfilter(Ytl,h,'replicate');
%             case {'linear','trapz'}
%                 Ytl = imfilter(Ytl,h,'replicate');
%         end
%         Yt(l,i,:,:,:) = Ytl;
%     end
% end

s = size(Yt);
Yt = reshape(Yt,[s(1)*s(2),s(3:end)]);
parfor i=1:size(Yt,1)
    Yti = squeeze(Yt(i,:,:,:));
    switch filterType
        case {'box','mean','average'}
            % filterSize = h;
            % padding = getcharin('padding',varargin,'replicate');
            % normalizationfactor = getcharin('normalizationfactor',varargin,1/prod(filterSize));
            % Yti = imboxfilt3(Yti,filterSize,'Padding',padding,'NormalizationFactor',normalizationfactor);
            Yti = imfilter(Yti,h,'replicate');
        case {'linear','trapz'}
            Yti = imfilter(Yti,h,'replicate');
    end
    Yt(i,:,:,:) = Yti;
end
Yt = reshape(Yt,s);

end
