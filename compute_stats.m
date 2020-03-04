function [mQ,stdQ,RQ,IQ] = compute_stats(Q,dt)

p = size(Q,4)-1;
N = size(Q,1);
n = size(Q,2);
mQ = mean(Q,1);
Qc = Q - repmat(mQ,N,1,1,1);
stdQ = zeros(n*2,p+1);
for t=0:p
    Qt = Q(:,:,:,t+1);
    stdQ(:,t+1) = std(Qt(:,:));
    % Qct = Qc(:,:,:,t+1);
    % stdQ(:,t+1) = sqrt(1/(N-1)*sum(Quct(:,:).^2));
end
RQ = zeros(n*2,n*2,p+1,p+1);
for t=0:p
    Qct = Qc(:,:,:,t+1);
    Qct = Qct(:,:)./repmat(stdQ(:,t+1)',[N,1]);
    for tt=0:p
        Qctt = Qc(:,:,:,tt+1);
        Qctt = Qctt(:,:)./repmat(stdQ(:,tt+1)',[N,1]);
        RQ(:,:,t+1,tt+1) = 1/(N-1)*Qct(:,:)'*Qctt(:,:);
    end
end
IQ = zeros(n*2,n*2,p+1);
for t=0:p-1
    fun = zeros(n*2,n*2,p-t+1);
    for tt=0:p-t
        fun(:,:,tt+1) = RQ(:,:,tt+t+1,tt+1);
    end
    IQt = 1/((p-t)*dt)*trapz((0:p-t)*dt,fun,3);
    IQ(:,:,t+1) = IQt;
end
IQ = reshape(IQ,n,2,n,2,p+1);

end
