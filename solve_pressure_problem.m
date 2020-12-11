function pressure = solve_pressure_problem(B,Rho,Gradinvrho,GradxN,GradyN,GradzN,LaplacianN)
% function pressure = solve_pressure_problem(B,Rho,Gradinvrho,GradxN,GradyN,GradzN,LaplacianN)

% [m,N] = size(B);
[m,N] = size(Rho);
pressure = zeros(m,N);
if ~verLessThan('matlab','9.2') % introduced in R2017a
    q = parallel.pool.DataQueue;
    afterEach(q, @nUpdateProgressBar);
else
    q = [];
end
textprogressbar('Solving pressure problem: ');
i = 0;
parfor l=1:N
    if ~verLessThan('matlab','9.2') % introduced in R2017a
        send(q,l);
    end
    b = B(:,l);
    rho = Rho(:,l);
    gradinvrho = Gradinvrho(:,:,l);
    if verLessThan('matlab','9.1') % compatibility (<R2016b)
        A = bsxfun(@times,gradinvrho(1,:),GradxN) + bsxfun(@times,gradinvrho(2,:),GradyN) + bsxfun(@times,gradinvrho(3,:),GradzN);
        A = A + bsxfun(@rdivide,LaplacianN,rho);
    else
        A = gradinvrho(1,:).*GradxN + gradinvrho(2,:).*GradyN + gradinvrho(3,:).*GradzN;
        A = A + LaplacianN./rho;
    end
    % A = repmat(gradinvrho(1,:),[m,1]).*GradxN + repmat(gradinvrho(2,:),[m,1]).*GradyN + repmat(gradinvrho(3,:),[m,1]).*GradzN;
    % A = A + LaplacianN./repmat(rho,[1,m]);
    p = A\b;
    pressure(:,l) = p;
end
textprogressbar(' done');

function nUpdateProgressBar(~)
i = i+1;
textprogressbar(i/N*100,sprintf('(%d/%d)',i,N));
end

end