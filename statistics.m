clc
clearvars
close all

tol = 1e-12;

% for g=2.^(4:7)
for g=2.^4
    gridname = ['Grid' num2str(g)];
    load(fullfile(gridname,['data' num2str(g) '.mat']),'Y');
    
    mY = mean(Y);
    [U,S,V,err] = svdtruncate(Y,tol,1);
    
    %save(fullfile(gridname,['res' num2str(g) '.mat']),'mY');
end