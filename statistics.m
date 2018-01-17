clc
clearvars
close all

tol = 1e-12;

for g=2.^(4:7)
    load(fullfile('Data',['data' num2str(g) '.mat']),'Y');
    
    mY = mean(Y);
    [U,S,V,err] = svdtruncate(Y,tol,1);
end