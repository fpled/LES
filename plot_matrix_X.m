figure('Name','Matrix')
clf
imagesc(X)
colorbar
axis image
set(gca,'FontSize',fontsize)
xlabel('$\beta$','Interpreter',interpreter)
ylabel('$\ell$','Interpreter',interpreter)
title('Matrix $[X]_{\ell,\beta}$ for $\ell=1,\dots,N$ and $\beta=1,\dots,Q$','Interpreter',interpreter)
mysaveas(gridpathname,'matrix_X',formats,renderer);
mymatlab2tikz(gridpathname,'matrix_X.tex');