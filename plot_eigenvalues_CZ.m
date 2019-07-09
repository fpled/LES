function plot_eigenvalues_CZ(times,dt,g,gridpathname,fontsize)

if g<2^7
    load(fullfile(gridpathname,'solution.mat'),'sig','R');
end

figure('Name','Evolution of eigenvalues')
clf
leg = cell(1,length(times));
c = 1;
for t=times
    c = c+1;
    if g<2^7
        Rt = R(t+1);
        semilogy(1:Rt,sig(1:Rt,t+1).^2,'LineStyle','-','Color',getfacecolor(c),'LineWidth',1);
    else
        load(fullfile(gridpathname,['PCA_t' num2str(t) '.mat']),'sigt','Rt');
        semilogy(1:Rt,sigt(:).^2,'LineStyle','-','Color',getfacecolor(c),'LineWidth',1);
    end
    leg{c} = ['t = ' num2str(t*100*dt) ' s'];
    hold on
end
hold off
grid on
box on
set(gca,'FontSize',fontsize)
xlabel('$\alpha$','Interpreter','latex')
ylabel('$\lambda_{\alpha}(t)$','Interpreter','latex')
legend(leg{:},'Location','NorthEastOutside')

end
