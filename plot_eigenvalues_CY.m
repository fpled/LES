function plot_eigenvalues_CY(times,dt,nCols,g,gridpathname,fontsize)

if g<2^7
    load(fullfile(gridpathname,'solution.mat'),'Sig','R');
end

figure('Name','Evolution of eigenvalues')
clf
leg = cell(1,length(times));
hdl = zeros(1,length(times));
color = distinguishable_colors(length(times));
c = 0;
for t=times
    c = c+1;
    if g<2^7
        Rt = R(t+1);
        hdl(c) = semilogy(1:Rt,Sig(1:Rt,t+1).^2,'LineStyle','-','Color',color(c,:),'LineWidth',1);
    else
        load(fullfile(gridpathname,['PCA_t' num2str(t) '.mat']),'Sigt','Rt');
        hdl(c) = semilogy(1:Rt,Sigt(1:Rt).^2,'LineStyle','-','Color',color(c,:),'LineWidth',1);
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
gridLegend(hdl,nCols,leg,'Location','BestOutside');
% legend(leg{:},'Location','NorthEastOutside')

end
