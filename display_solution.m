
tf = p*dt;
T = TIMEMODEL(0,tf,p);

Mscal = final(M,DDL('C'));

ampl = 0;
% ampl = getsize(M)/max(abs(mut))/5;
az_z = -37.5; % azimuth for the default 3D view with vertical z-axis
el_z = 30; % elevation for the default 3D view with vertical z-axis
az = az_z-90; % azimuth for the default 3D view with vertical y-axis
el = -el_z; % elevation for the default 3D view with vertical y-axis

ti = 11;

%% Mean
mU = TIMEMATRIX(mU,T);
mC = TIMEMATRIX(mC,T);
mtauTime = TIMEMATRIX(mtauTime,T);
mdivtauConv = TIMEMATRIX(mdivtauConv,T);
mdivtauDiff = TIMEMATRIX(mdivtauDiff,T);
mtauSurf = TIMEMATRIX(mtauSurf,T);
mtauInterf = TIMEMATRIX(mtauInterf,T);

for i=1:3
    evolSolution(M,mU,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_mean_u' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
    evolSolution(M,mtauTime,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_mean_tauTime' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
    evolSolution(M,mdivtauConv,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_mean_divtauConv' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
    evolSolution(M,mdivtauDiff,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_mean_divtauDiff' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
    evolSolution(M,mtauSurf,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_mean_tauSurf' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
end
evolSolution(Mscal,mC,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename','evol_mean_C','pathname',gridpathname,'FrameRate',framerate,...
    'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
evolSolution(Mscal,mtauInterf,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename','evol_mean_tauInterf','pathname',gridpathname,'FrameRate',framerate,...
    'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);

for i=1:3
    figure('Name',['Mean of velocity u_' num2str(i)])
    clf
    plot_sol(M,getmatrixatstep(mU,ti+1),'displ',i,'ampl',ampl);
    title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['mean_u' num2str(i) '_t' num2str(ti*100)],formats,renderer);
    
    figure('Name',['Mean of tauTime_' num2str(i)])
    clf
    plot_sol(M,getmatrixatstep(mtauTime,ti+1),'displ',i,'ampl',ampl);
    title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['mean_tauTime' num2str(i) '_t' num2str(ti*100)],formats,renderer);
    
    figure('Name',['Mean of div(tauConv)_' num2str(i)])
    clf
    plot_sol(M,getmatrixatstep(mdivtauConv,ti+1),'displ',i,'ampl',ampl);
    title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['mean_divtauConv' num2str(i) '_t' num2str(ti*100)],formats,renderer);
    
    figure('Name',['Mean of div(tauDiff)_' num2str(i)])
    clf
    plot_sol(M,getmatrixatstep(mdivtauDiff,ti+1),'displ',i,'ampl',ampl);
    title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['mean_divtauDiff' num2str(i) '_t' num2str(ti*100)],formats,renderer);
    
    figure('Name',['Mean of tauSurf ' num2str(i)])
    clf
    plot_sol(M,getmatrixatstep(mtauSurf,ti+1),'displ',i,'ampl',ampl);
    title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['mean_tauSurf' num2str(i) '_t' num2str(ti*100)],formats,renderer);
end

figure('Name','Mean of indicator function C')
clf
plot_sol(Mscal,getmatrixatstep(mC,ti+1));
title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
colormap(cmap)
colorbar
axis on
box on
view(az,el)
camup([0 1 0])
set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
mysaveas(gridpathname,['mean_C_t' num2str(ti*100)],formats,renderer);

figure('Name','Mean of tauInterf')
clf
plot_sol(Mscal,getmatrixatstep(mtauInterf,ti+1));
title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
colormap(cmap)
colorbar
axis on
box on
view(az,el)
camup([0 1 0])
set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
mysaveas(gridpathname,['mean_tauInterf_t' num2str(ti*100)],formats,renderer);

%% Variance
vU = TIMEMATRIX(vU,T);
vC = TIMEMATRIX(vC,T);
vtauTime = TIMEMATRIX(vtauTime,T);
vdivtauConv = TIMEMATRIX(vdivtauConv,T);
vdivtauDiff = TIMEMATRIX(vdivtauDiff,T);
vtauSurf = TIMEMATRIX(vtauSurf,T);
vtauInterf = TIMEMATRIX(vtauInterf,T);

for i=1:3
    evolSolution(M,vU,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_var_u' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
    evolSolution(M,vtauTime,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_var_tauTime' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
    evolSolution(M,vdivtauConv,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_var_divtauConv' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
    evolSolution(M,vdivtauDiff,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_var_divtauDiff' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
    evolSolution(M,vtauSurf,'displ',i,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename',['evol_var_tauSurf' num2str(i)],'pathname',gridpathname,'FrameRate',framerate,...
        'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
end
evolSolution(Mscal,vC,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename','evol_var_C','pathname',gridpathname,'FrameRate',framerate,...
    'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);
evolSolution(Mscal,vtauInterf,'colormap',cmap,'view',[az,el],'camup',[0 1 0],'filename','evol_var_tauInterf','pathname',gridpathname,'FrameRate',framerate,...
    'axison',true,'boxon',true,'boxstylefull',true,'noxtick',true,'noytick',true,'noztick',true);

for i=1:3
    figure('Name',['Variance of velocity u_' num2str(i)])
    clf
    plot_sol(M,getmatrixatstep(vU,ti+1),'displ',i,'ampl',ampl);
    title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['var_u' num2str(i) '_t' num2str(ti*100)],formats,renderer);
    
    figure('Name',['Variance of tauTime_' num2str(i)])
    clf
    plot_sol(M,getmatrixatstep(vtauTime,ti+1),'displ',i,'ampl',ampl);
    title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['var_tauTime' num2str(i) '_t' num2str(ti*100)],formats,renderer);
    
    figure('Name',['Variance of div(tauConv)_' num2str(i)])
    clf
    plot_sol(M,getmatrixatstep(vdivtauConv,ti+1),'displ',i,'ampl',ampl);
    title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['var_divtauConv' num2str(i) '_t' num2str(ti*100)],formats,renderer);
    
    figure('Name',['Variance of div(tauDiff)_' num2str(i)])
    clf
    plot_sol(M,getmatrixatstep(vdivtauDiff,ti+1),'displ',i,'ampl',ampl);
    title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['var_divtauDiff' num2str(i) '_t' num2str(ti*100)],formats,renderer);
    
    figure('Name',['Variance of tauSurf_' num2str(i)])
    clf
    plot_sol(M,getmatrixatstep(vtauSurf,ti+1),'displ',i,'ampl',ampl);
    title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
    colormap(cmap)
    colorbar
    axis on
    box on
    view(az,el)
    camup([0 1 0])
    set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
    mysaveas(gridpathname,['var_tauSurf' num2str(i) '_t' num2str(ti*100)],formats,renderer);
end

figure('Name','Variance of indicator function C')
clf
plot_sol(Mscal,getmatrixatstep(vC,ti+1));
title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
colormap(cmap)
colorbar
axis on
box on
view(az,el)
camup([0 1 0])
set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
mysaveas(gridname,['var_C_t' num2str(ti*100)],formats,renderer);

figure('Name','Variance of tauInterf')
clf
plot_sol(Mscal,getmatrixatstep(vtauInterf,ti+1));
title(['time ' num2str(ti*dt,'%.2f') ' s'],'FontSize',fontsize)
colormap(cmap)
colorbar
axis on
box on
view(az,el)
camup([0 1 0])
set(gca,'FontSize',fontsize,'BoxStyle','full','XTick',[],'YTick',[],'ZTick',[])
mysaveas(gridname,['var_tauInterf_t' num2str(ti*100)],formats,renderer);
