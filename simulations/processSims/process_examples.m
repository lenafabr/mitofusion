filename = '../examples/example_cog.snap.out';

% for each snapshot, we save information on each mitochondrion in the
% domain (here, cluster refers to a mitochondrion)
% read in positions of all clusters (positions)
% health content of all clusters (proteins)
% time points (timevals)
% number of clusters in each snapshot (nclust)
% id numbers of all the clusters (ids)
% states (-1=retro, 1 = antero, 0 = stopped)
%%
[positions, proteins,timevals,nregions,nclust,~,states,idtags] = readSnap(filename);

%% Watch animation of the simulated  mitochondria
cmat = BBVYWcolormap(100);


idmax = max(cellfun(@(x) max(x), idtags));
yvals = rand(idmax,1);

ct = 0;
nreg = nregions(1);

for fc = 1:10:length(positions)
    fc
    %%
    nclustcur = nclust(fc);
    pos = positions{fc};
    prots = proteins{fc};
    id = idtags{fc};
    state = states{fc};
    
    walkind = find(abs(state)==1);
    stopind = find(state==0);
    
    scatter(pos(walkind),yvals(id(walkind)),50,prots(walkind),'LineWidth',3)
    colormap(cmat(1:90,:));
    cb = colorbar;
    cb.LineWidth=2;
    cb.TickLabelInterpreter = 'latex';
    caxis([0,1])
    
    hold all
  %  quiver(pos(walkind),yvals(id(walkind)),state(walkind)*0.02,zeros(size(walkind)),'LineWidth',2,'MaxHeadSize',10,'AutoScale','off','Color','m')
        
    scatter(pos(stopind),yvals(id(stopind)),50,prots(stopind),'filled')
    
    rectangle('Position',[0 -0.1 1 1.2],'EdgeColor','k','LineWidth',3)
      
     set(gcf,'Color','w','Position',[-1557         392        1328         185])
    set(gca,'Visible','off','defaultTextInterpreter','latex','FontSize',16)
    title('Space Shuttle model. $M=100, n=5, \widehat{p}_s = 0.4, f_s = 0.5, \widehat{k}_d = 0.6$')
    hold off
    set(findall(gca, 'type', 'text'), 'visible', 'on')    
    
    set(gcf,'Position',[100 100 539*2 32*2])
    set(gca,'Position',[0.01 0.01 0.98 0.98])
    drawnow    
end
