%% Analyze between rats...
% this script is largly responcible for producing the figures
% in the analysis.

load("Mats/btRats.mat");

plt = 1;

if plt % plot connectivity matricies for one rat
    h = figure('Position', [700, 100, 700, 225]); hold on;
    subplot(1,3,1);
    imagesc(btRats.phase(1).CIJ{1})
    title('PRE');
    ax = gca;
    set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
    ax.FontSize = 14;
    
    subplot(1,3,2);
    imagesc(btRats.phase(2).CIJ{1})
    title({'Connectivity Matricies'; 'Maze'})
    ax = gca;
    set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
    ax.FontSize = 14;
    
    subplot(1,3,3);
    imagesc(btRats.phase(3).CIJ{1})
    title('POST');
    ax = gca;
    set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
    ax.FontSize = 14;
    
    
    figName = "Figs/CIJ/achCIJ";
    savefig(h, char(figName+".fig"));
    saveas(h,  char(figName+".png"));
end

if plt % Average similarity increases?
    h = figure; hold on;
    cVec = [[0 1 0]; [1 0 1]; [0 1 1]; [1 0 0]];
    for zz = 1:size(btRats.mSim,1)
        hold on;
        scatter(1:3,btRats.mSim(zz,:),70,cVec(zz,:),'filled');
        line(1:3,btRats.mSim(zz,:),'Color',cVec(zz,:))
        %plot(btRats.degs(zz,:)', 'LineWidth',3);
        
        ylabel('Average Similarity Matrix Values');
        set(gca,'XTick',[0 1 2 3 4])
        names = {'PRE';'MAZE';'POST'};
        set(gca,'xtick',[1:3],'xticklabel',names)
        xlim([.5 3.5]);
        ax = gca; ax.FontSize = 14;
    end
    
    figName = "Figs/avgSim";
    savefig(h, char(figName+".fig"));
    saveas(h,  char(figName+".png"));
end

if plt % Average degree increases?
    h = figure; hold on;
    cVec = [[0 1 0]; [1 0 1]; [0 1 1]; [1 0 0]];
    for zz = 1:size(btRats.degs,1)
        hold on;
        err = btRats.degsSEM(zz,:);
        errorbar(1:3, btRats.degs(zz,:), err,'k');
        scatter(1:3,btRats.degs(zz,:),70,cVec(zz,:),'filled');
        line(1:3,btRats.degs(zz,:),'Color',cVec(zz,:));
        
        %plot(btRats.degs(zz,:)', 'LineWidth',3);
        
        ylabel('Average Degree');
        set(gca,'XTick',[0 1 2 3 4])
        names = {'PRE';'MAZE';'POST'};
        set(gca,'xtick',[1:3],'xticklabel',names)
        xlim([.5 3.5]); ylim([0 210]);
        ax = gca; ax.FontSize = 14;
    end
    
    figName = "Figs/degDist/avgDeg";
    savefig(h, char(figName+".fig"));
    saveas(h,  char(figName+".png"));
end

if plt % plot community matricies for one rat
    h = figure('Position', [700, 700, 225, 700]); hold on;
    subplot(3,1,1);
    imagesc(btRats.phase(1).CIJCM{1})
    title({'Communities'; 'PRE'});
    ax = gca;
    set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
    ax.FontSize = 14;
    
    subplot(3,1,2);
    imagesc(btRats.phase(2).CIJCM{1})
    title({'Maze'})
    ax = gca;
    set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
    ax.FontSize = 14;
    
    subplot(3,1,3);
    imagesc(btRats.phase(3).CIJCM{1})
    title('POST');
    ax = gca;
    set(gca,'YTickLabel',[]); set(gca,'XTickLabel',[]);
    ax.FontSize = 14;
    
    
    figName = "Figs/CIJ/achCIJ_commV";
    savefig(h, char(figName+".fig"));
    saveas(h,  char(figName+".png"));
end


if plt % Average number of communities?
    h = figure; hold on;
    cVec = [[0 1 0]; [1 0 1]; [0 1 1]; [1 0 0]];
    for zz = 1:size(btRats.meanci,1)
        hold on; % Really low error
        err = btRats.meanciSEM(zz,:);
        errorbar(1:3, btRats.meanci(zz,:), err,'k');
        scatter(1:3,btRats.meanci(zz,:),70,cVec(zz,:),'filled');
        line(1:3,btRats.meanci(zz,:),'Color',cVec(zz,:))
        
        ylabel('Average Number of Communities');
        set(gca,'XTick',[0 1 2 3 4])
        names = {'PRE';'MAZE';'POST'};
        set(gca,'xtick',[1:3],'xticklabel',names)
        xlim([.5 3.5]); ylim([1 37])
        ax = gca; ax.FontSize = 14;
    end
    
    figName = "Figs/avgCi";
    savefig(h, char(figName+".fig"));
    saveas(h,  char(figName+".png"));
end

if plt % Average q?
    h = figure; hold on;
    cVec = [[0 1 0]; [1 0 1]; [0 1 1]; [1 0 0]];
    for zz = 1:size(btRats.meanq,1)
        hold on; % Really low error
        err = btRats.meanqSEM(zz,:);
        errorbar(1:3, btRats.meanq(zz,:), err,'k');
        scatter(1:3,btRats.meanq(zz,:),70,cVec(zz,:),'filled');
        line(1:3,btRats.meanq(zz,:),'Color',cVec(zz,:))
        
        ylabel('Average Q Value');
        set(gca,'XTick',[0 1 2 3 4])
        names = {'PRE';'MAZE';'POST'};
        set(gca,'xtick',[1:3],'xticklabel',names)
        xlim([.5 3.5]);
        ax = gca; ax.FontSize = 14;
    end
    
    figName = "Figs/avgQ";
    savefig(h, char(figName+".fig"));
    saveas(h,  char(figName+".png"));
end




if plt % plot degree distribution
    cVec = [[0 1 0]; [1 0 1]; [0 1 1]; [1 0 0]];
    h = figure;
    subplot(3,1,1); hold on;
    cellfun(@plot, btRats.phase(1).deg); title({' Degree Distribution'; 'PRE'})
    ylabel('Degree');
    ax = gca; ax.FontSize = 14;
    
    subplot(3,1,2); hold on;
    cellfun(@plot,btRats.phase(2).deg); title('MAZE')
    ylabel('Degree');
    ax = gca; ax.FontSize = 14;
    
    subplot(3,1,3); hold on;
    cellfun(@plot,btRats.phase(3).deg); title('POST')
    ylabel('Degree'); xlabel('Node');
    ax = gca; ax.FontSize = 14;
    
    figName = "Figs/degDist/degs";
    savefig(h, char(figName+".fig"));
    saveas(h,  char(figName+".png"));
    
end


%% Visuals testing (depriciated)
if 0
    CIJ_threshold = threshold_proportional(CIJ,.08);
    
    [u,v] = find(triu(CIJ,1));
    G = graph(u,v);
    h = plot(G,'Layout','force','Iterations',100);
    
    pos = [h.XData(:),h.YData(:),h.ZData(:)];  % extract coordinates
    [ex,ey,ez] = adjacency_plot_und(CIJ_threshold,pos);
    
    % Replot with communities
    plot3(ex,ey,ez,'k'); hold on;
    scatter3(pos(:,1),pos(:,2),pos(:,3),10,ciconsensus,'filled'); % plot each brain regions as a circle with size proportional to strength and color their community
    hold off; axis image; view([0,90]);
end


%%

if 0
    % Grab null model
    % The lack of correlations means that either the model
    % was wrong, or that the matricies are very particular...
    np = ratMats.null.pre;
    np(np<0) = 0;
    np(isnan(np)) = 0;
    mean2(np)
    [density,numNodes,numEdges] = density_dir(np);
    [id,od,deg] = degrees_dir(np);
    
    % To remove diag:
    %JIC(find(eye(size(JIC)))) = 0;
    
    % subtract null
    %figure; imagesc(CIJ - JIC)
    
end

%%

disp('fin')






