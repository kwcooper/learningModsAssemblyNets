
load("Mats/ratMats_Achilles_11012013.mat")
ratName = "Achilles_11012013";

degs = zeros(4,3);
tic
for rat = 1:4
    %load('NoveltySessInfoMatFiles/Achilles_10252013_sessInfo.mat')
    
    if rat == 1
        load("Mats/ratMats_Achilles_11012013.mat")
        ratName = "Achilles_11012013";
    end
    
    if rat == 2
        load("Mats/ratMats_Buddy_06272013.mat")
        ratName = "Buddy_06272013";
    end
    
    if rat == 3
        load("Mats/ratMats_Cicero_09012014.mat")
        ratName = "Cicero_09012014";
    end
    
    if rat == 4
        load("Mats/ratMats_Gatsby_08022013.mat")
        ratName = "Gatsby_08022013";
    end
    
    disp(strjoin(["Phase:", ratName]))
    
    % For each phase, run network analysis
    netStats = struct;
    fields = {'pre','maz','pst'};
    for i = 1:3
        disp(strjoin(["Running", fields{i}]))
        keyboard; 
        
        % Grab matrix
        CIJ = squeeze(ratMats.(fields{i}));
        if 0
        h = figure; imagesc(CIJ); title("Origional"); 
        figName = "Figs/CIJ/" + ratName + "_" + fields{i} + "_CIJ";
        savefig(h, char(figName+".fig"));
        saveas(h,  char(figName+".png"));
        end
        
        CIJ(CIJ<0) = 0;     % remove negitive weights.
        CIJ(isnan(CIJ)) = 0;
        
        %% Grab basic network stats
        
        symHuh = issymmetric(CIJ);
        uniqueHuh = unique(CIJ);
        [density,numNodes,numEdges] = density_dir(CIJ);
        
        % compute the shortest path matrix
        CIJ_recip = 1./CIJ;
        D_wei = distance_wei(CIJ_recip);
        
        % calculate path length and efficiency
        % no diag or inf path = 0,0
        [pth,eff] = charpath(D_wei,0,0);
        
        % Degree Calculations
        [id,od,deg] = degrees_dir(CIJ);
        
        % Generate report
        disp(strjoin(["Is symmetric:                   ", num2str(symHuh)]));
        disp(strjoin(["Is unique:                      ", num2str(length(uniqueHuh))]));
        disp(strjoin(["Number of Nodes:                ", num2str(numNodes)]));
        disp(strjoin(["Number of Edges:                ", num2str(numEdges)]));
        disp(strjoin(["Density:                        ", num2str(density)]));
        disp(strjoin(["Average Path Length (weighted): ", num2str(pth)]));
        disp(strjoin(["Efficiency:                     ", num2str(eff)]));
        disp(strjoin(["Average Degree:                 ", num2str(mean(deg))]));
        
        
        
        %% community detection
        numreps = 20;
        % find communities
        ci = zeros(length(ratMats.clusters),numreps);
        q = zeros(1,numreps);
        for irep = 1:numreps
            [ci(:,irep),q(irep)] = community_louvain(CIJ);
        end
        
        ci = fcn_relabel_partitions(ci); % relabel communities
        
        disp(strjoin(["mean num communities:           ", num2str(mean(max(ci)))]));
        disp(strjoin(["mean q:                         ", num2str(mean(q))]));
        
        % histogram of community labels
        %h = hist(ci,1:max(ci))
        
        % calculate a co-assignment matrix
        ag = agreement(ci);
        
        % consensus clustering
        ciconsensus = fcn_consensus_communities(ci,10);
        
        if 1 % run plot blocks again
            [gx,gy,idx] = fcn_plot_blocks(ciconsensus);
            
            % visualize matrix
            h = figure; imagesc(ag(idx,idx)); hold on;
            plot(gx,gy,'w'); hold off;
            title([ratName + " communities: " + fields{i}]);
            figName = "Figs/comm/" + ratName + "_" + fields{i} + "_comm";
            savefig(h, char(figName+".fig"));
            saveas(h,  char(figName+".png"));
        end
        
        %% Pack it all up
        netStats.name = ratName;
        netStats.(fields{i}).density   = density;
        netStats.(fields{i}).numNodes  = numNodes;
        netStats.(fields{i}).numEdges  = numEdges;
        netStats.(fields{i}).pth       = pth;
        netStats.(fields{i}).eff       = eff;
        netStats.(fields{i}).meanDeg   = mean(deg);
        netStats.(fields{i}).numComm   = max(ci);
        netStats.(fields{i}).ci        = ci;
        netStats.(fields{i}).ciconsens = ciconsensus;
        netStats.(fields{i}).q         = q;
        netStats.(fields{i}).CIJ       = CIJ;
        %netStats.(fields{i}).          =    ;
        
        disp(' ')
        degs(rat,i) = mean(deg);
    end
    
    saveName = "Mats/netStats_" + ratName + ".mat";
    save(saveName,"netStats");
end
t=toc; % Takes around 15s
disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))




%% Visuals testing
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



if 0
% Average degree increases
tst = [1,3,5; 2,4,8; 3,4,5; 2,2,3];

figure; hold on; 
for zz = 1:size(degs,1)
  %scatter(1:3,tst(zz,:),50,'filled');
  plot(degs(zz,:)');

end

end


















