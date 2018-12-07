% Network analysis on the spiking connectivity

% Define inter-rat struct
btRats.degs = zeros(4,3);
btRats.deg = {};
btRats.meanq = zeros(4,3);
plt = 0;

tic
% for each rat, run analysis
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
    
    disp(" ")
    disp(strjoin(["RAT:", ratName]))
    
    
    % For each phase, run the network analysis
    netStats = struct;
    fields = {'pre','maz','pst'};
    for i = 1:3% size(fields,2)
        disp(strjoin(["Phase:", fields{i}]))
        %keyboard; 
        
        % Grab matrix
        CIJ = squeeze(ratMats.(fields{i}));
        if plt
        h = figure; imagesc(CIJ); 
        title([ratName + " CIJ: " + fields{i}]); 
        figName = "Figs/CIJ/" + ratName + "_" + fields{i} + "_CIJ";
        savefig(h, char(figName+".fig"));
        saveas(h,  char(figName+".png"));
        end
        
        % remove negitive weights and NaN's.
        CIJ(CIJ<0) = 0;     
        CIJ(isnan(CIJ)) = 0;
        
        %% Compute basic network stats
        
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
        
        
        %% community detection
        numreps = 20;
        % find communities
        ci = zeros(length(ratMats.clusters),numreps);
        q = zeros(1,numreps);
        for irep = 1:numreps
            [ci(:,irep),q(irep)] = community_louvain(CIJ);
        end
        
        ci = fcn_relabel_partitions(ci); % relabel communities
        
        % histogram of community labels
        %h = hist(ci,1:max(ci))
        
        % calculate a co-assignment matrix
        ag = agreement(ci);
        
        % consensus clustering
        ciconsensus = fcn_consensus_communities(ci,10);
        [gx,gy,idx] = fcn_plot_blocks(ciconsensus);
        CIJ_comm = ag(idx,idx); 
        
        if plt % run plot blocks again
            
            
            % visualize matrix
            h = figure; imagesc(ag(idx,idx)); hold on;
            plot(gx,gy,'w'); hold off;
            title([ratName + " communities: " + fields{i}]);
            
            figName = "Figs/comm/" + ratName + "_" + fields{i} + "_comm";
            savefig(h, char(figName+".fig"));
            saveas(h,  char(figName+".png"));
        end
        
        %% Generate report
        disp(strjoin(["Is symmetric:                   ", num2str(symHuh)]));
        disp(strjoin(["Is unique:                      ", num2str(length(uniqueHuh))]));
        disp(strjoin(["Number of Nodes:                ", num2str(numNodes)]));
        disp(strjoin(["Number of Edges:                ", num2str(numEdges)]));
        disp(strjoin(["Density:                        ", num2str(density)]));
        disp(strjoin(["Avg Mat Sim:                    ", num2str(mean2(CIJ))]));
        disp(strjoin(["Average Path Length (weighted): ", num2str(pth)]));
        disp(strjoin(["Efficiency:                     ", num2str(eff)]));
        disp(strjoin(["Average Degree:                 ", num2str(mean(deg))]));
        disp(strjoin(["mean num communities:           ", num2str(mean(max(ci)))]));
        disp(strjoin(["mean q:                         ", num2str(mean(q))]));
        
        %% Pack it all up
        
        % Within rat
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
        
        % Between rat
        btRats.degs(rat,i)             = mean(deg);
        btRats.mSim(rat,i)             = mean2(CIJ);
        btRats.meanq(rat,i)            = mean(q);
        btRats.meanci(rat,i)           = mean(max(ci));
        btRats.phase(i).deg{rat}       = sort(deg, 'descend');
        btRats.phase(i).CIJ{rat}       = CIJ;
        btRats.phase(i).CIJCM{rat}     = CIJ_comm;
        btRats.phase(i).ci{rat}        = mean(max(ci));
        btRats.phase(i).q{rat}         = mean(q);
        disp(' ')
    end
    
    % save it for later
    saveName = "Mats/netStats_" + ratName + ".mat";
    save(saveName,"netStats");
end
t=toc; % Takes around 44s
disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
    
% save it for later
save("Mats/btRats.mat","btRats");


%% Analyze between rats...

if 0 % plot connectivity matricies for one rat
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

if 0 % Average similarity increases?
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

if 0 % Average degree increases?
    h = figure; hold on;
    cVec = [[0 1 0]; [1 0 1]; [0 1 1]; [1 0 0]];
    for zz = 1:size(btRats.degs,1)
        hold on;
        scatter(1:3,btRats.degs(zz,:),70,cVec(zz,:),'filled');
        line(1:3,btRats.degs(zz,:),'Color',cVec(zz,:))
        %plot(btRats.degs(zz,:)', 'LineWidth',3);
        
        ylabel('Average Degree');
        set(gca,'XTick',[0 1 2 3 4])
        names = {'PRE';'MAZE';'POST'};
        set(gca,'xtick',[1:3],'xticklabel',names)
        xlim([.5 3.5]);
        ax = gca; ax.FontSize = 14; 
    end
    
    figName = "Figs/degDist/avgDeg";
    savefig(h, char(figName+".fig"));
    saveas(h,  char(figName+".png"));
end

if 0 % Average number of communities?
    h = figure; hold on;
    cVec = [[0 1 0]; [1 0 1]; [0 1 1]; [1 0 0]];
    for zz = 1:size(btRats.meanci,1)
        hold on;
        scatter(1:3,btRats.meanci(zz,:),70,cVec(zz,:),'filled');
        line(1:3,btRats.meanci(zz,:),'Color',cVec(zz,:))
        %plot(btRats.degs(zz,:)', 'LineWidth',3);
        
        ylabel('Average Number of Communities');
        set(gca,'XTick',[0 1 2 3 4])
        names = {'PRE';'MAZE';'POST'};
        set(gca,'xtick',[1:3],'xticklabel',names)
        xlim([.5 3.5]);
        ax = gca; ax.FontSize = 14; 
    end
    
    figName = "Figs/avgCi";
    savefig(h, char(figName+".fig"));
    saveas(h,  char(figName+".png"));
end

if 0 % Average q?
    h = figure; hold on;
    cVec = [[0 1 0]; [1 0 1]; [0 1 1]; [1 0 0]];
    for zz = 1:size(btRats.meanq,1)
        hold on;
        scatter(1:3,btRats.meanq(zz,:),70,cVec(zz,:),'filled');
        line(1:3,btRats.meanq(zz,:),'Color',cVec(zz,:))
        %plot(btRats.degs(zz,:)', 'LineWidth',3);
        
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




if 0 % plot degree distribution
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
% 1. e? 2. 

%%
if 0 
N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,N*N*T+2*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+2*omega*N*(T-1);
B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
[S,Q] = genlouvain(B);
Q = Q/twomu
S = reshape(S,N,T);

end
% small omega .01

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
















