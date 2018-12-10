% Network analysis on the spiking connectivity

% Define inter-rat struct
btRats.degs = zeros(4,3);
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
        btRats.degsSEM(rat,i)          = std(deg)/sqrt(length(deg));
        btRats.mSim(rat,i)             = mean2(CIJ);
        btRats.q(rat,i,:)              = q;
        btRats.meanq(rat,i)            = mean(q);
        btRats.meanqSEM(rat,i)         = std(q)/sqrt(length(q));
        btRats.meanci(rat,i)           = mean(max(ci));
        btRats.meanciSEM(rat,i)        = std(max(ci))/sqrt(length(max(ci)));
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
disp('fin')
t=toc; % Takes around 44s
disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
    
% save it for later
save("Mats/btRats.mat","btRats"); 










