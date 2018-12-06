
load("Mats/ratMats_Achilles_11012013.mat")


CIJ = squeeze(ratMats.pre);
figure; imagesc(CIJ); title("Origional");
% remove negitive weights. 
CIJ(CIJ<0) = 0;
figure; imagesc(CIJ); title("No Neg");

if 0, imagesc(squeeze(ratMats.pre)), end

%% Grab basic network stats

symHuh = issymmetric(CIJ);
uniqueHuh = unique(CIJ);
[density,numNodes,numEdges] = density_dir(CIJ);

% compute the shortest path matrix
CIJ_recip = 1./CIJ;            
D_wei = distance_wei(CIJ_recip); 

% calculate path length and efficiency
[pth,eff] = charpath(D_wei);    

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
numreps = 5;

ci = zeros(length(ratMats.clusters),numreps);
q = zeros(1,numreps);
for irep = 1:numreps
    [ci(:,irep),q(irep)] = community_louvain(CIJ);
end

ci = fcn_relabel_partitions(ci); % relabel communities

% total number of comminities
%max(ci)

% generate histogram of community labels
%h = hist(ci,1:max(ci))

% What is the q value?

% calculate agreement matrix
d = agreement(ci);

% consensus clustering 
numconsensus = 10;
ciconsensus = fcn_consensus_communities(ci,numconsensus);

if 1 % run plot blocks again
[gx,gy,idx] = fcn_plot_blocks(ciconsensus);

% visualize matrix
imagesc(d(idx,idx)); hold on;
figure; plot(gx,gy,'w'); hold off;
end










  









