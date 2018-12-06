%% community detection
% Detecting communities is an important part of network analysis. It helps
% with pattern detection, dimensions reduction, and identifies hubs and
% functionally related groups of nodes. In this demo, we'll use modularity
% maximization to detect communities in the coactivation matrix dataset
% from the BCT.

%%
% Let's start with some familiar lines of code.

% clears workspace, closes figures
clear all; close all;
 
% add BCT to MATLAB path
addpath(genpath('/Users/richardbetzel/Documents/MATLAB+PSY-P457/'));

% load coactivation matrix
load('Coactivation_matrix.mat');

%% a quick visualization
% Before we do any analysis, let's take a look at the connectivity matrix.
%%
figure;                         % make new figure
imagesc(Coactivation_matrix);   % visualize weighted connectivity matrix

%% detecting communities
% Let's use the community detection code |community_louvain| that's part of
% the BCT to detect communities in the coactivation matrix. This code
% generates two outputs (you can type "help community_louvain" in the 
% command line to learn more).

% run louvain algorithm
[ci,q] = community_louvain(Coactivation_matrix);

%%
% The first output, |ci| is a vector of node's community assignments. it
% has dimensions equal to the number of nodes in the network and each
% element's value is an integer that tells us the community to which the
% corresponding node was assigned. If I wanted to know the community
% assignment of node 10, for instance, I would write:

% community assignment of node 10
ci(10)

%%
% Alternatively, if I wanted the list of nodes assigned to community 2, for
% instance, I would write:

% get all nodes assigned to community 2
find(ci == 2)

%%
% If I wanted to know the total number of communities, I would
% write something like:

% biggest community label
max(ci)

%%
% If I wanted to know the number of nodes in each community, I would write:

% generate histogram of community labels
h = hist(ci,1:max(ci))

%%
% It's worth mentioning that the label a community is given is fairly
% arbitrary. What matters is the nodes assigned to that community. For
% example, if I reassigned all nodes in community 2 and gave them the new
% label of community 100000001, the composition of the community doesn't
% change. It's still the same set of nodes.

%%
% The second output, |q|, is the modularity value associated with this 
% particular partition (set of communities). It is a quality function, and
% larger values generally indicate better partitions of a network into
% communities.

% how big is q
q

%% visualizing communities
% I've included a function called |fcn_plot_blocks| for visualizing
% communities. Basically, it reorders nodes according to some community
% labels and helps visualize the "community blocks." It works like this:

% run plot blocks
[gx,gy,idx] = fcn_plot_blocks(ci);

%%
% The first two outputs are for plotting a grid around communities and the
% third output is new ordering of nodes. Let's look at our communities

% visualize matrix
imagesc(Coactivation_matrix(idx,idx))

% add a hold to current axis, because we're going to add extra components
% to our figure
hold on;

% plot the blocks
plot(gx,gy,'r');
hold off;

%%
% Nice. This helps us see the communities more clearly.

%% stochasticity of louvain algorithm
% The algorithm we used to detect communities is called "non-deterministic"
% or "stochastic." What this means is that it can return different
% estimates of communities if I run it twice in a row. Let's test this by
% running the algorithm 5 times.

% define number of repetitions
numreps = 5;

% number of nodes
numnodes = length(ci);

% preallocate array
ci = zeros(numnodes,numreps);
q = zeros(1,numreps);

% create a loop
for irep = 1:numreps
    [ci(:,irep),q(irep)] = community_louvain(Coactivation_matrix);
end

% an extra step that relabels communities
ci = fcn_relabel_partitions(ci);

%%
% This is the first time we've used a loop in this class. Basically, a loop
% is a way of repeating a particular command a fixed number of times. Here,
% we want to run the Louvain algorithm |numreps| times. We can talk more
% about loops or you can read more here: <https://www.mathworks.com/help/matlab/ref/for.html> 
%%
% Let's look at the outputs of the algorithm.

% show list of q values
q

%%
% These q values aren't identical, which means we're probably getting
% different partitions for each repetition. Let's look at those.

% visualize communities
imagesc(ci); colorbar; colormap(jet(max(ci(:)) + 1));

%%
% They are not identical.
%%
% A useful way of visualizing variability across communities is by
% calculate a co-assignment matrix (or agreement matrix).

% calculate agreement matrix
d = agreement(ci);

%%
% Each cell in the agreement matrix counts the number of times that the
% corresponding pair of nodes were assigned to the same community. Let's
% visualize it.

% plot consensus matrix
imagesc(d);

%% resolving variability
% We can resolve this variability by computing "consensus communities." To
% do this, we use the code |fcn+consensus_communities| (you'll need to
% download code from <http://netwiki.amath.unc.edu/GenLouvain/GenLouvain>
% before you can use this function). When you download it, you'll get a
% directory called "GenLouvain-master". Put that directory right next to
% where your BCT directory is stored in your "MATLAB+PSY-P457" folder. If
% you do this, it'll work just fine.

%%
% Basically, what this function does is it takes a group of partitions and
% tries to generate an "average" or "consensus" partition. It does this by
% suppressing communities that are inconsistent across the group of
% partitions, while amplifying communities that are consistent. As input,
% you need to give this function your list of partitions, |ci| in this
% case, and the number of times you'll run the consensus clustering
% algorithm |numconsensus|.

% number of consensus runs
numconsensus = 10;

% run consensus clustering code
ciconsensus = fcn_consensus_communities(ci,numconsensus);

%%
% Let's visualize the agreement matrix, but ordered by the new communities.

% run plot blocks again
[gx,gy,idx] = fcn_plot_blocks(ciconsensus);

% visualize matrix
imagesc(d(idx,idx));
hold on;

% plot the blocks
plot(gx,gy,'w');
hold off;