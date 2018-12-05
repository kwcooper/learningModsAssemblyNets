% for playing with buszaki data

%load('NoveltySessInfoMatFiles/Achilles_10252013_sessInfo.mat')
load('NoveltySessInfoMatFiles/Achilles_11012013_sessInfo.mat')

% unpack the data
spikeTimes = sessInfo.Spikes.SpikeTimes;
spikeIDs = sessInfo.Spikes.SpikeIDs;

% split into phases
PRE = sessInfo.Epochs.PREEpoch;
MAZ = sessInfo.Epochs.MazeEpoch;
PST = sessInfo.Epochs.POSTEpoch;

% 
preSpikeTimes = spikeTimes(floor((PRE(1)))+1 : (floor(PRE(2))));
preSpikeIDs = spikeIDs(floor((PRE(1)))+1 : (floor(PRE(2))));
idTst = preSpikeIDs(1:10);
tmeTst = preSpikeTimes(1:10);


%% Unpack and sort the spikes by neuron
if 0 
% sort each spike time by cells (rows)
clust = unique(preSpikeIDs)';
spikes = zeros(length(clust),length(preSpikeTimes));
i = 1; % counter to keep track of rows
for clu = clust
    
  spkTms = preSpikeTimes(find(preSpikeIDs == clu));
  tmsInd = arrayfun(@(x)find(preSpikeTimes==x,1),spkTms); % grab indx's of times
  spikes(i,tmsInd) = spkTms';
  
  i = i + 1;
end
% figure; imagesc(spikes) % look at the initial raster
end

%% sort the spikes accounting for time 
clust = unique(preSpikeIDs)';
spikesT = zeros(length(clust),round((preSpikeTimes(end)*1000)));
i = 1; % counter to keep track of rows
for clu = clust
  % Grab the spike times from the specified cluster
  spkTms = preSpikeTimes(find(preSpikeIDs == clu));
  
  % for each time, update the spikes matrix
  % where each col is one us
  for tme = spkTms'
   spikesT(i,floor(tme*10000)) = 1; % Need to test for rounding errors
  end 

  i = i+1;
end

if 0 % View the raster plot
    figure; imagesc(-1 * spikesT); title('Spikes');
    xlabel('us'); ylabel('neuron'); colormap('gray');
end

numSpikes = sum(spikesT,2); % find the number of spikes for each neuron


%% Find the similarity between spike trains

% Victor & Purpura spike time distance with a cost of 1
if 0 % Very time consuming (Discontinued)
    vectLen = 1000;
    for n = 1:size(spikesT,1)
        disp(n)
      for nn = 1:size(spikesT,1)
        vpstdMat(n,nn) = spkd(tmeTst(find(idTst == clust(n))),tmeTst(find(idTst == clust(nn))),1);
      end
    end
    figure; imagesc(vpstdMat);
end

if 0 % vanRossum method (Discontinued)
    rsd_tc = .001;
    cosalpha = .5;
    MVRD = vanRossumMN(spikesT(n,1:vectLen),spikesT(nn,1:vectLen),rsd_tc,cosalpha);
end

% Convolution method
% Using a kernel, K, expand the width of each spike. 
% todo: define gaussian kernel...
kSize = 10000;
k = ones(1,kSize)/kSize;
spikesConv = zeros(size(spikesT,1),size(spikesT,2)+kSize-1);
for i = 1:size(spikesT,1)
    spikesConv(i,:) = conv(k,spikesT(i,:));
end
if 1, plot(spikesConv(1,:)), end

% Build correlation matrix
corrMat = corr(spikesConv');
imagesc(corrMat); title([num2str(kSize), ' kernel']);


% Test for variability
spikesRand = zeros(size(spikesT));
for i = 1:104
idx = sum(spikesT(i,:));
r = randperm(size(spikesT,2),idx);
spikesRand(i,r) = 1;
end

cRand = corr(spikesConvRand');
figure; imagesc(cRand)
plot(spikesRand(i,:))






