% for playing with buszaki data

%load('NoveltySessInfoMatFiles/Achilles_10252013_sessInfo.mat')
load('NoveltySessInfoMatFiles/Achilles_11012013_sessInfo.mat')

% unpack the data
spikeTimes = sessInfo.Spikes.SpikeTimes;
spikeIDs = sessInfo.Spikes.SpikeIDs;

PRE = sessInfo.Epochs.PREEpoch;
MAZ = sessInfo.Epochs.MazeEpoch;
PST = sessInfo.Epochs.POSTEpoch;


% 
preSpikeTimes = spikeTimes(floor((PRE(1)))+1 : (floor(PRE(2))));
preSpikeIDs = spikeIDs(floor((PRE(1)))+1 : (floor(PRE(2))));
idTst = preSpikeIDs(1:10);
tmeTst = preSpikeTimes(1:10);


%%
neurons = unique(idTst);    % find each spike cluster
find(idTst == 1032)         % grab spikes associated with the cluster
tmeTst(find(idTst == 1032)) % Grab the times for each neuron


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
figure; imagesc(spikes)



%%
% sort each spike by time
clust = unique(preSpikeIDs)';
spikesT = zeros(length(clust),round(max(preSpikeTimes)*1000));
i = 1; % counter to keep track of rows
for clu = clust
  % Grab the spike times from the specified cluster
  spkTms = preSpikeTimes(find(preSpikeIDs == clu));
  
  % for each time, update the spikes matrix
  % where each col is one ms
  for tme = spkTms'
   spikesT(i,round(tme*1000)) = 1; % Need to test for rounding errors
  end 

  i = i+1;
end

figure; imagesc(-1 * spikesT); title('Spikes');
xlabel('ms'); ylabel('neuron'); colormap('gray');

% find the number of spikes for each neuron
numSpikes = sum(spikesT,2);


%% Find the correlation between spike trains

%This needs serious attention

% Matlab xcorr method for windowing (not working)
if 0
    for n = 1:size(spikesT,1)
      for nn = 1:size(spikesT,1)
        xcoorMat(n,nn) = xcorr(spikesT(n,:),spikesT(nn,:));
      end
    end
end

% uses the Victor & Purpura spike time distance with a cost of 1
% Very time consuming
vectLen = 1000;
for n = 1:size(spikesT,1)
    disp(n)
  for nn = 1:size(spikesT,1)
    vpstdMat(n,nn) = spkd(spikesT(n,1:vectLen),spikesT(nn,1:vectLen),1);
  end
end
figure; imagesc(vpstdMat);

rsd_tc = .001;
cosalpha = .5;
MVRD = vanRossumMN(spikesT(n,1:vectLen),spikesT(nn,1:vectLen),rsd_tc,cosalpha);

% run -o- the mill corr
coorMat = corr(spikesT(:,:)');

% let's just look at positive correlations
coorMat(coorMat<0) = 0;
figure; imagesc(coorMat);

% i know neuron 21 and 61 have a high correlation
figure; imagesc(spikesT([21 61],:));














function [spikeMat] = sortSpikes()


end



% BONEYARD
% A = A( ismember( A, B ) )


% [sharedVals,idxsIntoA] = intersect(A,B);
% [~,Bsorting] = sort(B);
% A(idxsIntoA(Bsorting))




