
% Flags for plotting
viewCorr = 1;
viewNull = 1;
runNull  = 1;
viewRast = 0;

% Choose dataset
tic
for i = 1:4
    %load('NoveltySessInfoMatFiles/Achilles_10252013_sessInfo.mat')
    disp('Loading Data...');
    if i == 1
        load("NoveltySessInfoMatFiles/Achilles_11012013_sessInfo.mat")
        ratName = "Achilles_11012013";
    end
    
    if i == 2
        load("NoveltySessInfoMatFiles/Buddy_06272013_sessInfo.mat")
        ratName = "Buddy_06272013";
    end
    
    if i == 3
        load("NoveltySessInfoMatFiles/Cicero_09012014_sessInfo.mat")
        ratName = "Cicero_09012014";
    end
    
    if i == 4
        load("NoveltySessInfoMatFiles/Gatsby_08022013_sessInfo.mat")
        ratName = "Gatsby_08022013";
    end
    
    disp(strjoin(["running", ratName]))
    
    % unpack the data and split into phases
    spikeTimes = sessInfo.Spikes.SpikeTimes;
    spikeIDs = sessInfo.Spikes.SpikeIDs;
    clust = unique(spikeIDs)';
    
    PRE = sessInfo.Epochs.PREEpoch;
    MAZ = sessInfo.Epochs.MazeEpoch;
    PST = sessInfo.Epochs.POSTEpoch;
    phases = [PRE, MAZ, PST];
    
    %%
    labels = {"PRE"; "MAZ"; "PST"};
    ratsBuffer = zeros(3, size(clust,2), size(clust,2));
    ratsNullBuffer = zeros(3, size(clust,2), size(clust,2));
    count = 1;
    for i = [1,3,5]
        disp(strjoin([labels(count), "phase..."]));
        spikeTimesSeg = spikeTimes(floor((phases(i)))+1 : (floor(phases(i+1))));
        spikeIDsSeg = spikeIDs(floor((phases(i)))+1 : (floor(phases(i+1))));
        % Need to update this ^
        
        
        
        %% Unpack and sort the spikes by neuron
        if 0
            % sort each spike time by cells (rows)
            spikes = zeros(length(clust),length(spikeTimesSeg));
            i = 1; % counter to keep track of rows
            for clu = clust
                
                spkTms = spikeTimesSeg(find(spikeIDsSeg == clu));
                tmsInd = arrayfun(@(x)find(spikeTimesSeg==x,1),spkTms); % grab indx's of times
                spikes(i,tmsInd) = spkTms';
                
                i = i + 1;
            end
            % figure; imagesc(spikes) % look at the initial raster
        end
        
        
        %% sort the spikes accounting for time
        disp("Sorting Spikes...");
        spikesT = zeros(length(clust),round((spikeTimesSeg(end)*1000)));
        i = 1; % counter to keep track of rows
        for clu = clust
            % Grab the spike times from the specified cluster
            spkTms = spikeTimesSeg(find(spikeIDsSeg == clu));
            
            % for each time, update the spikes matrix
            % where each col is one us
            for tme = spkTms'
                spikesT(i,floor(tme*10000)) = 1; % Need to test for rounding errors
            end
            i = i+1;
        end
        
        if viewRast % View the raster plot
            figure; imagesc(-1 * spikesT); title("Spikes");
            xlabel("us"); ylabel("neuron"); colormap("gray");
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
        % use step or gaussian kernel, needs fuzzy logic toolbox
        % step method uses neither. Almost identical results
        % https://www.researchgate.net/post/Kernel_functions_to_convolve_spike_events
        disp("Running Convolution...");
        kSize = 10000;
        %k = ones(1,kSize)/kSize; % for step method uncomment
        k = gaussmf(linspace(0,100,kSize),[kSize/4 kSize/2]);
        spikesConv = zeros(size(spikesT,1),size(spikesT,2)+kSize-1);
        
        for i = 1:size(spikesT,1)
            spikesConv(i,:) = conv(k,spikesT(i,:));
        end
        if 0, figure; plot(spikesConv(1,:)), end
        
        % Build correlation matrix
        disp("Building Correlation Matrix...");
        corrMat = corr(spikesConv');
        if viewCorr, figure; imagesc(corrMat); end
        title([labels(count), " kernel size = ", num2str(kSize),]);
        
        ratsBuffer(count,:,:) = corrMat;
        
        %% run null model
        if runNull
            
            disp("Running Null...")
            % computed by randomly permuting the same number
            % of spikes and reconvolving
            % Gaussian Kernel is necessary
            % Almost 0 correlation...
            spikesRand = zeros(size(spikesT));
            for i = 1:size(spikesRand,1)
                idx = sum(spikesT(i,:));
                r = randperm(size(spikesT,2),idx);
                spikesRand(i,r) = 1;
            end
            
            kSize = 10000;
            k = gaussmf(linspace(0,100,kSize),[kSize/4 kSize/2]);
            spikesConv = zeros(size(spikesRand,1),size(spikesRand,2)+kSize-1);
            for i = 1:size(spikesT,1)
                spikesConv(i,:) = conv(k,spikesRand(i,:));
            end
            
            cRand = corr(spikesRand');
            %figure; imagesc(cRand^8) % couldn't see much until the 8th...
            if viewNull, figure; imagesc(cRand); end
            title(["Null: ", labels(count), " kernel size = ", num2str(kSize),]);
        end
        
        ratsNullBuffer(count,:,:) = cRand;
        
        count = count + 1;
        disp(' ')
    end
    
    % unpack, store, and save
    ratMats = struct;
    ratMats.pre        = squeeze(ratsBuffer(1,:,:));
    ratMats.maz        = squeeze(ratsBuffer(2,:,:));
    ratMats.pst        = squeeze(ratsBuffer(3,:,:));
    ratMats.null.pre   = squeeze(ratsNullBuffer(1,:,:));
    ratMats.null.maz   = squeeze(ratsNullBuffer(2,:,:));
    ratMats.null.pst   = squeeze(ratsNullBuffer(3,:,:));
    ratMats.numSpikes  = numSpikes;
    ratMats.clusters   = clust;
    ratMats.spikeTimes = spikeTimes;
    ratMats.spikeIDs   = spikeIDs;
    ratMats.phases     = phases;
    ratMats.name       = ratName;
    
    saveName = "Mats/ratMats_" + ratName + ".mat";
    save(saveName,"ratMats");
    toc
end
t=toc; % Takes around 5m
disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))

disp('fin')






