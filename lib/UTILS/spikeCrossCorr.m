load toytrains.mat

%t1 = faketrain
%t2 = faketrain


t1 = spikesA;
t2 = spikesC;


maxwin = 1; %set the maximum window to 1 sec.

c=1;
for i=1:length(t1)
    for j=1:length(t2)
     
        if abs(t2(j) - t1(i)) < maxwin
        raw(c) = t2(j) - t1(i);
        c = c+1;
    
        end
        
    end    
    
end


[n,xout] = hist(raw,1000);

bar(xout,n)

spkd(t1,spikesC,1)
spkd(t1,spikesB,1)
spkd(spikesB,spikesB,1)

