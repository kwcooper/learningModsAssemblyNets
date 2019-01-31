function [count,stats] = fcn_spikecondprob(u,v,window,nrand)
% fcn_spikecondprob     calculate conditional probability for fc study
%
%   [count,stats] = fcn_spikecondprob(spikes,window) takes a cell x time
%   input matrix and a window (measured in time points and NOT seconds),
%   computes the conditional probability that given a spike in cell i at
%   time t, you'll observe at least one spike in cell j within a window to
%   t + window. This measurement is biased by spike rate, so you compare
%   against null model in which spike times are randomized but spike rate
%   is preserved.
%   
%   fcn_spikecond()
%   inputs:     u,          cell list
%               v,          spike time list        
%               window,     measured in number of time points
%               nrand,      number of randomizations
%
%   outputs:    count,      count matrix -- # of times you see spike in j
%                           given spike in i.
%               stats,      includes mean, standard deviation, z-score, and
%                           p-values for conditional spikes.
%
%   Rick Betzel, Indiana University, 2019
%   kwc 190128
%%
vu = unique(v);
t = length(vu);
n = max(u);
%%
count = zeros(n);
for i = 1:t
    x = u(v == vu(i));
    y = u(v > vu(i) & v <= (vu(i) + window));
    h = hist(y,1:n) > 0;
    count(x,:) = count(x,:) + h;
end
%%
c = zeros(n,n,nrand);
p = zeros(n,n,2,nrand);
for irand = 1:nrand
    vr = v(randperm(length(v)));
    countr = zeros(n);
    for i = 1:t
        x = u(vr == vu(i));
        y = u(vr > vu(i) & vr <= (vu(i) + window));
        h = hist(y,1:n) > 0;
        countr(x,:) = countr(x,:) + h;
    end
    c(:,:,irand) = countr;
    p(:,:,1,irand) = count > countr;
    p(:,:,2,irand) = count < countr;
end
%%
stats.mu = nanmean(c,3);
stats.sd = nanstd(c,[],3);
stats.z = (count - stats.mu)./stats.sd;
stats.nrand = stats;
stats.pval = 1 - nanmean(p,4);