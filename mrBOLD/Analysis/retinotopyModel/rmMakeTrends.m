function [t, nt, dcid] = rmMakeTrends(params,stimSet,verbose)
% rmMakeTrends - make trends to add to GLM fit
%
% [t, nt, dcid] = rmMakeTrends(params, [verbose]);
%
% params are retinotopy parameters (see rmDefineParameters).
%
% number of trends is determined by params.analysis.nDCT:
% cos([0:0.5:nDCT]).
%
% 2006/02 SOD: wrote it.
% 2007/08 SOD: vectorized and changed logic.
% 2020/12 ISK: futher developed to account fo cross validation runs.

if ~exist('verbose','var') || isempty(verbose), verbose = prefsVerboseCheck; end
if notDefined('stimSet')
    stimSet = 1:length(params.stim);
end
% preparation
tf = [params.stim(stimSet).nFrames]./[params.stim(stimSet).nUniqueRep];
ndct   = [params.stim(stimSet).nDCT].*2+1;
t      = zeros(sum(tf),max(sum(ndct),1));
start1 = [0 cumsum(tf)];
start2 = [0 cumsum(ndct)];

% make them separately for every scan
dcid= zeros(1,numel(stimSet));
for n = 1:numel(stimSet)
    % stimulus length 
    tc = linspace(0,2*pi,tf(n))';

    % trends for one scan
    t(start1(n)+1:start1(n+1),start2(n)+1:start2(n+1)) = cos(tc*(0:0.5:params.stim(n).nDCT));
    
    % keep track of dc components
    dcid(n) = start2(n)+1;
end

% number of trends
nt = size(t,2);
if verbose,
    fprintf('[%s]: Removing %d trends from combined data.\n', mfilename, nt);
end


return;
