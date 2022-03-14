function model = rmGridFit_spatiotemporal(model,prediction,data,params,t)
% rmGridFit_spatiotemporal - core of one non-linear (exponential) Gaussian fit
% model = rmGridFit_oneGaussian(model,prediction,data,params);
%
% model

% input check 
if nargin < 4,
    error('Not enough arguments');
end
% cellfun(@(x) x(1,:,:), prediction)
% nchan = length(unique(params.analysis.combineNeuralChan));
nchan = getChanNumber(params);

trends         = t.trends;
t_id           = t.dcid+nchan;

% we compute mean rss but we need sum rss (old convention)
model.rss=single(model.rss./(size(prediction,1)-size(trends,2)+1));  
startfind = zeros(1,size(data,2));   
%-----------------------------------
%--- fit different receptive fields profiles
%--- another loop --- and a slow one too
%-----------------------------------
tic; progress = 0;
% 
% if params.useGPU
%     prediction = gpuArray(prediction);
%     data = gpuArray(data);
%     trends = gpuArray(trends);
% end

%%
% % % normalize channels
% if strcmp(params.analysis.pRFmodel{1}, 'st')
% if size(prediction,3) == 2
%     for ii = 1:size(prediction,2) % 720       82608           2
%         maxS = max(max(prediction(:,ii,1)));
%         maxT = max(max(prediction(:,ii,2)));
%         normTs(ii) = maxS / maxT;
%         prediction(:,ii,2) = prediction(:,ii,2) * normTs(ii);
%     end
% end
% end




%%

warning('off', 'MATLAB:lscov:RankDefDesignMat')
for n=1:numel(params.analysis.x0)
%     if mod(n,100) == 0
%         disp(n)
%         esttime = 1;
%     end
        %-----------------------------------
    % progress monitor (10 dots) and time indicator
    %-----------------------------------
    if floor(n./numel(params.analysis.x0).*10)>progress
        if progress==0
            % print out estimated time left
            esttime = toc.*10;
            if floor(esttime./3600)>0
                fprintf(1,'[%s]:Estimated processing time: %d hours.\t(%s)\n',...
                    mfilename, ceil(esttime./3600), datestr(now));
            else
                fprintf(1, '[%s]:Estimated processing time: %d minutes.\t(%s)\n',...
                    mfilename, ceil(esttime./60), datestr(now));
            end
            fprintf(1,'[%s]:Grid (x,y,sigma) fit:',mfilename);drawnow;
        end
        % progress monitor
        fprintf(1,'.');drawnow;
        progress = progress + 1;
    end

    %-----------------------------------
    %--- now apply glm to fit RF
    %-----------------------------------
    % minimum RSS fit
%     X    = [prediction(:,n,1) prediction(:,n,2) trends];
    if iscell(prediction)
        eachPrediction = cell2mat(cellfun(@(x) x(:,n,:), prediction, 'UniformOutput', false));
        X    = [squeeze(eachPrediction) trends];
    else
        X    = [squeeze(prediction(:,n,:)) trends];
    end

    %     X    = [squeeze(prediction(:,n,:)) trends];
    
%     Maxnorm = max(squeeze(prediction(:,n,:)));

    % This line takes up 30% of the time
    % lscov takes as long as the pinv method but provides the rss as well...
    [b,~,rss]    = lscov(gather(X),gather(data)); 
    
    
%   b  = (ones([size(X,2) size(data,2)]));
%    rss  = ones(1,size(data,2));
%     b   = X \ data;
%     pred= X * b; %ones(1,size(data,2))
%     rss = calccod(pred,data);
%     rs  = response - data;
%     rss = sum(rs,1) ;

% rssdata = sum(data.^2);

    % Compute RSS only for positive fits. The basic problem is
    % that if you have two complementary locations, you
    % could fit with a postive beta on the one that drives the signal or a
    % negative beta on the portion of the visual field that never sees the
    % stimulus. This would produce the same prediction. We don't like that
    %nkeep   = b(1,:)<0; % Now we only set the negative fits to inf.
    
    % To save time limit the rss computation to those we care about.
    % This line is takes up 60% of the time.... (replaced by lscov)
    if nchan == 2
%         nkeep1 = ~all([b(1,:); b(2,:)] > 0);
%         nkeep2 = ~all([b(1,:); b(2,:)] < 20);
%         nkeep =  (nkeep1 + nkeep2) > 1;
        c = [b(1,:); b(2,:)] ;
        nkeep = ~(all(c>0)  &  all(c<20));
    elseif nchan == 1
        c = b(1,:);
        nkeep = ~((c>0)  & (c<20));
    end
    rss(nkeep) = inf('single');
%     rss(nkeep) = 0;

    %-----------------------------------
    %--- store data with lower rss
    %-----------------------------------
%     minRssIndex = rss > model.rss;
%     findStop =  startfind > minRssIndex * n;
%     startfind(findStop) = n;

    minRssIndex = rss < model.rss;
    findStop =  startfind < minRssIndex * n;
    startfind(findStop) = n;

    % now update
    model.x0(minRssIndex)       = params.analysis.x0(n);
    model.y0(minRssIndex)       = params.analysis.y0(n);
    model.s(minRssIndex)        = params.analysis.sigmaMajor(n);
    model.s_major(minRssIndex)  = params.analysis.sigmaMajor(n);
    model.s_minor(minRssIndex)  = params.analysis.sigmaMajor(n);
    model.s_theta(minRssIndex)  = params.analysis.theta(n);
    model.exponent(minRssIndex) = params.analysis.exponent(n);
    model.rss(minRssIndex)      = rss(minRssIndex);
    model.b([1:nchan t_id],minRssIndex) = b(:,minRssIndex);

end

failedidx = startfind==0;
startfind(failedidx) = 1; % just give index of one, b/c betas are zero anyways

if iscell(prediction)
    pred_X = cell2mat(cellfun(@(x) x(:,startfind,:), prediction, 'UniformOutput', false));
else
    pred_X = prediction(:,startfind,:);
end

model.pred_X = pred_X;

% figure()
% subplot(3,1,1)
% % plot(pred_X(:,1,1)*model.b(1,1)); hold on
% % plot(pred_X(:,1,2)*model.b(2,1)); hold on
% plot(pred_X(:,12,1)*model.b(1,12)+pred_X(:,12,2)*model.b(2,12)) ; hold on
% plot(data(:,12),'k')
% ylim([-5 10])
% % % % 
% % subplot(3,1,2)
% % plot(pred_X(:,2,1)*model.b(1,2)); hold on
% % plot(pred_X(:,2,2)*model.b(2,2)); hold on
% plot(pred_X(:,2,1)*model.b(1,2)+pred_X(:,2,2)*model.b(2,2)) ; hold on
% 
% % % plot(data(:,2),'k')
% % % ylim([-5 10])
% % % 
% % % subplot(3,1,3)
% % % % plot(pred_X(:,3,1)*model.b(1,3)); hold on
% % % % plot(pred_X(:,3,2)*model.b(2,3)); hold on
% % % plot(pred_X(:,3,1)*model.b(1,3)+pred_X(:,3,2)*model.b(2,3)) ; hold on
% % % plot(data(:,3),'k')
% % % ylim([-5 10])



%warning('on', 'MATLAB:lscov:RankDefDesignMat')
% model.n = minRssIndex;
% Under some conditions, the grid fit never returns an acceptable fit, For
% example for onegaussian fits with data driven DC component, when the DC
% is artificially high. In this case some of the rss values remain Inf,
% which fails to interpolate and compute correct variance explained values.
% So we check it here and reset any Inf (bad fits) to rawrss, so the
% variance explained will be 0.
model.rss(model.rss==Inf)=model.rawrss(model.rss==Inf);

% Correct lscov. It returns the mean rss. To maintain compatibility with the
% sum rss this function expects, we have to multiply by the divisor. See
% the lscov docs for details.
model.rss=single(model.rss.*(size(prediction,1)-(size(trends,2)+1)));  

% X    = [prediction(:,toto,1) prediction(:,toto,2) trends];
% 
% model.X   = X;

% end time monitor
et  = toc;
if floor(esttime/3600)>0,
    fprintf(1,'Done[%d hours].\t(%s)\n', ceil(et/3600), datestr(now));
else
    fprintf(1,'Done[%d minutes].\t(%s)\n', ceil(et/60), datestr(now));
end;
drawnow;
return;


