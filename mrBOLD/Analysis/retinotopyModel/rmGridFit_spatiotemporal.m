function model = rmGridFit_spatiotemporal(model,prediction,data,params,t)
% rmGridFit_oneGaussianNonlinear - core of one non-linear (exponential) Gaussian fit
% CSS or compressive spatial summation model
%
% model = rmGridFit_oneGaussian(model,prediction,data,params);
%
% 2014/02 JW: duplicated from rmGridFit_oneGaussian, then exponent added to
% model

% input check 
if nargin < 4,
    error('Not enough arguments');
end

nchan = size(prediction,3);

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



%%
% % % normalize channels
% size(tmodel.run_preds,2)
if strcmp(params.analysis.pRFmodel{1}, 'st')
if size(prediction,3) == 2
    for ii = 1:size(prediction,2) % 720       82608           2
        maxS = max(max(prediction(:,ii,1)));
        maxT = max(max(prediction(:,ii,2)));
        normTs(ii) = maxS / maxT;
        prediction(:,ii,2) = prediction(:,ii,2) * normTs(ii);
    end
end
end




%%

warning('off', 'MATLAB:lscov:RankDefDesignMat')
for n=1:numel(params.analysis.x0),
    %-----------------------------------
    % progress monitor (10 dots) and time indicator
    %-----------------------------------
    if floor(n./numel(params.analysis.x0).*10)>progress,
        if progress==0,
            % print out estimated time left
            esttime = toc.*10;
            if floor(esttime./3600)>0,
                fprintf(1,'[%s]:Estimated processing time: %d hours.\t(%s)\n',...
                    mfilename, ceil(esttime./3600), datestr(now));
            else
                fprintf(1, '[%s]:Estimated processing time: %d minutes.\t(%s)\n',...
                    mfilename, ceil(esttime./60), datestr(now));
            end;
            fprintf(1,'[%s]:Grid (x,y,sigma) fit:',mfilename);drawnow;
        end;
        % progress monitor
        fprintf(1,'.');drawnow;
        progress = progress + 1;
    end;

    %-----------------------------------
    %--- now apply glm to fit RF
    %-----------------------------------
    % minimum RSS fit
%     X    = [prediction(:,n,1) prediction(:,n,2) trends];
    X    = [squeeze(prediction(:,n,:)) trends];
    
    Maxnorm = max(squeeze(prediction(:,n,:)));

    % This line takes up 30% of the time
    % lscov takes as long as the pinv method but provides the rss as well...
    [b,~,rss]    = lscov(X,data); 
    
    % Compute RSS only for positive fits. The basic problem is
    % that if you have two complementary locations, you
    % could fit with a postive beta on the one that drives the signal or a
    % negative beta on the portion of the visual field that never sees the
    % stimulus. This would produce the same prediction. We don't like that
    %nkeep   = b(1,:)<0; % Now we only set the negative fits to inf.
    
    % To save time limit the rss computation to those we care about.
    % This line is takes up 60% of the time.... (replaced by lscov)
    if nchan == 2
        nkeep = ~all([b(1,:); b(2,:)] > 0);
    elseif nchan == 1
        nkeep = b(1,:) < 0;
    end
    rss(nkeep) = inf('single');
    
    %-----------------------------------
    %--- store data with lower rss
    %-----------------------------------
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
    
end;
failedidx = startfind==0;
startfind(failedidx) = 1; % just give index of one, b/c betas are zero anyways
pred_X = prediction(:,startfind,:);
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


