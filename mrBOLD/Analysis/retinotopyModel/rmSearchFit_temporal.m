function model = rmSearchFit_temporal(model, data, params, wProcess, t)
% rmSearchFit_temporal - wrapper for 'fine' one Gaussian fit
%
% model = rmSearchFit_oneGaussianNonlinear(model, data, params, wProcess, t);
%
% 2008/01 SOD: split of from rmSearchFit.
% 2010/02 SOD: cleanup.
% 2015/02 JW: branched from rmSearchFit_oneGaussian, now includes
%             non-linear model (see Kay et al, 2013, on Compressive
%             Spatial Summation)
% 2020/04 KIS: added spatiotemporal

% fminsearch options
searchOptions = params.analysis.fmins.options;
expandRange   = params.analysis.fmins.expandRange;
 
% convert to double just in case
params.analysis.X = double(params.analysis.X);
params.analysis.Y = double(params.analysis.Y);


% % % [cst] parameters % % %% % %% % %% % %% % %% % %% % %% % %% % %
if isfile(params.analysis.predFile)
    disp('***st predfile exists --- loading...')
    params.analysis.predFile
    load(params.analysis.predFile);
else
    error("does not have Grid")
end

% if size(tmodel.chan_preds,1) == 3
%     tmodel.run_preds=cat(1,tmodel.run_preds,tmodel.run_preds);
% end

% % %% % %% % %% % %% % %% % %% % %% % %% % %% % %% % %% % %% % %% % %


%params.analysis.allstimimages_unconvolved = double(params.analysis.allstimimages_unconvolved);
data = double(data);

% get starting upper and lower range and reset TolFun 
% (raw rss computation (similar to norm) and TolFun adjustments)
model.s = model.s_major;
% [range, TolFun] = rmSearchFit_range(params,model,data);
% [range, TolFun] = rmSearchFit_range_temporal(params,model,data);
[range, TolFun] =  rmSearchFit_range_temporal_linear(params,model,data);
% [range TolFun] = rmSearchFit_range(params,model,data);
% [range2 TolFun2] = rmSearchFit_range(params,model,data);
% range.start(3,:) = range.lower(3,:);

% amount of negative fits
nNegFit  = 0;
vethresh = params.analysis.fmins.vethresh;
trends   = t.trends;
t_id     = t.dcid+tmodel.num_channels;


%-----------------------------------
% Go for each voxel
%-----------------------------------
progress = 0;tic;
for ii = 1:numel(wProcess),
    if progress==0,
        esttime = toc.*10;
        if floor(esttime./3600)>0
            fprintf(1,'[%s]:Estimated processing time: %d voxels: %d hours.\n',...
                mfilename,numel(wProcess),ceil(esttime./3600));
        else
            fprintf(1,'[%s]:Estimated processing time: %d voxels: %d minutes.\n',...
                mfilename,numel(wProcess),ceil(esttime./60));
        end
        fprintf(1,'[%s]:Nonlinear optimization (x,y,sigma, exponent):',mfilename);
        progress = progress + 1;
    end

    % progress monitor (10 dots)
    if floor(ii./numel(wProcess)*10)>progress,
        % print out estimated time left
        fprintf(1,'.');drawnow;
        progress = progress + 1;
    end

    % volume index
    vi = wProcess(ii);
    vData = data(:,ii);
    
    % reset tolFun: Precision of evaluation function. 
    % We define RMS improvement relative to the initial raw 'no-fit' data
    % RMS. So, 1 means stop if there is less than 1% improvement on the fit:
    % searchOptions = optimset(searchOptions,'tolFun',optimget(params.analysis.fmins.options,'tolFun')./100.*rawrss);
    % optimset and optimget are a little slow so:
    searchOptions.TolFun = TolFun(ii);

    % actual fitting routine
    if searchOptions.MaxIter>0
%         outParams = ...
%             fmincon(@(x) rmModelSearchFit_temporal(x,vData,...
%             params.analysis.X,...
%             params.analysis.Y,...
%             params.analysis.allstimimages_unconvolved,...
%             params.analysis.Hrf,... 
%             params.analysis.scans,...
%             trends),...
%             range.start(:,vi),[],[],[],[],range.lower(:,vi),range.upper(:,vi),...
%             @(x) distanceCon(x,range.start(:,vi),range.step(:,vi).*expandRange),searchOptions); 
%            
        outParams = ...
            fmincon(@(x) rmModelSearchFit_temporal(x,vData,...
            params.analysis.X,...
            params.analysis.Y,...
            tmodel,...
            tmodel.irfs.hrf,...
            params.analysis.scans,...
            trends),...
            range.start(:,vi),[],[],[],[],range.lower(:,vi),range.upper(:,vi),...
            @(x) distanceCon(x,range.start(:,vi),range.step(:,vi).*expandRange),searchOptions);
        
        %cstmodel.tmodel.irfs.hrf

    else
        outParams = range.start(:,vi);
    end

    %[ outParams range.lower(:,vi) range.start(:,vi) range.upper(:,vi)]

    % make RF, prediction and get rss,b
    Xv = params.analysis.X-outParams(1);
    Yv = params.analysis.Y-outParams(2);
    n  = outParams(4);
    rf = exp( (Yv.*Yv + Xv.*Xv) ./ (-2.*(outParams(3).^2)) );
    
    
    %%%%% [cst]  params
    stim = tmodel.chan_preds;
    nChan = tmodel.num_channels;
    hrf = tmodel.irfs.hrf;
    % chan == 1
    for cc = 1:nChan
        
        if size(tmodel.chan_preds,1) == 3 % for abc
            pred = cat(1, (stim{1,cc}*rf).^n, (stim{2,cc}*rf).^n ,(stim{3,cc}*rf).^n);
%             pred = cat(1,pred,pred);
        else
            pred = (stim{cc}*rf).^n;
        end
        
        pred = {double(pred)};
        
        pred_hrf = cellfun(@(X, Y) convolve_vecs(X, Y, tmodel.fs, 1 / tmodel.tr), ...
            pred, hrf, 'uni', false);
        pred_hrf = cellfun(@transpose,pred_hrf,'UniformOutput',false);
        pred_hrf=cell2mat(pred_hrf)';
        
        % normalization
        if cc ==2
            pred_hrf = pred_hrf*tmodel.normT;
        end
        
        prediction{cc} = pred_hrf;
        
    end
    if nChan == 2
        preds = cell2mat(prediction);
        normTs = max(preds(:,1))/max(preds(:,2));
        prediction{2} = preds(:,2) *normTs;
    end

    X  = [cell2mat(prediction) trends];
    b    = pinv(X)*vData;
    rss  = norm(vData-X*b).^2;

    % store results only if the first beta is positive, somehow fmincon
    % outputs negative fits. If the fit is negative keep old (grid) fit. We
    % do adjust the rss, so it won't be accidentally counted as a 'good'
    % fit. 
    if nChan == 1
        if b(1)>0
            model.x0(vi)         = outParams(1);
            model.y0(vi)         = outParams(2);
            model.s(vi)          = outParams(3);
            model.s_major(vi)    = outParams(3);
            model.s_minor(vi)    = outParams(3);
            model.s_theta(vi)    = 0;
            model.exponent(vi)   = outParams(4);
            model.rss(vi)        = rss;
            model.b([1 t_id],vi) = b;
            model.pred_X(:,vi,1) = prediction{1};

        else
            % change the percent variance explained to be just under the
            % current vethresh. So it counts as a 'coarse'-fit but can still be
            % included in later 'fine'-fits
            model.rss(vi)  = (1-max((vethresh-0.01),0)).*model.rawrss(vi);
            nNegFit = nNegFit + 1;
        end
    elseif nChan == 2
        if b(1)>0 && b(2)>0
            model.x0(vi)         = outParams(1);
            model.y0(vi)         = outParams(2);
            model.s(vi)          = outParams(3);
            model.s_major(vi)    = outParams(3);
            model.s_minor(vi)    = outParams(3);
            model.s_theta(vi)    = 0;
            model.exponent(vi)   = outParams(4);
            model.rss(vi)        = rss;
            model.b([1 2 t_id],vi) = b;
            model.pred_X(:,vi,1) = prediction{1};
            model.pred_X(:,vi,2) = prediction{2};

        else
            % change the percent variance explained to be just under the
            % current vethresh. So it counts as a 'coarse'-fit but can still be
            % included in later 'fine'-fits
            model.rss(vi)  = (1-max((vethresh-0.01),0)).*model.rawrss(vi);
            nNegFit = nNegFit + 1;
        end
    end
        
end

% end time monitor
et  = toc;
if floor(et/3600)>0,
    fprintf(1,'Done [%d hours].\n',ceil(et/3600));
else
    fprintf(1,'Done [%d minutes].\n',ceil(et/60));
end;
fprintf(1,'[%s]:Removed negative fits: %d (%.1f%%).\n',...
    mfilename,nNegFit,nNegFit./numel(wProcess).*100);
return;



%-----------------------------------
% make sure that the pRF can only be moved "step" away from original
% position "startParams" - for the one Gaussian model
function [C, Ceq]=distanceCon(x,startParams,step)
Ceq = [];
dist = x([1 2])-startParams([1 2]);
C = norm(dist) - step;
return;
%-----------------------------------

