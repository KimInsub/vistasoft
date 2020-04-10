function view = rmGridFit(view,params)
% rmGridFit - fit retinotopic model along grid (coarse stage)
%
% view=rmGridFit(view,params);
%
% Brute force fitting of predictions based upon premade receptive fields
% (rmDefineParams, x, y, sigma)
%
% Output is saved in structure model, which should be accessed through
% rmSet and rmGet.
%
%
% 2005/12 SOD: wrote it.
% 2006/12 SOD: converted calculations to single precision. This
% speeds things up considerably and is kinder on memory use.
% 2007/03 SOD: incorporated true coarse to fine search and trimmed code
% considerably.
% 2008/01 SOD: split of actual fitting from rmGridFit. This allows a more
% general use of this code with several fitting procedure options.

if notDefined('view'),   error('Need view struct'); end;
if notDefined('params'), error('Need params'); end;


%-----------------------------------
%--- For speed we do our computations in single precision.
%--- But we output in double (for compatibility).
%-----------------------------------
params.analysis.x0               = single(params.analysis.x0);
params.analysis.y0               = single(params.analysis.y0);
params.analysis.sigmaMajor       = single(params.analysis.sigmaMajor);
params.analysis.sigmaMinor       = single(params.analysis.sigmaMinor);
params.analysis.theta            = single(params.analysis.theta);
params.analysis.exponent         = single(params.analysis.exponent);
params.analysis.X                = single(params.analysis.X);
params.analysis.Y                = single(params.analysis.Y);
params.analysis.allstimimages    = single(params.analysis.allstimimages);
params.analysis.sigmaRatio       = single(params.analysis.sigmaRatio);
params.analysis.sigmaRatioInfVal = single(params.analysis.sigmaRatioInfVal);
params.analysis.sigmaRatioMaxVal = single(params.analysis.sigmaRatioMaxVal);


% Accessing the trends needs to move inside the slice loop. This is because
% the size of trends can change inside the loop, causing errors upon the
% second iteration of the loop. In particular the size of trends can be
% changed due to temporal decimation, making the size incommensurate with
% the size of 'data', unless we remake the trends on each loop.
% %-----------------------------------
% %--- make trends to fit with the model (discrete cosine set)
% %-----------------------------------
% [trends, ntrends, dcid]  = rmMakeTrends(params);
% trends = single(trends);


% we can also specify other nuisance factors that should be removed, e.g.
% large fixation dot changes.
%if isfield(params.analysis,'allnuisance')
%    trends = [trends single(rmDecimate(params.analysis.allnuisance,params.analysis.coarseDecimate))];
%end

%-----------------------------------
%--- now loop over slices
%--- but initiate stuff first
%-----------------------------------
switch lower(params.wData),
    case {'fig','roi'},
        loopSlices = 1;
    otherwise,
        loopSlices = 1:params.analysis.nSlices;
end;
nSlices = length(loopSlices);


%-----------------------------------
%--- make all predictions first
%-----------------------------------
n = numel(params.analysis.x0);  % Is this the position?

% What is s?  Something about the grid?  Why is this always 1000?
s = [[1:ceil(n./1000):n-2] n+1]; %#ok<NBRAK>

% if we have a nonlinear model, then we cannot pre-convolve the stimulus
% with the hRF. instead we make predictions with the unconvolved images and
% then convolve with the hRF afterwards
if ~checkfields(params, 'analysis', 'nonlinear') || ~params.analysis.nonlinear
   
    %�for a lineaer model, use the pre-convolved stimulus images
    % the decimate function is not Matlab 2013b
    allstimimages = rmDecimate(params.analysis.allstimimages,...
        params.analysis.coarseDecimate);
    
    prediction = zeros(size(allstimimages,1),n,'single');
    fprintf(1,'[%s]:Making %d model samples:',mfilename,n);
    drawnow;tic;
    for n=1:numel(s)-1,
        % make rfs
        rf   = rfGaussian2d(params.analysis.X, params.analysis.Y,...
            params.analysis.sigmaMajor(s(n):s(n+1)-1), ...
            params.analysis.sigmaMinor(s(n):s(n+1)-1), ...
            params.analysis.theta(s(n):s(n+1)-1), ...
            params.analysis.x0(s(n):s(n+1)-1), ...
            params.analysis.y0(s(n):s(n+1)-1));
        % convolve with stimulus
        pred = allstimimages*rf;
        
        % store
        prediction(:,s(n):s(n+1)-1) = pred;
        if ismember(n, round((1:10)/10* numel(s)-1)), % every 10% draw a dot
            fprintf(1,'.');drawnow;
        end
    end;
    
    clear n s rf pred;
    fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now));
    drawnow;
    
% cst model    
elseif strcmp(params.analysis.pRFmodel{1}, 'cst')
    
    predictionFile = ['./cst_seq-' params.analysis.stimseq, ...
        '-tm-' params.analysis.temporaltype, ...
        '_prediction.mat'];
    if ~isfile(predictionFile)
        
        allstimimages = params.analysis.allstimimages_unconvolved;
%         save('allstimimages.mat','allstimimages','-v7.3')
        fprintf(1,'[%s]:Making %d model samples:',mfilename,n);
        drawnow;tic;
        
        %%%% convert back to ones %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% make it as cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        allstimimages(find(allstimimages)) = 1;
        cellimage = num2cell(allstimimages,1)';
        
        %%%% set temporal model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fit_exps = {'Exp2'};
        % make this so we can parse as an inputs
%         temp_type = '1ch-glm';
        temp_type = params.analysis.temporaltype;
        dohrf = 2;   
        
     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        
        mn = length(cellimage);  % Is this the position?
        ms = [[1:ceil(mn./1000):mn-2] mn+1]; %#ok<NBRAK>
        tic
        fprintf('Generating irf for %s model...\n', temp_type)
        for mn = 1:numel(ms)-1
            
            stimirf_base=[];
            stimirf_t_s=[];
            stimirf_t_t=[];
            
            tmodel = stModel(temp_type, fit_exps,'default');
            tmodel.stim = cellimage(ms(mn):ms(mn+1)-1);
            
       
            
            [tmodel.onsets, tmodel.offsets, dur] = cellfun(@cst_codestim, ...
                tmodel.stim, 'uni', false);
            
% % % % %             % gap implementation 
% % % % %             sfiles = tmodel.stimfiles; [nruns_max, nsess] = size(sfiles);
% % % % %             empty_cells = cellfun(@isempty, sfiles);
% % % % %             fs = tmodel.fs; gd = tmodel.gap_dur;
% % % % %             
% % % % %             % dummy stuff for now
% % % % %             nruns_max =size(tmodel.onsets,1);
% % % % %             cat_list = {'all'}; tt = repmat(cat_list,1,length(tmodel.offsets{1}));
% % % % %             
% % % % %             im=[];
% % % % %             for a= 1:size(tmodel.offsets,1)
% % % % %                 im{a,1} = tt;
% % % % %             end
% % % % %             for cc = 1:length(cat_list)   
% % % % %                 cat_idxs = cellfun(@(X) find(strcmp(cat_list(cc), X)), ...
% % % % %                     im, 'uni', false);
% % % % %                 on_idxs = cellfun(@(X, Y) round(fs * X(Y)), ...
% % % % %                     tmodel.onsets, cat_idxs, 'uni', false); on_idxs(empty_cells) = {1};
% % % % %                 off_idxs = cellfun(@(X, Y) round(fs * (X(Y) - gd / 2)), ...
% % % % %                     tmodel.offsets, cat_idxs, 'uni', false); off_idxs(empty_cells) = {1};
% % % % %                 % find frame indices of gap offset times
% % % % %                 goff_idxs = cellfun(@(X, Y) round(fs * (X(Y) + gd / 2)), ...
% % % % %                     tmodel.offsets, cat_idxs, 'uni', false); goff_idxs(empty_cells) = {1};
% % % % %                 % compile all frame indices during stimuli and gaps
% % % % %                 stim_idxs = cellfun(@code_stim_idxs, ...
% % % % %                     on_idxs, off_idxs, 'uni', false);
% % % % %                 gap_idxs = cellfun(@code_stim_idxs, ...
% % % % %                     off_idxs, goff_idxs, 'uni', false);
% % % % %                 % code stimulus as a step function with gaps at offsets
% % % % %                 cc_x = repmat({cc}, nruns_max, nsess);
% % % % %                 stims = cellfun(@code_stim_vec, tmodel.stim, stim_idxs, ...
% % % % %                     cc_x, repmat({1}, nruns_max, nsess), 'uni', false);
% % % % %                 stims = cellfun(@code_stim_vec, tmodel.stim, gap_idxs, ...
% % % % %                     cc_x, repmat({0}, nruns_max, nsess), 'uni', false);
% % % % %             end
% % % % %             stims(empty_cells) = {[]}; tmodel.stim = stims;
                        
            
            %%%% make IRFs
            tmodel = cst_pred_runs(tmodel,dohrf);
            
            
            
            %%%% transpose and sum up two channels 
            stimirf_base = cellfun(@transpose,tmodel.pixel_preds,'UniformOutput',false);
            
            if contains(temp_type,'2ch')
                
                stimirf_t_s  = cell2mat(cellfun(@(X) X(1,:), stimirf_base, 'uni', false))'; %  sustained
                stimirf_t_t  = cell2mat(cellfun(@(X) X(2,:), stimirf_base, 'uni', false))'; % transient
                stimirf_base = cell2mat(cellfun(@sum, stimirf_base, 'uni', false))';        % sum
                
                stimirf_chan_s(:,ms(mn):ms(mn+1)-1) = stimirf_t_s;
                stimirf_chan_t(:,ms(mn):ms(mn+1)-1) = stimirf_t_t;                
            else
                stimirf_base=cell2mat(stimirf_base)'; %        30000     X   1412
            end
            
            stimirf(:,ms(mn):ms(mn+1)-1) = stimirf_base;
            
            if ismember(mn, round((1:10)/10* numel(ms)-1)), % every 10% draw a dot
                fprintf(1,'(irf)');drawnow;
            end
            

        end
        
        if tmodel.num_channels == 2
            tmodel.chan_preds{1} = stimirf_chan_s;
            tmodel.chan_preds{2} = stimirf_chan_t;
        else
            tmodel.chan_preds{1} = stimirf_base;
        end
        
        
        % temporal channel normalization
        % need to think about ways of normalizing in the future.
        temporal_channel_normalization=false;
        if temporal_channel_normalization
             tmodel.normT = max(max(stimirf_chan_s)) /  max(max(stimirf_chan_t));
        else
            tmodel.normT = 1;
        end
        
        
%         rf   = rfGaussian2d(params.analysis.X, params.analysis.Y,...
%                 params.analysis.sigmaMajor(1), ...
%                 params.analysis.sigmaMinor(1), ...
%                 params.analysis.theta(1), ...
%                 params.analysis.x0(1), ...
%                 params.analysis.y0(1));
%             
%             for cc=1:tmodel.num_channels
%                 pred = tmodel.chan_preds{cc}*rf;
%                 
%                 % apply css
%                 pred = bsxfun(@power, pred, 0.5);
%                 pred = double(pred);
%                 
%                 % apply hrf
%                 pred_cell = {pred};
%                 
%                 % use mrvista HRF
%                 % hrf = params.analysis.Hrf;
%                 hrf = tmodel.irfs.hrf;
%                 curhrf = repmat(hrf, 1, 1);
%                 pred_hrf = cellfun(@(X, Y) convolve_vecs(X, Y, tmodel.fs, 1 / tmodel.tr), ...
%                     pred_cell, curhrf, 'uni', false);
%                 pred_hrf = cellfun(@transpose,pred_hrf,'UniformOutput',false);
%                 pred_hrf=cell2mat(pred_hrf)';
%                 
%                 
%                 % store
%                 prediction(:,cc) = pred_hrf;
%                 %             prediction{n} = pred_hrf;
%                 
%             end
%             maxvals = max(prediction)
%             normTs = maxvals(1) / maxvals(2)
%         %             prediction{n} = pred_hrf;
%         

% % %         % plot check
% % %         
% % % plot(tmodel.chan_preds{1})
% % % plot(tmodel.chan_preds{2})
% % % 
% % % xlim([2000 6000])
% % % ttmodel = pred_runs(tmodel);

        
        toc
        clear stimirf_t stimirf_chan_s stimirf_chan_t stimirf_t_s stimirf_t_t
   
                
        prediction = zeros(size(stimirf,1)/tmodel.fs,n,tmodel.num_channels,'single');
        % loop over grid
        tic
        for n=1:numel(s)-1
            
            fprintf('Generating predictors for %s model...%d/%d \n', temp_type,n,numel(s)-1)

            % make rfs
            rf   = rfGaussian2d(params.analysis.X, params.analysis.Y,...
                params.analysis.sigmaMajor(s(n):s(n+1)-1), ...
                params.analysis.sigmaMinor(s(n):s(n+1)-1), ...
                params.analysis.theta(s(n):s(n+1)-1), ...
                params.analysis.x0(s(n):s(n+1)-1), ...
                params.analysis.y0(s(n):s(n+1)-1));
            
            % convolve rf with stimulus for each t-channel
            for cc=1:tmodel.num_channels
                pred = tmodel.chan_preds{cc}*rf;
                
                % apply css
                pred = bsxfun(@power, pred, params.analysis.exponent(s(n):s(n+1)-1)');
                pred = double(pred);
                
                % apply hrf
                pred_cell = num2cell(pred,1)';
                npixel_max = size(s(n):s(n+1)-1,2);
                
                % use mrvista HRF
                %   params.stim(n).images = filter(params.analysis.Hrf{n}, 1, params.stim(n).images'); % images: pixels by time (so images': time x pixels)
%                 hrf = params.analysis.Hrf;
                hrf = tmodel.irfs.hrf;
                curhrf = repmat(hrf, npixel_max, 1);
                pred_hrf = cellfun(@(X, Y) convolve_vecs(X, Y, tmodel.fs, 1 / tmodel.tr), ...
                    pred_cell, curhrf, 'uni', false);
                pred_hrf = cellfun(@transpose,pred_hrf,'UniformOutput',false);
                pred_hrf=cell2mat(pred_hrf)';
                
                if cc ==2
                    pred_hrf = pred_hrf*tmodel.normT;
                end    
                
                % store
                prediction(:,s(n):s(n+1)-1,cc) = pred_hrf;
                %             prediction{n} = pred_hrf;
                
            end
            if ismember(n, round((1:10)/10* numel(s)-1)), % every 10% draw a dot
                fprintf(1,'(=.=)');drawnow;
            end
            
        end
        
        tmodel.run_preds = prediction;
        params.analysis.temporal.tmodel = tmodel;
        save(predictionFile, 'tmodel', '-v7.3');
        
        clear n s rf pred pred_hrf pred_cell;
        fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now));
        drawnow;
        
    else
        load(predictionFile)
    end
    
 

    
else
    allstimimages = params.analysis.allstimimages_unconvolved;
    prediction = zeros(size(allstimimages,1),n,'single');
    fprintf(1,'[%s]:Making %d model samples:',mfilename,n);
    drawnow;tic;
    for n=1:numel(s)-1,
        % make rfs
        rf   = rfGaussian2d(params.analysis.X, params.analysis.Y,...
            params.analysis.sigmaMajor(s(n):s(n+1)-1), ...
            params.analysis.sigmaMinor(s(n):s(n+1)-1), ...
            params.analysis.theta(s(n):s(n+1)-1), ...
            params.analysis.x0(s(n):s(n+1)-1), ...
            params.analysis.y0(s(n):s(n+1)-1));
        % convolve with stimulus
        pred = allstimimages*rf;
        
        % store
        prediction(:,s(n):s(n+1)-1) = pred;
        if ismember(n, round((1:10)/10* numel(s)-1)), % every 10% draw a dot
            fprintf(1,'.');drawnow;
        end
    end;
    
    scans = params.analysis.scan_number;
    
    % for nonlinear model, do the hRF convolution after the prediction
    if checkfields(params, 'analysis', 'nonlinear') && params.analysis.nonlinear
        prediction = bsxfun(@power, prediction, params.analysis.exponent');
        for scan = 1:numel(params.stim)
            inds = scans == scan;
            hrf = params.analysis.Hrf{scan};
            prediction(inds,:) = filter(hrf, 1, prediction(inds,:));
        end
    end
    
    % decimate predictions after convolution instead of before
    fprintf(1,'[%s]:Decimating data\n',mfilename);
    prediction = rmDecimate(prediction, params.analysis.coarseDecimate);

    clear n s rf pred;
    fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now));
    drawnow;
    
end

% go loop over slices
for slice=loopSlices,
    %-----------------------------------
    % Place datasets behind each other. This is a rather crude way of
    % stimultaneously fitting both. Due to this we cannot
    % prewhiten (we could zeropad/let the trends deal with this/not care).
    %-----------------------------------
    
    % [CST] for each channel level usage
    switch lower(params.analysis.pRFmodel{1})
        case {'cst'}
            files = dir(sprintf('synBOLD*seq*%s*%s.mat',...
                params.analysis.stimseq, params.analysis.temporaltype));
            load(files.name)
            for vv = 1:size(synBOLD{1,1},1)
               chandata(:,:,vv) =  cell2mat(cellfun(@(X) X(vv,:), synBOLD,'UniformOutput', false)')';
            end

            % CURRENTLY< CST ONLY WORKS with SYNTH PC
            params.analysis.calcPC=0;
    end
    
    [data, params] = rmLoadData(view, params, slice,...
        params.analysis.coarseToFine);
    

    
    
    
    % for speed convert to single and remove NaNs (should not be
    % there anyway!), 
    % We could remove NaN from data and put back later, so computations are
    % even faster. Actually, this would make the code a whole lot more
    % complicated. An easier way to do this is to create an ROI that only
    % includes non NaN data.
    
    % remove trends from data so they do not count in the percent variance
    % explained calculation later.
    data(isnan(data)) = 0;
    data = single(data);
    
    %-----------------------------------
    %--- make trends to fit with the model (discrete cosine set)
    %-----------------------------------
    [trends, ntrends, dcid] = rmMakeTrends(params);
    trends = single(trends);


    % The "trends" are cosine functions placed into the columns of a
    % matrix.  The betas for these are solved by the matrix multiplication
    % done here.  
    trendBetas = pinv(trends)*data;
    %if isfield(params.analysis,'allnuisance')
    %    trendBetas(ntrends+1:end) = 0;
    %    ntrends = ntrends + size(params.analysis.allnuisance,2);
    %end
    
    %%%%% don't do this for single pulse
    if size(data,1) > 30
        data = data - trends*trendBetas;
    end
    
    % reset DC component by specific data-period (if requested)
    if params.analysis.dc.datadriven
        [data, trendBetas] = rmEstimateDC(data,trendBetas,params,trends,dcid);
    end
    % decimate (if requested)

    fprintf('[%s]: Decimating data:...\n', mfilename)
    data   = rmDecimate(data,params.analysis.coarseDecimate);
    trends = rmDecimate(trends,params.analysis.coarseDecimate);
  
    % compute rss raw data for variance computation later
    rssdata = sum(data.^2);

    %-----------------------------------
    % initiate stuff on first loop
    %-----------------------------------
    if slice == 1,
        fprintf(1,'[%s]:Number of voxels: %d.\n',mfilename,size(data,2));drawnow;
        model = initiateModel(params, nSlices, size(data,2), ntrends);
        
        % These sets are necessary because the params struct is not saved.  Hence,
        % we need to pull out the critical values from params and put them
        % specifically in each model{mm}.
        
        % If there are roi data, get the coordinates.
        if strcmp(params.wData,'roi');
            for mm = 1:numel(model),
                model{mm} = rmSet(model{mm},'roiCoords',rmGet(params,'roiCoords'));
                model{mm} = rmSet(model{mm},'roiIndex',rmGet(params,'roiIndex'));
                model{mm} = rmSet(model{mm},'roiName',rmGet(params,'roiName'));
            end
        end
        
        % For all cases, put in number of data points. 
        for mm = 1:numel(model),
            model{mm} = rmSet(model{mm},'npoints',size(data,1));
        end
    end

    %-----------------------------------
    % Extract  the data from a slice and put it in a
    % temporary structure that will be modified throughout.
    %-----------------------------------
    s = rmSliceGet(model,slice);

    % initiateModel fills the rss-field with Infs. We reset them here
    % to a more data-driven maximum value of sum(data.^2)
    for n=1:numel(s),
        s{n}.rawrss       = rssdata;
    end;

    %-----------------------------------
    % At this point, we have a slice of data and we will fit different
    % pRF models.  This is the part of the code that is the slowest.
    %-----------------------------------
    if params.analysis.dc.datadriven
        t.trends = [];
        t.dcid   = [];
    else
        t.trends = trends(:,dcid);
        t.dcid   = dcid;       
    end

    % Run the grid fit estimate for this slice, s{}.
    for n=1:numel(params.analysis.pRFmodel)
        switch lower(params.analysis.pRFmodel{n}),
            case {'onegaussian','one gaussian','default'}
                s{n}=rmGridFit_oneGaussian(s{n},prediction,data,params,t);
                
            case {'1dgaussian','1d gaussian','1d'}
                s{n}=rmGridFit_oneGaussian(s{n},prediction,data,params,t);
                
            case {'onegaussiann','one gaussian with nuisance factors'}
                t.trends = trends;
                t.dcid   = 1:size(trends,2);
                s{n}=rmGridFit_oneGaussian(s{n},prediction,data,params,t);
                
            case {'onegaussianlinked','one gaussian linked to neighbors'}
                t.trends = trends;
                t.dcid   = 1:size(trends,2);
                s{n}=rmGridFit_oneGaussianLink(s{n},prediction,data,params,t,view);
            
            case {'addgaussian','add one gaussian'}
                [residuals s{n}] = rmComputeResiduals(view,params,s{n},slice,[true params.analysis.coarseDecimate>1]);
                t.dcid = t.dcid + 1;
                s{n}=rmGridFit_oneGaussian(s{n},prediction,residuals,params,t);
                trendBetas = zeros(size(trendBetas));
                
            case {'addgaussianlinked','add one gaussian linked to neighbors'}
                [residuals s{n}] = rmComputeResiduals(view,params,s{n},slice,...
                    [params.analysis.coarseToFine params.analysis.coarseDecimate>1]);
                t.dcid = t.dcid + 1;
                s{n}=rmGridFit_oneGaussianLink(s{n},prediction,residuals,params,t,view);
                trendBetas = zeros(size(trendBetas));
                
            case {'onegaussianunsigned','one gaussian unsigned','unsigned'}
                s{n}=rmGridFit_oneGaussianUnsigned(s{n},prediction,data,params,t);

            case {'twogaussianstog','tog','two gaussians on the same position'}
                s{n}=rmGridFit_twoGaussiansToG(s{n},prediction,data,params,t);

            case {'twogaussiansdog','dog','difference of gaussians'}
                s{n}=rmGridFit_twoGaussiansDoG(s{n},prediction,data,params,t);

            case {'two1dgaussiansdog','1ddog','1d dog','1d difference of gaussians'}
                s{n}=rmGridFit_twoGaussiansDoG(s{n},prediction,data,params,t);

            case {'twogaussiansdogfixed','dogf','difference of gaussians fixed'}
                % remake predictions using sigmaRatio and betaRatio
                ii = numel(params.analysis.x0);
                m = [[1:ceil(ii./1000):ii-2] ii+1]; %#ok<NBRAK>
                allstimimages = rmDecimate(params.analysis.allstimimages,...
                    params.analysis.coarseDecimate);
                % compute sigma 2
                sigma2.major = params.analysis.sigmaMajor .* params.analysis.sigmaRatioFixedValue(1);
                sigma2.minor = params.analysis.sigmaMajor .* params.analysis.sigmaRatioFixedValue(1);
                sigma2.theta = params.analysis.theta;
                
                fprintf(1,'[%s]:Making %d surround:',mfilename,n);
                drawnow;tic;
                for ii=1:numel(m)-1,
                    % make rfs
                    rf   = rfGaussian2d(params.analysis.X, params.analysis.Y,...
                        sigma2.major(m(ii):m(ii+1)-1), ...
                        sigma2.minor(m(ii):m(ii+1)-1), ...
                        sigma2.theta(m(ii):m(ii+1)-1), ...
                        params.analysis.x0(m(ii):m(ii+1)-1), ...
                        params.analysis.y0(m(ii):m(ii+1)-1));
                    % convolve with stimulus
                    pred = allstimimages*rf;
                    
                    % store
                    prediction(:,m(ii):m(ii+1)-1) = prediction(:,m(ii):m(ii+1)-1) + params.analysis.betaRatio(1).*pred;
                    fprintf(1,'.');drawnow;
                end;
                clear rf pred;
                fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now));
                drawnow;
                
                s{n}=rmGridFit_oneGaussian(s{n},prediction,data,params,t);

            case {'twogaussiansdogbetafixed','dogbf','difference of gaussians beta fixed'}
                s{n}=rmGridFit_twoGaussiansDoGbetafixed(s{n},prediction,data,params,t);
                
            case {'twogaussiansposonly','two gaussians','two prfs'}
                s{n}=rmGridFit_twoGaussiansPosOnly(s{n},prediction,data,params,t);
                
            case {'twogaussiansmirror','two gaussians mirrored','mirror'}
                s{n}=rmGridFit_twoGaussiansMirror(s{n},params.analysis.mirror,data,params,t);
                 
			case {'shiftedgaussians','two shifted gaussians'}
				s{n}=rmGridFit_shiftedGaussians(s{n},params.analysis.pRFshift,data,params,t);


            case {'oneovalgaussian','one oval gaussian','one oval gaussian without theta'}
                s{n}=rmGridFit_oneOvalGaussian(s{n},prediction,data,params,t);
                
            case {'css' 'onegaussiannonlinear', 'onegaussianexponent' }
                s{n}=rmGridFit_oneGaussianNonlinear(s{n},prediction,data,params,t);

            case {'cst'}
                prediction = tmodel.run_preds ;
%                 sus=rmGridFit_oneGaussianNonlinear(s{1},prediction(:,:,1),data,params,t);
%                 tran=rmGridFit_oneGaussianNonlinear(s{1},prediction(:,:,2),data,params,t);
%                 s{n}=rmGridFit_oneGaussianNonlinear(s{n},prediction(:,:,1),data,params,t);

fitFile = ['./cst_seq-' params.analysis.stimseq, ...
        '-tm-' params.analysis.temporaltype, ...
        '_fit.mat']
save(fitFile,'s','prediction','data','params','t','-v7.3')


% spatial temporal change the name
                s{n}=rmGridFit_spatiotemporal(s{n},prediction,data,params,t);
                
% figure()
% % a=prediction(:,:,1);
% % plot(a(:,1:1000))
% figure()
% b=prediction(:,:,2);
% figure()
% plot(data)

% figure()
% a=prediction(:,s{1}.n,1);
% b=prediction(:,s{1}.n,2);
% 
% plot(a); hold on
% plot(b); hold on
% plot(a+b); hold on
% 
% figure()
% b=prediction(:,:,2);
% figure()
% plot(data)
% 
% 166447

%          rss: 112.8717
%       rawrss: 116.9028

% 2d PRF
%          rss: 10.0340
%       rawrss: 56.4629
% x0: 14.3644
% y0: 2.0393
% s: 10
% s_major: 10
% s_minor: 10

                % 2d css
%             rss: 10.0340
%          rawrss: 56.4629
                 
            otherwise
                fprintf('[%s]:Unknown pRF model: %s: IGNORED!',mfilename,params.analysis.pRFmodel{n});
        end
    end

    %-----------------------------------
    % now put back the trends to the fits
    %-----------------------------------
      
    for mm=1:numel(s)
        nB = size(s{mm}.b,1);
        s{mm}.b(nB-ntrends+1:end,:) = s{mm}.b(nB-ntrends+1:end,:)+trendBetas;
    end
    
    %-----------------------------------
    % now we put back the temporary data from that slice
    %-----------------------------------
    model = rmSliceSet(model,s,slice);  
end;


%-----------------------------------
% recreate complete model if we used coarse sampling
%-----------------------------------
if params.analysis.coarseToFine,
    model = rmInterpolate(view, model, params);
end;

%-----------------------------------
% save and return output (if run interactively)
%-----------------------------------
rmFile = rmSave(view,model,params,1,'gFit');
view = viewSet(view,'rmFile',rmFile);


% that's it
return;
%-----------------------------------


%-----------------------------------
function model = initiateModel(params,d1,d2,nt)
% make the model struct with rmSet
fillwithzeros       = zeros(d1,d2);
fillwithinfs        = ones(d1,d2).*Inf;
   

        
% add a small number to sigmas because a pRF with 0 sigma does not exist
smallnumber   = 0.001;

% initiate all models
model = cell(numel(params.analysis.pRFmodel),1);
for n=1:numel(params.analysis.pRFmodel),
    % minimum for each model
    model{n} = rmSet;
    model{n} = rmSet(model{n},'x'   ,fillwithzeros);
    model{n} = rmSet(model{n},'y'   ,fillwithzeros);
    model{n} = rmSet(model{n},'s'   ,fillwithzeros+smallnumber);
    model{n} = rmSet(model{n},'rawrss',fillwithzeros);
    model{n} = rmSet(model{n},'rss' ,fillwithinfs);
    model{n} = rmSet(model{n},'df'  ,0);
    model{n} = rmSet(model{n},'ntrends',nt);
    % store hrf too since it is part of the model
    % fix me: we need to store all HRFs for each stimuli.
    % We could just store the entire params.stim struct.
    model{n} = rmSet(model{n},'whrf'     ,params.stim(1).hrfType);
    model{n} = rmSet(model{n},'hrfparams',params.stim(1).hrfParams);
    model{n} = rmSet(model{n},'hrfmax'   ,params.analysis.HrfMaxResponse);

    %--- model specific
    % These model description names are important, they later guide the
    % refine stage.
    switch lower(params.analysis.pRFmodel{n}),
        case {'onegaussian','one gaussian','default','standard',...
              'onegaussiann','one gaussian with nuisance factors'}
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+1));
            model{n} = rmSet(model{n},'desc','2D pRF fit (x,y,sigma, positive only)');
            
        case {'1dgaussian','1d gaussian','1d'}
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+1));
            model{n} = rmSet(model{n},'desc','1D pRF fit (x,sigma, positive only)');
            
        case {'oneovalgaussian','one oval gaussian','oval'}
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+1));
            model{n} = rmSet(model{n},'desc','oval 2D pRF fit (x,y,sigma_major,sigma_minor,theta)');
            
        case {'oneovalgaussianwithouttheta','one oval gaussian without theta'}
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+1));
            model{n} = rmSet(model{n},'desc','radial oval 2D pRF fit (x,y,sigma_major,sigma_minor)');
            
        case {'onegaussianunsigned','one gaussian unsigned','unsigned'}
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+1));
            model{n} = rmSet(model{n},'desc','unsigned 2D pRF fit (x,y,sigma)');
            
        case {'onegaussianlinked','one gaussian linked to neighbors'}
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+1));
            model{n} = rmSet(model{n},'desc','Linked 2D pRF fit (x,y,sigma)');

        case {'twogaussiansdog','dog','difference of gaussians'}
            model{n} = rmSet(model{n},'s2'  ,fillwithzeros+smallnumber);
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+2));
            model{n} = rmSet(model{n},'desc','Difference 2D pRF fit (x,y,sigma,sigma2, center=positive)');

        case {'two1dgaussiansdog','1ddog','1d difference of gaussians'}
            model{n} = rmSet(model{n},'s2'  ,fillwithzeros+smallnumber);
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+2));
            model{n} = rmSet(model{n},'desc','Difference 1D pRF fit (x,sigma, sigma2, center=positive)');
            
        case {'twogaussiansdogfixed','dogf','difference of gaussians fixed'}
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+1));
            model{n} = rmSet(model{n},'desc','Difference 2D pRF fit fixed (x,y,sigma,sigma2, center=positive)');
            
        case {'twogaussiansdogbetafixed','dogbf','difference of gaussians beta fixed'}
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+1));
            model{n} = rmSet(model{n},'desc','Difference 2D pRF fit beta fixed (x,y,sigma,sigma2, center=positive)');
            
        case {'twogaussianstog','tog','two gaussians on the same position'}
            model{n} = rmSet(model{n},'s2'  ,fillwithzeros+smallnumber);
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+2));
            model{n} = rmSet(model{n},'desc','Double 2D pRF fit (x,y,sigma,sigma2, center=positive)');

        case {'twogaussiansposonly','two gaussians','two prfs'}
            model{n} = rmSet(model{n},'x2'  ,fillwithzeros);
            model{n} = rmSet(model{n},'y2'  ,fillwithzeros);
            model{n} = rmSet(model{n},'s2'  ,fillwithzeros+smallnumber);
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+2));
            model{n} = rmSet(model{n},'desc','Two independent 2D pRF fit (2*(x,y,sigma, positive only))');
 
        case {'twogaussiansmirror','two gaussians mirrored','mirror'}
            model{n} = rmSet(model{n},'x2'  ,fillwithzeros);
            model{n} = rmSet(model{n},'y2'  ,fillwithzeros);
            model{n} = rmSet(model{n},'s2'  ,fillwithzeros+smallnumber);
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+1));
            model{n} = rmSet(model{n},'desc','Mirrored 2D pRF fit (2*(x,y,sigma, positive only))');
 

		case {'shiftedgaussians','two shifted gaussians'}
			model{n} = rmSet(model{n},'x2'  ,fillwithzeros);
			model{n} = rmSet(model{n},'y2'  ,fillwithzeros);
			model{n} = rmSet(model{n},'s2'  ,fillwithzeros+smallnumber);
			model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+1));
			model{n} = rmSet(model{n},'desc','Shifted 2D pRF fit (2*(x,y,sigma, positive only))');

        case {'addgaussian','sequentialgaussians','two sequential gaussians'}
            model{n} = rmSet(model{n},'x2'  ,fillwithzeros);
            model{n} = rmSet(model{n},'y2'  ,fillwithzeros);
            model{n} = rmSet(model{n},'s2'  ,fillwithzeros+smallnumber);
            model{n} = rmSet(model{n},'rss2',fillwithinfs);
            model{n} = rmSet(model{n},'rawrss2',fillwithinfs);
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+2));
            model{n} = rmSet(model{n},'desc','Sequential 2D pRF fit (2*(x,y,sigma, positive only))');
            
        case {'addgaussianlinked','sequentialgaussianslinked','two sequential gaussians linked to neighbors'}
            model{n} = rmSet(model{n},'x2'  ,fillwithzeros);
            model{n} = rmSet(model{n},'y2'  ,fillwithzeros);
            model{n} = rmSet(model{n},'s2'  ,fillwithzeros+smallnumber);
            model{n} = rmSet(model{n},'rss2',fillwithinfs);
            model{n} = rmSet(model{n},'rawrss2',fillwithinfs);
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+2));
            model{n} = rmSet(model{n},'desc','Linked sequential 2D pRF fit (2*(x,y,sigma, positive only))');

         case {'css' 'onegaussiannonlinear', 'onegaussianexponent'}
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+1));
            model{n} = rmSet(model{n},'desc','2D nonlinear pRF fit (x,y,sigma,exponent, positive only)');
            model{n} = rmSet(model{n},'exponent', fillwithzeros+1);
              
         case {'cssboxcar' 'onegaussiannonlinearboxcar', 'onegaussianexponentboxcar'}
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+2));
            model{n} = rmSet(model{n},'desc','2D nonlinear pRF fit with boxcar (x,y,sigma,exponent, positive only)');
            model{n} = rmSet(model{n},'exponent', fillwithzeros+1);
            
        case {'cst'}
            model{n} = rmSet(model{n},'b'   ,zeros(d1,d2,nt+2));
            model{n} = rmSet(model{n},'desc','2D CSS nonlinear spatiotemporal pRF fit');
            model{n} = rmSet(model{n},'exponent', fillwithzeros+1);
            
        otherwise
            fprintf('Unknown pRF model: %s: IGNORED!',mfilename,params.analysis.pRFmodel{n})
    end

end;

return;
%-----------------------------------
    