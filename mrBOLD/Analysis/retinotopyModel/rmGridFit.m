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
% 2020/12 IK: Crossvalidation & spatiotemporal model added

if notDefined('view'),   error('Need view struct'); end
if notDefined('params'), error('Need params'); end 


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
switch lower(params.wData)
    case {'fig','roi'}
        loopSlices = 1;
    otherwise
        loopSlices = 1:params.analysis.nSlices;
end
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
    for n=1:numel(s)-1
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
        if ismember(n, round((1:10)/10* numel(s)-1)) % every 10% draw a dot
            fprintf(1,'.');drawnow;
        end
    end
    
    clear n s rf pred;
    fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now));
    drawnow;
    
    % nonlinear spatio-temproal model prediction-grid creation
elseif strcmp(params.analysis.pRFmodel{1}, 'st')
    
    % cache grid -- as it takes very very long time to generate
    % if we already created the gird, simply just load it!
    % stimGrid, prediction, grid are the same thing
    if exist(params.analysis.predFile, 'file')
        disp('***st predfile exists --- loading...')
        disp(params.analysis.predFile)
        load(params.analysis.predFile);
    else
        stimGrid = rmGridstPred(params);
    end
    % assign each grid to params
    for es = 1:length(params.stim)
        params.stim(es).prediction = stimGrid(es).prediction;
    end
    
    prediction=cat(1,stimGrid.prediction);
    clear stimGrid;
    % other nonlinear models (ex) CSS
else
    allstimimages = params.analysis.allstimimages_unconvolved;
    prediction = zeros(size(allstimimages,1),n,'single');
    fprintf(1,'[%s]:Making %d model samples:',mfilename,n);
    drawnow;tic;
    for n=1:numel(s)-1
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
        if ismember(n, round((1:10)/10* numel(s)-1)) % every 10% draw a dot
            fprintf(1,'.');drawnow;
        end
    end
    
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

% % % normalize channels
% size(tmodel.run_preds,2)
% if strcmp(params.analysis.pRFmodel{1}, 'st')
% if size(params.stim(1).prediction,3) == 2
%     for ii = 1:size(params.stim(1).prediction,2) % 720       82608           2
%         maxS = max(max(params.stim(1).prediction(:,ii,1)));
%         maxT = max(max(params.stim(1).prediction(:,ii,2)));
%         normTs(ii) = maxS / maxT;
% %         tmodel.run_preds(:,ii,2) = tmodel.run_preds(:,ii,2) * normTs(ii);
%     end
% end
% end

% go loop over slices
for slice=loopSlices
    %-----------------------------------
    % Place datasets behind each other. This is a rather crude way of
    % stimultaneously fitting both. Due to this we cannot
    % prewhiten (we could zeropad/let the trends deal with this/not care).
    %-----------------------------------
 
    
    % Prepare data structure for cross-validation
    if params.analysis.cv == 1
        
        % reset the prediction before..
        % we are doing cross validation here
        prediction = [];
        
        % perfom k-fold (3 fold in current example) cross valiation
        rng('default') % For reproducibility
        
        nStim = length(params.stim);
        cv_split = cvpartition(nStim,'KFold',3);
        numFolds = cv_split.NumTestSets;
        
        % create: data{fold}, prediction{fold}
        % train_grid (predictions) => Time X GRID X Channel
        % data => Time X Voxel
        % prediction and data should match in their time domain
        for fold = 1:numFolds
            trainSet = find(cv_split.training(fold));
            testSet  =  find(cv_split.test(fold));
            
            % concat and load  traning     data---
            [traindata{fold},params] = rmLoadData(view, params, slice, ...
                params.analysis.coarseToFine, [], trainSet);
           % [data, params] = rmLoadData(view, params, slice,...
            % concat and load  testing     data---
            [testdata{fold},params] = rmLoadData(view, params, slice, ...
                params.analysis.coarseToFine, [], testSet);

            % gather the prediction  data---
            train_grid = cell(length(trainSet),1);
            [train_grid{:}] = params.stim(trainSet).prediction;
            train_grid = cell2mat(train_grid);
            
            test_grid = cell(length(testSet),1);
            [test_grid{:}] = params.stim(testSet).prediction;
            test_grid = cell2mat(test_grid);

            % [IK] edited to account for cross validation
            [train_trend, train_ntrends, train_dcid] = rmMakeTrends(params,trainSet);
            [test_trend, test_ntrend, test_dcid] = rmMakeTrends(params,testSet);

            
            % save cv information to a strcuct
            df(fold).info = cv_split;
            df(fold).numFolds = numFolds;
            df(fold).train_set = trainSet;
            df(fold).test_set = testSet;
            df(fold).train_data = traindata{fold};
            df(fold).test_data = testdata{fold};
            df(fold).train_grid = train_grid;
            df(fold).test_grid = test_grid;

            df(fold).train_trend = train_trend;
            df(fold).test_trend = test_trend;
            
            df(fold).train_ntrends = train_ntrends;
            df(fold).test_ntrend = test_ntrend;
            
            df(fold).train_dcid = train_dcid;
            df(fold).test_dcid = test_dcid;

        end
        
        % remove variables that suck up storage resources
        clear train_grid; 
        params.stim = rmfield( params.stim , 'prediction' ) ;
%         params.stim.images_org = [];
%         params.analysis.allstimimages = [];
%         params.analysis.allstimimages_unconvolved = [];
        
        params.stim = rmfield( params.stim , 'images_org' ) ;
        params.analysis = rmfield( params.analysis , 'allstimimages' ) ;
        params.analysis = rmfield( params.analysis , 'allstimimages_unconvolved' ) ;

    else
        [data, params] = rmLoadData(view, params, slice,...
            params.analysis.coarseToFine);
        df.train_set = '';
        df.test_set = '';

    end
    
    % SOLVE & Tidy and store data into 'df' structure
    if params.analysis.cv == 1
        numFolds = df(1).numFolds;
        for fold = 1:numFolds
            
            % assign variables according to each fold
            data = df(fold).train_data;
            prediction = df(fold).train_grid;
            trends = df(fold).train_trend; trends = single(trends);
            ntrends = df(fold).train_ntrends;
            dcid = df(fold).train_dcid;
            
            % solve GRID! - for train_data set
            model = rmGridSolve(params,data,prediction,trends,ntrends,dcid,slice,nSlices);
            
            if params.analysis.coarseToFine
                model = rmInterpolate(view, model, params);
            end
            
            % save and clean df output
            df(fold).x0            = model{1}.x0;
            df(fold).y0            = model{1}.y0;
            df(fold).sigma         = model{1}.sigma;
            df(fold).exponent      = model{1}.exponent;
            df(fold).beta          = model{1}.beta;
            df(fold).train_pred    = model{1}.pred_X;
            
            switch lower(params.wData)
                case {'roi'}
                    df(fold).roi_name      = params.roi.name;
                    df(fold).coords        = params.roi.coords;
                    df(fold).coordsIndex   = params.roi.coordsIndex;
                otherwise
                    df(fold).roi_name      = 'all';
                    df(fold).coords        = [];
                    df(fold).coordsIndex   = [];
            end
            

        end
    else
        %-----------------------------------
        %--- make trends to fit with the model (discrete cosine set)
        %-----------------------------------
        [trends, ntrends, dcid] = rmMakeTrends(params);
        trends = single(trends);

        model = rmGridSolve(params,data,prediction,trends,ntrends,dcid,slice,nSlices);
        
    end
    
    
    
end


%-----------------------------------
% recreate complete model if we used coarse sampling
%-----------------------------------
if params.analysis.coarseToFine
    model = rmInterpolate(view, model, params);
end

%-----------------------------------
% get Timecouse results
%-----------------------------------
% % load no smooting data;
% [data, params] = rmLoadData(view, params, slice,...
%     params.analysis.coarseToFine);
% data(isnan(data)) = 0;
% rawdata = single(data);
%  
% [trends, ntrends, dcid] = rmMakeTrends(params);
% trends = single(trends);
% trendBetas = pinv(trends)*data;
% if params.analysis.doDetrend
%     data = data - trends*trendBetas;
%     
% end
% 
% rawdata   = rmDecimate(rawdata,params.analysis.coarseDecimate);
% data   = rmDecimate(data,params.analysis.coarseDecimate);
% trends = rmDecimate(trends,params.analysis.coarseDecimate);
% 
% if ~params.analysis.doDetrend 
%     trends = zeros(size(trends));
% end

% tc.rawdata = rawdata;
% tc.data = data;
% tc.trends = trends;
% tc.predictions = prediction;
% tc.beta = model{1}.beta;
% 
% for mm = 1:numel(s)
%     
%     pp = s{mm}.pred_X;
%     
%     if size(model{1}.beta,2) == 1
%         bb = squeeze(tc.beta(mm,:,:));
%         bb = num2cell(bb,1);
%        
%     else
%         bb = squeeze(tc.beta(mm,:,:));
%         bb = num2cell(bb,2); bb = cellfun(@transpose,bb,'UniformOutput',false);
%     end
%     
%     
%     predictors = [];
%     for i = 1:size(pp,2)
%         predictors{i} = [squeeze(pp(:,i,:))  tc.trends];
%     end
%     predictors = predictors';
%     tc.X{mm} = predictors;
%     tc.result_tc{mm} =s{mm}.pred_X;
%     tc.result_beta_tc{mm} = cell2mat(cellfun(@(x,y) x*y, predictors,bb, 'UniformOutput',false)');
%     
%     res = tc.rawdata - tc.result_beta_tc{mm};
%     res_var = sum(res .^ 2) ./ sum((tc.rawdata - mean(tc.rawdata)) .^ 2);
%     tc.varexp{mm} = 1 - res_var;
% 
% end
% 



%-----------------------------------
% save and return output (if run interactively)
%-----------------------------------
rmFile = rmSave(view,model,params,1,'gFit',df);
view = viewSet(view,'rmFile',rmFile);

 
% that's it
return;
%-----------------------------------

