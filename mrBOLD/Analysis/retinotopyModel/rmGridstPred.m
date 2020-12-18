function stimGrid = rmGridstPred(params)
%% get Variables

c = Constants.getTemporalParams.temporalParams;
temp_type = params.analysis.temporalModel;

for i = 1:length(c)
    if strcmp(c{i}.type, temp_type)
        idx = i;
    end
end

temporal_param = c{idx}.prm;
fs             = c{idx}.fs;
num_channels   = c{idx}.num_channels;
scale_max = @(x) x./max(x(:));

% define hrf
hrf = canonical_hrf(1 / fs, [5 14 28]);

% define TR
tr = Constants.getTemporalParams.tr; % seconds
tnorm = Constants.getTemporalParams.tnorm; % seconds

%% batch loop
% n = numel(params.analysis.x0);
% s = [[1:ceil(n./1000):n-2] n+1];
% 

%% apply spatiotermpoal filter
for es = 1:length(params.stim)
    
    % batch loop
    n = numel(params.analysis.x0);
    s = [[1:ceil(n./1000):n-2] n+1];

    allstimimages = params.stim(es).images_unconvolved;
    t = 0.001 : 0.001 : size(allstimimages,2)/fs;
    rsp =st_tModel(temp_type,temporal_param, allstimimages, t); % rsp = time (ms) x space
    fprintf(1,'[%s]:Making %d model samples:',mfilename,n);
    prediction = zeros(size(rsp{1},1)/fs,n,num_channels,'single'); % prediction = time (s) x space
    % loop over grid
    tic
    for n=1:numel(s)-1
        % make rfs
        rf   = rfGaussian2d(params.analysis.X, params.analysis.Y,...
            params.analysis.sigmaMajor(s(n):s(n+1)-1), ...
            params.analysis.sigmaMinor(s(n):s(n+1)-1), ...
            params.analysis.theta(s(n):s(n+1)-1), ...
            params.analysis.x0(s(n):s(n+1)-1), ...
            params.analysis.y0(s(n):s(n+1)-1));
        
        % convolve rf with stimulus for each t-channel
        for cc=1:num_channels
            pred = rsp{cc}*rf;
            % apply css if needed
            pred = bsxfun(@power, pred, params.analysis.exponent(s(n):s(n+1)-1)');
            pred = double(pred);
            
            % apply hrf
            pred_cell = num2cell(pred,1)';
            npixel_max = size(s(n):s(n+1)-1,2);
            
            % use mrvista HRF
            %   params.stim(n).images = filter(params.analysis.Hrf{n}, 1, params.stim(n).images'); % images: pixels by time (so images': time x pixels)
            %                 hrf = params.analysis.Hrf;
            curhrf = repmat({hrf'}, npixel_max, 1);
            pred_hrf = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 /tr), ...
                pred_cell, curhrf, 'uni', false);
            pred_hrf = cellfun(@transpose,pred_hrf,'UniformOutput',false);
            pred_hrf=cell2mat(pred_hrf)';
            
            
            %             
            if cc ==2
                pred_hrf = pred_hrf*tnorm;
            end
            
            % store
            prediction(:,s(n):s(n+1)-1,cc) = pred_hrf;
            %             prediction{n} = pred_hrf;
            
        end
        
        
        if ismember(n, round((1:10)/10* numel(s)-1)) % every 10% draw a dot
            fprintf('Generating predictors for %s model...%d/%d (stimulus = %d) \n', ...
                params.analysis.temporalModel,n,numel(s)-1, es); drawnow;
        end
        stimGrid(es).prediction = prediction;

    end
    
    
    clear n s rf pred pred_hrf pred_cell;
    fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now));
    drawnow;
end
    
save(params.analysis.predFile, 'stimGrid', '-v7.3');


end





% 
% 
% if isfile(params.analysis.predFile)
%     disp('***st predfile exists --- loading...')
%     params.analysis.predFile
%     load(params.analysis.predFile);
% elseif contains(params.analysis.predFile,'abc')
%     a =strrep(params.analysis.predFile,'abc','a');
%     load(a);
%     chan_preds1 = tmodel.chan_preds;
%     preds{1} = tmodel.run_preds;
%     
%     b =strrep(params.analysis.predFile,'abc','b');
%     load(b);
%     chan_preds2 = tmodel.chan_preds;
%     preds{2} = tmodel.run_preds;
%     
%     c =strrep(params.analysis.predFile,'abc','c');
%     load(c);
%     chan_preds3 = tmodel.chan_preds;
%     
%     preds{3} = tmodel.run_preds;
%     
%     tmodel.run_preds = cat(1,preds{1},preds{2}, preds{3});
%     tmodel.chan_preds=cat(1,chan_preds1,chan_preds2,chan_preds3);
%     save(params.analysis.predFile, 'tmodel', '-v7.3');
% else
%     fprintf(1,'[%s]:Making %d model samples:',mfilename,n);
%     tmodel = st_createIRF_grid(params);
%     
%     
%     prediction = zeros(size(tmodel.chan_preds{1},1)/tmodel.fs,n,tmodel.num_channels,'single');
%     % loop over grid
%     tic
%     for n=1:numel(s)-1
%         
%         
%         % make rfs
%         rf   = rfGaussian2d(params.analysis.X, params.analysis.Y,...
%             params.analysis.sigmaMajor(s(n):s(n+1)-1), ...
%             params.analysis.sigmaMinor(s(n):s(n+1)-1), ...
%             params.analysis.theta(s(n):s(n+1)-1), ...
%             params.analysis.x0(s(n):s(n+1)-1), ...
%             params.analysis.y0(s(n):s(n+1)-1));
%         
%         % convolve rf with stimulus for each t-channel
%         for cc=1:tmodel.num_channels
%             pred = tmodel.chan_preds{cc}*rf;
%             
%             % apply css if needed
%             pred = bsxfun(@power, pred, params.analysis.exponent(s(n):s(n+1)-1)');
%             pred = double(pred);
%             
%             % apply hrf
%             pred_cell = num2cell(pred,1)';
%             npixel_max = size(s(n):s(n+1)-1,2);
%             
%             % use mrvista HRF
%             %   params.stim(n).images = filter(params.analysis.Hrf{n}, 1, params.stim(n).images'); % images: pixels by time (so images': time x pixels)
%             %                 hrf = params.analysis.Hrf;
%             hrf = tmodel.irfs.hrf;
%             curhrf = repmat(hrf, npixel_max, 1);
%             pred_hrf = cellfun(@(X, Y) convolve_vecs(X, Y, tmodel.fs, 1 / tmodel.tr), ...
%                 pred_cell, curhrf, 'uni', false);
%             pred_hrf = cellfun(@transpose,pred_hrf,'UniformOutput',false);
%             pred_hrf=cell2mat(pred_hrf)';
%             
%             if cc ==2
%                 pred_hrf = pred_hrf*tmodel.normT;
%             end
%             
%             % store
%             prediction(:,s(n):s(n+1)-1,cc) = pred_hrf;
%             %             prediction{n} = pred_hrf;
%             
%         end
%         if ismember(n, round((1:10)/10* numel(s)-1)), % every 10% draw a dot
%             fprintf('Generating predictors for %s model...%d/%d \n', ...
%                 params.analysis.temporalModel,n,numel(s)-1)
%         end
%         
%     end
%     
%     tmodel.run_preds = prediction;
%     %         params.analysis.temporal.tmodel = tmodel;
%     save(params.analysis.predFile, 'tmodel', '-v7.3');
%     
%     clear n s rf pred pred_hrf pred_cell;
%     fprintf(1, 'Done[%d min].\t(%s)\n', round(toc/60), datestr(now));
%     drawnow;
% end
% 
