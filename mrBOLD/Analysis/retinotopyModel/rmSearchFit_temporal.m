function model = rmSearchFit_temporal(model, data, params, wProcess, t,trainSet)
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
fprec = 1e-6;
fmin_options = optimoptions('fmincon', 'Display', 'iter', ...
    'StepTolerance', fprec, 'UseParallel', true);

expandRange   = params.analysis.fmins.expandRange;

% convert to double just in case
params.analysis.X = double(params.analysis.X);
params.analysis.Y = double(params.analysis.Y);

% tNorm?
do_stig_t_Norm = 0;

% get Stimulus info
stim = cell(length(trainSet),1);
[stim{:}] = params.stim(trainSet).images_unconvolved;
stim = cellfun(@transpose,stim,'UniformOutput',false);
stim = cellfun(@double,stim,'UniformOutput',false);

%% get Temporal params
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

hrf = canonical_hrf(1 / fs, [5 14 28]);
tr = Constants.getTemporalParams.tr; % seconds
tnorm = Constants.getTemporalParams.tnorm; % seconds

% asign temproal params to model
model.tau_s      = ones(size(model.x0))*temporal_param(1);
model.tau_ae     = ones(size(model.x0))*temporal_param(2);
model.Lp         = ones(size(model.x0))*temporal_param(3);
model.Kp         = ones(size(model.x0))*temporal_param(4);
model.Kn         = ones(size(model.x0))*temporal_param(5);
model.weight     = ones(size(model.x0))*temporal_param(6);

% tau_ae 
temporal_param(2) = temporal_param(2)/10000;
% remove shift param
temporal_param(end) = [];

%%
%params.analysis.allstimimages_unconvolved = double(params.analysis.allstimimages_unconvolved);
data = double(data);

% get starting upper and lower range and reset TolFun
% (raw rss computation (similar to norm) and TolFun adjustments)
model.s = model.s_major;
% [range, TolFun] = rmSearchFit_range(params,model,data);
% [range, TolFun] = rmSearchFit_range_temporal(params,model,data);
[range, TolFun] =  rmSearchFit_range_temporal_linear(params,model,data);


% amount of negative fits
nNegFit  = 0;
vethresh = params.analysis.fmins.vethresh;
trends   = t.trends;
t_id     = t.dcid+num_channels;

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
    vData=num2cell(reshape(vData,size(stim{1},1)/fs,[]),1)';
%     vData = vData(trainSet);
    
    % create objective function
     [comp_ws, conv_nb, pred_bs, obj_fun]  = st_obj_fun_2ch_exp_sig(stim,vData,params.analysis.X,params.analysis.Y);
    
    
    % reset tolFun: Precision of evaluation function.
    % We define RMS improvement relative to the initial raw 'no-fit' data
    % RMS. So, 1 means stop if there is less than 1% improvement on the fit:
    % searchOptions = optimset(searchOptions,'tolFun',optimget(params.analysis.fmins.options,'tolFun')./100.*rawrss);
    % optimset and optimget are a little slow so:
%     searchOptions.TolFun = TolFun(ii);
    
    % actual fitting routine
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
        %         outParams = ...
        %             fmincon(@(x) rmModelSearchFit_temporal(x,vData,...
        %             params.analysis.X,...
        %             params.analysis.Y,...
        %             tmodel,...
        %             tmodel.irfs.hrf,...
        %             params.analysis.scans,...
        %             trends),...
        %             range.start(:,vi),[],[],[],[],range.lower(:,vi),range.upper(:,vi),...
        %             @(x) distanceCon(x,range.start(:,vi),range.step(:,vi).*expandRange),searchOptions);
        %
        %
        
        init_p = [range.start(1:3,vi)' temporal_param];
        if do_stig_t_Norm  == 1
            lower_p = [range.lower(1:3,vi)' 4  1 .01 .1 .1 0.5];
            upper_p = [range.upper(1:3,vi)' 20 4 .5   6  6 0.5];
        else
            lower_p = [range.lower(1:3,vi)' 4  1 .01 .1 .1 0];
            upper_p = [range.upper(1:3,vi)' 20 4 .5   6  6 1];
        end
        
%         outParams = init_p;
        outParams = fmincon(obj_fun, init_p, [], [], [], [], ...
            lower_p, upper_p, ...
            @(x) distanceCon(x,range.start(1:3,vi),range.step(:,vi).*expandRange), fmin_options);
%         
%         outParams = fmincon(obj_fun, init_p, [], [], [], [], ...
%             lower_p, upper_p, ...
%             @(x) distanceCon(x,range.start(1:3,vi),range.step(:,vi).*expandRange), searchOptions);

        

%         outParams2 = fmincon(obj_fun, init_p, [], [], [], [], ...
%             lower_p, upper_p, ...
%             [], fmin_options); 
%         
%         outParams3 = fmincon(obj_fun, init_p, [], [], [], [], ...
%             lower_p, upper_p, ...
%             @(x) distanceCon(x,range.start(1:3,vi),range.step(:,vi).*expandRange), fmin_options);
%         
%         outParams4 = fmincon(obj_fun, init_p, [], [], [], [], ...
%             lower_p, upper_p, ...
%             [], searchOptions);
        
   
    
    %%
%     preds = conv_nb(stim,outParams(1),outParams(2),outParams(3),outParams(4),outParams(5),outParams(6),outParams(7),outParams(8),outParams(9));
%     preds = cell2mat(preds);
%     
%     X  = [preds trends];
%     b    = pinv(X)*cell2mat(vData);
%     rss  = norm(cell2mat(vData)-X*b).^2;

    %% for debugging
% %     % make RF, prediction and get rss,b
    Xv = params.analysis.X-outParams(1);
    Yv = params.analysis.Y-outParams(2);
%     n  = outParams(4);
    rf = exp( (Yv.*Yv + Xv.*Xv) ./ (-2.*(outParams(3).^2)) );
    
    rfstim=cellfun(@(x,y) x*y, stim,repmat({rf},size(stim,1),1),'UniformOutput',false);
    rfstim=cellfun(@transpose, rfstim,'UniformOutput',false);
    samplet = 0.001 : 0.001 : size(rfstim{1},2)/fs;

    outParams(5) =    outParams(5) * 10000;
    rsp = cellfun(@(ty, tp, st, t) st_tModel(ty,tp,  st, t), ...
    repmat({temp_type},size(stim,1),1), repmat({[outParams(4:end) 0]},size(stim,1),1), ...
    rfstim, repmat({samplet},size(stim,1),1), 'UniformOutput', false);

    rsp = vertcat(rsp{:});
    curhrf = repmat({hrf'}, size(rsp));
    preds = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 /tr), ...
        rsp, curhrf, 'uni', false);
    preds = cell2mat(preds);

    if num_channels == 2
        preds(:,3) = []; % remove the added version
        if do_stig_t_Norm == 1
            normTs = max(preds(:,1))/max(preds(:,2));
            preds(:,2)= preds(:,2) * normTs;            
        end
    end
    
    X  = [preds trends];
    b    = pinv(X)*cell2mat(vData);
    rss  = norm(cell2mat(vData)-X*b).^2;
    
    %%  
    % store results only if the first beta is positive, somehow fmincon
    % outputs negative fits. If the fit is negative keep old (grid) fit. We
    % do adjust the rss, so it won't be accidentally counted as a 'good'
    % fit.
    if num_channels == 1
        if b(1)>0
            model.x0(vi)         = outParams(1);
            model.y0(vi)         = outParams(2);
            model.s(vi)          = outParams(3);
            model.s_major(vi)    = outParams(3);
            model.s_minor(vi)    = outParams(3);
            model.s_theta(vi)    = 0;
            model.exponent(vi)   = 1;
%             model.exponent(vi)   = outParams(4);
            model.rss(vi)        = rss;
            model.b([1 t_id],vi) = b;
            model.pred_X(:,vi,1) = preds(:,1);
            
        else
            % change the percent variance explained to be just under the
            % current vethresh. So it counts as a 'coarse'-fit but can still be
            % included in later 'fine'-fits
            model.rss(vi)  = (1-max((vethresh-0.01),0)).*model.rawrss(vi);
            nNegFit = nNegFit + 1;
        end
    elseif num_channels == 2
        if b(1)>0 && b(2)>0
            model.x0(vi)         = outParams(1);
            model.y0(vi)         = outParams(2);
            model.s(vi)          = outParams(3);
            model.s_major(vi)    = outParams(3);
            model.s_minor(vi)    = outParams(3);
            model.s_theta(vi)    = 0;
            model.exponent(vi)   = 1;
            model.tau_s(vi)      = outParams(4);
            model.tau_ae(vi)     = outParams(5)*10000;
            model.Lp(vi)         = outParams(6);
            model.Kp(vi)         = outParams(7);
            model.Kn(vi)         = outParams(8);
            model.weight(vi)     = outParams(9);

%             model.exponent(vi)   = outParams(4);
            model.rss(vi)        = rss;
            model.b([1 2 t_id],vi) = b;
            model.pred_X(:,vi,1) = preds(:,1);
            model.pred_X(:,vi,2) = preds(:,2);
            
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

