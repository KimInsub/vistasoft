function e = rmModelSearchFit_temporal(p,Y,Xv,Yv,tmodel, hrf, scan_num, t)
% rmModelSearchFit_temporal - actual fit function of rmSearchFit
%
% error = rmModelSearchFit(p,Y,trends,Xgrid,YGrid,stimulusMatrix);
%
% Basic barebones fit of a single time-series. Error is returned in
% percentage: 100% is RSS of unfitted time-series. This way we can quantify
% the improvement of the fit independend of the variation in the raw
% time-series.
%
% 2006/06 SOD: wrote it.
% 2006/12 SOD: modifications for fmincon, this is litterally called >10000
% times so we cut every corner possible. 
% 2010/02 SOD: evaluated lscov this did not improve performance here (using
% profiler)

% fprintf(1,'[cst] search optimization .... \n');

% make RF (taken from rfGaussian2d)
Xv = Xv - p(1);   % positive x0 moves center right
Yv = Yv - p(2);   % positive y0 moves center up
RF = exp( (Yv.*Yv + Xv.*Xv) ./ (-2.*(p(3).^2)) );


%%%%% [cst]  params
stim = tmodel.chan_preds;
nChan = tmodel.num_channels;

% chan == 1
for cc = 1:nChan
    
    if size(tmodel.chan_preds,1) == 3 % for abc        
        pred = cat(1, (stim{1,cc}*RF).^p(4), (stim{2,cc}*RF).^p(4) ,(stim{3,cc}*RF).^p(4));
%         pred = cat(1,pred,pred);
    else
        
        pred = (stim{cc}*RF).^p(4); % for CSS

    end

%     pred = (stim{cc}*RF);

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
%     X    = [prediction(:,n,1) prediction(:,n,2) trends];
if nChan == 2
   preds = cell2mat(prediction);
   normTs = max(preds(:,1))/max(preds(:,2));
   prediction{2} = preds(:,2) *normTs;
end
X = [cell2mat(prediction) t];

% fit - inlining pinv
%b = pinv(X)*Y; 
[U,S,V] = svd(X,0);

s = diag(S); 
tol = numel(X) * eps(max(s));
r = sum(s > tol);
if (r == 0)
    pinvX = zeros(size(X'));
else
    s = diag(ones(r,1)./s(1:r));
    pinvX = V(:,1:r)*s*U(:,1:r)';
end
b = pinvX*Y;

% do for both negative and positive fits
e = norm(Y - X*b);
% 
% compute residual sum of squares (e)
% e = norm(Y - X*abs(b));

% if nChan == 1
%     if b(1)>0
%         e = norm(Y - X*b);
%     else
%         e = norm(Y).*(1+sum(abs(b(1))));
%     end
%     
% elseif nChan == 2
%     
%     if b(1)>0  &&  b(2)>0
%         e = norm(Y - X*b);
%     else
%         e = norm(Y).*(1+sum(abs(b(1))));
%     end
% end
return;



 