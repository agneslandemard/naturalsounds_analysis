%%% Predict human data based on ferret components
% parameters of analysis
clear I;
I.target = 'human';
% if I.target  = 'human' -> fig S9
% if I.target = 'ferret' -> fig S10

%% Load and format components

clear LP;
LP.experiment = 'natural';
LP.name = version_name;
FC = LoadComponents(LP);

% load human components
load([additional_path 'Human/HumanComponentsForPrediction.mat'],'human_R','human_R1','human_R2');

% load ferret components
load([analysis_path LP.experiment '_' LP.name '/FerretComponentsForPrediction.mat'],'ferret_R','ferret_R1','ferret_R2');

X = load([data_path FC.hemis{1} FC.data_suffix],'param');
P = X.param;

% keep only sounds that are common to both species
ns = find(isnan(nanmean(human_R,[2 3])));
SoundstoKeep = setdiff(1:P.snd.Nsound,ns);
ModelstoKeep = [1 5]; % models to consider

switch I.target
    case 'human'
        target_R = human_R(SoundstoKeep,ModelstoKeep,:);
        target_R1 = human_R1(SoundstoKeep,ModelstoKeep,:);
        target_R2 = human_R2(SoundstoKeep,ModelstoKeep,:);
        
        predictor_R = ferret_R(SoundstoKeep,ModelstoKeep,:);
        predictor_R1 = ferret_R1(SoundstoKeep,ModelstoKeep,:);
        predictor_R2 = ferret_R2(SoundstoKeep,ModelstoKeep,:);
        
        CompOrder = 1:size(target_R,3);
    case 'ferret'
        target_R = ferret_R(SoundstoKeep,ModelstoKeep,:);  
        target_R1 = ferret_R1(SoundstoKeep,ModelstoKeep,:);
        target_R2 = ferret_R2(SoundstoKeep,ModelstoKeep,:);
   
        predictor_R = human_R(SoundstoKeep,ModelstoKeep,:);
        predictor_R1 = human_R1(SoundstoKeep,ModelstoKeep,:);
        predictor_R2 = human_R2(SoundstoKeep,ModelstoKeep,:);
        
         CompOrder = [7 5 3 2 6 1 4 8];
        % can be modified
    otherwise
        error('Target must be human or ferret not %s', I.target);
end

%% Quantify variance in target for synthetic and natural-synthetic parts of components response

switch I.target
    case 'human'
        load([additional_path 'Human/HumanData.mat'],'SrecoFull_human_allSubjects','xr','yr','xl','yl','n_subj');
        S = reshape(SrecoFull_human_allSubjects,FC.n_stims_per_cond,FC.n_conds,xr*yr+xl*yl,n_subj,2);
        S = reshape(nanmean(S(SoundstoKeep,ModelstoKeep,:,13:17,:),[4 5]),length(SoundstoKeep)*length(ModelstoKeep),xr*yr+xl*yl);
        W = pinv(reshape(human_R(SoundstoKeep,ModelstoKeep,:),length(SoundstoKeep)*length(ModelstoKeep),size(human_R,3)))*S; 
        RMS_w = sqrt(nanmean(W.^2,2));

    case 'ferret'
        RMS_w = sqrt(nanmean(FC.W_ica_train.^2,2));
        
end

target_R1_n = nan(size(target_R1));
target_R2_n = nan(size(target_R2));
for i = 1 : size(target_R,3)
    target_R1_n(:,:,i) = target_R1(:,:,i).*RMS_w(i);
    target_R2_n(:,:,i) = target_R2(:,:,i).*RMS_w(i);
end

%% Prediction of target components

[n_stims_per_cond, n_conds, n_target_components, ~] = size(target_R);
[~, ~, n_predictor_components,~] = size(predictor_R);
assert(size(target_R,1) == size(predictor_R,1))
assert(size(target_R,2) == size(predictor_R,2))

Y = reshape(target_R, [n_stims_per_cond*n_conds, n_target_components]);
X = reshape(predictor_R, [n_stims_per_cond*n_conds, n_predictor_components]);

% folds used for cross-validation
ResetRandStream2(1);
n_folds = 9;
folds = subdivide(n_stims_per_cond, n_folds)';
folds = repmat(folds, 1, n_conds);
folds = folds(randperm(n_stims_per_cond),:);
folds_unwrap = reshape(folds, [n_stims_per_cond*n_conds, 1]);

% prediction
Yh = regress_predictions_from_3way_crossval(X, Y, ...
    'test_folds', folds_unwrap, 'train_folds', folds_unwrap, ...
    'demean_feats',true , 'std_feats', true, 'method', 'ridge','K', 2.^(-100:100));

% split out conditions
target_prediction = reshape(Yh, [n_stims_per_cond, n_conds, n_target_components]);

%% Plot measured vs predicted component responses, for synthetic and difference natural - synthetic

figh = figure;
n_rows = 2;
n_cols = n_target_components;

set(figh, 'Position', [0, 0, 200*n_cols, 200*n_rows]);
for ic = 1:n_target_components
    
    i = CompOrder(ic);
    
    % Plot measured vs predicted response to synthetic
    subplot(n_rows, n_cols, ic);
    hold on;
    X = cat(2,[target_prediction(:,end,i); target_R(:,end,i)],[target_prediction(:,1,i)-target_prediction(:,end,i);target_R(:,1,i)-target_R(:,end,i)]);
    maxrange = max(diff([min(X); max(X)]));
    bounds = min([target_prediction(:,end,i); target_R(:,end,i)])+[0 maxrange];
    bounds = bounds + [-1 1]*diff(bounds)*0.1;
    plot(bounds, bounds, 'r--','LineWidth', 2);
    scatter(target_prediction(:,end,i), target_R(:,end,i),[],P.plt.SoundColors(SoundstoKeep,:),'filled')
    xlabel('Predicted');
    ylabel('Measured');
    title([upper(I.target(1)) num2str(i)]);
    xlim(bounds);
    ylim(bounds);
    axis equal tight
    if bounds(1)<0
        vline(0,'k--')
        hline(0,'k--')
    end
    
    
    % Plot measured vs predicted difference between natural and synthetic
    subplot(n_rows, n_cols, ic+n_cols);
    hold on;
    bounds=min([target_prediction(:,1,i)-target_prediction(:,end,i); target_R(:,1,i)-target_R(:,end,i)])+[0 maxrange];
    bounds = bounds + [-1 1]*diff(bounds)*0.1;
    plot(bounds, bounds, 'r--','LineWidth', 2);
    scatter(target_prediction(:,1,i)-target_prediction(:,end,i), target_R(:,1,i)-target_R(:,end,i),[],P.plt.SoundColors(SoundstoKeep,:),'filled')
    xlabel('Predicted');
    ylabel('Measured');
    xlim(bounds);
    ylim(bounds);
    axis equal tight
    
    if bounds(1)<0
        vline(0,'k--')
        hline(0,'k--')
    end
end


%% Prediction using independent splits, for quantification

Nbs = 1000; % number of bootstraps
[n_stims_per_cond, n_conds, n_target_components, ~] = size(target_R1);
bstraps = randi(n_stims_per_cond,n_stims_per_cond,Nbs);

target_variance = nan(n_target_components,2,Nbs);
predicted_variance = nan(n_target_components,2,Nbs,2);
 

splits = {target_R1,target_R2;
    predictor_R1,predictor_R2};

ResetRandStream2(1);
nc_nse_all = nan(n_target_components,2,2);

% Check test-retest of components to decide if noise-correction
trt = nan(n_target_components,1);
for ic = 1:n_target_components
    trt(ic,1)= NSE(squeeze(target_R1(:,end,ic)),squeeze(target_R2(:,end,ic)),1)  ;  
    trt(ic,2)= NSE(squeeze(target_R1(:,1,ic)-target_R1(:,end,ic)),squeeze(target_R2(:,1,ic)-target_R2(:,end,ic)),1)  ; 
end 


for spl = 1:2
    
    target_R1_tmp=splits{1,1};
    predictor_R_tmp=splits{2,spl};
    
    [~, ~, n_predictor_components,~] = size(predictor_R_tmp);
    assert(size(target_R1_tmp,1)==size(predictor_R_tmp,1))
    assert(size(target_R1_tmp,2)==size(predictor_R_tmp,2))
    
    Y = reshape(target_R1_tmp, [n_stims_per_cond*n_conds, n_target_components]);
    X = reshape(predictor_R_tmp, [n_stims_per_cond*n_conds, n_predictor_components]);
    
    
    % folds used for cross-validation
    
    n_folds = 9;
    folds = subdivide(n_stims_per_cond, n_folds)';
    folds = repmat(folds, 1, n_conds);
    folds = folds(randperm(n_stims_per_cond),:);
    folds_unwrap = reshape(folds, [n_stims_per_cond*n_conds, 1]);
    
    % prediction
    Yh = regress_predictions_from_3way_crossval(X, Y, ...
        'test_folds', folds_unwrap, 'train_folds', folds_unwrap, ...
        'demean_feats', true, 'std_feats', true, 'method','ridge','K', 2.^(-100:100));
    
    % split out conditions
    target_prediction1 = reshape(Yh, [n_stims_per_cond, n_conds, n_target_components]);
    
    % split 2
    target_R2_tmp=splits{1,2};
    predictor_R_tmp=splits{2,setdiff(1:2,spl)};
    
    [n_stims_per_cond, n_conds, n_target_components, ~] = size(target_R2);
    [~, ~, n_predictor_components,~] = size(predictor_R_tmp);
    assert(size(target_R2_tmp,1)==size(predictor_R_tmp,1))
    assert(size(target_R2_tmp,2)==size(predictor_R_tmp,2))
    
    Y = reshape(target_R2_tmp, [n_stims_per_cond*n_conds, n_target_components]);
    X = reshape(predictor_R_tmp, [n_stims_per_cond*n_conds, n_predictor_components]);
    
    % prediction
    Yh = regress_predictions_from_3way_crossval(X, Y, ...
        'test_folds', folds_unwrap, 'train_folds', folds_unwrap, ...
        'demean_feats', true, 'std_feats', true, 'method', 'ridge','K', 2.^(-100:100));
    
    % split out conditions
    target_prediction2 = reshape(Yh, [n_stims_per_cond, n_conds, n_target_components]);
    
    for bs = 1:Nbs
        
        for i = 1 : n_target_components
            
            % natural
            target_variance(i,1,bs) = 0.25*(var(target_R1_n(bstraps(:,bs),end,i)+target_R2_n(bstraps(:,bs),end,i))-var(target_R1_n(bstraps(:,bs),end,i)-target_R2_n(bstraps(:,bs),end,i)));
            % difference natural -synthetic
            target_variance(i,2,bs) = 0.25*(var((target_R1_n(bstraps(:,bs),1,i)-target_R1_n(bstraps(:,bs),end,i))+(target_R2_n(bstraps(:,bs),1,i)-target_R2_n(bstraps(:,bs),end,i)))-var((target_R1_n(bstraps(:,bs),1,i)-target_R1_n(bstraps(:,bs),end,i))-(target_R2_n(bstraps(:,bs),1,i)-target_R2_n(bstraps(:,bs),end,i))));
            
        end
        % Noise corrected NSE
        nse_tmp = nan(n_target_components,2);
        for ic = 1:n_target_components
            
            if trt(ic,1) < 0.4
                % only mm=shared
                nse_tmp(ic,1) = NSE_noise_corrected(squeeze(target_prediction1(bstraps(:,bs),end,ic)),squeeze(target_prediction2(bstraps(:,bs),end,ic)),...
                    squeeze(target_R1_tmp(bstraps(:,bs),end,ic)),squeeze(target_R2_tmp(bstraps(:,bs),end,ic)),1);
                
            else
                
                nse_tmp(ic,1) = NSE(0.5*(squeeze(target_prediction1(bstraps(:,bs),end,ic))+squeeze(target_prediction2(bstraps(:,bs),end,ic))),...
                    0.5*(squeeze(target_R1_tmp(bstraps(:,bs),end,ic))+squeeze(target_R2_tmp(bstraps(:,bs),end,ic))),1);
            end
            
            
            if trt(ic,2) < 0.4
                nse_tmp(ic,2) = NSE_noise_corrected(squeeze(target_prediction1(bstraps(:,bs),1,ic)-target_prediction1(bstraps(:,bs),end,ic)),squeeze(target_prediction2(bstraps(:,bs),1,ic)-target_prediction2(bstraps(:,bs),end,ic)),...
                    squeeze(target_R1_tmp(bstraps(:,bs),1,ic)-target_R1_tmp(bstraps(:,bs),end,ic)),squeeze(target_R2_tmp(bstraps(:,bs),1,ic)-target_R2_tmp(bstraps(:,bs),end,ic)),1);
                
            else
                nse_tmp(ic,2) = NSE(0.5*(squeeze(target_prediction2(bstraps(:,bs),1,ic)-target_prediction2(bstraps(:,bs),end,ic))+squeeze(target_prediction1(bstraps(:,bs),1,ic)-target_prediction1(bstraps(:,bs),end,ic))),...
                    0.5*(squeeze(target_R1_tmp(bstraps(:,bs),1,ic)-target_R1_tmp(bstraps(:,bs),end,ic))+squeeze(target_R2_tmp(bstraps(:,bs),1,ic)-target_R2_tmp(bstraps(:,bs),end,ic))),1);
                
            end
        end
        nc_nse_all = sign(1-nse_tmp).*(1-nse_tmp).^2;
     
        predicted_variance(:,:,bs,spl) = nc_nse_all.*target_variance(:,:,bs);
        
    end
    
end

%% Plot summary of variance explained for synthetic and natural-synthetic
% parts of components' responses.

target_variance_n = target_variance./sum(mat2vec(snm(target_variance,3)));
predicted_variance_n = snm(predicted_variance,4)./sum(mat2vec(snm(target_variance,3)));

xpos = @(i,j) (i-1)*3+j;

error_target = cat(3,abs(prctile(target_variance_n,2.5,3)-nanmedian(target_variance_n,[3])),...
    prctile(target_variance_n,97.5,3)-nanmedian(target_variance_n,[3]));

error_predicted = cat(3,abs(prctile(predicted_variance_n,2.5,3)-nanmedian(predicted_variance_n,[3])),...
    prctile(predicted_variance_n,97.5,3)-nanmedian(predicted_variance_n,[3]));

figure('Position',[440 430 1079 368]);
subplot(1,3,1:2)

for i = 1 : n_target_components
    ic = CompOrder(i);
    hold all
    
    % Variance of synthetic
    b = bar(xpos(i,1),nanmedian(target_variance_n(ic,1,:),3));
    b.EdgeColor = 'k';
    b.FaceColor = [1 1 1];
    errorbar(xpos(i,1),nanmedian(target_variance_n(ic,1,:),[1 3]),error_target(ic,1,1),error_target(ic,1,2),'k')%,'horizontal')
    
    % Predicted variance
    c = bar(xpos(i,2),nanmedian(predicted_variance_n(ic,1,:),3));
    c.FaceColor = [1 1 1]*0.8;
    errorbar(xpos(i,2),nanmedian(predicted_variance_n(ic,1,:),3),error_predicted(ic,1,1),error_predicted(ic,1,2),'k')%,'horizontal')
    
    % Variance of difference natural - synthetic
    b = bar(xpos(i,3*(n_target_components+1)),nanmedian(target_variance_n(ic,2,:),3));
    b.EdgeColor = 'k';
    b.FaceColor = [1 1 1];
    errorbar(xpos(i,3*(n_target_components+1)),nanmedian(target_variance_n(ic,2,:),[1 3]),error_target(ic,2,1),error_target(ic,2,2),'k')%,'horizontal')
    
    % Predicted variance
    d = bar(xpos(i,3*(n_target_components+1)+1),nanmedian(predicted_variance_n(ic,2,:),3));
    d.FaceColor = [1 1 1]*0.3;
    errorbar(xpos(i,3*(n_target_components+1)+1),nanmedian(predicted_variance_n(ic,2,:),[1 3]),error_predicted(ic,2,1),error_predicted(ic,2,2),'k')
    
    
end
sgtitle(['Predict ' I.target])

xticks([xpos(1:n_target_components,1) xpos(1:n_target_components,3*(n_target_components+1))])
xticklabels(repmat(1:n_target_components,1,2))
xlabel([I.target ' components'])
ylabel('Variance (normalized)')


% same but averaged across components

target_variance_n = target_variance./sum(snm(target_variance,[1 3]),2);
predicted_variance_n = snm(predicted_variance,4)./sum(snm(target_variance,[1 3]),2);

xpos = @(i,j) (i-1)*3+j;

error_target = cat(2,abs(prctile(snm(target_variance_n,1),2.5,2)-nanmedian(snm(target_variance_n,1),2)),...
    abs(prctile(snm(target_variance_n,1),97.5,2)-nanmedian(snm(target_variance_n,1),2)));

error_predicted = cat(2,abs(prctile(snm(predicted_variance_n,1),2.5,2)-nanmedian(snm(predicted_variance_n,1),2)),...
    prctile(snm(predicted_variance_n,1),97.5,2)-nanmedian(snm(predicted_variance_n,1),2));


xpost = @(i,j) (i-1)*3+j;
subplot(1,3,3)
hold all
b = bar(1,nanmedian(nanmean(target_variance_n(:,1,:),1),[3]));
b.EdgeColor = 'k';
b.FaceColor = [1 1 1];
errorbar(1,nanmedian(nanmean(target_variance_n(:,1,:),1),[3]),error_target(1,1),error_target(1,2),'k')

c = bar(2,nanmedian(nanmean(predicted_variance_n(:,1,:),1),[ 3]));
c.FaceColor = [1 1 1]*0.8;
errorbar(2,nanmedian(nanmean(predicted_variance_n(:,1,:),1),[ 3]),error_predicted(1,1),error_predicted(1,2),'k')

b = bar(4,nanmedian(nanmean(target_variance_n(:,2,:),1),[3]));
b.EdgeColor = 'k';
b.FaceColor = [1 1 1];
errorbar(4,nanmedian(nanmean(target_variance_n(:,2,:),1),[3]),error_target(2,1),error_target(2,2),'k')

c = bar(5,nanmedian(nanmean(predicted_variance_n(:,2,:),1), 3));
c.FaceColor = [1 1 1]*0.3;
errorbar(5,nanmedian(nanmean(predicted_variance_n(:,2,:),1), 3),error_predicted(2,1),error_predicted(2,2),'k')

xticks([])
ylabel('Variance (normalized)')
xlabel('Average across components')
legend([b,c,d],'Variance of each part of IC','mm variance explained','orig-mm variance explained','Location','NorthOutside');

   