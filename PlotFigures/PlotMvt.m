%% Plot Mvt for different categories of sounds

experiment = 'vocalization';

% Choose feat = sounds or sessions, will be the points displayed and the
% statistics will be over it too
feature = 'sounds';
display_big_categories = 1;

% Break sounds into major categories
switch experiment
    case 'natural'
        
        animals = {'A','T','C'};
        grouped_cats={'Ferrets','Speech','Music','Others'};
        
    case 'vocalization'
        
        animals = {'A','T'};        
        if display_big_categories
            grouped_cats = {'Ferrets','Speech','Music'};
        else
            grouped_cats = {'Fights','Kit','Kits','Speech','Music'};
        end
end
n_categories = length(grouped_cats);
    
    
figure('Position', [257 619 779 245]);
clear ax
for an = 1:length(animals)
    
    % load movement for this animal
    load([additional_path '/Movement/Movement_' upper(experiment(1))  animals{an} '.mat'],'RawMvt','param');
    
    % Remove empty sessions (with no video data)
    EmptySessions = find(isnan(snm(RawMvt,[1 2 3])));
    RawMvt(:,:,:,EmptySessions) = [];
    
    % Normalize by std across all sounds
    NormMvt = squeeze(nanmean(RawMvt,3)./std(nanmean(RawMvt,3),[],[1 2]));
    
    switch feature
        case 'sounds'
            dim = 3;
        case 'sessions'
            dim = 1;
    end
    
    p = nan(n_categories,1);
    groups = cell(n_categories,1);
    xpos = @(cat) 1.5+3*(cat-1);
    clear stats;
    ax(an) = subplot(1,length(animals),an);
    hold all
    for cat = 1:n_categories
        sds = SelectSounds(grouped_cats{cat},param);
        clr = mean(param.plt.SoundColors(sds,:),1);
        
        % display single points
        switch feature
            case 'sounds'
                Npoints = length(sds);
            case 'sessions'
                Npoints = size(NormMvt,3);
        end
        jitterw = 0.2;
        xpos1 = repmat(xpos(cat)-0.5,Npoints,1)+(rand(Npoints,1)-0.5)*jitterw;
        xpos2 = repmat(xpos(cat)+0.5,Npoints,1)+(rand(Npoints,1)-0.5)*jitterw;
        line([xpos1 xpos2]',[squeeze(nanmean(NormMvt(sds,1,:),dim)) squeeze(nanmean(NormMvt(sds,2,:),dim))]','Color',0.8*[1 1 1],'LineWidth',1)
        scatter(xpos1,squeeze(nanmean(NormMvt(sds,1,:),dim)),[],clr,'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
        scatter(xpos2,squeeze(nanmean(NormMvt(sds,2,:),dim)),[],clr,'MarkerEdgeAlpha',0.6);
        
        % Plot distribution of evoked movement for natural and synthetic sounds
        bplot(mat2vec(nanmean(NormMvt(sds,1,:),dim)),xpos(cat)-0.5,'color',0.5*clr,'linewidth',1);
        bplot(mat2vec(nanmean(NormMvt(sds,2,:),dim)),xpos(cat)+0.5,'color',0.5*clr,'linewidth',1);
        
        % stats nat vs synth
        [p(cat),~,stats.(grouped_cats{cat})] = signrank(mat2vec(nanmean(NormMvt(sds,1,:),dim)), mat2vec(nanmean(NormMvt(sds,2,:),dim)));
        groups{cat} = {xpos(cat)-0.5, xpos(cat)+0.5};
       
    end
    ylabel('Magnitude of Movement (A.U)')
    xticks(xpos(1:n_categories));
    xticklabels(grouped_cats)
    title(animals{an})
    
    % Stats between natural sounds of different categories
    combs = nchoosek(1:n_categories,2);
    p2 = nan(nchoosek(n_categories,2),1);
    groups2 = cell(nchoosek(n_categories,2),1);
    for cb = 1 : nchoosek(n_categories,2)
        sds1 = SelectSounds(grouped_cats{combs(cb,1)},param);
        sds2 = SelectSounds(grouped_cats{combs(cb,2)},param);
        [p2(cb),~,stats] = ranksum(mat2vec(nanmean(NormMvt(sds1,1,:),dim)), mat2vec(nanmean(NormMvt(sds2,1,:),dim)));
        groups2{cb} = xpos(combs(cb,:))-0.5;
    end
    
    % Correct for multiple comparison
    groups_all = [groups; groups2];
    [~,~,corrp_all] = fdr_bh([p; p2]);
    ns = find(corrp_all > 0.05);
    corrp_all(ns) = [];
    groups_all(ns) = [];
    
    % plot
    sigstar(groups_all, corrp_all)
end
linkaxes(ax)