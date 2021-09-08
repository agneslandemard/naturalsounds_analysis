function [corr_pvals_indiv,corr_pvals_bigCats,stats_cats] = PlotNSEBoxplots(D,sbj,SoundType,method)
% Boxplot display and quantification across sound categories and subjects

signifWhat = 'bigCats' ;% bigCats just merges speech and music together for stats
cond = D.param.snd.Nmodels; % will look at difference between natural and this model

% define colors for categories
colors.Ferrets = [0.949    0.604       0.72];
colors.Speech = [0.1882    0.7490    0.6196];
colors.Music = [0.0549    0.3020    0.5843 ];
try
    colors.Others = snm(D.param.plt.SoundColors(SelectSounds('Others',D.param),:),1);
catch
end
    
f1 = figure('Position', [440 279 length(sbj.Subjects)*160 489]);

ratioall = cell(length(sbj.Subjects),1);
ratioall_allsbj = [];
ratioctrl = cell(length(sbj.Subjects),1);
xpos = @(h,s) s+(h-1)*(length(SoundType)+1);

m = 0;
pvals_allCats=nan(length(sbj.Subjects)+1,length(SoundType),length(SoundType));
pvals_indiv=nan(length(sbj.Subjects)+1,length(SoundType));
pvals_bigCats=nan(length(sbj.Subjects)+1,2);

sbj.Subject_list = [sbj.Subject_list 'GroupMean'];

for f = 1 : length(sbj.Subjects)
    Resp = D.(sbj.Subjects{f});

    AvResp = snm(Resp(D.param.snd.idxSound(:,[1 cond]),:,:),[1 3]);
    
    X_orig1 = snm(Resp(D.param.snd.idxSound(:,1),:,1),3);
    X_mm1 = snm(Resp(D.param.snd.idxSound(:,cond),:,1),3);
    
    X_orig2 = snm(Resp(D.param.snd.idxSound(:,1),:,2),3);
    X_mm2 = snm(Resp(D.param.snd.idxSound(:,cond),:,2),3);
  
    switch method
        case 'ncNSE_bycat'
            % for nse by category, voxel-wise, noise-corrected
            dim = 1;
            a1 = 0.5*(X_orig1.^2+X_orig2.^2-(X_orig1-X_orig2).^2);
            b1 = 0.5*(X_mm1.^2+X_mm2.^2-(X_mm1-X_mm2).^2);
            c = 0.25*(X_orig1.*X_mm1+X_orig1.*X_mm2+X_orig2.*X_mm1+X_orig2.*X_mm2);
            
            a2 = 0.5*(nanmean(X_orig1.^2,dim)+nanmean(X_orig2.^2,dim)-nanmean((X_orig1-X_orig2).^2,dim));
            b2 = 0.5*(nanmean(X_mm1.^2,dim)+nanmean(X_mm2.^2,dim)-nanmean((X_mm1-X_mm2).^2,dim));
            d = 0.5*(nanmean(X_orig1,dim)+nanmean(X_orig2,dim));
            e = 0.5*(nanmean(X_mm1,dim)+nanmean(X_mm2,dim));
            
            ratioall{f} = nanmedian((a1+b1-2*c)./(a2+b2-2*d.*e),2);
           
        case 'ncNSE_bycat_adapted'

            dim = 1;
            a1 = 0.5*(X_orig1.^2+X_orig2.^2-(X_orig1-X_orig2).^2);
            b1 = 0.5*(X_mm1.^2+X_mm2.^2-(X_mm1-X_mm2).^2);
            c = 0.25*(X_orig1.*X_mm1+X_orig1.*X_mm2+X_orig2.*X_mm1+X_orig2.*X_mm2);
            
            a2 = 0.5*(nanmean(X_orig1.^2,dim)+nanmean(X_orig2.^2,dim)-nanmean((X_orig1-X_orig2).^2,dim));
            b2 = 0.5*(nanmean(X_mm1.^2,dim)+nanmean(X_mm2.^2,dim)-nanmean((X_mm1-X_mm2).^2,dim));
            d = 0.5*(nanmean(X_orig1,dim)+nanmean(X_orig2,dim));
            e = 0.5*(nanmean(X_mm1,dim)+nanmean(X_mm2,dim));
            
            mes_all = (a1+b1-2*c)./(a2+b2-2*d.*e);
            
            only1rep = isnan(nanmean(X_mm2,1));
            
            a1 = 0.5*(X_orig1.^2+X_orig2.^2-(X_orig1-X_orig2).^2);
            b1 = 0.5*(2*X_mm1.^2-(X_orig1-X_orig2).^2);
            c = 0.5*(X_orig1.*X_mm1+X_orig2.*X_mm1);
            
            a2 = 0.5*(nanmean(X_orig1.^2,dim)+nanmean(X_orig2.^2,dim)-nanmean((X_orig1-X_orig2).^2,dim));
            b2 = 0.5*(2*nanmean(X_mm1.^2,dim)-nanmean((X_orig1-X_orig2).^2,dim));
            d = 0.5*(nanmean(X_orig1,dim)+nanmean(X_orig2,dim));
            e = nanmean(X_mm1,dim);
            
            tmp = (a1+b1-2*c)./(a2+b2-2*d.*e);
            mes_all(:,only1rep) = tmp(:,only1rep);
         
            ratioall{f} = nanmedian(mes_all,2);
            
        case 'RawDiff'
            ratioall{f} = 0.5*(nanmean((X_orig1-X_mm2).^2,2)+nanmean((X_orig2-X_mm1).^2,2));
            ratioctrl{f} = 0.5*(nanmean((X_orig1-X_orig2).^2,2)+nanmean((X_mm1-X_mm2).^2,2));
            
        case 'mNSE'
            ratioall{f} = 0.5*(NSE_bis(X_orig1,X_mm2,2,1)+NSE_bis(X_orig2,X_mm1,2,1));
            ratioctrl{f} = 0.5*(NSE_bis(X_orig1,X_orig2,2,1)+NSE_bis(X_mm2,X_mm1,2,1));
           
        case 'mncNSE'
            X_mm1 = X_mm1-AvResp;
            X_mm2 = X_mm2-AvResp;
            X_orig1 = X_orig1-AvResp;
            X_orig2 = X_orig2-AvResp;
            
            if contains(sbj.Subjects{f},'allSubjects')
                [n_stims,n_vx] = size(X_orig1);
                X_orig1bySubj = reshape(X_orig1,n_stims,n_vx/n_subj,n_subj);
                X_orig2bySubj = reshape(X_orig2,n_stims,n_vx/n_subj,n_subj);
                
                X_mm1bySubj = reshape(X_mm1,n_stims,n_vx/n_subj,n_subj);
                
                X_mm2bySubj = reshape(X_mm2,n_stims,n_vx/n_subj,n_subj);
                
                ratiobysubj = nan(n_stims,n_subj);
                for sb = 1 :n_subj
                    ratiobysubj(:,sb) = NSE_bis_noise_corrected(X_orig1bySubj(:,:,sb),X_orig2bySubj(:,:,sb),X_mm1bySubj(:,:,sb),X_mm2bySubj(:,:,sb),2,1);
                end
                ratioall{f} = nanmean(ratiobysubj,2);
                
            else
                ratioall{f} = NSE_bis_noise_corrected(X_orig1,X_orig2,X_mm1,X_mm2,2,1);
                
            end
            
        case 'ncNSE'
            ratioall{f} = NSE_noise_corrected(X_orig1,X_orig2,X_mm1,X_mm2,2);
        
        case 'wTestRetest'
            denom = 0.5*(nanmean(sqrt(nansum((X_mm1 - nanmean(X_mm1,1)).^2,2)),1) + nanmean(sqrt(nansum((X_orig1 - nanmean(X_orig1,1)).^2,2)),1)+...
                nanmean(sqrt(nansum((X_mm2 - nanmean(X_mm2,1)).^2,2)),1) + nanmean(sqrt(nansum((X_orig2 - nanmean(X_orig2,1)).^2,2)),1));
            
            numall = 0.5*(sqrt(nansum((X_mm1 - X_orig2).^2,2))+sqrt(nansum((X_mm2 - X_orig1).^2,2)));
            testretest = [sqrt(nansum((X_mm1 - X_mm2).^2,2)); sqrt(nansum((X_orig1 - X_orig2).^2,2))];
            
            
            ratioall{f} = numall/denom;
            ratioctrl{f} = testretest/denom;
    end
    
    ratioall_allsbj = cat(2,ratioall_allsbj,ratioall{f});
    
    
    figure(f1)
    hold all
    for sd = 1:length(SoundType)
        if any(~isnan(ratioall{f}(SelectSounds(SoundType{sd},D.param))))
            
            if ~isempty(ratioctrl{f})
                bplot(ratioctrl{f}(SelectSounds(SoundType{sd},D.param)),xpos(f,sd),'color',0.8*[1 1 1])
            end
            
            bplot(ratioall{f}(SelectSounds(SoundType{sd},D.param)),xpos(f,sd),'color',colors.(SoundType{sd}));
            hold on
            xpos1 = repmat(xpos(f,sd),length(SelectSounds(SoundType{sd},D.param)),1)+(rand(length(SelectSounds(SoundType{sd},D.param)),1)-0.5)*0.3;
    
            scatter(xpos1,ratioall{f}(SelectSounds(SoundType{sd},D.param)),[],colors.(SoundType{sd}),'filled','MarkerEdgeAlpha',0.4,'MarkerFaceAlpha',0.4)
                 
            for sd2 = 1:sd-1
                if any(~isnan(ratioall{f}(SelectSounds(SoundType{sd2},D.param))))
                    [pvals_allCats(f,sd,sd2)] = ranksum(ratioall{f}(SelectSounds(SoundType{sd},D.param)),ratioall{f}(SelectSounds(SoundType{sd2},D.param)));
                    
                end
            end
            pvals_indiv(f,sd) = signrank(ratioall{f}(SelectSounds(SoundType{sd},D.param)));
        end
        
      
    end
     if any(~isnan(ratioall{f}(SelectSounds('Ferrets',D.param)))) && any(~isnan(ratioall{f}(SelectSounds('SpeechMusic',D.param))))
         [pvals_bigCats(f,1),~,stats] = ranksum(ratioall{f}(SelectSounds('Ferrets',D.param)),ratioall{f}(SelectSounds('SpeechMusic',D.param)));
         stats_cats.(sbj.Subject_list{f}).medians(1)=nanmedian(ratioall{f}(SelectSounds('Ferrets',D.param)));
         stats_cats.(sbj.Subject_list{f}).medians(2)=nanmedian(ratioall{f}(SelectSounds('SpeechMusic',D.param)));
         stats_cats.(sbj.Subject_list{f}).zval=stats.zval;
         stats_cats.(sbj.Subject_list{f}).signedrank=stats.ranksum;
                
     end
    
     try
        pvals_bigCats(f,2) = ranksum(ratioall{f}(SelectSounds('SpeechMusic',D.param)),ratioall{f}(SelectSounds('Others',D.param)));
     catch
     end
    m = max(m,max(ratioall{f}));
end

ff = length(sbj.Subjects)+1;

% also plot average across subjects
ratiomean = nanmean(ratioall_allsbj,2);

if any(~isnan(ratiomean(SelectSounds('Ferrets',D.param)))) && any(~isnan(ratiomean(SelectSounds('SpeechMusic',D.param))))
    [pvals_bigCats(ff,1),~,stats] = ranksum(ratiomean(SelectSounds('Ferrets',D.param)),ratiomean(SelectSounds('SpeechMusic',D.param)));
    stats_cats.(sbj.Subject_list{ff}).medians(1)=nanmedian(ratiomean(SelectSounds('Ferrets',D.param)));
    stats_cats.(sbj.Subject_list{ff}).medians(2)=nanmedian(ratiomean(SelectSounds('SpeechMusic',D.param)));
    stats_cats.(sbj.Subject_list{ff}).zval=stats.zval;
    stats_cats.(sbj.Subject_list{ff}).signedrank=stats.ranksum;
end

try
    pvals_bigCats(ff,2) = ranksum(ratiomean(SelectSounds('SpeechMusic',D.param)),ratiomean(SelectSounds('Others',D.param)));
catch
end


for sd = 1:length(SoundType)
    if any(~isnan(ratiomean(SelectSounds(SoundType{sd},D.param))))
        
        
        bplot(ratiomean(SelectSounds(SoundType{sd},D.param)),xpos(ff,sd),'color',colors.(SoundType{sd}));
        hold on
        xpos1 = repmat(xpos(ff,sd),length(SelectSounds(SoundType{sd},D.param)),1)+(rand(length(SelectSounds(SoundType{sd},D.param)),1)-0.5)*0.3;
        
        scatter(xpos1,ratiomean(SelectSounds(SoundType{sd},D.param)),[],colors.(SoundType{sd}),'filled','MarkerEdgeAlpha',0.4,'MarkerFaceAlpha',0.4)
        
        for sd2 = 1:sd-1
            if any(~isnan(ratiomean(SelectSounds(SoundType{sd2},D.param))))
                pvals_allCats(ff,sd,sd2) = ranksum(ratiomean(SelectSounds(SoundType{sd},D.param)),ratiomean(SelectSounds(SoundType{sd2},D.param)));
            end
        end
        [pvals_indiv(ff,sd)] = signrank(ratiomean(SelectSounds(SoundType{sd},D.param)));
      
    end
    
    
end
if any(~isnan(ratiomean(SelectSounds('Ferrets',D.param)))) && any(~isnan(ratiomean(SelectSounds('SpeechMusic',D.param))))
    [pvals_bigCats(ff,1)] = ranksum(ratiomean(SelectSounds('Ferrets',D.param)),ratiomean(SelectSounds('SpeechMusic',D.param)));
  
end

try
    pvals_bigCats(ff,2) = ranksum(ratiomean(SelectSounds('SpeechMusic',D.param)),ratiomean(SelectSounds('Others',D.param)));
catch
end

m = max(m,max(ratiomean));


% Correct for multiple comparison
corr_pvals_bigCats = nan(size(pvals_bigCats));
[~,~,corr_pvals_bigCats(~isnan(pvals_bigCats(:)))] = fdr_bh(pvals_bigCats(~isnan(pvals_bigCats(:))));
corr_pvals_bigCats = reshape(corr_pvals_bigCats,size(pvals_bigCats));

corr_pvals_allCats = nan(size(pvals_allCats));
[~,~,corr_pvals_allCats(~isnan(pvals_allCats(:)))] = fdr_bh(pvals_allCats(~isnan(pvals_allCats(:))));
corr_pvals_allCats = reshape(corr_pvals_allCats,size(pvals_allCats));


corr_pvals_indiv = nan(size(pvals_indiv));
[~,~,corr_pvals_indiv(~isnan(pvals_indiv(:)))] = fdr_bh(pvals_indiv(~isnan(pvals_indiv(:))));
corr_pvals_indiv = reshape(corr_pvals_indiv,size(pvals_indiv));

for f = 1:length(sbj.Subjects)+1
    
    stats_cats.(sbj.Subject_list{f}).pval = squeeze(corr_pvals_bigCats(f,:,:));
    
    switch signifWhat
        case 'bigCats'
            SpMIdx = 0.5*(find(strcmp(SoundType,'Speech'))+find(strcmp(SoundType,'Music')));
            if corr_pvals_bigCats(f,1)<0.05
                sigstar({[xpos(f,find(strcmp(SoundType,'Ferrets'))),xpos(f,SpMIdx)]},corr_pvals_bigCats(f,1));
            end
            if corr_pvals_bigCats(f,2)<0.05
                sigstar({[xpos(f,SpMIdx),xpos(f,find(strcmp(SoundType,'Others')))]},corr_pvals_bigCats(f,2));
            end
        case 'allCats'
            for sd = 1:length(SoundType)
                for sd2 = 1:sd-1
                    if corr_pvals(f,sd,sd2)<0.05
                        sigstar({[xpos(f,sd),xpos(f,sd2)]},corr_pvals_allCats(f,sd,sd2));
                    end
                end
            end
    end
end



xlim([0 xpos(f,length(SoundType))+2]);
if contains(method,'nc')
    hline(0,'k-')
    ylim([-0.1 1.2*m])
  
end 
yticks(0:0.5:1)

name = regexprep(method,'_',' ');

title(name)
xticks(xpos(1:length(sbj.Subjects)+1,2))
xticklabels(sbj.Subject_list)
ylabel('Distance natural - synthetic')



