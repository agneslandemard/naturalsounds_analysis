function [mes_from_pac, mes_all, mult_sbj_error ] = NSE_disFromPAC(RespFull,SubjIdx,param,dist,sounds_list,analysis)
% RespFull should be in the format Stims * Voxels(=1 if no
% distinction) * 2 reps
% SubjIdx is a vector of lenght nb of voxels indicating which voxels belong
% to which subject/hemi
% dist is struct with fields "distances" and "distance_to_pac" (of
% size(Voxels)

if nargin < 6
    analysis = 'nse_across_sounds';
end
if isempty(SubjIdx)
    SubjIdx = ones(size(RespFull,2),1);
end

cond = param.snd.Nmodels; % will compare natural sounds to this model (spectrotemporal)
n_subj = max(SubjIdx);
n_sounds = length(sounds_list);

X_orig1 = RespFull(param.snd.idxSound(:,1),:,1);
X_orig2 = RespFull(param.snd.idxSound(:,1),:,2);
X_mm1 = RespFull(param.snd.idxSound(:,cond),:,1);
X_mm2 = RespFull(param.snd.idxSound(:,cond),:,2);

switch analysis
    
    case 'diff_resp'
        mes_all = snm(RespFull(param.snd.idxSound(:,1),:,:)-RespFull(param.snd.idxSound(:,cond),:,:),3)./nanstd(RespFull(param.snd.idxSound(:,[1 cond]),:,:),[],[1 3]);
        
    case 'av_resp'
        mes_all = snm(RespFull(param.snd.idxSound(:,1),:,:),3);
        
    case 'rel'
        mes_all = RespFull(param.snd.idxSound(:,1),:,1)-RespFull(param.snd.idxSound(:,1),:,2);
        
    case 'nse_across_sounds'
        dim = 1;
        
        a1 = 0.5*(X_orig1.^2+X_orig2.^2-(X_orig1-X_orig2).^2);
        b1 = 0.5*(X_mm1.^2+X_mm2.^2-(X_mm1-X_mm2).^2);
        c = 0.25*(X_orig1.*X_mm1+X_orig1.*X_mm2+X_orig2.*X_mm1+X_orig2.*X_mm2);
        
        a2 = 0.5*(nanmean(X_orig1.^2,dim)+nanmean(X_orig2.^2,dim)-nanmean((X_orig1-X_orig2).^2,dim));
        b2 = 0.5*(nanmean(X_mm1.^2,dim)+nanmean(X_mm2.^2,dim)-nanmean((X_mm1-X_mm2).^2,dim));
        d = 0.5*(nanmean(X_orig1,dim)+nanmean(X_orig2,dim));
        e = 0.5*(nanmean(X_mm1,dim)+nanmean(X_mm2,dim));
        
        mes_all = (a1+b1-2*c)./(a2+b2-2*d.*e);
        
    case 'nse_across_sounds_adapted'
        dim = 1;
        a1 = 0.5*(X_orig1.^2+X_orig2.^2-(X_orig1-X_orig2).^2);
        b1 = 0.5*(X_mm1.^2+X_mm2.^2-(X_mm1-X_mm2).^2);
        c = 0.25*(X_orig1.*X_mm1+X_orig1.*X_mm2+X_orig2.*X_mm1+X_orig2.*X_mm2);
        
        a2 = 0.5*(nanmean(X_orig1.^2,dim)+nanmean(X_orig2.^2,dim)-nanmean((X_orig1-X_orig2).^2,dim));
        b2 = 0.5*(nanmean(X_mm1.^2,dim)+nanmean(X_mm2.^2,dim)-nanmean((X_mm1-X_mm2).^2,dim));
        d = 0.5*(nanmean(X_orig1,dim)+nanmean(X_orig2,dim));
        e = 0.5*(nanmean(X_mm1,dim)+nanmean(X_mm2,dim));
        
        mes_all = (a1+b1-2*c)./(a2+b2-2*d.*e);
        
        % when only one rep, use only natural sounds to estimate the noise
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
        

    case 'nse_across_sounds_bycat'
        mes_all = nan(size(X_orig1));
        dim = 1;
        
        for sd = 1 : n_sounds
            SdsIdx = SelectSounds(sounds_list{sd}, param);
            
            a1 = 0.5*(X_orig1(SdsIdx,:).^2+X_orig2(SdsIdx,:).^2-(X_orig1(SdsIdx,:)-X_orig2(SdsIdx,:)).^2);
            b1 = 0.5*(X_mm1(SdsIdx,:).^2+X_mm2(SdsIdx,:).^2-(X_mm1(SdsIdx,:)-X_mm2(SdsIdx,:)).^2);
            c = 0.25*(X_orig1(SdsIdx,:).*X_mm1(SdsIdx,:)+X_orig1(SdsIdx,:).*X_mm2(SdsIdx,:)+X_orig2(SdsIdx,:).*X_mm1(SdsIdx,:)+X_orig2(SdsIdx,:).*X_mm2(SdsIdx,:));
            
            a2 = 0.5*(nanmean(X_orig1(SdsIdx,:).^2,dim)+nanmean(X_orig2(SdsIdx,:).^2,dim)-nanmean((X_orig1(SdsIdx,:)-X_orig2(SdsIdx,:)).^2,dim));
            b2 = 0.5*(nanmean(X_mm1(SdsIdx,:).^2,dim)+nanmean(X_mm2(SdsIdx,:).^2,dim)-nanmean((X_mm1(SdsIdx,:)-X_mm2(SdsIdx,:)).^2,dim));
            d = 0.5*(nanmean(X_orig1(SdsIdx,:),dim)+nanmean(X_orig2(SdsIdx,:),dim));
            e = 0.5*(nanmean(X_mm1(SdsIdx,:),dim)+nanmean(X_mm2(SdsIdx,:),dim));
            
            mes_all(SdsIdx,:) = (a1+b1-2*c)./(a2+b2-2*d.*e);
        end
end


if ~isempty(dist)
    n_bins = length(dist.distances);
    mes_from_pac = nan(n_sounds,n_bins,n_subj,2);
    mes_save = nan(size(mes_all,1),n_bins,n_subj);
    mult_sbj_error = nan(n_sounds,n_bins);
    
   
    for sd = 1 : n_sounds
        SdsIdx = SelectSounds(sounds_list{sd}, param);
        for sb = 1: n_subj
            xi = SubjIdx == sb;
            mes_subj = mes_all(:,xi);
            
            for d = 1 : n_bins
                if ~isempty(find(dist.distance_to_pac(xi) == dist.distances(d),1))
                    % first median across voxels and then mean across sounds
                    Resp = mes_subj(SdsIdx,dist.distance_to_pac(xi) == dist.distances(d));
                    mes_from_pac(sd,d,sb,1) = squeeze(nanmean(nanmedian(Resp,2),1));
                    mes_from_pac(sd,d,sb,2) = squeeze(nanstd(nanmedian(Resp,2),[],1)./sqrt(length(find(SdsIdx))));
                    mes_save(SdsIdx,d,sb) = nanmedian(Resp,2);
                    
                end
            end
        end
         mult_sbj_error(sd,:) = nanstd(nanmean(mes_save(SdsIdx,:,:),3),1)./sqrt(length(find(SdsIdx)));
    end
   
else
    mes_from_pac = [];
end

end