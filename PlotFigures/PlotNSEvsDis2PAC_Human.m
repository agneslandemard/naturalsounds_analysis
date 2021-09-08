%%% Plots all NSE vs distance to PAC curves

% - for all sounds, without noise-correction, with test-retest displayed
% - for all sounds, with noise-correction
% - by category, noise-corrected
% - NSE for each sound to quantify differences across categories

%% NSE nat vs synth + NSE test-retest as function of distance to PAC

h = load([additional_path 'Coordinates/distances_to_pac_human.mat']);

n_conds = D.param.snd.Nmodels;
HD = reshape(SrecoFull_human_allSubjects,200,xl*yl+xr*yr,n_subj,2);
local_nse_h = nan(n_conds,length(h.distance_humans),n_subj);
local_nse_trt_h = nan(n_conds,length(h.distance_humans),n_subj);

% select subjects for each model remove those without tonotopy
stk = cell(n_conds, 1);
for md = 2 : n_conds - 1
    stk{md} = setdiff(1:17,[4 9 10 12]); % subjects to keep for submodels
end
stk{n_conds} = setdiff(2:17,[4 9 10 12]); % subjects to keep for full model (need to remove a subject present in both paradigm)

figure('Position',[680 829 606 149]);
for md = 2 : n_conds
  
    for s = 1:n_subj
            if s <= 12 % subjects from paradigm I, estimate noise using natural sounds only

                nse_single_subj= 0.5*(NSE(HD(D.param.snd.idxSound(:,1),:,s,1),HD(D.param.snd.idxSound(:,md),:,s,1),1)+...
                    NSE(HD(D.param.snd.idxSound(:,1),:,s,2),HD(D.param.snd.idxSound(:,md),:,s,1),1));
                
                nse_trt_single_subj = NSE(HD(D.param.snd.idxSound(:,1),:,s,1),HD(D.param.snd.idxSound(:,1),:,s,2),1);
                
                if all(isnan(nse_single_subj))
                    % if no data for this model, don't show test-retest
                    % curve for this subject
                    nse_trt_single_subj = nan(size(nse_trt_single_subj));
                end
                
            else
                nse_single_subj= 0.5*(NSE(HD(D.param.snd.idxSound(:,1),:,s,1),HD(D.param.snd.idxSound(:,md),:,s,2),1)+...
                    NSE(HD(D.param.snd.idxSound(:,1),:,s,2),HD(D.param.snd.idxSound(:,md),:,s,1),1));
                
                nse_trt_single_subj = NSE(HD(D.param.snd.idxSound(:,1),:,s,1),HD(D.param.snd.idxSound(:,1),:,s,2),1);
                
                if all(isnan(nse_single_subj))
                      % if no data for this model, don't show test-retest
                    % curve for this subject
                    nse_trt_single_subj = nan(size(nse_trt_single_subj));
                end
            end
            
            % remove unreliable voxels
            NRV = nse_trt_single_subj > 0.4;
            nse_trt_single_subj(NRV) = nan;
            nse_single_subj(NRV) = nan;
            
            for d = 1 : length(h.distance_humans)
                local_nse_h(md,d,s) = nanmedian(nse_single_subj(h.distance_to_pac_human(:,s) == h.distance_humans(d)));
                local_nse_trt_h(md,d,s) = nanmedian(nse_trt_single_subj(h.distance_to_pac_human(:,s) == h.distance_humans(d)));
            end
    end
    
    subplot(1,n_conds-1,md-1)
    hold all
    plot(h.distance_humans(1:end-2),squeeze(local_nse_h(md,1:end-2,stk{md})),'color',[232 191 88]./255,'LineWidth',0.5)
    plot(h.distance_humans(1:end-2),snm(local_nse_h(md,1:end-2,stk{md}),3),'color',[173 133 43]./255,'LineWidth',1)
    
    plot(h.distance_humans(1:end-2),squeeze(local_nse_trt_h(md,1:end-2,stk{md})),'--','color',[232 191 88]./255,'LineWidth',0.5)
    plot(h.distance_humans(1:end-2),snm(local_nse_trt_h(md,1:end-2,stk{md}),3),'--','color',[173 133 43]./255,'LineWidth',1)
    

    ylim([-0.05 1.1])
    yticks(0:0.2:1)
    xticks(h.distance_humans(1:2:end))
    xlim([0 h.distance_humans(end-2)])
    
     if md == 2
        xlabel('Distance from center of PAC (mm)')
        ylabel('NSE')
     end
end


%% With noise-correction

n_conds = D.param.snd.Nmodels;
n_vxs = xl*yl+xr*yr;
HD = reshape(SrecoFull_human_allSubjects,200,n_vxs,n_subj,2);

% select subjects for each model remove those without tonotopy
stk = cell(n_conds, 1);
for md = 2 : n_conds - 1
    stk{md} = setdiff(1:17,[4 9 10 12]); % subjects to keep for submodels
end
stk{n_conds} = setdiff(2:17,[4 9 10 12]); % subjects to keep for full model (need to remove a subject present in both paradigm)

h = load([additional_path 'Coordinates/distances_to_pac_human.mat']);

local_nse_h = nan(n_conds,length(h.distance_humans),n_subj);
for md = 2 : n_conds
    
    for s = 1:n_subj
        or = HD(D.param.snd.idxSound(:,1),:,s,:);
        mm = HD(D.param.snd.idxSound(:,md),:,s,:);
        
        if s <= 12
            % for subjects from paradigm I, noise-correct using natural
            % sounds only (not enough repetition for synthetic sounds)
            disp('using same noise')
            nse_single_subj = NSE_noise_corrected_adapted(or(:,:,:,1),or(:,:,:,2),mm(:,:,:,1),1);
            
        else
            % for subjects from paradigm II, noise-correct using natural
            % and synthetic sounds
            disp('using different noise')
            nse_single_subj = NSE_noise_corrected(or(:,:,:,1),or(:,:,:,2),mm(:,:,:,1),mm(:,:,:,2),1);
            
        end
        
        % remove unreliable voxels
        nse_trt_single_subj = NSE(or(:,:,:,1),or(:,:,:,2),1);
        NRV = nse_trt_single_subj > 0.4;
        nse_single_subj(NRV)=nan;
        
        for d = 1 : length(h.distance_humans)
            local_nse_h(md,d,s) = nanmedian(nse_single_subj(h.distance_to_pac_human(:,s) == h.distance_humans(d)));
        end
        
    end
end

figure('Position',[680 829 606 149]);
stop = 2;
for md = 2 : n_conds
    axt(md) = subplot(1,n_conds-1,md-1);
    hold all
    plot(h.distance_humans(1:end-stop),squeeze(local_nse_h(md,1:end-stop,stk{md})),'color',[232 191 88]./255)
    plot(h.distance_humans(1:end-stop),snm(local_nse_h(md,1:end-stop,stk{md}),3),'color',[173 133 43]./255,'LineWidth',2)
    
    ylim([0 1])
    yticks(0:0.2:1)
    xticks(h.distance_humans(1:2:end))
     if md == 2
        xlabel('Distance from center of PAC (mm)')
        ylabel('NSE (noise-corrected)')
     end

end
linkaxes(axt)

%% NSE by category as function of distance to PAC

SoundType = {'Speech','Music','Others'};
plot_type = 'bounded';

dist = load([additional_path 'Coordinates/distances_to_pac_human.mat']);
dist.distance_to_pac = dist.distance_to_pac_human(:);
dist.distances = dist.distance_humans(1:end-2);

colors = [0.1882    0.7490    0.6196 ;...
    0.0549    0.3020    0.5843 ;...
    snm(D.param.plt.SoundColors(SelectSounds('Others',D.param),:),1)];

SubjIdx = mat2vec(repmat(1:n_subj,[xl*yl+xr*yr 1]));

% keep only reliable voxels
HD = reshape(SrecoFull_human_allSubjects,200,xl*yl+xr*yr,n_subj,2);
nse_trt = NSE(HD(D.param.snd.idxSound(:,1),:,:,1),HD(D.param.snd.idxSound(:,1),:,:,2),1);
NRV = nse_trt > 0.4;
tmpS = SrecoFull_human_allSubjects;
tmpS(:,NRV(:),:)=nan;

% Compute NSE
[nse_by_dis,~,err] = NSE_disFromPAC(tmpS,SubjIdx,D.param,dist,SoundType,'nse_across_sounds_adapted');

figure; 
hold all
stk = setdiff(1:17,[4 9 10 12]); % remove subjects without tonotopy
for sd = 1:length(SoundType)
    switch plot_type
        case 'normal'
            plot(dist.distances+2.5,snm(nse_by_dis(sd,:,stk ,1),3),'color',colors(sd,:),'LineWIdth',2);
            plot(dist.distances+2.5,squeeze(nse_by_dis(sd,:,stk ,1)),'color',colors(sd,:),'LineWIdth',0.3);
            
        case 'bounded'
            boundedline(dist.distances+2.5,snm(nse_by_dis(sd,:,stk ,1),3),err(sd,:),'cmap',colors(sd,:),'alpha');
    end
end
xlabel('Distance from center of PAC (mm)')
title('Human')
yticks(0:0.5:1)
ylim([-0.1 1.4])
xticks(dist.distances(1:2:end))

%% Quantify NSE across sounds of different categories

SoundType = {'Speech','Music','Others'};
method = 'ncNSE_bycat_adapted';

nse_trt = NSE(SrecoFull_human_allSubjects(D.param.snd.idxSound(:,1),:,1),SrecoFull_human_allSubjects(D.param.snd.idxSound(:,1),:,2),1);
NRV = nse_trt > 0.4;
SrecoFull_human_allSubjects_rmNRV = SrecoFull_human_allSubjects;
SrecoFull_human_allSubjects_rmNRV(:,NRV,:)=nan;
SrecoFull_human_allSubjects_rmNRV = reshape(SrecoFull_human_allSubjects_rmNRV,200,xl*yl+xr*yr,n_subj,2);

clear sbj
sbj.Subjects = {};
sbj.Subject_list = {};
stk = setdiff(1:17,[4 9 10 12]);
for sb = 1 : length(stk)
    D.(['h' num2str(sb)]) = squeeze(SrecoFull_human_allSubjects_rmNRV(:,:,stk(sb),:));
    sbj.Subjects = cat(2,sbj.Subjects,['h' num2str(sb)]);
    sbj.Subject_list = cat(2,sbj.Subject_list,['h' num2str(sb)]);
end

% Plots boxplot for each human subject, and in the end the average across subjects 
[corr_pvals_indiv,corr_pvals_bigcats] = PlotNSEBoxplots(D,sbj,SoundType,method);

