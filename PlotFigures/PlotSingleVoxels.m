%%% Plot responses for single voxels 
% Fig 1C and 2A-B

%% Load recomposed data

clear P
P.experiment = 'natural';
P.n_ics = NbICs;
P.name = version_name;

D = LoadDenoisedData(P);

%% Load Raw data

DP = D;
DP.fullRaw = 1; % can be set to 0 to get data after CCA but before DSS
B = LoadRawData(DP);

%% Plot raw and denoised timecourse for chosen single voxel (Fig 1)
ExampleVoxel = 35056;
datatype_list = {'raw','denoised'};

onset = D.param.exp.PreStimSilence;
n_conds = D.param.snd.Nmodels;

figure;
for d = 1:length(datatype_list)
    
    data = datatype_list{d};
    
    switch data
         case 'raw'
            DD = B;
        case 'denoised'
            DD = D;
    end
    
    R = snm(DD.Sreco(:,ExampleVoxel),2);
    R_tc = snm(DD.SrecoTC(:,:,ExampleVoxel),3);
    bounds = [min(R(:)), max(R(:))];
    bounds = bounds + 0.1*[-1 1]*diff(bounds);
    [n_stims_per_cond, n_conds] = size(D.param.snd.idxSound);
    R_tc = reshape(R_tc,20,n_stims_per_cond, n_conds);
    
    sound_colors = reshape(D.param.plt.SoundColors',[3,n_stims_per_cond, n_conds]);
    
    subplot(1,length(datatype_list), d)
    hold all
    patch([onset onset+10 onset+10 onset],[-0.3 -0.3 1 1],[1 1 1].*0.8,'FaceAlpha',0.6,'EdgeColor','none');
    for k = 1:n_stims_per_cond
        plot(R_tc(:,k,1), '-', 'LineWidth', 0.5, 'Color', sound_colors(:, k, 1));
    end
    
    plot(mean(mean(R_tc,3),2), 'k-', 'LineWidth', 2);
    axis tight
    ylabel('Denoised %CBV')
    yticks(0:0.4:max(mat2vec(R_tc(:,:,1))))
    yticklabels(100.*(0:0.4:max(mat2vec(R_tc(:,:,1)))))
    xlabel('Time (s)')
    ylim([-0.3 1])
    
    xticks(onset+(-5:5:10))
    xticklabels(-5:5:10)
    title(datatype_list{d})
end

%% Plot figure 2A-D: test-retest and nat vs synth for two single voxels

figure;
 
% primary and non-primary voxel in NTL, slices 16 and 9 respectively
% their position on the slice can be viewed at the end of this script

ExampleVoxels = [39672 35056];

for vx = 1:length(ExampleVoxels)
    
    
    [n_stims_per_cond, n_conds] = size(D.param.snd.idxSound);
    
    sound_colors = reshape(D.param.plt.SoundColors',[3,n_stims_per_cond, n_conds]);
    R_rp = reshape(snm(D.SrecoFull(:,ExampleVoxels(vx),:),2),n_stims_per_cond, n_conds,2);
    bounds = [min(mat2vec(R_rp(:,[1 n_conds],:))), max(mat2vec(R_rp(:,[1 n_conds],:)))];
    bounds = bounds + 0.1*[-1 1]*diff(bounds);
    
    % Natural vs synthetic
    subplot(length(ExampleVoxels), 2, (vx-1)*2+1)
    hold on
    newR = nan(n_stims_per_cond,2);
    for k = 1:n_stims_per_cond
        rdrp1 = randi(2); rdrp2 = setdiff(1:2,rdrp1);
        
        %rdrp1 = 1; rdrp2 = 2;
        scatter(R_rp(k,1,rdrp1), R_rp(k,end,rdrp2),[],sound_colors(:, k, 1)','filled');
        
        xlim(bounds); ylim(bounds);
        xlabel('Natural Sounds (%dCBV)');
        ylabel('spectrotemporal')
        newR(k,:) = [R_rp(k,1,rdrp1); R_rp(k,end,rdrp2)];
    end
    
    plot(bounds, bounds, 'r--', 'LineWidth', 2)
    xticks(round(bounds(2),1)-2*round(diff(bounds)/3,1):round(diff(bounds)/3,1):round(bounds(2),1))
    yticks(round(bounds(2),1)-2*round(diff(bounds)/3,1):round(diff(bounds)/3,1):round(bounds(2),1))
    xticklabels(xticks.*100)
    yticklabels(yticks.*100)
    
    title(['NSE = ' num2str(NSE(newR(:,1), newR(:,2)),2)])% ', nc =  ' num2str(nse,2)])
    axis equal tight
    
    % Test-retest
    subplot(length(ExampleVoxels), 2, (vx-1)*2+2)
    if n_conds == 5
        symbols = {'o','','','', '+'};
    else
        symbols = {'o', '+'};
    end
    hold all
    
   newR = nan(n_stims_per_cond,n_conds,2);
    for k = 1:n_stims_per_cond
        for j = 1:n_conds-1:n_conds
            rdrp1 = randi(2); rdrp2 = setdiff(1:2,rdrp1);
            
            cl = sound_colors(:,k,j);
            scatter(R_rp(k,j,rdrp1),R_rp(k,j,rdrp2),[],cl(:,:)',symbols{j});
            newR(k,j,:) = [R_rp(k,j,rdrp1) R_rp(k,j,rdrp2)];
        end
       
    end
    xlim(bounds); ylim(bounds);
    plot(bounds, bounds, 'r--')
    xticks(round(bounds(2),1)-2*round(diff(bounds)/3,1):round(diff(bounds)/3,1):round(bounds(2),1))
    yticks(round(bounds(2),1)-2*round(diff(bounds)/3,1):round(diff(bounds)/3,1):round(bounds(2),1))
    xticklabels(xticks.*100)
    yticklabels(yticks.*100)
    
    title(['NSE = ' num2str(NSE(mat2vec(newR(:,[1 n_conds],1)),mat2vec(newR(:,[1 n_conds],2))),2)])
    axis equal tight
    xlabel('test')
    ylabel('retest')

end

%% Check position of chosen voxel

ExampleVoxel = 35056;
type = 'locate';
infos2.voxel = ExampleVoxel;
infos2.plot = 1;
vx2 = LocateVoxel(type,infos2,D);

%% Choose other voxel

type = 'choose';
infos2.hemi = 'NTL';
infos2.slice = 5;
ExampleVoxel = LocateVoxel(type,infos2,D);


%% Plot human voxels, denoised

load([additional_path '/Human/HumanData.mat'],'SrecoFull_human_allSubjects','n_subj','xl','yl','xr','yr')

tmp =  reshape(SrecoFull_human_allSubjects,200,xl*yl+xr*yr,n_subj,2);

ex_subj = 13;
D_LH = reshape(tmp(:,1:xl*yl,ex_subj,:),200,xl,yl,2);
ex_primary_voxel = reshape(squeeze(D_LH(:,21, 25, :)),40,5,2);
ex_non_primary_voxel = reshape(squeeze(D_LH(:,18, 41, :)),40,5,2);
clear cat
Resp = cat(4,ex_primary_voxel,ex_non_primary_voxel);

figure('Position', [440 464 458 334])
for vx = 1:2

    [n_stims_per_cond, n_conds] = size(D.param.snd.idxSound);
    sound_colors = reshape(D.param.plt.SoundColors',[3,n_stims_per_cond, n_conds]);
    R_rp = Resp(:,:,:,vx);
    bounds = [min(mat2vec(R_rp(:,[1 n_conds],:))), max(mat2vec(R_rp(:,[1 n_conds],:)))];
    bounds = bounds + 0.1*[-1 1]*diff(bounds);
    
    % Natural vs synthetic
    subplot(2, 2, (vx-1)*2+1)
    hold on
    newR = nan(n_stims_per_cond,2);
    for k = 1:n_stims_per_cond
        rdrp1 = randi(2); rdrp2 = setdiff(1:2,rdrp1);
        scatter(R_rp(k,1,rdrp1), R_rp(k,end,rdrp2),[],sound_colors(:, k, 1)','filled'); 
        xlim(bounds); ylim(bounds);
        xlabel('Natural Sounds (%dCBV)');
        ylabel('spectrotemporal')
        newR(k,:) = [R_rp(k,1,rdrp1); R_rp(k,end,rdrp2)];
    end
    
    plot(bounds, bounds, 'r--', 'LineWidth', 2)
    
    title(['NSE = ' num2str(NSE(newR(:,1), newR(:,2)),2)])% ', nc =  ' num2str(nse,2)])
    axis equal tight
    
    % Test-retest
    subplot(2, 2, (vx-1)*2+2)
    if n_conds == 5
        symbols = {'o','','','', '+'};
    else
        symbols = {'o', '+'};
    end
    hold all
    
   newR = nan(n_stims_per_cond,n_conds,2);
    for k = 1:n_stims_per_cond
        for j = 1:n_conds-1:n_conds
            rdrp1 = randi(2); rdrp2 = setdiff(1:2,rdrp1);
            cl = sound_colors(:,k,j);
            scatter(R_rp(k,j,rdrp1),R_rp(k,j,rdrp2),[],cl(:,:)',symbols{j});
            newR(k,j,:) = [R_rp(k,j,rdrp1) R_rp(k,j,rdrp2)];
        end
       
    end
    xlim(bounds); ylim(bounds);
    plot(bounds, bounds, 'r--')
    nsev = 0.5*(NSE(mat2vec(newR(:,[1],1)),mat2vec(newR(:,[1],2)))+NSE(mat2vec(newR(:,[5],1)),mat2vec(newR(:,[5],2))));
    title(['NSE = ' num2str(nsev,2)])
    axis equal tight
    xlabel('test')
    ylabel('retest')

end
