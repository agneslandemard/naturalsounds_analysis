%%%  Plot results for experiment II (Fig 4 & Fig S11)

% Load data from experiment II
P.experiment = 'vocalization';
P.n_ics = 8;
P.name = 'rmyVersion';

V = LoadVenoisedVata(P);

A_Idx = find(contains(V.hemis,'A'));
T_Idx = find(contains(V.hemis,'T'));
C_Idx = find(contains(V.hemis,'C'));

V.SrecoA = V.Sreco(:,ismember(V.si,A_Idx));
V.SrecoT = V.Sreco(:,ismember(V.si,T_Idx));
V.SrecoC = V.Sreco(:,ismember(V.si,C_Idx));

V.SrecoFullA = V.SrecoFull(:,ismember(V.si,A_Idx),:);
V.SrecoFullT = V.SrecoFull(:,ismember(V.si,T_Idx),:);
V.SrecoFullC = V.SrecoFull(:,ismember(V.si,C_Idx),:);

load([additional_path '/Coordinates/distances_to_pac_' P.experiment '.mat'])

%% plot NSE Map (panel V)

hemis_to_plot = 1:V.n_hemis;

n_conds = V.param.snd.Nmodels; % will compare natural sounds to this model (spectrotemporal)

% define color scale
clear Color
Color.ColorAxis = [0 1];
Color.cm = cmap_from_name('lightblue-to-yellow1');

figure('Position',[136 354 1305 415]);
n_cols = 5;
nsm = NSE_noise_corrected(V.SrecoFull(V.param.snd.idxSound(:,1),:,1),V.SrecoFull(V.param.snd.idxSound(:,1),:,2),...
    V.SrecoFull(V.param.snd.idxSound(:,n_conds),:,1),V.SrecoFull(V.param.snd.idxSound(:,n_conds),:,2),1);

for h = 1: length(hemis_to_plot)
    hemi = hemis_to_plot(h);
    
    xi = V.si == hemi;
    nse_single_subj = nsm(xi);
    X = load([data_path V.hemis{hemi} V.data_suffix '.mat'],'param');
    Q = X.param;
    
    subplot(length(hemis_to_plot),n_cols,1+(h-1)*n_cols)
    PlotTopView(nse_single_subj,Q,Color);
    title([V.hemis{hemi} ', mean= ' num2str(mean(nse_single_subj),2) ])
    
end

 %% Plot difference maps (panel E)

SoundType = {'Ferrets','Speech','Music','Others'};
m = 2;
clear Color
Color.ColorAxis = [-m m];
Color.cm = cmap_from_name('cbrewer-blue-red');

denom_by_pix = nanstd(V.Sreco(V.param.snd.idxSound(:,[1 n_conds]),:),[],1);
diff_all = V.Sreco(V.param.snd.idxSound(:,1),:,1) - V.Sreco(V.param.snd.idxSound(:,n_conds),:,1);
    
for sd = 1:length(SoundType)

    SdsIdx = SelectSounds(SoundType{sd},V.param);
    
    for h = 1: length(hemis_to_plot)
        hemi = hemis_to_plot(h);
        
        % pick out one subject
        xi = V.si == hemi;
        diff_single_subj = nanmean(diff_all(SdsIdx,xi),1)./denom_by_pix(xi);
        
        X = load([data_path V.hemis{hemi} V.data_suffix '.mat'],'param');
        Q = X.param;
        
        subplot(length(hemis_to_plot),n_cols,1+sd+(h-1)*n_cols)
        PlotTopView(diff_single_subj,Q,Color);
        title([SoundType{sd} ' ' V.hemis{hemi}])
    end
end

%% Plot NSE as function of distance to PAC, by category (panel G)

plot_type = 'bounded';
SoundType = {'Ferrets','Speech','Music'};

% Ferret A
dist = load([additional_path '/Coordinates/distances_to_pac_' P.experiment '.mat']);
dist.distance_to_pac = dist.distance_to_pac(ismember(V.si,A_Idx));
nse_by_dis = NSE_disFromPAC(V.SrecoFullA,V.si(ismember(V.si,A_Idx)),V.param,dist,SoundType);

colors = [0.949    0.604	0.722  ;...
        0.1882    0.7490    0.6196 ;...
        0.0549    0.3020    0.5843 ];
    
figure('Position',[489 462 500 336]);
axf1 = subplot(1,3,1);
hold all
for sd = 1:length(SoundType)
   switch plot_type
        case 'normal'
            plot((dist.distances+2.5)./10,snm(nse_by_dis(sd,:,:,1),3),'color',colors(sd,:),'LineWIdth',2);
            plot((dist.distances+2.5)./10,squeeze(nse_by_dis(sd,:,:,1)),'color',colors(sd,:),'LineWIdth',0.3);
        case 'bounded'
            boundedline((dist.distances+2.5)./10,snm(nse_by_dis(sd,:,:,1),3),snm(nse_by_dis(sd,:,:,2),3),'cmap',colors(sd,:),'alpha');
   end
end
xlabel('Distance to PAC (mm)')
title('Ferret A')
yticks(0:0.2:0.4)
xlim([0 4.5])
yticks(0:0.5:1)
ylim([0 1])
xticks(dist.distances(1:2:end)./10)

% Ferret T
dist = load([additional_path '/Coordinates/distances_to_pac_' P.experiment '.mat']);
dist.distance_to_pac = dist.distance_to_pac(ismember(V.si,T_Idx));
nse_by_dis = NSE_disFromPAC(V.SrecoFullT,V.si(ismember(V.si,T_Idx))-1,V.param,dist,SoundType);

colors = [0.949    0.604	0.722  ;...
        0.1882    0.7490    0.6196 ;...
        0.0549    0.3020    0.5843 ];

axf2 = subplot(1,3,2);
hold all
for sd = 1:length(SoundType)
   switch plot_type
        case 'normal'
            plot((dist.distances+2.5)./10,snm(nse_by_dis(sd,:,:,1),3),'color',colors(sd,:),'LineWIdth',2);
            plot((dist.distances+2.5)./10,squeeze(nse_by_dis(sd,:,:,1)),'color',colors(sd,:),'LineWIdth',0.3);
        case 'bounded'
            boundedline((dist.distances+2.5)./10,snm(nse_by_dis(sd,:,:,1),3),snm(nse_by_dis(sd,:,:,2),3),'cmap',colors(sd,:),'alpha');
   end
end
xlabel('Distance to PAC (mm)')
title('Ferret T')
%yticks(0:0.2:0.4)
xlim([0 4.5])

yticks(0:0.5:1)
ylim([0 1])
xticks(dist.distances(1:2:end)./10)

% Ferret C
dist = load([additional_path '/Coordinates/distances_to_pac_' P.experiment '.mat']);
dist.distance_to_pac = dist.distance_to_pac(ismember(V.si,C_Idx));
nse_by_dis = NSE_disFromPAC(V.SrecoFullT,V.si(ismember(V.si,C_Idx))-1,V.param,dist,SoundType);

colors = [0.949    0.604	0.722  ;...
        0.1882    0.7490    0.6196 ;...
        0.0549    0.3020    0.5843 ];

axf2 = subplot(1,3,3);
hold all
for sd = 1:length(SoundType)
   switch plot_type
        case 'normal'
            plot((dist.distances+2.5)./10,snm(nse_by_dis(sd,:,:,1),3),'color',colors(sd,:),'LineWIdth',2);
            plot((dist.distances+2.5)./10,squeeze(nse_by_dis(sd,:,:,1)),'color',colors(sd,:),'LineWIdth',0.3);
        case 'bounded'
            boundedline((dist.distances+2.5)./10,snm(nse_by_dis(sd,:,:,1),3),snm(nse_by_dis(sd,:,:,2),3),'cmap',colors(sd,:),'alpha');
   end
end
xlabel('Distance to PAC (mm)')
title('Ferret C')
xlim([0 4.5])

yticks(0:0.5:1)
ylim([0 1])
xticks(dist.distances(1:2:end)./10)
 
%%
method = 'ncNSE_bycat';
SoundType = {'Ferrets','Speech','Music'};
sbj.Subjects = {'SrecoFullA','SrecoFullT','SrecoFullC'};
sbj.Subject_list = {'A','T','C'};
 
[corr_pvals_indiv,corr_pvals_bigCats,stats] = PlotNSEBoxplots(V,sbj,SoundType,method);
yticks(-0.2:0.2:1)
yticklabels(-0.2:0.2:1)
ylim([-0.2 1])

%%

clear nse_by_dis
dist = load([additional_path '/Coordinates/distances_to_pac_' P.experiment '.mat']);
dist_tmp = dist;
figure('Position', [680 806 171 172]);
sdcat = 1;
hold all
for f = 1:length(sbj.Subject_list)
    F_Idx = find(contains(V.hemis,sbj.Subject_list{f}));
    dist_tmp.distance_to_pac = dist.distance_to_pac(ismember(V.si,F_Idx));
    
    nse_by_dis(:,:,:,f) = snm(NSE_disFromPAC(V.(['SrecoFull' sbj.Subject_list{f}]),V.si(ismember(V.si,F_Idx))-min(F_Idx)+1,V.param,dist_tmp,SoundType),3);
    plot((dist.distances+2.5)./10,squeeze(nse_by_dis(sdcat,:,1,f)),'color',mean(V.param.plt.SoundColors(SelectSounds('Ferrets',V.param),:),1))
 
end

plot((dist.distances+2.5)./10,snm(nse_by_dis(sdcat,:,1,:),4),'color',mean(V.param.plt.SoundColors(SelectSounds('Ferrets',V.param),:),1),'LineWidth',2)

xlim([0 4.5])

yticks(0:0.2:1)
ylim([-0.1 1])
xticks(dist.distances(1:2:end)./10)
xlabel('Distance to PAC')
ylabel('NSE(noise-corrected)')
title('vocaliizations')