%% Compare raw and denoised responses

P.experiment = 'natural';
P.n_ics = NbICs;
P.name = version_name;

% Load denoised data
D = LoadDenoisedData(P);

% Load raw data
B = LoadRawData(D);

resp_den = D.SrecoFull;
resp_raw = B.SrecoFull;


%% S1A : Plot scatter

ScatterRawDenoised(resp_den,resp_raw,'corr')

%% S1B : Plot map

hemi = 4;
%hemi = 1;
clear Color
Color.cm = cmap_from_name('cbrewer-reds');

figure;
% load hemi info
xi = D.si == hemi;
X = load([data_path D.hemis{hemi} D.data_suffix '.mat'],'param');
Q = X.param;

% Correlation denoised / raw
subplot(121)
diff_all = nan(1,size(resp_den ,2));
for pix = 1:size(diff_all ,2)
    diff_all(pix) = (1/2)*(corr(resp_den(:,pix,1),resp_raw(:,pix,2),'rows','pairwise')+corr(resp_den(:,pix,2),resp_raw(:,pix,1),'rows','pairwise'));
end
toplot = diff_all;
diff_single_subj = toplot(xi);
M = quantile(toplot,0.95);
Color.ColorAxis = [0 M];
PlotTopView(diff_single_subj,Q,Color);
title('Corr den/raw')

% Ceiling : sqrt of test-retest raw
subplot(122)
diff_all = nan(1,size(resp_den ,2));
for pix = 1:size(diff_all ,2)
    diff_all(pix) = corr(resp_raw(:,pix,1),resp_raw(:,pix,2),'rows','pairwise');
end
toplot = sqrt(diff_all);
diff_single_subj = toplot(xi);
PlotTopView(diff_single_subj,Q,Color);
title('Ceiling')

sgtitle([D.hemis{hemi},'Corr raw and reco, cb: ' num2str(M,2)])
    
%% S1C : correlation original/denoised by nb of components DSS
% from PlotCVResults.m, need to compute this before

load([D.data_path '/CV_NbComponents.mat'])

figure('Position', [440 408 226 390])
hold all
CV_corr_med = squeeze(nanmedian(snm(CV_corr,3),2));
CI = [prctile(CV_corr_med,2.5,2) prctile(CV_corr_med,97.5,2)];
l = boundedline(n_components,nanmean(CV_corr_med,2),abs(CI-nanmean(CV_corr_med,2)),'cmap',[0 0 0]);
l.LineWidth = 1;

% Ceiling
baseCorr_med = squeeze(nanmedian(real(sqrt(baseCorr)),1));
CI = [prctile(baseCorr_med,2.5) prctile(baseCorr_med,97.5)];
boundedline(n_components,repmat(nanmean(baseCorr_med),length(n_components),1),repmat(abs(CI-nanmean(baseCorr_med)),length(n_components),1),'cmap',0.6.*[1 1 1])

xlabel('#Components')
ylabel('Correlation (mean across voxels)')
vline(n_components(NbICs),'r:')

%% S11B : exemple slice before/after removing out-of-cortex components
% here for vocalizations

P.experiment = 'vocalization';
P.n_ics = NbICs;
P.name = version_name;

% Load denoised data
DV = LoadDenoisedData(P);

% Load raw data, after denoising CCA
BV = LoadRawData(DV);

% Load raw data, before denoising CCA
DP = DV;
DP.experiment = P.experiment;
DP.fullRaw = 1;
B0V = LoadRawData(DP);
clear DP

% Plot difference nat - mm , before and after CCA correction
% for an example slice
hemi = 1;
slice = 6;
cm = cmap_from_name('cbrewer-blue-red');

figure('Position',[440 515 245 283]);

X = load([data_path D.hemis{hemi} D.data_suffix '.mat'],'param');

DiffNatMM_before = snm(B0V.Sreco(D.param.snd.idxSound(:,1),:)-B0V.Sreco(D.param.snd.idxSound(:,end),:),1);
DiffNatMM_before = Pixs2Mat(DiffNatMM_before(D.si == hemi)',X.param.msk);

DiffNatMM_after = snm(BV.Sreco(D.param.snd.idxSound(:,1),:)-BV.Sreco(D.param.snd.idxSound(:,end),:),1);
DiffNatMM_after = Pixs2Mat(DiffNatMM_after(D.si == hemi)',X.param.msk);

m = quantile(mat2vec(DiffNatMM_before(:,:,slice)),.95);

imagesc([DiffNatMM_before(:,:,slice); DiffNatMM_after(:,:,slice)],'AlphaData',~isnan([DiffNatMM_before(:,:,slice); DiffNatMM_after(:,:,slice)]))
axis equal tight
colormap(cm); caxis([-m m]);
set(gca,'XTick',[],'YTick',[])
ylabel('after / / before')
title([D.hemis{hemi} ', slice ' num2str(slice)])

