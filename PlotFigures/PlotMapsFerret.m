%% Plot ferret maps

% to vary time window
P.RespWindow = 3:11;

% Load data
P.experiment = 'natural';
P.n_ics = NbICs;
P.name = version_name;
D = LoadDenoisedData(P);

%% NSE maps

% define color scale
clear Color
Color.ColorAxis = [0 1];
Color.cm = cmap_from_name('lightblue-to-yellow1');

% load distances to pac
load([additional_path '/Coordinates/distances_to_pac_' P.experiment '.mat'])

figure;
n_conds = D.param.snd.Nmodels;
hemis_to_plot = 1:D.n_hemis;

for md = 2 : n_conds
    
    % NSE between selected model and natural
    nsm = NSE_noise_corrected(D.SrecoFull(D.param.snd.idxSound(:,1),:,1),D.SrecoFull(D.param.snd.idxSound(:,1),:,2),...
        D.SrecoFull(D.param.snd.idxSound(:,md),:,1),D.SrecoFull(D.param.snd.idxSound(:,md),:,2),1);
    
    for j = 1:length(hemis_to_plot)
        
        hemi = hemis_to_plot(j);
        xi = D.si == hemi;
        nse_single_subj = nsm(xi);
        X = load([data_path  D.hemis{hemi} D.data_suffix '.mat'],'param');
        Q = X.param;
        
        subplot(length(hemis_to_plot),n_conds-1,md-1+(j-1)*(n_conds-1))
        PlotTopView(nse_single_subj,Q,Color);
        title(regexprep(D.param.snd.models{md},'synth_',''))
        
    end
end

%% Difference Maps

SoundType = {'Ferrets','Speech','Music','Others'};
n_cols = length(SoundType);
m = 2;
clear Color
Color.ColorAxis = [-m m];
Color.cm = cmap_from_name('cbrewer-blue-red');

denom_by_pix = nanstd(D.Sreco(D.param.snd.idxSound(:,[1 n_conds]),:),[],1);
diff_all = D.Sreco(D.param.snd.idxSound(:,1),:,1) - D.Sreco(D.param.snd.idxSound(:,n_conds),:,1);

figure;   
for sd = 1:length(SoundType)

    SdsIdx = SelectSounds(SoundType{sd},D.param);
    
    for h = 1: length(hemis_to_plot)
        hemi = hemis_to_plot(h);

        xi = D.si == hemi;
        diff_single_subj = nanmean(diff_all(SdsIdx,xi),1)./denom_by_pix(xi);
        
        X = load([data_path D.hemis{hemi} D.data_suffix '.mat'],'param');
        Q = X.param;
        
        subplot(length(hemis_to_plot),n_cols,sd+(h-1)*n_cols)
        PlotTopView(diff_single_subj,Q,Color);
        title([SoundType{sd} ' ' D.hemis{hemi}])
    end
end


