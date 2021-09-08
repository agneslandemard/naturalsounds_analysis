%% FULL DENOISING PIPELINE

%%% Define experiment type and version name
experiment = 'vocalization'; % natural or vocalization
version_name = [experiment '_myVersion' ] ;

%%% Where to load and save data
save_path = [analysis_path version_name '/'];

if ~ exist(save_path,'dir')
    mkdir(save_path)
end

% Experiment parameters
switch experiment
    case 'natural' % experiment I
        hemis = {'NAL', 'NAR', 'NTR', 'NTL'};
        n_stims = 200;
        animals = {'A','T'};
        
    case 'vocalization' % experiment II
        hemis = {'VAL', 'VTL', 'VTR', 'VCL'};
        n_stims = 120;
        animals = {'A','T','C'};

end
n_tps = 20;
baseline_tps = 1:7;
n_reps = 4;
n_hemis = length(hemis);
n_animals = length(animals);
        
% load file defining folds for cross-validation
load([additional_path 'Leaveout/leaveout_all_' experiment '.mat'],'leaveout_all');
n_runs = max(leaveout_all);

%% Define denoising parameters

% part I
cca_denoise = 1; % choose whether or not to perform CCA.
% If not, parameters below are irrelevant.
cca_params.recenterIn = 1;
cca_params.recenterOut = 1;
cca_params.pc2Keep = 250;
cca_params.cc2Remove = 20;
cca_params.what2Remove = 'out';
cca_params.baseline_tps = baseline_tps;

normalize_by_trial = 1;

% part II
n_components = [1:10 15 20]; % values to explore for cross-validation
finalK = 8; % final value chosen to save denoised data and perform ICA
run_ica = 1; % choose whether or not to perform ICA on denoised data
% If not, parameters below are irrelevant.
ica_params.n_random_inits = 10;
ica_params.plt_figures = false;
ica_params.randseed = 1;

%% DENOISING PART I : CCA WITH OUT OF ROI

fprintf('Denoising part I: CCA with out-of-ROI signal \n'); drawnow;

D_corr = cell(1,n_hemis);
slice_refs = cell(1,n_hemis);

n_total_voxels = 0;
clear cat

for j = 1 : n_hemis
    
    sliceList = dir(fullfile(data_path, [hemis{j} '_*.mat']));
    idx = arrayfun(@(x)(~isempty(regexp(x.name,[hemis{j} '_\d+'],'once'))),sliceList);
    sliceList = sliceList(idx);
    
    % Sort list according to last number of the file name
    lastDigit = arrayfun(@(x)((regexp(x.name,'_(\d+).mat','tokens'))),sliceList);
    lastDigit = cellfun(@(x)(str2double(x{1})),lastDigit);
    [~,idx] = sort(lastDigit);
    sliceList = sliceList(idx);
    
    n_slices = length(sliceList);
    
    D_corr_hemi = nan(n_tps, n_stims, n_reps, 0);
    SliceRef = nan(0,1);
    
    for ii = 1 : n_slices
        
        % Load data for this slice
        load(fullfile(sliceList(ii).folder,sliceList(ii).name),'Din','Dout','Mask_slice');
        % Each D is a matrix time x stim x repetition x pixels
        
        if cca_denoise
            disp(['Correcting slice #' num2str(ii)])
            if exist('Mask_slice','var')
                cca_params.Mask_slice = Mask_slice; %this will lead to display of CC weights
            end
            Din_corr = CCAcorrection(Din,Dout,cca_params);
            % the idea is to do a CCA between what's in and what's outside the
            % mask. It suppposes that there's no neural signal directly related to sound
            % outside of the mask
            
        else
            Din_corr = Din;
        end
        
        % save data for this slice, and keep track of voxel indices
        SliceRef = cat(1,SliceRef,ii*ones(size(Din,4),1));
        D_corr_hemi = cat(4,D_corr_hemi,Din_corr);
        
    end
    
    slice_refs{j} = SliceRef;
    D_corr{j} = D_corr_hemi;
    n_total_voxels = n_total_voxels + length(SliceRef);
    
end

%% Normalize by baseline
disp('Normalize hemi by baseline')
D_norm = cell(1, n_hemis);
for j = 1 : n_hemis
    if normalize_by_trial
        
        baseline = nanmean(D_corr{j}(baseline_tps,:,:,:),1);
        
        % subtract baseline separately for each presentation of each sound
        D_norm{j} = bsxfun(@minus, D_corr{j}, baseline);
        
        % and then divide by the mean across time, repetitions, and stimuli
        D_norm{j} = bsxfun(@times, D_norm{j}, 1./nanmean(nanmean(baseline,3),2)); % global mean
    else
        
        baseline = nanmean(D_corr{j}(baseline_tps,:,:,:),[1 2 3]);
        D_norm{j} = bsxfun(@minus, D_corr{j}, baseline);
        D_norm{j} = bsxfun(@times, D_norm{j}, 1./baseline); % global mean
        
    end
    
end
save([save_path '/D_norm.mat'],'D_norm','leaveout_all','slice_refs','n_total_voxels','hemis','-v7.3');

%% DENOISING PART II : DSS analysis
for an = 1:n_animals
    
    hemis_to_keep  = contains(hemis,animals{an});
    tmp_hemis = hemis(hemis_to_keep);
    tmp_n_hemis = length(tmp_hemis);
    tmp_D_norm = D_norm(hemis_to_keep);
    
    tmp_slice_refs = slice_refs(hemis_to_keep);
    
    for run = 1:n_runs
        
        fprintf('DSS analysis\n'); drawnow;
        
        % define left out sounds
        leaveout = leaveout_all == run;
        
        % Compute whitened data matrix
        fprintf('Whitening\n'); drawnow;
        D_white = cell(1, tmp_n_hemis);
        for j = 1:tmp_n_hemis
            
            [n_tps, n_train_stims, n_reps, n_voxels] = size(tmp_D_norm{j}(:,~leaveout,:,:));
            
            D_white{j} = nan([n_tps*n_train_stims, n_reps, n_voxels]);
            
            for i = 1:n_reps
                
                n_slices = max(tmp_slice_refs{j});
                
                for ii = 1:n_slices
                    
                    vox_in_slice = find(tmp_slice_refs{j}==ii);
                    
                    % select one rep and one slice
                    % concatenate time and stimulus dimension
                    Z = reshape(tmp_D_norm{j}(:,~leaveout,i,vox_in_slice), [n_tps*n_train_stims, length(vox_in_slice)]);
                    
                    % remove missing tps/voxels
                    
                    missing_tps = any(isnan(Z),2);
                    missing_voxels = any(isnan(Z),1);
                    
                    Z = Z(~missing_tps,~missing_voxels);
                    
                    % demean rows (i.e. pattern across voxels for single timepoint)
                    % if not already
                    repmean = mean(Z,2);
                    Z = bsxfun(@minus, Z, repmean);
                    
                    % whiten
                    [U,~,V] = svd(Z, 'econ');
                    
                    Z_white = U*V';
                    clear U V;
                    
                    % add back to matrix
                    D_white{j}(~missing_tps,i,vox_in_slice(~missing_voxels)) = Z_white;
                    clear vox_in_slice Z_white;
                    
                end
                
            end
        end
        
        % Compute DSS
        US_all = [];
        US_nobias = [];
        for j = 1:tmp_n_hemis
            
            % average across reps
            Xm = squeeze(nanmean(D_white{j},2));
            
            % reshape to separate out time and stimuli
            n_voxels = size(D_white{j},3);
            Xmr = reshape(Xm, [n_tps, n_train_stims, n_voxels]);
            
            % select timepoints not in the baseline and rewrap stimuli and time
            Xbias = Xmr(setdiff(1:n_tps, baseline_tps),:,:);
            Xbias = reshape(Xbias, [(n_tps-length(baseline_tps))*n_train_stims, n_voxels]);
            
            % find nan voxels
            nn_vx = find(~isnan(mean(Xbias,1)));
            
            % PCA
            [U,S,V] = svd(Xbias(:,nn_vx), 'econ');
            US_all = cat(2, US_all, U*S);
            
            % compute PCs for the original data with the baseline period
            US_nobias = cat(2, US_nobias, Xm(:,nn_vx)*V);
            
        end
        clear Xmr X Xbias;
        
        
        % PCA on the concatenated PCs from all hemispheres
        [~,~,V_dss] = svd(US_all, 'econ');
        U_dss = US_nobias * V_dss;
        
        % Save inferred components for this run
        save([save_path '/' animals{an} '/Run_' num2str(run) '_components.mat'],'U_dss','leaveout')
        
    end
    
    % Format data and prepare two splits of data (average across even or odd reps)
    
    Dodd = [];
    Deven = [];
    si = [];
    for j = 1:length(tmp_D_norm)
        Dtmp = tmp_D_norm{j};
        Dodd = cat(3, Dodd, squeeze(nanmean(Dtmp(:,:,1:2:end,:),3))); % odd reps
        Deven = cat(3, Deven, squeeze(nanmean(Dtmp(:,:,2:2:end,:),3))); % even reps
        si = cat(2, si, j*ones(1,size(Dtmp,4)));
 
    end
    Dall = cat(4, Dodd, Deven);
    clear Dodd Deven
    
    n_total_voxels = size(Dall,3);
    
    % Save recomposed data for chosen nb of ICs
    if exist('finalK','var')
        % number of components to keep
        K = finalK;
        D_denoised = nan(n_tps,n_stims,n_total_voxels,2);
        
        for run = 1 : n_runs
            load([save_path '/Run_' num2str(run) '_components.mat'],'U_dss','leaveout')
            
            n_train_stims = length(find(~leaveout));
            n_test_stims = length(find(leaveout));
            
            Dtrain = nanmean(Dall(:,~leaveout,:,:),4);
            Dtestodd = Dall(:,leaveout,:,1);
            Dtesteven = Dall(:,leaveout,:,2);
            
            Dtrain = reshape(Dtrain, n_tps*n_train_stims, n_total_voxels);
            Dtestodd = reshape(Dtestodd, n_tps*n_test_stims, n_total_voxels);
            Dtesteven = reshape(Dtesteven, n_tps*n_test_stims, n_total_voxels);
            
            W_dss = pinv(U_dss(:,1:K))*Dtrain;
            W_dss_pinv = pinv(Dtrain)*U_dss(:,1:K);
            
            disp(['Selected top ' num2str(K) ' components.'])
            
            Rtestodd = Dtestodd * W_dss_pinv;
            Rtesteven = Dtesteven * W_dss_pinv;
            
            D_denoised(:,leaveout,:,1)=reshape(Rtestodd * W_dss,n_tps,n_test_stims, n_total_voxels);
            D_denoised(:,leaveout,:,2)=reshape(Rtesteven * W_dss,n_tps,n_test_stims, n_total_voxels);
            
        end
        
        % save data recomposed with this amount of components
        save([save_path '/' animals{an} '/D_denoised_K' num2str(K) '.mat'],'D_denoised','si','-v7.3');
        
    end
end


%% Load data from all animals
Dodd = [];
Deven = [];
si = [];  % hemi index
ai = []; % animal index
for j = 1:length(D_norm)
    Dtmp = D_norm{j};
    Dodd = cat(3, Dodd, squeeze(nanmean(Dtmp(:,:,1:2:end,:),3))); % odd reps
    Deven = cat(3, Deven, squeeze(nanmean(Dtmp(:,:,2:2:end,:),3))); % even reps
    si = cat(2, si, j*ones(1,size(Dtmp,4)));
    a = find(strcmp(hemis{j}(2),animals));
    ai = cat(2, ai, a*ones(1,size(Dtmp,4)));

end
Dall = cat(4, Dodd, Deven);
clear Dodd Deven

n_total_voxels = size(Dall,3);

%% Cross validate nb of components

RespWindow = baseline_tps(end)+(3:11);

% Test-retest correlation
baseCorr = nan(n_total_voxels,1);
for px = 1:n_total_voxels
    baseCorr(px) = corr(squeeze(nanmean(Dall(RespWindow,:,px,1),1))',squeeze(nanmean(Dall(RespWindow,:,px,2),1))','rows','pairwise');
end

CV_corr = nan(length(n_components),n_total_voxels,2);

for kk = 1:length(n_components)
    
    % number of components to keep
    K = n_components(kk);
    D_denoised = nan(size(Dall));
    
    for an = 1:length(animals)
        for run = 1 : n_runs
            
            load([save_path '/' animals{an} '/Run_' num2str(run) '_components.mat'],'U_dss','leaveout')
            n_vxs = length(find(ai == an));
            n_train_stims = length(find(~leaveout));
            n_test_stims = length(find(leaveout));
            Dtrain = nanmean(Dall(:,~leaveout,ai == an,:),4);
            Dtestodd = Dall(:,leaveout,ai == an,1);
            Dtesteven = Dall(:,leaveout,ai == an,2);
            
            Dtrain = reshape(Dtrain, n_tps*n_train_stims, n_vxs);
            Dtestodd = reshape(Dtestodd, n_tps*n_test_stims, n_vxs);
            Dtesteven = reshape(Dtesteven, n_tps*n_test_stims, n_vxs);
            
            miss_tps = any(isnan(Dtrain),2);
            
            W_dss = pinv(U_dss(~miss_tps,1:K))*Dtrain(~miss_tps,:);
            W_dss_pinv = pinv(Dtrain(~miss_tps,:))*U_dss(~miss_tps,1:K);
            
            D_dss = U_dss(:,1:K) * W_dss;
            disp(['Selected top ' num2str(K) ' components.'])
            
            
            Rtestodd = Dtestodd * W_dss_pinv;
            Rtesteven = Dtesteven* W_dss_pinv;
            D_denoised(:,leaveout,ai == an,1) = reshape(Rtestodd * W_dss,n_tps,n_test_stims, n_vxs);
            D_denoised(:,leaveout,ai == an,2) = reshape(Rtesteven * W_dss,n_tps,n_test_stims, n_vxs);
        end
    end
    
    for px = 1 :n_total_voxels
        CV_corr(kk,px,1) = corr(squeeze(nanmean(Dall(RespWindow,:,px,1),1))',squeeze(nanmean(D_denoised(RespWindow,:,px,2),1))','rows','pairwise');
        CV_corr(kk,px,2) = corr(squeeze(nanmean(Dall(RespWindow,:,px,2),1))',squeeze(nanmean(D_denoised(RespWindow,:,px,1),1))','rows','pairwise');
    end
    
    
end

figure;
plot(n_components,nanmedian(nanmean(CV_corr,3),2),'k','LineWidth',1)
x = get(gca,'xlim');
hold on;
plot(x,repmat(sqrt(nanmedian(baseCorr)),2,1),'r-');
xlabel('# Components')
ylabel('Median correlation across voxels between raw and denoised')
save([save_path '/CV_NbComponents.mat'],'n_components','CV_corr','baseCorr')

%% Run ICA on denoised data

n_ics = finalK;  % number of ICs

if run_ica
    
    % load denoised data
    D_denoised = nan(size(Dall));
    for an = 1:n_animals
        data_an = load([save_path '/' animals{an} '/D_denoised_K' num2str(n_ics) '.mat'],'D_denoised');
        D_denoised(:,:,ai == an,:) = data_an.D_denoised;
    end
    
    Dodd = reshape(Dall(:,:,:,1), n_tps*n_stims, n_total_voxels);
    Deven = reshape(Dall(:,:,:,2), n_tps*n_stims, n_total_voxels);
    Dav = reshape(nanmean(Dall,4), n_tps*n_stims, n_total_voxels);
    
    D_denoised_tmp = reshape(nanmean(D_denoised,4),n_tps*n_stims, n_total_voxels);
    
    % Run ICA
    [R_ica_train_tc, W_ica_train] = ...
        nonparametric_ica(D_denoised_tmp,n_ics, ica_params.n_random_inits, ica_params.plt_figures, ica_params.randseed);
    
    
    % orient so that response during the stimulus period is positive
    R_ica_train_tc = reshape(R_ica_train_tc, [n_tps, n_stims, n_ics]);
    RespWindow = ceil([size(R_ica_train_tc,1)/2  size(R_ica_train_tc,1)]); % Response window, here considered second half of trial.
    sgn = sign(mean(mean(R_ica_train_tc(RespWindow,:,:),1), 2));
    R_ica_train_tc = bsxfun(@times, R_ica_train_tc, sgn);
    W_ica_train = bsxfun(@times, W_ica_train, squeeze(sgn));
    R_ica_train_tc = reshape(R_ica_train_tc, [n_tps*n_stims,  n_ics]);
    W_ica_train_inverse = pinv(Dav) * R_ica_train_tc;
    
    % Get response to test stimuli
    
    R_ica_test_odd_tc = Dodd * W_ica_train_inverse;
    R_ica_test_even_tc = Deven * W_ica_train_inverse;
    
    R_ica_test_odd_tc = reshape(R_ica_test_odd_tc, [n_tps,n_stims, n_ics]);
    R_ica_test_even_tc = reshape(R_ica_test_even_tc, [n_tps,n_stims,  n_ics]);
    R_ica_train_tc = reshape(R_ica_train_tc, [n_tps,n_stims,  n_ics]);
    
    % save results
    save([save_path '/ICA_results' num2str(n_ics) '.mat'], ...
        'R_ica_train_tc', 'W_ica_train','si', 'R_ica_test_odd_tc', 'R_ica_test_even_tc', 'W_ica_train_inverse');
    
end


