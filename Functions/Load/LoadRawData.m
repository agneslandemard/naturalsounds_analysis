function B = LoadRawData(P)

if ~ isfield(P,'normalize_by_trial') 
    P.normalize_by_trial = 1;
end
% P should be the output of LoadRecomposedData, and will load associated
% raw data
global data_path

if isfield(P,'fullRaw') && P.fullRaw
    % data before any denoising
    tmp.D_norm = cell(1, P.n_hemis);
    for j = 1 : P.n_hemis
        n_slices = length(dir([data_path P.hemis{j} '_*.mat']))-1;
        D_hemi = [];
        for ii = 1 : n_slices
            load([data_path P.hemis{j} '_' num2str(ii) '.mat'],'Din');
            D_hemi = cat(4,D_hemi,Din);
        end
        
        % Reshape D_norm (data before DSS/ICA) to fit reconstructed format
        if P.normalize_by_trial
            baseline = nanmean(D_hemi(1:P.param.exp.PreStimSilence,:,:,:),1);
            tmp.D_norm{j} = bsxfun(@minus, D_hemi, baseline);
            tmp.D_norm{j} = bsxfun(@times, tmp.D_norm{j}, 1./nanmean(baseline,[2 3]));
        else
            baseline = nanmean(D_hemi(1:P.param.exp.PreStimSilence,:,:,:),[1 2 3]);
            tmp.D_norm{j} = bsxfun(@minus, D_hemi, baseline);
            tmp.D_norm{j} = bsxfun(@times, tmp.D_norm{j}, 1./baseline);
       
        end
    end
    
else
    % data after denoising part I (CCA) but before part II (DSS)
    tmp = load([P.data_path '/D_norm.mat']);  
end

[n_tps, n_stims ,n_voxels, n_reps] = size(P.SrecoFullTC);

B.Sreco = nan(n_stims, n_voxels);
B.SrecoTC = nan(n_tps, n_stims, n_voxels);
B.SrecoFull = nan(n_stims, n_voxels, n_reps);
B.SrecoFullTC = nan(n_tps, n_stims, n_voxels, n_reps);
B.SrecoAllReps = nan(n_stims, n_voxels, P.param.exp.nTrials);
B.n_hemis = length(tmp.D_norm);
B.data_path = P.data_path;

for h = 1:B.n_hemis
    
    Dhemi = tmp.D_norm{h};
%     if P.normalize_by_trial
%         baseline = nanmean(Dhemi(1:P.param.exp.PreStimSilence,:,:,:),1);
%         Dhemi = bsxfun(@minus, Dhemi, baseline);  
%     end
    
    B.Sreco(:,P.si == h) = snm(Dhemi(P.param.exp.PreStimSilence+P.RespWindow,:,:,:),[1 3]);
    
    B.SrecoTC(:,:,P.si == h) = snm(Dhemi,3);
    
    B.SrecoFull(:,P.si == h,1) = snm(Dhemi(P.param.exp.PreStimSilence+P.RespWindow,:,[1 3],:),[1 3]);
    B.SrecoFull(:,P.si == h,2) = snm(Dhemi(P.param.exp.PreStimSilence+P.RespWindow,:,[2 4],:),[1 3]);
    
    B.SrecoFullTC(:,:,P.si == h,1) = snm(Dhemi(:,:,[1 3],:),3);
    B.SrecoFullTC(:,:,P.si == h,2) = snm(Dhemi(:,:,[2 4],:),3);
    
    B.SrecoAllReps(:,P.si == h,:) = permute(snm(Dhemi(P.param.exp.PreStimSilence+P.RespWindow,:,:,:),1),[1 3 2]);
   
end

end