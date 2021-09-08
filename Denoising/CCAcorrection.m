function Din_corr = CCAcorrection(Din, Dout, cca_params)
% Performs CCA between regions "in" and "out" of ROI to denoise "in" by
% removing shared signals (considered as artefacts), as described in:
% Distinct higher-order representations of natural sounds in human and ferret
% auditory cortex (Landemard A, Bimbard C, Demené C, Shamma S,
% Norman-Haigneré S, Boubenec Y)
%
%%% INPUTS
%  - Din ("in" region) and Dout ("out" region) should  have the following format: 
%    n_timepoints x n_stims x n_repetitions x n_voxels
%  - cca_params specify parameters for the analysis
%    Example : 
%   cca_params.recenter = 1;
%   cca_params.pc2Keep = 250;
%   cca_params.cc2Remove = 20;
%   cca_params.baseline_tps = 1:7;
%
%%% OUTPUT
% - Din_corr : denoised version of Din

% If in and out are not specified, both will be equally centered or
% not.
if or(~isfield(cca_params,'recenterIn'),~isfield(cca_params,'recenterOut'))
    disp(['Recenter both ' cca_params.recenter])
    cca_params.recenterIn = cca_params.recenter;
    cca_params.recenterOut = cca_params.recenter;
end

[n_tps, n_stims, n_reps, n_voxels_in] = size(Din);
n_voxels_out = size(Dout,4);

% get 'in' data
if cca_params.recenterIn
    Din_r = Din-nanmean(Din,3);
else
    Din_r = Din;
end

% get 'out' data
if cca_params.recenterOut
    Dout_r = Dout-nanmean(Dout,3);
else
    Dout_r = Dout;
end
muOutBaseline = snm(Dout(cca_params.baseline_tps,:,:,:),[1 2 3]);

% Reshape data in 2D Matrix
Din = reshape(Din, [n_tps*n_stims*n_reps, n_voxels_in]);
Din_r = reshape(Din_r, [n_tps*n_stims*n_reps,n_voxels_in]);
Dout = reshape(Dout, [n_tps*n_stims*n_reps,n_voxels_out]);
Dout_r = reshape(Dout_r, [n_tps*n_stims*n_reps,n_voxels_out]);

% Remove missing time points if any
tmptsNaN = any(isnan(Din),2);
if any(tmptsNaN)
    disp('Removing missing time points')
    Din = Din(~tmptsNaN,:);
    Din_r = Din_r(~tmptsNaN,:);
    Dout = Dout(~tmptsNaN,:);
    Dout_r = Dout_r(~tmptsNaN,:);
end

% Remove missing voxels in out if any
voxNaN = any(isnan(Dout));
if any(voxNaN)
    disp('Removing NaN voxels in out')
    Dout = Dout(:,~voxNaN);
    Dout_r = Dout_r(:,~voxNaN);
    muOutBaseline = muOutBaseline(~voxNaN);
end

% gets SVDs
muIn = mean(Din_r); % should be zero if recentered
muOut = mean(Dout_r); % should be zero if recentered
[Uin,~,Vin] = svd(Din_r-muIn, 'econ');
[Uout,~,Vout] = svd(Dout_r-muOut, 'econ');

% get reduced, whitened data
Din_w = Uin(:,1:cca_params.pc2Keep)*Vin(:,1:cca_params.pc2Keep)';
Dout_w = Uout(:,1:cca_params.pc2Keep)*Vout(:,1:cca_params.pc2Keep)';

% perform cca
Dcat = [Din_w,Dout_w];
[Ucc_r,~,~] = svd(Dcat,'econ');

% get projection of the non-centered 'out' data
B = pinv(Dout_r-muOut)*Ucc_r; % get weights for 'out' with centered data

Ucc_out = (Dout-muOutBaseline')*B; % reproject non-centered data to get CC timecourses
W_out = pinv(Ucc_r(:,1:cca_params.cc2Remove))*Dout_r;

% denoise original data
W_in = pinv(Ucc_r(:,1:cca_params.cc2Remove))*Din_r;
Din_corr = Din - Ucc_out(:,1:cca_params.cc2Remove)*W_in;

if isfield(cca_params,'Mask_slice')
    
    prm.ManualMask = cca_params.Mask_slice == 0;
    W2 = Pixs2Mat(W_out',prm);
    
    prm.ManualMask = cca_params.Mask_slice == 1;
    W3 = Pixs2Mat(W_in',prm);
    % show weights of removed CCs for both in and out
    imtool3D(nansum(cat(4,W2,W3),4))
    
end

% reshape it
inROI_corr = nan(n_tps*n_stims*n_reps,n_voxels_in);
inROI_corr(~tmptsNaN,:) = Din_corr;
inROI_corr = reshape(inROI_corr,[n_tps,n_stims,n_reps,n_voxels_in]);

Din_corr = inROI_corr;
end