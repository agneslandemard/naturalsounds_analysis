%% GUIDE on how to get original organisation of voxels based on linearized version.
% Two examples are provided for raw data (no previous steps needed) or for
% denoised volumes (runs after FullDenoising.m)
% Note that in both cases, only voxels inside the cortex ROI are
% considered. The rest of the image will be filled with NaNs.
% Voxels in 1D will be replaced by 2 or 3D: image dimension 1 (cortical depth, from
% dorso-lateral to ventro-medial), image dimension 2 (from dorso-medial to
% ventro-latera), slices.

% This process is also done by function Pixs2Mat, and the reverse by Mat2Pixs

%% For raw data 

hemi_id = 'NTR';
slice = 1; 
load([data_path hemi_id '_' num2str(slice) '.mat'] ,'Din');
load([data_path hemi_id '_Anat&Param.mat'] ,'Anat','param');
        
Din = permute(Din,[4 1 2 3]); 
s = size(Din);

% select cortex ROI for this slice
slice_roi = param.msk.ManualMask(:,:,slice);
[x,y] = size(slice_roi); % get image dimensions
voxels_in_roi = find(slice_roi(:)); 

% Reorganize voxels in ROI according to original organisation
Din_orig = nan([x*y s(2:end)]);
Din_orig(voxels_in_roi,:) = Din(:,:);
Din_orig = reshape(Din_orig,[x y s(2:end)]);

% Plot side to side full anatomical image and anatomical version of 
% newly-formatted matrix (square root of mean across all timepoints).
% Note that this one only contains voxels inside the cortex ROI.
figure;
imagesc([Anat(:,:,slice) sqrt(mean(Din_orig(:,:,:),3))])
axis equal tight; colormap hot
xlabel('dorso-medial --> ventro-lateral')
ylabel('ventro-medial --> dorso-lateral (in depth)')

%% For denoised data

% load denoised data
P.experiment = 'natural';
P.version = 'myVersion';
P.n_comps = 8;
D = LoadDenoisedData(P);

% example for one hemisphere
hemi = 1;
Dc = permute(D.Sreco(:,:,si == hemi),[3 1 2]); 
s = size(Dc);

% load appropriate mask
load([data_path hemis{hemi} '_Anat&Param.mat'],'param');

% select cortex ROI for all slices
all_slices_roi = param.msk.ManualMask;
[x,y,n_slices] = size(all_slices_roi); % get image dimensions
voxels_in_roi = find(all_slices_roi(:)); 

Dc_or = nan([x*y*n_slices s(2:end)]);
Dc_or(voxels_in_roi,:) = Dc(:,:);
Dc_or = reshape(Dc_or,[x y n_slices s(2:end)]);
% now 3 first dimensions of Dc_or are x (image dimension 1), y (image dimension 2), slices. 