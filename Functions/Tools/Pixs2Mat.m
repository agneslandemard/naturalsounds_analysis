function [IN_mat]=Pixs2Mat(OUT_mat,param)
% Projects back first dim of input matrix on 2 (or 3)D voxel space
% Reverse operation is Mat2Pixs;

s = size(OUT_mat);
ms = size(param.ManualMask);
param.PixInds = find(param.ManualMask(:));

IN_mat = nan([prod(ms) s(2:end)]);

IN_mat(param.PixInds,:) = OUT_mat(:,:);

IN_mat = reshape(IN_mat,[ms s(2:end)]);

end