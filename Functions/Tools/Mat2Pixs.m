function [OUT_mat]=Mat2Pixs(IN_mat,param)
% Converts the two first dimensions of input matrix into one, with only
% pixels inside of the mask in param.ManualMask 
% Reverse operation is Pixs2Mat

sz=length(size(param.ManualMask));

s=size(IN_mat); 
if length(s)>sz
    OUT_mat=reshape(IN_mat,[prod(s(1:sz)) s(sz+1:end)]);
else
    OUT_mat=reshape(IN_mat,[prod(s(1:sz)) 1]);
end

param.PixInds=find(param.ManualMask(:));
OUT_mat=OUT_mat(param.PixInds,:);
if length(s)>sz
    OUT_mat=reshape(OUT_mat,[length(param.PixInds) s(sz+1:end)]);
else
    OUT_mat=reshape(OUT_mat,[length(param.PixInds) 1]);
end

end

