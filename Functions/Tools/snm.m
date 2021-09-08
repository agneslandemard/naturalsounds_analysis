function y = snm(x,dims)
% just nanmean across chosen dimensions and squeeze result

    if nargin < 2
        dims = 1;
    end
    
    y = squeeze(nanmean(x,dims));
end