function nse = NSE_noise_corrected_adapted(x1,x2,y1,dim)
% case when only one rep for y 

if nargin < 4
    dim = 1;
    x1 = x1(:);
    x2 = x2(:);
    y1 = y1(:);
    
end
    
    a = 0.5 * (nanmean(x1.^2,dim) + nanmean(x2.^2,dim) - nanmean((x1-x2).^2,dim));
    b = nanmean(y1.^2,dim) - 0.5*nanmean((x1-x2).^2,dim);
    c = 0.5 * (nanmean(x1.*y1,dim)+nanmean(x2.*y1,dim));
    d = 0.5 * (nanmean(x1,dim)+nanmean(x2,dim));
    e = nanmean(y1,dim);
   
    nse = (a+b-2*c)./(a+b-2*d.*e);

end