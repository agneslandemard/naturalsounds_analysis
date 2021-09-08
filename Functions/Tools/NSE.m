function nse = NSE(x,y,dim)

if nargin == 3 
    nse=nanmean((x-y).^2,dim)./(nanmean(x.^2,dim)+nanmean(y.^2,dim)-2*nanmean(x,dim).*nanmean(y,dim));
else
    x(isnan(x)==1)=[];
    y(isnan(y)==1)=[];
    x=x(:);
    y=y(:);
    if ~isequal(size(x),size(y))
        error('Vectors should be the same size')
    else
        nse=nanmean((x-y).^2)/(nanmean(x.^2)+nanmean(y.^2)-2*nanmean(x)*nanmean(y));
    end
end

end