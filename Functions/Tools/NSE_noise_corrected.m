function [nse,num,denom]= NSE_noise_corrected(x1,x2,y1,y2,dim)

if nargin<5
    dim=1;
    x1=x1(:);
    x2=x2(:);
    y1=y1(:);
    y2=y2(:);
end
    
    a=0.5*(nanmean(x1.^2,dim)+nanmean(x2.^2,dim)-nanmean((x1-x2).^2,dim));
    b=0.5*(nanmean(y1.^2,dim)+nanmean(y2.^2,dim)-nanmean((y1-y2).^2,dim));
    c=0.25*(nanmean(x1.*y1,dim)+nanmean(x1.*y2,dim)+nanmean(x2.*y1,dim)+nanmean(x2.*y2,dim));
    d=0.5*(nanmean(x1,dim)+nanmean(x2,dim));
    e=0.5*(nanmean(y1,dim)+nanmean(y2,dim));
   
    nse=(a+b-2*c)./(a+b-2*d.*e);
%     disp(a)
%     disp(b)
%     disp(['num = ' num2str((a+b-2*c))])
%     disp(['denom = ' num2str((a+b-2*d*e))])
    num = (a+b-2*c);
    denom = (a+b-2*d.*e);

end