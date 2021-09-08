%% Compute similarity between raw and denoised

function ScatterRawDenoised(resp_den,resp_raw,metric)
% resp_den and resp_raw should be formatted stims * voxels * reps(should be 2)
% metric should be 'corr' or 'nse'

if nargin < 3 
    metric = 'corr';
end

switch metric
    case 'corr'
        % correlation
        TR_raw = nan(1,size(resp_den ,2));
        TR_denraw = nan(1,size(resp_den ,2));
        TR_den = nan(1,size(resp_den ,2));
        for pix = 1:size(resp_den ,2)
            TR_denraw(pix) = (1/2)*(corr(resp_den(:,pix,1),resp_raw(:,pix,2),'rows','pairwise')+corr(resp_den(:,pix,2),resp_raw(:,pix,1),'rows','pairwise'));
            TR_den(pix) = corr(resp_den(:,pix,1),resp_den(:,pix,2),'rows','pairwise');
            TR_raw(pix) = corr(resp_raw(:,pix,1),resp_raw(:,pix,2),'rows','pairwise');
        end
    case 'nse'
        % same with nse
        TR_denraw = (1/2)*(NSE(resp_den(:,:,1),resp_raw(:,:,2),1)+NSE(resp_den(:,:,2),resp_raw(:,:,1),1));
        TR_den = NSE(resp_den(:,:,1),resp_den(:,:,2),1);
        TR_raw = NSE(resp_raw(:,:,1),resp_raw(:,:,2),1);
end

%% plot scatter reliability denoised / original

n_bins = 100;

figure('Position',[440 498 652 300]);

limits = [min([TR_denraw TR_den TR_raw]) max([TR_denraw TR_den TR_raw])];

subplot(121)
cm = cbrewer('seq','Greens',100);
cm = [cm(5:100,:); repmat(cm(100,:),30,1)];

hold all
binscatter(TR_raw,TR_den,n_bins)
ylim(limits)
xlim(limits)
colormap(gca,cm)
axis equal tight
l = line(limits,limits);
l.Color = 'k';
l.LineStyle = '--';
l.LineWidth = .5;
c = caxis;

vline(0,'k:')
hline(0,'k:')

% Indicate median across voxels
scatter(median(TR_raw),median(TR_den),60,'r*')
colorbar

xlabel('Split half reliability (Pearson Correlation)')
ylabel('Denoised reliability (Pearson Correlation)')
title('Denoised two splits')
disp(['Median across voxels for denoised two splits: ' num2str(median(TR_den(:)),2)])


subplot(122)
hold all
cm = cbrewer('seq','Blues',100);
cm = [cm(5:100,:); repmat(cm(100,:),30,1)];
binscatter(TR_raw,TR_denraw,n_bins)
ylim(limits)
xlim(limits)

colormap(gca,cm)
axis equal tight
l = line(limits,limits);
l.Color = 'k';
l.LineStyle = '--';
l.LineWidth = .5;
caxis(c);
colorbar

% Indicate median across voxels
scatter(median(TR_raw),median(TR_denraw),60,'r*')

xlabel('Split half reliability (Pearson Correlation)')
ylabel('Denoised reliability (Pearson Correlation)')
title('Denoised one split')
disp(['Median across voxels for denoised one split: ' num2str(median(TR_denraw),2)])


% Draw ceiling (sqrt of raw test-retest)
disp(['Median across voxels for ceiling: ' num2str(median(real(sqrt(TR_raw))),2)])
[~,edges] = histcounts(TR_raw);
bound = real(sqrt(edges));
bound(bound == 0) = nan;
plot(edges(:),bound,'k-','LineWidth',1);

vline(0,'k:')
hline(0,'k:')

end
