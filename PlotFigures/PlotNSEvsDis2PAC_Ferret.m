%%%  Revisions_eLife
%% Compute NSE maps and NSE as function of distance to PAC (Fig 2E&F)

% to vary time window
P.RespWindow = 3:11;

% Load data
P.experiment = 'natural';
P.n_ics = NbICs;
P.name = version_name;
D = LoadDenoisedData(P);

%% NSE nat vs synth + NSE test-retest (without noise-correction)

load([additional_path '/Coordinates/distances_to_pac_' P.experiment '.mat'])

figure('Position',[680 829 606 149]);
n_conds = D.param.snd.Nmodels;
local_nse = nan(n_conds,length(distances),D.n_hemis);
local_nse_trt = nan(n_conds,length(distances),D.n_hemis);

% group by ferret
si_byf = nan(size(D.si));
si_byf(D.si <= 2)=1;
si_byf(D.si >= 3)=2;

for md = 2 : n_conds
    
    % NSE between selected model and natural 
    nsm = 0.5*(NSE(D.SrecoFull(D.param.snd.idxSound(:,1),:,1),D.SrecoFull(D.param.snd.idxSound(:,md),:,2),1) + ...
        NSE(D.SrecoFull(D.param.snd.idxSound(:,1),:,2),D.SrecoFull(D.param.snd.idxSound(:,md),:,1),1));
    
    nsm_trt = NSE(D.SrecoFull(D.param.snd.idxSound(:,[1]),:,1),D.SrecoFull(D.param.snd.idxSound(:,[1]),:,2),1);
   
    for j = 1 : 2
        % pick out one subject
       % xi = D.si == j;
        xi = si_byf == j;
        nse_single_subj = nsm(xi);
        nse_trt_single_subj = nsm_trt(xi);
      
       for d = 1 : length(distances)
            local_nse(md,d,j) = nanmedian(nse_single_subj(distance_to_pac(xi)==distances(d)));
            local_nse_trt(md,d,j) = nanmedian(nse_trt_single_subj(distance_to_pac(xi)==distances(d)));
       
        end
        
    end

    subplot(1,n_conds-1,md-1)
    % ferrets
    hold all
    plot(distances(2:end)./10,squeeze(local_nse(md,1:end-1,:)),'color',0.6.*[1 1 1])
    plot(distances(2:end)./10,snm(local_nse(md,1:end-1,:),3),'k','LineWidth',2)
    
    plot(distances(2:end)./10,squeeze(local_nse_trt(md,1:end-1,:)),'--','color',0.6.*[1 1 1])
    plot(distances(2:end)./10,snm(local_nse_trt(md,1:end-1,:),3),'k--','LineWidth',2)
    

    if md == 2
        xlabel('Distance from center of PAC (mm)')
        ylabel('NSE')
    end
    
    ylim([0 1])
    yticks(0:0.2:1)
    xticks(distances(1:2:end-1)./10)
    xlim([0 distances(end-1)./10])
    
     
end

%% Noise-corrected NSE

figure('Position',[680 829 606 149]);
n_conds = D.param.snd.Nmodels;
local_nse = nan(n_conds,length(distances),D.n_hemis);
local_nse_h = nan(n_conds,length(distances),D.n_hemis);

for md = 2 : n_conds
    
    % NSE between selected model and natural 
    nsm = NSE_noise_corrected(D.SrecoFull(D.param.snd.idxSound(:,[1]),:,1),D.SrecoFull(D.param.snd.idxSound(:,[1]),:,2),...
        D.SrecoFull(D.param.snd.idxSound(:,md),:,1),D.SrecoFull(D.param.snd.idxSound(:,md),:,2),1);
   
    for j = 1 : 2
        % pick out one subject
        xi = si_byf == j;
        nse_single_subj = nsm(xi);
       
       for d = 1 : length(distances)
            local_nse(md,d,j) = nanmedian(nse_single_subj(distance_to_pac(xi)==distances(d)));
            
        end
        
    end

    subplot(1,n_conds-1,md-1)
    % ferrets
    hold all
    plot(distances(2:end)./10,squeeze(local_nse(md,1:end-1,:)),'color',0.6.*[1 1 1])
    plot(distances(2:end)./10,snm(local_nse(md,1:end-1,:),3),'k','LineWidth',2)
   
    if md == 2
        xlabel('Distance from center of PAC (mm)')
        ylabel('NSE')
    end
   
    ylim([0 1.1])
    yticks(0:0.2:1)
    xticks(distances(1:2:end-1)./10)
    
    xticklabels(distances(1:2:end-1)./10)
    xlim([0 distances(end-1)./10])
    

end
