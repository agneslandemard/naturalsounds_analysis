%% Load ferret data for both experiments
clear P1
P1.experiment = 'natural';
P1.n_ics = NbICs;
P1.name = version_name;
D = LoadDenoisedData(P1);

clear P2
P2.experiment = 'vocalization';
P2.n_ics = NbICs;
P2.name = version_name;
V = LoadDenoisedData(P2);

%% Compute ferret slopes

distf = load([additional_path 'Coordinates/distances_to_pac_natural.mat']);
distf.distances = distf.distances(1:end-1);

% Ferret subjects
animals1 = {'A','T','C'};
SubjIdF = nan(size(D.si));
for h = 1:D.n_hemis
    f = find(strcmp(animals1,D.hemis{h}(2)));
    SubjIdF(D.si==h)=f;
end

% Compute slopes for all sounds, in exp I
SoundType = {'All'};

nse_by_dis_all_f = NSE_disFromPAC(D.SrecoFull,SubjIdF,D.param,distf,SoundType);
slopes_all_f = distf.distances'\squeeze(nse_by_dis_all_f(:,:,:,1));

% Compute slopes for categories of sounds, in exp II

SubjIdF2 = nan(size(V.si));
for h = 1:V.n_hemis
    f = find(strcmp(animals1,V.hemis{h}(2)));
    SubjIdF2(V.si==h)=f;
end

SoundType = {'Speech','Ferrets'};

distf2 = load([additional_path 'Coordinates/distances_to_pac_vocalization.mat']);
distf2.distances = distf2.distances(1:end-1);
nse_by_dis = NSE_disFromPAC(V.SrecoFull,SubjIdF2,V.param,distf2,SoundType);

slopes_speech_f = distf2.distances'\squeeze(nse_by_dis(1,:,:,1));
slopes_vocs_f = distf2.distances'\squeeze(nse_by_dis(2,:,:,1));

%% Compute human slopes

dist = load([additional_path 'Coordinates/distances_to_pac_human.mat']);
HD = reshape(SrecoFull_human_allSubjects,200,xl*yl+xr*yr,n_subj,2);

% remove unreliable voxels
nse_trt = NSE(HD(D.param.snd.idxSound(:,1),:,:,1),HD(D.param.snd.idxSound(:,1),:,:,2),1);
NRV = nse_trt > 0.4;
tmpS = SrecoFull_human_allSubjects;
tmpS(:,NRV(:),:) = nan;

SubjId = nan(xl*yl+xr*yr,n_subj);
stk = setdiff(2:17,[4 9 10 12]); % remove subjects without tonotopy
for s = 1:n_subj
    if ismember(s,stk)
        SubjId(:,s)=s;
    end
end
dist.distance_to_pac = dist.distance_to_pac_human(:);
dist.distances = dist.distance_humans(1:end-2);
SoundType = {'Speech'};

nse_by_dis_s = NSE_disFromPAC(tmpS,SubjId(:),D.param,dist,SoundType,'nse_across_sounds_adapted');
slopes_speech_h = dist.distances'\squeeze(nse_by_dis_s(:,:,:,1));

SoundType = {'All'};

nse_by_dis = NSE_disFromPAC(tmpS,SubjId(:),D.param,dist,SoundType,'nse_across_sounds_adapted');
slopes_all_h = dist.distances'\squeeze(nse_by_dis(:,:,:,1));

%% Plot slope values

figure('Position',[680 837 755 141])
hold all
% human slopes
plot(slopes_all_h.*[1 ;1], [1 2], 'color',[232 191 88]./255,'LineWidth',2)
plot(slopes_speech_h.*[1 ;1],[3 4] , 'color',[232 191 88]./255,'LineWidth',2)

% exp I ferret slopes
animals1 = {'A','T'};
flines1 = {'--','-'};
for f = 1:length(animals1)
    plot(slopes_all_f(f).*[1 ;1], [1 2], 'color',0.6.*[1 1 1],'LineWidth',2,'LineStyle',flines1{f})
end  

% exp II ferret slopes
animals2 = {'A','T','C'};
flines2 = {'--','-',':'};
for f = 1:length(animals2)
    plot(slopes_vocs_f(f).*[1 ;1], [3 4], 'color',0.*[1 1 1],'LineWidth',2,'LineStyle',flines2{f})
    plot(slopes_speech_f(f).*[1 ;1], [3 4], 'color',0.3.*[1 1 1],'LineWidth',2,'LineStyle',flines2{f})
end

xlabel('Slopes (∆NSE/mm)')
yticks([1.5 3.5])
yticklabels({'Exp I','Exp II'})


%% Statistical testing
% for exp I, all sounds
clear stats
for f = 1 : 2
    [p(f),~ ,stats(f)] = signtest(slopes_all_h - slopes_all_f(f));
end

% for exp II, speech 
clear stats2
for f = 1:3
   [p(f,1),~, stats2(f,1)] = signtest(slopes_speech_h - slopes_speech_f(f));
   [p(f,2),~, stats2(f,2)] = signtest(slopes_speech_h - slopes_vocs_f(f));
end




