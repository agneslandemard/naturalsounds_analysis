function param = LoadGlobals_NSD(param,plt)
% Defines category and plotting parameters
if nargin < 2
    plt = 0;
end 

switch param.exp.Exp
    
    case 'AllMM' % experiment I
        param.snd.Categories = {'Ferret vocalizations','English/french speech','Foreign speech','Vocal music',...
            'Instrumental music','Human non vocal', 'Animal non vocal','Non-speech vocal','Mechanical sounds','Environmental sounds'};
        param.snd.CategoryLabels = [1 1 1 1 10 8 6 10 6 2 5 5 4 5 6 5 6 9 5 9 5 6 2 6 8 5 3 3 3 3 3 3 6 5 10 10 7 9 4 7];
        colors = [255 103 218; 0 103 78; 48 191 158; 0 152 211; 14 77 149; 201 51 30; 255 103 131; 123 65 149; 255 139 21; 118 118 118]/255;
        param.snd.SpeechMusicLabels = logical((param.snd.CategoryLabels>=2).*(param.snd.CategoryLabels<=5));
        param.snd.FerretLabels = logical(param.snd.CategoryLabels==1);
        param.snd.SpeechLabels=logical((param.snd.CategoryLabels>=2).*(param.snd.CategoryLabels<=3));
        param.snd.MusicLabels=logical((param.snd.CategoryLabels>=4).*(param.snd.CategoryLabels<=5));
        
    case 'VocalisationsSpM' % experiment II
        
        param.snd.Categories = {'Fights','Kit','Kits','Speech','Music'};
        param.snd.CategoryLabels = [5*ones(1,2) ones(1,5) 2*ones(1,17) 3*ones(1,8) 5*ones(1,12) 4*ones(1,14) 5*ones(1,2) ];
        
        clrs = cbrewer('seq','RdPu',3);
        colors = [clrs(end:-1:1,:);0.1882    0.7490    0.6196 ; 0.0549    0.3020    0.5843];
        
        param.snd.SpeechMusicLabels = logical(param.snd.CategoryLabels>=4);
        param.snd.FerretLabels = logical(param.snd.CategoryLabels<=3);
        param.snd.SpeechLabels = logical(param.snd.CategoryLabels==4);
        param.snd.MusicLabels = logical(param.snd.CategoryLabels==5);
        
end

for k = 1:param.snd.Nmodels
    param.snd.idxSound(:,k) = (1:param.snd.Nsound)+(k-1)*param.snd.Nsound;
end

if isfield(param.snd,'CategoryLabels')
    color_list = nan(param.snd.Nsound,3);
    for i = 1:param.snd.Nsound
        color_list(i,:) = colors(param.snd.CategoryLabels(i),:);
    end
    
    param.plt.SoundColors = repmat(color_list,[param.snd.Nmodels,1]);
end

param.plt.CT = cbrewer('div', 'RdBu', 64);
param.plt.CT = param.plt.CT(end:-1:1,:);


% If specified, plots legend for category colors
if plt == 1
     [~,ord]=sort(param.snd.CategoryLabels,'descend');
     
     figure; 
     scatter(ones(param.snd.Nsound,1),1:param.snd.Nsound,[],param.plt.SoundColors(param.snd.idxSound(ord,1),:),'filled')
     text(1.2*ones(param.snd.Nsound,1),1:param.snd.Nsound,regexprep(param.snd.SoundList(ord),'_',' '),'filled')
     
     hold all
     for k = 1:length(param.snd.Categories)
         scatter(1,0,[],colors(k,:),'filled','DisplayName',param.snd.Categories{k})
         hold on
     end
     legend('show','location','northwest')
end
