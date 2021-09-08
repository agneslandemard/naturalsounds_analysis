%%% Plot Human Maps:
% - NSE nat vs synth
% - Diff nat - synth

% get grid properties for plotting
load([additional_path '/Human/HumanGrid.mat'],'grid_info')

%% NSE maps

% colormap properties
clear Color
Color.ColorAxis = [0 1];
Color.cm = cmap_from_name('lightblue-to-yellow1');

toPlot = grid_info;

hemi = 2; %1=RH, 2=LH

stk = 1:12; % subjects from paradigm I only

Stmp = SrecoFull_human_allSubjects;
HD = reshape(Stmp,200,xl*yl+xr*yr,n_subj,2);

% 1) average across subjects
Split1 = snm(HD(:,:,stk,1),3);
Split2 = snm(HD(:,:,stk,2),3);
nse_trt = NSE(Split1(D.param.snd.idxSound(:,1),:) ,Split2(D.param.snd.idxSound(:,1),:),1);

for md = 2:D.param.snd.Nmodels
    
    % 2) compute NSE
    diff_all = NSE_noise_corrected_adapted(Split1(D.param.snd.idxSound(:,1),:),...
        Split2(D.param.snd.idxSound(:,1),:),Split1(D.param.snd.idxSound(:,md),:),1);
    
    % 3) voxel selection
    diff_all(nse_trt > 0.4)=nan;
    
    % prepare plotting
    D_LH = reshape(diff_all(1:xl*yl),xl,yl);
    D_RH = reshape(diff_all(xl*yl+1:end),xr,yr);
    
    toPlot.grid_data{1} = D_RH;
    toPlot.grid_data{2} = D_LH;
    
    surftoPlot = grid2surface(toPlot);
 
    plot_fsaverage_1D_overlay(surftoPlot(:,hemi),'lh',Color.cm ,Color.ColorAxis);
    set(gcf,'Position',[100+(md-1)*300 487 310 303])
    set(gcf, 'Renderer', 'opengl');
   
end   

%% Normalized difference nat vs synth
% colormap properties
clear Color
m = 2;
Color.ColorAxis = m*[-1 1];
Color.cm = cmap_from_name('cbrewer-blue-red');

% get properties for plotting
toPlot = grid_info;
stk = [2 8 11 13:17]; % subjects with multiple reps
SoundType = {'Speech','Music','Others'};
hemi = 2; %1=RH, 2=LH
plotWhat = 'diff_resp';

% 1) average across subjects
Stmp = SrecoFull_human_allSubjects;
Stmp = reshape(Stmp,200,xl*yl+xr*yr,n_subj,2);
Stmp  = snm(Stmp(:,:,stk,:),3);

% 2) Compute difference nat vs synth 
[~, diff_all_sds] = NSE_disFromPAC(Stmp,...
    [],D.param,[],SoundType,plotWhat);

% 3) Select voxels
nse_trt = NSE(Stmp(D.param.snd.idxSound(:,1),:,1),Stmp(D.param.snd.idxSound(:,1),:,2),1);
NRV = nse_trt > 0.4;
diff_all_sds(:,NRV) = nan;

for sd = 1:length(SoundType)
    
    SdsIdx = SelectSounds(SoundType{sd},D.param);
    
    % for selected sounds
    diff_all = snm(diff_all_sds(SdsIdx,:),1);
    
    % prepare plotting
    D_LH = reshape(diff_all(1:xl*yl),xl,yl);
    D_RH = reshape(diff_all(xl*yl+1:end),xr,yr);
    
    toPlot.grid_data{1} = D_RH;
    toPlot.grid_data{2} = D_LH;
    
    surftoPlot = grid2surface(toPlot);

    plot_fsaverage_1D_overlay(surftoPlot(:,hemi),'lh',Color.cm ,Color.ColorAxis);
    set(gcf,'Position',[100+(sd-1)*300 487 310 303])
    set(gcf, 'Renderer', 'opengl');

end
