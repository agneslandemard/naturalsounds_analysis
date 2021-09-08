%% Plot Human Components

%% Load data 

h = load([additional_path 'Human/HumanComponentsforPrediction.mat']);

P.exp.Exp = 'AllMM';
P = LoadGlobalsNSD(P);

sound_colors = permute(reshape(P.plt.SoundColors, [P.snd.Nsound, P.snd.Nmodels, 3]), [3, 1, 2]);
symbols = {'o','x', 'square', 'diamond', '+'};

% load cortical model for sounds
CM = load([additional_path 'CorticalModel/CM_natural.mat'],'mods','freqs','freqs_list','specmods_list','tempmods_list','stims','conditions');

n_ics = size(h.human_R,3);

%% Plot figure 3 B-E

figure;

for ic = 1:n_ics
    
    % Select component
    R = h.human_R(:,:,ic);
     
    % Responses to natural vs synthetic sounds
    subplot(n_ics,4,(ic-1)*4+1)
    bounds = [min(R(:)), max(R(:))];
    bounds = bounds + 0.1*[-1 1]*diff(bounds);
    hold on
    for k = 1:P.snd.Nsound
        scatter(R(k,1), R(k,P.snd.Nmodels),10,sound_colors(:, k, P.snd.Nmodels)','filled');
        xlim(bounds); ylim(bounds);
    end
    xlabel('Natural Sounds');
    ylabel('spectrotemporal');
    title(['NSE = ' num2str(NSE(R(:,1), R(:,P.snd.Nmodels)),2)])
    plot(bounds, bounds, 'r:', 'LineWidth', 0.3)
    axis equal tight
     
    
    % Test-retest 
    subplot(n_ics,4,(ic-1)*4+2)
    hold all
    for md = 1:P.snd.Nmodels-1:P.snd.Nmodels
        cl = sound_colors(:,:,md);
        scatter(mat2vec(human_R1(:,md,ic)),mat2vec(human_R2(:,md,ic)),10,cl(:,:)',symbols{md},'LineWidth',0.2);
    end
    mi = min([mat2vec(human_R1(:,:,ic));mat2vec(human_R2(:,:,ic))]);
    ma = max([mat2vec(human_R1(:,:,ic)); mat2vec(human_R2(:,:,ic))]);
    line([mi ma],[mi ma],'LineStyle',':','Color','k','LineWidth',0.3)
    title(['NSE = ' num2str(NSE(human_R1(:,:,ic),mat2vec(human_R2(:,:,ic))),2)])
    axis equal tight
    xlabel('test')
    ylabel('retest')
    
    % Correlation with frequency content
    subplot(n_ics,4,(ic-1)*4+3)
    n_freqbins = length(CM.freqs_list);
    r = corr(CM.freqs,R(:),'rows','pairwise');
    min_freq = 4; % remove lowest frequencies
    plot(r(min_freq:end), 'k-', 'LineWidth', 0.75);
    hold on;
    plot(zeros(n_freqbins-min_freq+1,1), 'r:', 'LineWidth', 0.3);
    bounds = [min(r(min_freq:end)) max(r(min_freq:end))];
    bounds = bounds.*(1/diff(bounds));
    ylim(bounds);
    xlim([1, n_freqbins-min_freq+1]);
    set(gca, 'XTick', 1:2:n_freqbins-min_freq+1, 'XTickLabel', round(CM.freqs_list(min_freq:2:n_freqbins)));
    xlabel('Frequency'); ylabel('Correlation');
    set(gca,'DataAspectRatio',[4 1 1])
    
    % Correlation with modulation content
    subplot(n_ics,4,(ic-1)*4+4)
    n_spec = length(CM.specmods_list);
    n_temp = length(CM.tempmods_list);
    r = corr(CM.mods(:,:),R(:),'rows','pairwise');
    r = reshape(r,n_spec,n_temp);
    corr_scale = cmap_from_name('cbrewer-blue-red');
    imagesc(flipud(r), 0.7*[-1, 1]);
    axis equal tight
    colormap(corr_scale);
    flip_yvals = fliplr(CM.specmods_list);
    set(gca, ...
        'YTick', 2:2:n_spec, ...
        'XTick', 2:2:n_temp,...
        'YTickLabel', flip_yvals(2:2:n_spec), ...
    'XTickLabel', CM.tempmods_list(2:2:n_temp));
    xlabel('Temp (Hz)'); ylabel('Spec (cyc/oct)');
 
    
end

