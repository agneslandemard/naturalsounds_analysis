function PlotTopView(Im,param,Color)

global additional_path

if nargin < 3
    Color.ColorAxis=[0 1];
end

if ~isfield(Color,'cm')
    Color.cm = cbrewer('div','Spectral',64,'pchip');
    Color.cm = Color.cm(end:-1:1,:);
end

% Load coordinates for top view
load([additional_path 'Coordinates/TopView_Coordinates_' param.exp.Animal '_' param.exp.Session '.mat'],'CoordinatesAll');

color_range = [-Inf linspace(Color.ColorAxis(1),Color.ColorAxis(end),size(Color.cm,1)-1) Inf];

% assign color to all values
ns_colors = zeros(length(Im),3);
for x = 1:length(Im)
    for n = 1:length(color_range)-1
        if Im(x)>=color_range(n) && Im(x)<color_range(n+1)
            ns_colors(x,:) = Color.cm(n,:);
            break;
        end
    end
end

% Compute top view
ns = Pixs2Mat(ns_colors, param.msk);

[n_slices, n_layers, n_points,~] = size(CoordinatesAll);
MapFromAbove = nan([n_slices, n_layers, n_points, 3]);
for ii = 1:n_slices
    for lay = 1:n_layers
        for pt = 1:n_points
            if any(squeeze(ns(CoordinatesAll(ii,lay,pt,1),CoordinatesAll(ii,lay,pt,2),ii,:))' ~= [0 0 0])
                MapFromAbove(ii,lay,pt,:) = ns(CoordinatesAll(ii,lay,pt,1),CoordinatesAll(ii,lay,pt,2),ii,:);
            end
        end
    end
end


% Plot top view
h = imagesc(rot90(squeeze(nanmean(MapFromAbove,2)),3));
h.AlphaData = rot90(~isnan(squeeze(nanmean(MapFromAbove,[2 4]))),3);
set(gca,'DataAspectRatio',[7 40 1])
set(gca,'XTick',[],'YTick',[])
colormap(Color.cm)

end
