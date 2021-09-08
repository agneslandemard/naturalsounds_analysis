function data = LoadDenoisedData(load_params)

global analysis_path data_path

% Define default parameters
if ~isfield(load_params,'n_comps')
    load_params.n_comps = 8;
end
if ~isfield(load_params,'name')
    load_params.name = 'myVersion';
end
if ~isfield(load_params,'by_ferret')
    load_params.by_ferret = 1;
end

proc_path = [analysis_path load_params.experiment '_' load_params.name ];
data = load([proc_path '/D_norm.mat'],'hemis');

if load_params.by_ferret
    
    data.SrecoFullTC = [];
    data.si = [];
    for h = 1 : length(data.hemis)
        ferret = data.hemis{h}(2);
        % Load and recompose data based on runs
        D = load([proc_path '/' ferret '/D_denoised_K' num2str(load_params.n_comps) '.mat']);
        hI = D.si == find(strcmp(D.hemis,data.hemis{h}));   
       
        data.SrecoFullTC = cat(3,data.SrecoFullTC, D.D_denoised(:,:,hI,:));
        data.si = cat(2,data.si, h*ones(1,length(find(hI)))); % hemisphere identifier
        
    end
   
else
    % Load hemispheres names
    
    
    % Load and recompose data based on runs
    D = load([proc_path '/D_denoised_K' num2str(load_params.n_comps) '.mat']);
    
    data.SrecoFullTC = D.D_denoised;
    data.si = D.si; % hemisphere identifier
end


% Define other parameters
data.experiment = load_params.experiment;
data.n_hemis = length(data.hemis);
data.data_suffix = '_Anat&Param';
X = load([data_path  data.hemis{1} data.data_suffix '.mat'],'param');
data.param = X.param; clear X
data.data_path = proc_path;


% Define response window and time-averaged response
if ~isfield(load_params,'RespWindow')
    load_params.RespWindow = data.param.proc.MeanInterval;
end
data.SrecoTC = snm(data.SrecoFullTC,4);
data.Sreco = snm(data.SrecoTC(data.param.exp.PreStimSilence+load_params.RespWindow,:,:),1);
data.SrecoFull = snm(data.SrecoFullTC(data.param.exp.PreStimSilence+load_params.RespWindow,:,:,:),1);
data.RespWindow = load_params.RespWindow;

end
