function [data,pltparam] =  LoadComponents(load_params)

global analysis_path data_path

if ~isfield(load_params,'n_ics')
    load_params.n_ics = 8;
end
comps_path = [analysis_path load_params.experiment '_' load_params.name];
data = load([comps_path '/ICA_results' num2str(load_params.n_ics) '.mat']);

tmp = load([comps_path '/D_norm.mat'],'hemis');
data.hemis = tmp.hemis;
switch load_params.experiment
    case 'natural'
        
        data.conditions = {'originals', 'cochlear', 'tempmod', 'specmod', 'spectempmod'};
        data.n_stims_per_cond = 40;
        data.n_conds = 5;
        
        pltparam.symbols = {'square','o', 'x', '+', 'v'};
        pltparam.fill=0;

    case 'vocalization'
        data.conditions = {'originals', 'spectempmod'};
        data.n_stims_per_cond = 60;
        data.n_conds = 2;
        
        pltparam.symbols = {'o','o'};
        pltparam.fill=1;   
end

data.experiment = load_params.experiment;
data.data_suffix='_Anat&Param';
data.n_stims = data.n_stims_per_cond*data.n_conds;
data.n_tps = size(data.R_ica_train_tc,1);
data.data_path = comps_path;
data.n_ics = load_params.n_ics;

data.R_ica_all_tc = reshape(data.R_ica_train_tc,data.n_tps,  data.n_stims_per_cond, data.n_conds,load_params.n_ics);

X = load([data_path data.hemis{1} data.data_suffix '.mat'],'param');
data.R_ica_all = reshape(snm(data.R_ica_all_tc(X.param.exp.PreStimSilence+X.param.proc.MeanInterval,:,:,:),1), data.n_stims_per_cond, data.n_conds,load_params.n_ics);
data.R_ica_train = reshape(snm(data.R_ica_train_tc(X.param.exp.PreStimSilence+X.param.proc.MeanInterval,:,:,:),1), size(data.R_ica_train_tc,2)/data.n_conds, data.n_conds,load_params.n_ics);
data.R_ica_test_even = reshape(snm(data.R_ica_test_even_tc(X.param.exp.PreStimSilence+X.param.proc.MeanInterval,:,:,:),1), size(data.R_ica_test_even_tc,2)/data.n_conds, data.n_conds,load_params.n_ics);
data.R_ica_test_odd = reshape(snm(data.R_ica_test_odd_tc(X.param.exp.PreStimSilence+X.param.proc.MeanInterval,:,:,:),1),size(data.R_ica_test_odd_tc,2)/data.n_conds, data.n_conds,load_params.n_ics);

end






