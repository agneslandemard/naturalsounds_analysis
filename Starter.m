%% Starter
% Source code for analysis described in:
% Distinct higher-order representations of natural sounds in human and ferret
% auditory cortex (2021)
% Landemard A, Bimbard C, Demené C, Shamma S, Norman-Haigneré S, Boubenec Y
% 2021, Sept. 8

%% Define paths 

global data_path additional_path analysis_path code_path

% TO ADAPT before starting
project_directory = 'naturalsoundsdata/'; % define data folder location
code_path = 'naturalsoundsanalysis/'; % define code folder location

data_path = [project_directory 'fUSData/'];
additional_path = [project_directory 'AdditionalData/'];
analysis_path = [project_directory 'Analysis/'];

addpath(genpath(code_path));
