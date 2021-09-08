%%% Plots main figures (and some sup) from: 
% Distinct higher-order representations of natural sounds in human and ferret
% auditory cortex (2021)
% Landemard A, Bimbard C, Demené C, Shamma S, Norman-Haigneré S, Boubenec Y
% 2021

% Paths need to be defined in Starter.m first
Starter 

% denoised version that will be used to plot figures
version_name = 'myVersion';
NbICs = 8;

%% Plot example single voxels
% Raw and denoised time-courses for example single voxels
% Nat vs Synth and test-retest time-averaged responses
% Viewing of voxels' spatial position

open PlotSingleVoxels % Fig 1C, Fig 2A-D

%% NSE natural vs synthetic for exp I, ferrets and humans 

% Top view maps, NSE and difference
open PlotMapsFerret % Fig 2E, Fig 4C
open PlotMapsHuman % Fig 2E, Fig 4C

% NSE as function of distance to center of PAC
open PlotNSEvsDis2PAC_Ferret % Fig 2F, S5
open PlotNSEvsDis2PAC_Human % Fig 2F, S5

% Quantification of NSE vs distance to center of PAC curve slopes and
% compare between ferrets and humans
open CompareFerretsHumansSlopes % Fig 2G, Fig 4F

%% Plot example components for exp I

open PlotFerretComponents % Fig 3, Fig S6

%% Plot vocalization results , exp I & II (Fig 4)

% for exp II (panels D,G)
% NSE maps
% Maps of difference nat - synth
% NSE as function to center of PAC, by category of sounds
open PlotVocalizationsResults % Fig 4A-E


%% Other analyses shown in supplementary figures

% Fig S11
% Quantifies movement as recorded by video camera for different sounds
open PlotMvt

% Fig S9 & S10
% Predict data across species
open CrossSpeciesPredictions

% Fig S8
% Plot human components
open PlotHumanComponents

% fig S2
% Quantify effect of denoising
open MeasureDenoisingEffect


%% Implementation of denoising

% Denoising procedure
% to implement full denoising procedure: 
% this will lead to processed data used for figures (as in
% natural_myVersion folder)
open FullDenoising

% Code implementing CCA correction :
open CCACorrection

% Simulation for CCA correction
open SimulationCCA

%% Notice on how to get voxels in their original spatial organisation 
% from raw or denoised data
open GetOriginalVoxels 
