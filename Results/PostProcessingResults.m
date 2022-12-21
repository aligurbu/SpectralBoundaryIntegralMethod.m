%% Post-processing the flow around the RBC motion
clear all; close all; clc;
addpath(genpath('../../SpectralBoundaryIntegralMethod.m'))
Starttime = tic;

%% Input the model and parameters for the analysis from Models folder
LoadElasticRBC_Shear_N16
% LoadElasticRBC_Parabolic_N16
% LoadMemViscosityRBC_Shear_N16
% LoadMemViscosityRBC_Parabolic_N16

%% Set up output file
fidCoord = fopen(['Coord_',name,'.dat'],'r');
fidMemFor = fopen(['MemFor_',name,'.dat'],'r');
fidSol = fopen(['Sol_',name,'.dat'],'r');
fidTime = fopen(['Time_',name,'.dat'],'r');
Time = fread(fidTime,'double');
NSTEPS = length(Time);

%%
timeStepIncrement = 1;
FRAMES = floor(linspace(1,NSTEPS,5));
TimesOfFrames = Time(FRAMES)/RefShearRate % in seconds

%% Choose a view angle
viewInd = [0 90];

%% Choose a background color for visualization
blackBackground = true; % if false then white background

%% Visualization settings
VisualizeSettings

%% Post-processing the flow around the RBC motion


%%
fclose('all');
RunTime = toc(Starttime)